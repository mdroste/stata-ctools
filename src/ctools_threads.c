/*
 * ctools_threads.c - Thread pool utilities for ctools
 *
 * Provides a persistent thread pool for efficient parallel execution, and
 * the thread-count policy shared by every command.
 */

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE  /* sched_getaffinity, CPU_ALLOC */
#endif
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <pwd.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <sched.h>
#endif
#endif

#include "ctools_threads.h"
#include "ctools_config.h"
#include "ctools_file.h"

/* ===========================================================================
   Unified Thread Management Implementation

   The default comes from ctools_thread_policy_get(); the threads() option
   and the session setting override it for one plugin call. The precedence
   is documented in ctools_config.h.
   =========================================================================== */

/* User-specified thread limit (0 = use default) */
static int g_ctools_max_threads = 0;

/* Resolved default, captured once from the OS and environment.
 * This is immune to pollution from other Stata plugins (reghdfe, gtools, etc.)
 * that may call omp_set_num_threads() and lower the OpenMP thread count. */
static int g_ctools_system_threads = 0;

static ctools_thread_policy g_thread_policy;
static int g_thread_policy_ready = 0;

#define CTOOLS_SETTINGS_HEADER "# ctools settings; change with: ctools, threads(#) permanently\n"

int ctools_parse_thread_value(const char *text, int list)
{
    long long value = 0;
    if (text == NULL) return 0;
    while (*text == ' ' || *text == '\t') text++;
    if (*text == '\0') return 0;
    if (*text < '0' || *text > '9') return -1;
    while (*text >= '0' && *text <= '9') {
        value = value * 10 + (*text++ - '0');
        if (value > INT_MAX) value = INT_MAX;  /* clamped later */
    }
    while (*text == ' ' || *text == '\t' || *text == '\r' || *text == '\n') text++;
    if (*text != '\0' && !(list && *text == ',')) return -1;
    return value > 0 ? (int)value : -1;
}

void ctools_thread_policy_resolve(ctools_thread_policy *p)
{
    int limit = p->available > 0 ? p->available : 1;
    if (p->omp_limit > 0 && p->omp_limit < limit) limit = p->omp_limit;
    if (limit > CTOOLS_MAX_THREADS) limit = CTOOLS_MAX_THREADS;
    p->limit = limit;

    if (p->env_ctools > 0) {
        p->threads = p->env_ctools;
        p->source = "CTOOLS_NUM_THREADS";
    } else if (p->saved > 0) {
        p->threads = p->saved;
        p->source = "saved setting";
    } else if (p->env_omp > 0) {
        p->threads = p->env_omp;
        p->source = "OMP_NUM_THREADS";
    } else {
        p->threads = CTOOLS_DEFAULT_THREAD_CAP;
        p->source = "automatic";
    }
    if (p->threads > limit) p->threads = limit;
}

#if defined(__linux__)
/* CPUs granted by one cgroup's quota (v2 cpu.max, v1 CFS), or 0 if none. */
static int cgroup_quota_cpus(const char *dir, int v2)
{
    char path[4200], word[32] = "";
    long long quota = 0, period = 0;
    FILE *fp;
    int ok;

    snprintf(path, sizeof(path), v2 ? "%s/cpu.max" : "%s/cpu.cfs_quota_us", dir);
    if ((fp = fopen(path, "r")) == NULL) return 0;
    ok = v2 ? fscanf(fp, "%31s %lld", word, &period) == 2 : fscanf(fp, "%lld", &quota) == 1;
    fclose(fp);
    if (!ok) return 0;
    if (v2) {
        if (strcmp(word, "max") == 0) return 0;
        quota = strtoll(word, NULL, 10);
    } else {
        snprintf(path, sizeof(path), "%s/cpu.cfs_period_us", dir);
        if ((fp = fopen(path, "r")) == NULL) return 0;
        ok = fscanf(fp, "%lld", &period) == 1;
        fclose(fp);
        if (!ok) return 0;
    }
    if (quota <= 0 || period <= 0) return 0;
    quota = quota / period + (quota % period != 0);
    return quota > INT_MAX ? INT_MAX : (int)quota;
}

/* Tightest quota from this process's cgroup up to the hierarchy root. */
static int cgroup_walk(const char *root, const char *rel, int v2)
{
    char dir[4096];
    size_t base = strlen(root);
    int best = 0;

    if (rel[0] != '/' || strstr(rel, "..") != NULL) rel = "";
    if ((size_t)snprintf(dir, sizeof(dir), "%s%s", root, rel) >= sizeof(dir)) return 0;
    for (;;) {
        size_t len = strlen(dir);
        while (len > base && dir[len - 1] == '/') dir[--len] = '\0';
        int n = cgroup_quota_cpus(dir, v2);
        if (n > 0 && (best == 0 || n < best)) best = n;
        if (len <= base) return best;
        *strrchr(dir + base, '/') = '\0';
    }
}

/* CPU quota of the cgroup that contains this process, or 0 if unlimited.
 * Containers (docker --cpus, Kubernetes limits) set one without changing
 * the affinity mask. */
static int cgroup_cpu_limit(void)
{
    char line[4096], v1[4096] = "", v2[4096] = "", controllers[4100];
    int have_v1 = 0, have_v2 = 0, best = 0;
    FILE *fp = fopen("/proc/self/cgroup", "r");

    if (fp == NULL) return 0;
    /* Lines are "hierarchy-id:controllers:path"; v2 is "0::path". */
    while (fgets(line, sizeof(line), fp) != NULL) {
        char *ctl = strchr(line, ':');
        char *path = ctl ? strchr(ctl + 1, ':') : NULL;
        if (path == NULL) continue;
        *ctl++ = '\0';
        *path++ = '\0';
        path[strcspn(path, "\r\n")] = '\0';
        if (strcmp(line, "0") == 0 && *ctl == '\0') {
            snprintf(v2, sizeof(v2), "%s", path);
            have_v2 = 1;
            continue;
        }
        snprintf(controllers, sizeof(controllers), ",%s,", ctl);
        if (strstr(controllers, ",cpu,") != NULL) {
            snprintf(v1, sizeof(v1), "%s", path);
            have_v1 = 1;
        }
    }
    fclose(fp);

    if (have_v2) best = cgroup_walk("/sys/fs/cgroup", v2, 1);
    if (best == 0 && have_v1) {
        best = cgroup_walk("/sys/fs/cgroup/cpu", v1, 0);
        if (best == 0) best = cgroup_walk("/sys/fs/cgroup/cpu,cpuacct", v1, 0);
    }
    return best;
}

/* CPUs in this thread's affinity mask (taskset, Slurm, numactl), or 0. */
static int linux_affinity_cpus(void)
{
#ifdef CPU_ALLOC_SIZE
    for (int ncpu = 1024; ncpu <= 65536; ncpu *= 2) {
        cpu_set_t *set = CPU_ALLOC(ncpu);
        size_t size = CPU_ALLOC_SIZE(ncpu);
        int n, error;
        if (set == NULL) return 0;
        CPU_ZERO_S(size, set);
        n = sched_getaffinity(0, size, set) == 0 ? CPU_COUNT_S(size, set) : -1;
        error = errno;
        CPU_FREE(set);
        if (n >= 0) return n;
        if (error != EINVAL) return 0;  /* EINVAL: mask larger than the set */
    }
#endif
    return 0;
}
#endif

/* Logical CPUs on the machine, and those this process may run on. */
static void detect_cpus(int *detected, int *available)
{
    int total = 0, usable = 0;
#if defined(_WIN32)
    DWORD_PTR process_mask = 0, system_mask = 0;
#ifdef ALL_PROCESSOR_GROUPS
    total = (int)GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
#else
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    total = (int)si.dwNumberOfProcessors;
#endif
    /* Both masks are zero when threads already span processor groups. */
    if (GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask)) {
        for (; process_mask != 0; process_mask &= process_mask - 1) usable++;
    }
#elif defined(__APPLE__)
    size_t len = sizeof(total);
    if (sysctlbyname("hw.logicalcpu_max", &total, &len, NULL, 0) != 0) total = 0;
    len = sizeof(usable);
    if (sysctlbyname("hw.activecpu", &usable, &len, NULL, 0) != 0) usable = 0;
#else
    total = (int)sysconf(_SC_NPROCESSORS_ONLN);
#if defined(__linux__)
    usable = linux_affinity_cpus();
    int quota = cgroup_cpu_limit();
    if (quota > 0 && (usable < 1 || quota < usable)) usable = quota;
#endif
#endif
    if (total < 1) total = omp_get_max_threads();  /* last resort fallback */
    if (usable < 1 || usable > total) usable = total;
    *detected = total;
    *available = usable;
}

/* Per-user settings file: %APPDATA%\ctools\ctools.conf on Windows,
 * ~/Library/Application Support/ctools/ctools.conf on macOS, and
 * $XDG_CONFIG_HOME/ctools/ctools.conf (default ~/.config) elsewhere.
 * It lives outside the ado tree so reinstalling ctools keeps it. */
static int thread_settings_path(char *out, size_t size)
{
    int n;
#if defined(_WIN32)
    wchar_t wide[2048];
    char base[4096];
    DWORD len = GetEnvironmentVariableW(L"APPDATA", wide, 2048);
    if (len == 0 || len >= 2048) return -1;
    if (!WideCharToMultiByte(CP_UTF8, 0, wide, -1, base, (int)sizeof(base), NULL, NULL)) return -1;
    n = snprintf(out, size, "%s\\ctools\\ctools.conf", base);
#else
    const char *home = getenv("HOME");
#if !defined(__APPLE__)
    const char *xdg = getenv("XDG_CONFIG_HOME");
    if (xdg != NULL && xdg[0] == '/') {
        n = snprintf(out, size, "%s/ctools/ctools.conf", xdg);
        return n > 0 && (size_t)n < size ? 0 : -1;
    }
#endif
    if (home == NULL || home[0] != '/') {
        struct passwd *pw = getpwuid(getuid());
        home = pw != NULL ? pw->pw_dir : NULL;
    }
    if (home == NULL || home[0] != '/') return -1;
#if defined(__APPLE__)
    n = snprintf(out, size, "%s/Library/Application Support/ctools/ctools.conf", home);
#else
    n = snprintf(out, size, "%s/.config/ctools/ctools.conf", home);
#endif
#endif
    return n > 0 && (size_t)n < size ? 0 : -1;
}

/* Value text of a "threads=" line, or NULL for any other line. */
static const char *thread_setting_value(const char *line)
{
    if (strncmp(line, "\xEF\xBB\xBF", 3) == 0) line += 3;  /* editors may add a BOM */
    while (*line == ' ' || *line == '\t') line++;
    if (strncmp(line, "threads", 7) != 0) return NULL;
    line += 7;
    while (*line == ' ' || *line == '\t') line++;
    return *line == '=' ? line + 1 : NULL;
}

/* Last threads= value in the settings file, encoded like
 * ctools_parse_thread_value. A missing file is not an error. */
static int read_thread_setting(const char *path)
{
    char line[1024];
    const char *text;
    int value = 0;
    FILE *fp = ctools_fopen(path, "r");

    if (fp == NULL) return 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if ((text = thread_setting_value(line)) != NULL) {
            value = ctools_parse_thread_value(text, 0);
        }
    }
    fclose(fp);
    return value;
}

#if defined(_WIN32)
#define CTOOLS_PATH_SEP '\\'
static int make_dir(const char *path)
{
    wchar_t *wide = ctools_wpath(path);
    BOOL ok;
    if (wide == NULL) return -1;
    ok = CreateDirectoryW(wide, NULL) || GetLastError() == ERROR_ALREADY_EXISTS;
    free(wide);
    if (!ok) errno = EACCES;
    return ok ? 0 : -1;
}
static void remove_empty_dir(const char *path)
{
    wchar_t *wide = ctools_wpath(path);
    if (wide != NULL) RemoveDirectoryW(wide);
    free(wide);
}
#else
#define CTOOLS_PATH_SEP '/'
static int make_dir(const char *path)
{
    return mkdir(path, 0700) == 0 || errno == EEXIST ? 0 : -1;
}
static void remove_empty_dir(const char *path)
{
    (void)rmdir(path);
}
#endif

/* Create any missing directories above a file path. */
static int make_parent_dirs(char *path)
{
    char *sep = strrchr(path, CTOOLS_PATH_SEP);
    int rc;
    if (sep == NULL || sep == path) return 0;
    *sep = '\0';
    rc = make_dir(path);
    if (rc != 0 && make_parent_dirs(path) == 0) rc = make_dir(path);
    *sep = CTOOLS_PATH_SEP;
    return rc;
}

const ctools_thread_policy *ctools_thread_policy_get(void)
{
    ctools_thread_policy *p = &g_thread_policy;
    if (!g_thread_policy_ready) {
        detect_cpus(&p->detected, &p->available);
        p->omp_limit = ctools_parse_thread_value(getenv("OMP_THREAD_LIMIT"), 0);
        p->env_ctools = ctools_parse_thread_value(getenv("CTOOLS_NUM_THREADS"), 0);
        p->env_omp = ctools_parse_thread_value(getenv("OMP_NUM_THREADS"), 1);
        if (thread_settings_path(p->settings_path, sizeof(p->settings_path)) == 0) {
            p->saved = read_thread_setting(p->settings_path);
        } else {
            p->settings_path[0] = '\0';
        }
        ctools_thread_policy_resolve(p);
        g_thread_policy_ready = 1;
    }
    return p;
}

int ctools_clamp_thread_request(int n)
{
    int limit = ctools_thread_policy_get()->limit;
    return n > limit ? limit : n;
}

int ctools_save_thread_setting(int n)
{
    ctools_thread_policy *p = &g_thread_policy;
    char path[4096], temp[4200], line[1024];
    const char *text;
    int keep = 0, ok;
    unsigned long pid;
    FILE *in, *out;

    if (n < 0) n = 0;
    if (n > CTOOLS_MAX_THREADS) n = CTOOLS_MAX_THREADS;
    (void)ctools_thread_policy_get();
    if (p->settings_path[0] == '\0') {
        errno = ENOENT;
        return -1;
    }
    snprintf(path, sizeof(path), "%s", p->settings_path);
    if (n > 0 && make_parent_dirs(path) != 0) return -1;
#if defined(_WIN32)
    pid = (unsigned long)GetCurrentProcessId();
#else
    pid = (unsigned long)getpid();
#endif

    /* Rewrite through a temporary file and rename it over the original, so
     * a concurrent reader sees the old or the new file, never a partial one.
     * Lines other than threads= are kept for settings added later. */
    in = ctools_fopen(path, "r");
    if (in == NULL && n == 0) {
        /* Missing is already removed; unreadable is still removable. */
        if (errno != ENOENT && ctools_remove(path) != 0 && errno != ENOENT) return -1;
        goto saved;
    }
    snprintf(temp, sizeof(temp), "%s.%lu.tmp", path, pid);
    out = ctools_fopen(temp, "w");
    if (out == NULL) {
        if (in != NULL) fclose(in);
        return -1;
    }
    fputs(CTOOLS_SETTINGS_HEADER, out);
    while (in != NULL && fgets(line, sizeof(line), in) != NULL) {
        if (thread_setting_value(line) != NULL || strcmp(line, CTOOLS_SETTINGS_HEADER) == 0) continue;
        for (text = line; *text == ' ' || *text == '\t'; text++) {}
        if (*text != '#' && *text != '\r' && *text != '\n' && *text != '\0') keep = 1;
        fputs(line, out);
    }
    if (in != NULL) fclose(in);
    if (n > 0) fprintf(out, "threads=%d\n", n);
    ok = !ferror(out);
    ok = fclose(out) == 0 && ok;
    if (ok && n == 0 && !keep) {
        /* Nothing else is configured: remove the file, and its ctools
         * directory when that is now empty. */
        ctools_remove(temp);
        ok = ctools_remove(path) == 0 || errno == ENOENT;
        if (ok) {
            char *sep = strrchr(path, CTOOLS_PATH_SEP);
            if (sep != NULL) {
                *sep = '\0';
                remove_empty_dir(path);
            }
        }
    } else if (ok) {
#if defined(_WIN32)
        ok = ctools_move_file(temp, path, 1) == 0;
#else
        ok = rename(temp, path) == 0;
#endif
    }
    if (!ok) {
        int error = errno;
        ctools_remove(temp);
        errno = error;
        return -1;
    }

saved:
    p->saved = n;
    ctools_thread_policy_resolve(p);
    g_ctools_system_threads = p->threads;
    return 0;
}

/*
 * Capture the resolved default thread count.
 * Safe to call multiple times (only resolves once).
 */
void ctools_init_thread_count(void)
{
    if (g_ctools_system_threads == 0) {
        g_ctools_system_threads = ctools_thread_policy_get()->threads;
    }
}

/*
 * Get the current maximum thread count for ctools operations.
 * Returns the user-set limit if specified, otherwise the default
 * resolved at init (immune to external OpenMP pollution).
 */
int ctools_get_max_threads(void)
{
    if (g_ctools_max_threads > 0) {
        return g_ctools_max_threads;
    }
    /* Use the captured default, immune to omp_set_num_threads()
     * calls by external plugins (reghdfe, gtools, etc.) */
    if (g_ctools_system_threads > 0) {
        return g_ctools_system_threads;
    }
    /* Not initialized (native harnesses): follow the OpenMP runtime. */
    return ctools_workspace_threads(omp_get_max_threads(), 0, 0);
}

/*
 * Set the maximum thread count for ctools operations.
 * Pass 0 to reset to the resolved default.
 * Pass n > 0 to set a specific limit. Callers that take the value from the
 * user clamp it with ctools_clamp_thread_request() first.
 */
void ctools_set_max_threads(int n)
{
    if (n < 0) n = 0;
    if (n > CTOOLS_MAX_THREADS) n = CTOOLS_MAX_THREADS;
    g_ctools_max_threads = n;

    /* Also set OpenMP thread limit if OpenMP is enabled and n > 0 */
#ifdef _OPENMP
    if (n > 0) {
        omp_set_num_threads(n);
    }
#endif
}

/*
 * Reset thread limit to the resolved default.
 * Also restores the OpenMP runtime thread count so that OpenMP regions
 * without a num_threads() clause use the same default.
 */
void ctools_reset_max_threads(void)
{
    g_ctools_max_threads = 0;

#ifdef _OPENMP
    /* Always restore the default — protects against both
     * our own threads() option and external plugin pollution */
    if (g_ctools_system_threads > 0) {
        omp_set_num_threads(g_ctools_system_threads);
    }
#endif
}

/* ===========================================================================
   Persistent Thread Pool Implementation
   =========================================================================== */

/*
 * Worker thread function for persistent pool.
 * Waits for work items and executes them until shutdown.
 */
static void *persistent_worker_thread(void *arg)
{
    ctools_persistent_pool *pool = (ctools_persistent_pool *)arg;

    while (1) {
        ctools_work_item *item = NULL;

        /* Wait for work */
        pthread_mutex_lock(&pool->queue_mutex);

        while (pool->queue_head == NULL && !pool->shutdown) {
            pthread_cond_wait(&pool->work_available, &pool->queue_mutex);
        }

        if (pool->shutdown && pool->queue_head == NULL) {
            pthread_mutex_unlock(&pool->queue_mutex);
            break;
        }

        /* Dequeue work item */
        item = pool->queue_head;
        if (item != NULL) {
            pool->queue_head = item->next;
            if (pool->queue_head == NULL) {
                pool->queue_tail = NULL;
            }
            pool->queue_size--;
            pool->active_workers++;
        }

        pthread_mutex_unlock(&pool->queue_mutex);

        /* Execute work item */
        if (item != NULL) {
            void *result = NULL;
            if (item->func != NULL) {
                result = item->func(item->arg);
            }
            free(item);

            /* Signal completion */
            pthread_mutex_lock(&pool->queue_mutex);
            pool->active_workers--;
            pool->pending_items--;
            if (result != NULL) {
                pool->has_error = 1;  /* Track failure (mutex held) */
            }
            if (pool->pending_items == 0) {
                pthread_cond_broadcast(&pool->work_complete);
            }
            pthread_mutex_unlock(&pool->queue_mutex);
        }
    }

    return NULL;
}

/*
 * Initialize a persistent thread pool.
 */
int ctools_persistent_pool_init(ctools_persistent_pool *pool, size_t num_workers)
{
    size_t i;

    if (pool == NULL || num_workers == 0 || num_workers > CTOOLS_MAX_THREADS) {
        return -1;
    }

    memset(pool, 0, sizeof(*pool));

    pool->workers = (pthread_t *)ctools_safe_malloc2(num_workers, sizeof(pthread_t));
    if (pool->workers == NULL) {
        return -1;
    }

    if (pthread_mutex_init(&pool->queue_mutex, NULL) != 0) {
        free(pool->workers);
        return -1;
    }

    if (pthread_cond_init(&pool->work_available, NULL) != 0) {
        pthread_mutex_destroy(&pool->queue_mutex);
        free(pool->workers);
        return -1;
    }

    if (pthread_cond_init(&pool->work_complete, NULL) != 0) {
        pthread_cond_destroy(&pool->work_available);
        pthread_mutex_destroy(&pool->queue_mutex);
        free(pool->workers);
        return -1;
    }

    pool->num_workers = num_workers;
    pool->shutdown = 0;
    pool->initialized = 1;

    pthread_attr_t attr;
    int have_attr = ctools_worker_attr_init(&attr) == 0;

    /* Create worker threads */
    for (i = 0; i < num_workers; i++) {
        if (pthread_create(&pool->workers[i], have_attr ? &attr : NULL,
                           persistent_worker_thread, pool) != 0) {
            /* Failed - shut down and clean up */
            if (have_attr) pthread_attr_destroy(&attr);
            pthread_mutex_lock(&pool->queue_mutex);
            pool->shutdown = 1;
            pthread_cond_broadcast(&pool->work_available);
            pthread_mutex_unlock(&pool->queue_mutex);

            for (size_t j = 0; j < i; j++) {
                pthread_join(pool->workers[j], NULL);
            }

            pthread_cond_destroy(&pool->work_complete);
            pthread_cond_destroy(&pool->work_available);
            pthread_mutex_destroy(&pool->queue_mutex);
            free(pool->workers);
            pool->workers = NULL;
            pool->num_workers = 0;
            pool->initialized = 0;
            return -1;
        }
    }

    if (have_attr) pthread_attr_destroy(&attr);
    return 0;
}

int ctools_worker_attr_init(pthread_attr_t *attr)
{
    size_t size = 0;
    if (pthread_attr_init(attr) != 0) return -1;
    /* Best effort: if the request is refused, the default size remains. */
    if (pthread_attr_getstacksize(attr, &size) != 0 || size < CTOOLS_WORKER_STACK_SIZE)
        (void)pthread_attr_setstacksize(attr, CTOOLS_WORKER_STACK_SIZE);
    return 0;
}

/*
 * Submit work to the persistent pool.
 */
int ctools_persistent_pool_submit(ctools_persistent_pool *pool,
                                   ctools_thread_func func, void *arg)
{
    ctools_work_item *item;

    if (pool == NULL || !pool->initialized || pool->shutdown) {
        return -1;
    }

    item = (ctools_work_item *)malloc(sizeof(ctools_work_item));
    if (item == NULL) {
        return -1;
    }

    item->func = func;
    item->arg = arg;
    item->next = NULL;

    pthread_mutex_lock(&pool->queue_mutex);

    /* Enqueue */
    if (pool->queue_tail != NULL) {
        pool->queue_tail->next = item;
    } else {
        pool->queue_head = item;
    }
    pool->queue_tail = item;
    pool->queue_size++;
    pool->pending_items++;

    pthread_cond_signal(&pool->work_available);
    pthread_mutex_unlock(&pool->queue_mutex);

    return 0;
}

/*
 * Submit multiple work items to the pool, all or none: every item is
 * allocated before any is queued, so on failure nothing has run and callers
 * may run the whole batch themselves.
 */
int ctools_persistent_pool_submit_batch(ctools_persistent_pool *pool,
                                         ctools_thread_func func,
                                         void *args, size_t count, size_t arg_size)
{
    size_t i;
    char *arg_ptr = (char *)args;
    ctools_work_item *head = NULL, *tail = NULL;

    if (pool == NULL || !pool->initialized || pool->shutdown) {
        return -1;
    }
    if (count == 0) {
        return 0;
    }

    for (i = 0; i < count; i++) {
        ctools_work_item *item = (ctools_work_item *)malloc(sizeof(ctools_work_item));
        if (item == NULL) {
            while (head != NULL) {
                ctools_work_item *next = head->next;
                free(head);
                head = next;
            }
            return -1;
        }
        item->func = func;
        item->arg = (args != NULL) ? (void *)(arg_ptr + i * arg_size) : NULL;
        item->next = NULL;
        if (tail != NULL) {
            tail->next = item;
        } else {
            head = item;
        }
        tail = item;
    }

    pthread_mutex_lock(&pool->queue_mutex);

    /* Enqueue the whole list */
    if (pool->queue_tail != NULL) {
        pool->queue_tail->next = head;
    } else {
        pool->queue_head = head;
    }
    pool->queue_tail = tail;
    pool->queue_size += count;
    pool->pending_items += count;

    pthread_cond_broadcast(&pool->work_available);
    pthread_mutex_unlock(&pool->queue_mutex);

    return 0;
}

/*
 * Wait for all submitted work to complete.
 * Returns 0 if all work succeeded, -1 if any work failed.
 */
int ctools_persistent_pool_wait(ctools_persistent_pool *pool)
{
    int had_error;

    if (pool == NULL || !pool->initialized) {
        return -1;
    }

    pthread_mutex_lock(&pool->queue_mutex);

    while (pool->pending_items > 0) {
        pthread_cond_wait(&pool->work_complete, &pool->queue_mutex);
    }

    /* Check and clear error flag (mutex held).
     * Error flag is set by workers when work fails, and cleared here after reading.
     * This ensures errors are reported exactly once per batch of work. */
    had_error = pool->has_error;
    pool->has_error = 0;

    pthread_mutex_unlock(&pool->queue_mutex);

    return had_error ? -1 : 0;
}

/*
 * Shut down and free the persistent pool.
 */
void ctools_persistent_pool_destroy(ctools_persistent_pool *pool)
{
    size_t i;
    ctools_work_item *item;

    if (pool == NULL || !pool->initialized) {
        return;
    }

    /* Signal shutdown */
    pthread_mutex_lock(&pool->queue_mutex);
    pool->shutdown = 1;
    pthread_cond_broadcast(&pool->work_available);
    pthread_mutex_unlock(&pool->queue_mutex);

    /* Wait for workers to finish */
    for (i = 0; i < pool->num_workers; i++) {
        pthread_join(pool->workers[i], NULL);
    }

    /* Free remaining queue items */
    while (pool->queue_head != NULL) {
        item = pool->queue_head;
        pool->queue_head = item->next;
        free(item);
    }

    /* Cleanup */
    pthread_cond_destroy(&pool->work_complete);
    pthread_cond_destroy(&pool->work_available);
    pthread_mutex_destroy(&pool->queue_mutex);
    free(pool->workers);

    pool->initialized = 0;
}

/* Global persistent pool (singleton) */
static ctools_persistent_pool g_global_pool;
static volatile int g_global_pool_initialized = 0;
static volatile int g_global_pool_lock = 0;

/*
 * Simple spinlock acquire/release for global pool protection.
 * Avoids PTHREAD_MUTEX_INITIALIZER which can crash on Windows DLL load.
 */
static void global_pool_lock_acquire(void)
{
#if defined(_WIN32)
    while (InterlockedCompareExchange((volatile long *)&g_global_pool_lock, 1, 0) != 0) {
        /* Spin - in practice this is rarely contended */
    }
#elif defined(__GNUC__) || defined(__clang__)
    while (__sync_lock_test_and_set(&g_global_pool_lock, 1)) {
        /* Spin */
    }
#else
    /* Fallback: no locking (unsafe for concurrent init, but rare) */
#endif
}

static void global_pool_lock_release(void)
{
#if defined(_WIN32)
    InterlockedExchange((volatile long *)&g_global_pool_lock, 0);
#elif defined(__GNUC__) || defined(__clang__)
    __sync_lock_release(&g_global_pool_lock);
#else
    g_global_pool_lock = 0;
#endif
}

/*
 * Get a global persistent pool (singleton).
 *
 * Uses double-checked locking with proper memory barriers to ensure
 * thread-safe lazy initialization on all architectures (including ARM).
 */
ctools_persistent_pool *ctools_get_global_pool(void)
{
    /* Fast path: already initialized
     * Need memory barrier after reading flag to ensure we see fully
     * initialized pool on weakly-ordered architectures (ARM, etc.) */
#if defined(_WIN32) && defined(_MSC_VER)
    if (g_global_pool_initialized) {
        _ReadBarrier();
        MemoryBarrier();
        return &g_global_pool;
    }
#elif defined(__GNUC__) || defined(__clang__)
    if (__atomic_load_n(&g_global_pool_initialized, __ATOMIC_ACQUIRE)) {
        return &g_global_pool;
    }
#else
    if (g_global_pool_initialized) {
        return &g_global_pool;
    }
#endif

    global_pool_lock_acquire();

    /* Double-check under lock */
    if (!g_global_pool_initialized) {
        int num_threads = ctools_get_max_threads();
        if (num_threads < 1) num_threads = 1;
        if (ctools_persistent_pool_init(&g_global_pool, (size_t)num_threads) == 0) {
            /* Memory barrier before setting flag ensures pool is fully
             * visible to other threads before they see initialized=1 */
#if defined(_WIN32) && defined(_MSC_VER)
            MemoryBarrier();
            _WriteBarrier();
            g_global_pool_initialized = 1;
#elif defined(__GNUC__) || defined(__clang__)
            __atomic_store_n(&g_global_pool_initialized, 1, __ATOMIC_RELEASE);
#else
            g_global_pool_initialized = 1;
#endif
        }
    }

    global_pool_lock_release();

    return g_global_pool_initialized ? &g_global_pool : NULL;
}

/*
 * Destroy the global persistent pool.
 */
void ctools_destroy_global_pool(void)
{
    global_pool_lock_acquire();

    if (g_global_pool_initialized) {
        ctools_persistent_pool_destroy(&g_global_pool);
        g_global_pool_initialized = 0;
    }

    global_pool_lock_release();
}
