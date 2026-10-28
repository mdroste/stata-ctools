/* Thread policy: parsing, precedence, CPU limits, cgroup quotas and the
 * per-user settings file. Modes:
 *   unit <tmpdir> <expected settings path>
 *   report    prints the policy resolved from this process's environment */
#undef NDEBUG  /* some compilers define it at -O1; these asserts are the test */
#include "ctools_threads.c"
#include <assert.h>

static ctools_thread_policy resolve(int available, int omp_limit, int env, int saved, int omp)
{
    ctools_thread_policy p;
    memset(&p, 0, sizeof(p));
    p.detected = available;
    p.available = available;
    p.omp_limit = omp_limit;
    p.env_ctools = env;
    p.saved = saved;
    p.env_omp = omp;
    ctools_thread_policy_resolve(&p);
    return p;
}

static void put(const char *path, const char *text)
{
    FILE *fp = ctools_fopen(path, "wb");
    assert(fp != NULL);
    fputs(text, fp);
    assert(fclose(fp) == 0);
}

/* Text mode: the file uses the platform's line endings. */
static void get(const char *path, char *text, size_t size)
{
    FILE *fp = ctools_fopen(path, "r");
    size_t n;
    assert(fp != NULL);
    n = fread(text, 1, size - 1, fp);
    text[n] = '\0';
    fclose(fp);
}

/* Whether the directory holding a file still exists. */
static int parent_exists(const char *path)
{
    char dir[4096];
    snprintf(dir, sizeof(dir), "%s", path);
    *strrchr(dir, CTOOLS_PATH_SEP) = '\0';
#if defined(_WIN32)
    wchar_t *wide = ctools_wpath(dir);
    DWORD attributes = wide ? GetFileAttributesW(wide) : INVALID_FILE_ATTRIBUTES;
    free(wide);
    return attributes != INVALID_FILE_ATTRIBUTES;
#else
    struct stat st;
    return stat(dir, &st) == 0;
#endif
}

static void parsing(void)
{
    assert(ctools_parse_thread_value(NULL, 0) == 0);
    assert(ctools_parse_thread_value("", 0) == 0);
    assert(ctools_parse_thread_value(" \t", 0) == 0);
    assert(ctools_parse_thread_value("8", 0) == 8);
    assert(ctools_parse_thread_value(" 8 \r\n", 0) == 8);
    assert(ctools_parse_thread_value("0", 0) == -1);
    assert(ctools_parse_thread_value("-2", 0) == -1);
    assert(ctools_parse_thread_value("+2", 0) == -1);
    assert(ctools_parse_thread_value("8x", 0) == -1);
    assert(ctools_parse_thread_value("8 9", 0) == -1);
    assert(ctools_parse_thread_value("8,4", 0) == -1);
    assert(ctools_parse_thread_value("8,4", 1) == 8);
    assert(ctools_parse_thread_value("99999999999999999999", 0) == INT_MAX);
}

static void precedence(void)
{
    ctools_thread_policy p;
    p = resolve(200, 0, 0, 0, 0);
    assert(p.threads == CTOOLS_DEFAULT_THREAD_CAP && p.limit == 200);
    assert(strcmp(p.source, "automatic") == 0);
    p = resolve(8, 0, 0, 0, 0);
    assert(p.threads == 8 && strcmp(p.source, "automatic") == 0);
    p = resolve(200, 0, 0, 0, 64);
    assert(p.threads == 64 && strcmp(p.source, "OMP_NUM_THREADS") == 0);
    p = resolve(200, 0, 0, 32, 64);
    assert(p.threads == 32 && strcmp(p.source, "saved setting") == 0);
    p = resolve(200, 0, 48, 32, 64);
    assert(p.threads == 48 && strcmp(p.source, "CTOOLS_NUM_THREADS") == 0);
    /* Every source is clamped to usable CPUs, OMP_THREAD_LIMIT and 256. */
    p = resolve(4, 0, 48, 0, 0);
    assert(p.threads == 4 && strcmp(p.source, "CTOOLS_NUM_THREADS") == 0);
    p = resolve(200, 2, 0, 32, 0);
    assert(p.limit == 2 && p.threads == 2);
    p = resolve(1000, 0, 999, 0, 0);
    assert(p.limit == CTOOLS_MAX_THREADS && p.threads == CTOOLS_MAX_THREADS);
    /* Invalid inputs are ignored. */
    p = resolve(200, -1, -1, -1, -1);
    assert(p.limit == 200 && p.threads == CTOOLS_DEFAULT_THREAD_CAP);
    assert(strcmp(p.source, "automatic") == 0);
    p = resolve(0, 0, 0, 0, 0);
    assert(p.limit == 1 && p.threads == 1);
}

#if defined(__linux__)
static void cgroups(const char *tmp)
{
    char root[4096], path[4200];
    snprintf(root, sizeof(root), "%s/cg", tmp);
    assert(mkdir(root, 0700) == 0);
    snprintf(path, sizeof(path), "%s/a", root);
    assert(mkdir(path, 0700) == 0);
    snprintf(path, sizeof(path), "%s/a/b", root);
    assert(mkdir(path, 0700) == 0);

    /* v2: the tightest quota on the path wins, rounded up to whole CPUs. */
    assert(cgroup_walk(root, "/a/b", 1) == 0);
    snprintf(path, sizeof(path), "%s/a/b/cpu.max", root);
    put(path, "max 100000\n");
    snprintf(path, sizeof(path), "%s/a/cpu.max", root);
    put(path, "400000 100000\n");
    snprintf(path, sizeof(path), "%s/cpu.max", root);
    put(path, "150000 100000\n");
    assert(cgroup_walk(root, "/a/b", 1) == 2);
    assert(cgroup_walk(root, "/a/b/", 1) == 2);
    put(path, "max 100000\n");
    assert(cgroup_walk(root, "/a/b", 1) == 4);
    assert(cgroup_walk(root, "/", 1) == 0);
    assert(cgroup_walk(root, "/../a", 1) == 0);  /* never leaves the root */
    assert(cgroup_walk(root, "/missing/x", 1) == 0);
    snprintf(path, sizeof(path), "%s/a/cpu.max", root);
    put(path, "garbage\n");
    assert(cgroup_walk(root, "/a/b", 1) == 0);

    /* v1: CFS quota over period; -1 means unlimited. */
    snprintf(path, sizeof(path), "%s/cpu.cfs_quota_us", root);
    put(path, "-1\n");
    snprintf(path, sizeof(path), "%s/cpu.cfs_period_us", root);
    put(path, "100000\n");
    snprintf(path, sizeof(path), "%s/a/cpu.cfs_quota_us", root);
    put(path, "250000\n");
    snprintf(path, sizeof(path), "%s/a/cpu.cfs_period_us", root);
    put(path, "100000\n");
    assert(cgroup_walk(root, "/a/b", 0) == 3);
    assert(cgroup_walk(root, "/", 0) == 0);
}
#endif

static void settings(const char *expected)
{
    const ctools_thread_policy *p = ctools_thread_policy_get();
    char text[1024];
    const char *path = p->settings_path;

    assert(strcmp(path, expected) == 0);
    assert(p->saved == 0);
    ctools_init_thread_count();

    /* Save creates missing directories and becomes this process's default. */
    assert(ctools_save_thread_setting(3) == 0);
    assert(p->saved == 3 && read_thread_setting(path) == 3);
    assert(strcmp(p->source, "saved setting") == 0 || p->env_ctools > 0);
    assert(ctools_get_max_threads() == p->threads);
    get(path, text, sizeof(text));
    assert(strcmp(text, CTOOLS_SETTINGS_HEADER "threads=3\n") == 0);

    /* Unknown keys and comments survive a rewrite; the last threads= wins. */
    put(path, "future = 1\n# note\nthreads = 6\nthreads=7\n");
    assert(read_thread_setting(path) == 7);
    assert(ctools_save_thread_setting(12) == 0);
    get(path, text, sizeof(text));
    assert(strcmp(text, CTOOLS_SETTINGS_HEADER "future = 1\n# note\nthreads=12\n") == 0);
    assert(ctools_save_thread_setting(0) == 0);
    get(path, text, sizeof(text));
    assert(strcmp(text, CTOOLS_SETTINGS_HEADER "future = 1\n# note\n") == 0);
    assert(read_thread_setting(path) == 0 && p->saved == 0);

    /* Hand-edited files: BOM and CRLF are accepted, bad values are invalid. */
    put(path, "\xEF\xBB\xBFthreads=5\r\n");
    assert(read_thread_setting(path) == 5);
    put(path, "threads=abc\n");
    assert(read_thread_setting(path) == -1);
    put(path, "threadsx=4\nthreads 4\n");
    assert(read_thread_setting(path) == 0);

    /* Removing the only setting removes the file and its now-empty
     * directory; removing twice is fine, and saving recreates both. */
    put(path, CTOOLS_SETTINGS_HEADER "threads=5\n");
    assert(ctools_save_thread_setting(0) == 0);
    assert(ctools_fopen(path, "r") == NULL && !parent_exists(path));
    assert(ctools_save_thread_setting(0) == 0);
    assert(ctools_save_thread_setting(-4) == 0);
    assert(ctools_save_thread_setting(100000) == 0 && read_thread_setting(path) == CTOOLS_MAX_THREADS);
    assert(ctools_save_thread_setting(0) == 0);
}

int main(int argc, char **argv)
{
    if (argc == 4 && strcmp(argv[1], "unit") == 0) {
        parsing();
        precedence();
#if defined(__linux__)
        cgroups(argv[2]);
#endif
        settings(argv[3]);
        puts("unit ok");
        return 0;
    }
    if (argc == 2 && strcmp(argv[1], "report") == 0) {
        const ctools_thread_policy *p = ctools_thread_policy_get();
        ctools_init_thread_count();
        ctools_reset_max_threads();
        printf("threads=%d source=%s limit=%d detected=%d available=%d saved=%d "
               "env=%d omp=%d omp_limit=%d max=%d clamp=%d\n",
               p->threads, p->source, p->limit, p->detected, p->available, p->saved,
               p->env_ctools, p->env_omp, p->omp_limit, ctools_get_max_threads(),
               ctools_clamp_thread_request(100000));
        return 0;
    }
    return 2;
}
