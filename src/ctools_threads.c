/*
 * ctools_threads.c - Thread pool utilities for ctools
 *
 * Provides a persistent thread pool for efficient parallel execution.
 */

#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "ctools_threads.h"
#include "ctools_config.h"

/* ===========================================================================
   Unified Thread Management Implementation

   Global thread limit that can be set via threads() option in Stata commands.
   When not set (0), defaults to omp_get_max_threads().
   =========================================================================== */

/* User-specified thread limit (0 = use default) */
static int g_ctools_max_threads = 0;

/*
 * Get the current maximum thread count for ctools operations.
 * Returns the user-set limit if specified, otherwise omp_get_max_threads().
 */
int ctools_get_max_threads(void)
{
    if (g_ctools_max_threads > 0) {
        return g_ctools_max_threads;
    }
    /* Default to runtime detection via OpenMP */
    return omp_get_max_threads();
}

/*
 * Set the maximum thread count for ctools operations.
 * Pass 0 to reset to default (omp_get_max_threads()).
 * Pass n > 0 to set a specific limit.
 */
void ctools_set_max_threads(int n)
{
    if (n < 0) n = 0;
    g_ctools_max_threads = n;

    /* Also set OpenMP thread limit if OpenMP is enabled and n > 0 */
#ifdef _OPENMP
    if (n > 0) {
        omp_set_num_threads(n);
    }
#endif
}

/*
 * Reset thread limit to default (runtime-detected).
 */
void ctools_reset_max_threads(void)
{
    g_ctools_max_threads = 0;
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

    if (pool == NULL || num_workers == 0) {
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

    /* Create worker threads */
    for (i = 0; i < num_workers; i++) {
        if (pthread_create(&pool->workers[i], NULL, persistent_worker_thread, pool) != 0) {
            /* Failed - shut down and clean up */
            pool->shutdown = 1;
            pthread_cond_broadcast(&pool->work_available);

            for (size_t j = 0; j < i; j++) {
                pthread_join(pool->workers[j], NULL);
            }

            pthread_cond_destroy(&pool->work_complete);
            pthread_cond_destroy(&pool->work_available);
            pthread_mutex_destroy(&pool->queue_mutex);
            free(pool->workers);
            pool->initialized = 0;
            return -1;
        }
    }

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
 * Submit multiple work items to the pool.
 */
int ctools_persistent_pool_submit_batch(ctools_persistent_pool *pool,
                                         ctools_thread_func func,
                                         void *args, size_t count, size_t arg_size)
{
    size_t i;
    char *arg_ptr = (char *)args;

    if (pool == NULL || !pool->initialized || pool->shutdown) {
        return -1;
    }

    for (i = 0; i < count; i++) {
        void *arg = (args != NULL) ? (void *)(arg_ptr + i * arg_size) : NULL;
        if (ctools_persistent_pool_submit(pool, func, arg) != 0) {
            /* Drain already-submitted items before returning to prevent
             * use-after-free if the caller frees args on error. */
            ctools_persistent_pool_wait(pool);
            return -1;
        }
    }

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
