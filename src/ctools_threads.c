/*
 * ctools_threads.c - Thread pool utilities for ctools
 *
 * Provides a simple interface for parallel execution patterns.
 * Includes both one-shot pools and persistent pools for efficiency.
 */

#include <stdlib.h>
#include <string.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "ctools_threads.h"
#include "ctools_config.h"

/*
 * Initialize a thread pool for N threads.
 */
int ctools_pool_init(ctools_thread_pool *pool, size_t count,
                     void *args, size_t arg_size)
{
    if (pool == NULL || count == 0) {
        return -1;
    }

    pool->threads = (pthread_t *)ctools_safe_malloc2(count, sizeof(pthread_t));
    if (pool->threads == NULL) {
        return -1;
    }

    pool->args = args;
    pool->count = count;
    pool->arg_size = arg_size;
    pool->allocated = 1;

    return 0;
}

/*
 * Launch all threads in the pool with the given function.
 */
int ctools_pool_launch(ctools_thread_pool *pool, ctools_thread_func func)
{
    size_t i;
    char *arg_ptr;

    if (pool == NULL || pool->threads == NULL || func == NULL) {
        return -1;
    }

    arg_ptr = (char *)pool->args;

    for (i = 0; i < pool->count; i++) {
        void *arg = (pool->args != NULL) ? (void *)(arg_ptr + i * pool->arg_size) : NULL;
        if (pthread_create(&pool->threads[i], NULL, func, arg) != 0) {
            /* Failed to create thread - join any already created */
            size_t j;
            for (j = 0; j < i; j++) {
                pthread_join(pool->threads[j], NULL);
            }
            return -1;
        }
    }

    return 0;
}

/*
 * Wait for all threads in the pool to complete.
 */
int ctools_pool_join(ctools_thread_pool *pool)
{
    size_t i;
    int all_success = 1;
    void *retval;

    if (pool == NULL || pool->threads == NULL) {
        return -1;
    }

    for (i = 0; i < pool->count; i++) {
        pthread_join(pool->threads[i], &retval);
        if (retval != NULL) {
            all_success = 0;
        }
    }

    return all_success ? 0 : -1;
}

/*
 * Free resources associated with the thread pool.
 */
void ctools_pool_free(ctools_thread_pool *pool)
{
    if (pool == NULL) {
        return;
    }

    if (pool->allocated && pool->threads != NULL) {
        free(pool->threads);
    }

    pool->threads = NULL;
    pool->args = NULL;
    pool->count = 0;
    pool->arg_size = 0;
    pool->allocated = 0;
}

/*
 * Convenience function: launch threads and wait for completion.
 */
int ctools_pool_run(ctools_thread_pool *pool, ctools_thread_func func)
{
    int rc;

    rc = ctools_pool_launch(pool, func);
    if (rc != 0) {
        return rc;
    }

    return ctools_pool_join(pool);
}

/*
 * One-shot parallel execution.
 */
int ctools_parallel_run(ctools_thread_func func, void *args,
                        size_t count, size_t arg_size)
{
    ctools_thread_pool pool;
    int rc;

    if (count == 0) {
        return 0;  /* Nothing to do */
    }

    rc = ctools_pool_init(&pool, count, args, arg_size);
    if (rc != 0) {
        return rc;
    }

    rc = ctools_pool_run(&pool, func);

    ctools_pool_free(&pool);

    return rc;
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
                pool->has_error = 1;  /* Track failure */
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

    /* Check and clear error flag */
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
 */
ctools_persistent_pool *ctools_get_global_pool(void)
{
    /* Fast path: already initialized */
    if (g_global_pool_initialized) {
        return &g_global_pool;
    }

    global_pool_lock_acquire();

    if (!g_global_pool_initialized) {
        if (ctools_persistent_pool_init(&g_global_pool, CTOOLS_IO_MAX_THREADS) == 0) {
            g_global_pool_initialized = 1;
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
