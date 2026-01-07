/*
 * ctools_threads.c - Thread pool utilities for ctools
 *
 * Provides a simple interface for parallel execution patterns.
 */

#include <stdlib.h>
#include <string.h>

#include "ctools_threads.h"

/*
 * Initialize a thread pool for N threads.
 */
int ctools_pool_init(ctools_thread_pool *pool, size_t count,
                     void *args, size_t arg_size)
{
    if (pool == NULL || count == 0) {
        return -1;
    }

    pool->threads = (pthread_t *)malloc(count * sizeof(pthread_t));
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
