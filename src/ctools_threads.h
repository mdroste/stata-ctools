/*
 * ctools_threads.h - Thread pool utilities for ctools
 *
 * Provides a simple interface for parallel execution patterns common
 * in ctools: launch N threads, wait for all to complete, check results.
 *
 * Part of the ctools suite for Stata.
 */

#ifndef CTOOLS_THREADS_H
#define CTOOLS_THREADS_H

#include <pthread.h>
#include <stdlib.h>
#include "ctools_types.h"

/*
 * Thread function signature.
 * Takes a void* argument, returns void*.
 * Return NULL for success, non-NULL for failure.
 */
typedef void *(*ctools_thread_func)(void *);

/*
 * Thread pool structure for managing parallel execution.
 */
typedef struct {
    pthread_t *threads;     /* Array of thread handles */
    void *args;             /* Array of thread arguments (caller-managed) */
    size_t count;           /* Number of threads */
    size_t arg_size;        /* Size of each argument struct */
    int allocated;          /* Whether we allocated the threads array */
} ctools_thread_pool;

/*
 * Initialize a thread pool for N threads.
 * Returns 0 on success, -1 on allocation failure.
 *
 * Parameters:
 *   pool     - Thread pool structure to initialize
 *   count    - Number of threads to manage
 *   args     - Pointer to array of argument structures (caller allocates)
 *   arg_size - Size of each argument structure
 *
 * Example:
 *   my_args_t args[8];
 *   ctools_thread_pool pool;
 *   ctools_pool_init(&pool, 8, args, sizeof(my_args_t));
 */
int ctools_pool_init(ctools_thread_pool *pool, size_t count,
                     void *args, size_t arg_size);

/*
 * Launch all threads in the pool with the given function.
 * Each thread receives &args[i] as its argument.
 * Returns 0 on success, -1 if any pthread_create fails.
 *
 * Example:
 *   ctools_pool_launch(&pool, my_thread_func);
 */
int ctools_pool_launch(ctools_thread_pool *pool, ctools_thread_func func);

/*
 * Wait for all threads in the pool to complete.
 * Returns 0 if all threads returned NULL (success).
 * Returns -1 if any thread returned non-NULL (failure).
 *
 * Example:
 *   int result = ctools_pool_join(&pool);
 */
int ctools_pool_join(ctools_thread_pool *pool);

/*
 * Free resources associated with the thread pool.
 * Does NOT free the args array (caller's responsibility).
 *
 * Example:
 *   ctools_pool_free(&pool);
 */
void ctools_pool_free(ctools_thread_pool *pool);

/*
 * Convenience function: launch threads and wait for completion.
 * Combines ctools_pool_launch and ctools_pool_join.
 * Returns 0 on success, -1 on failure.
 *
 * Example:
 *   int result = ctools_pool_run(&pool, my_thread_func);
 */
int ctools_pool_run(ctools_thread_pool *pool, ctools_thread_func func);

/*
 * One-shot parallel execution.
 * Allocates pool, launches threads, waits, and cleans up.
 * Returns 0 on success, -1 on failure.
 *
 * Parameters:
 *   func     - Thread function to execute
 *   args     - Array of argument structures
 *   count    - Number of threads/arguments
 *   arg_size - Size of each argument structure
 *
 * Example:
 *   my_args_t args[8];
 *   // ... fill in args ...
 *   int result = ctools_parallel_run(my_func, args, 8, sizeof(my_args_t));
 */
int ctools_parallel_run(ctools_thread_func func, void *args,
                        size_t count, size_t arg_size);

/*
 * Macro for simple parallel for-each pattern.
 * Executes func(args[i]) for i = 0..count-1 in parallel.
 *
 * Example:
 *   CTOOLS_PARALLEL_FOR(my_func, my_args, 8, sizeof(my_args_t), result);
 *   if (result != 0) { handle_error(); }
 */
#define CTOOLS_PARALLEL_FOR(func, args, count, arg_size, result) \
    do { \
        (result) = ctools_parallel_run((func), (args), (count), (arg_size)); \
    } while (0)

// Thread argument structure for variable I/O operations
// Used by ctools_data_io for parallel load/store processing
typedef struct {
    stata_variable *var;        /* Variable structure (in/out) */
    int var_idx;                /* 1-based Stata variable index */
    size_t obs1;                /* First observation (1-based) */
    size_t nobs;                /* Number of observations */
    int is_string;              /* Non-zero if string variable (load only) */
    int success;                /* 1 on success, 0 on failure (out) */
} ctools_var_io_args;

#endif /* CTOOLS_THREADS_H */
