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

/* ===========================================================================
   Persistent Thread Pool

   A reusable thread pool that keeps worker threads alive between tasks.
   Eliminates thread creation overhead for repeated parallel operations.
   =========================================================================== */

/*
 * Work item for the persistent pool.
 * The pool maintains a queue of these items for workers to process.
 */
typedef struct ctools_work_item {
    ctools_thread_func func;           /* Function to execute */
    void *arg;                          /* Argument to pass to function */
    struct ctools_work_item *next;      /* Next item in queue */
} ctools_work_item;

/*
 * Persistent thread pool structure.
 * Workers wait on a condition variable for new work.
 */
typedef struct {
    pthread_t *workers;                 /* Array of worker threads */
    size_t num_workers;                 /* Number of worker threads */

    /* Work queue */
    ctools_work_item *queue_head;       /* Head of work queue */
    ctools_work_item *queue_tail;       /* Tail of work queue */
    size_t queue_size;                  /* Number of items in queue */

    /* Synchronization */
    pthread_mutex_t queue_mutex;        /* Protects queue access */
    pthread_cond_t work_available;      /* Signaled when work is added */
    pthread_cond_t work_complete;       /* Signaled when all work is done */

    /* State */
    size_t active_workers;              /* Number of workers currently executing */
    size_t pending_items;               /* Items submitted but not yet complete */
    int shutdown;                       /* Set to 1 to terminate workers */
    int initialized;                    /* 1 if pool is initialized */
    int has_error;                      /* Set to 1 if any work item failed */
} ctools_persistent_pool;

/*
 * Initialize a persistent thread pool.
 * Creates num_workers threads that wait for work.
 * Returns 0 on success, -1 on failure.
 */
int ctools_persistent_pool_init(ctools_persistent_pool *pool, size_t num_workers);

/*
 * Submit work to the persistent pool.
 * The function will be called with the given argument by a worker thread.
 * Returns 0 on success, -1 on failure.
 */
int ctools_persistent_pool_submit(ctools_persistent_pool *pool,
                                   ctools_thread_func func, void *arg);

/*
 * Submit multiple work items to the pool.
 * Submits count items, where item i calls func with args + (i * arg_size).
 * Returns 0 on success, -1 on failure.
 */
int ctools_persistent_pool_submit_batch(ctools_persistent_pool *pool,
                                         ctools_thread_func func,
                                         void *args, size_t count, size_t arg_size);

/*
 * Wait for all submitted work to complete.
 * Returns 0 if all work succeeded, -1 if any work failed.
 */
int ctools_persistent_pool_wait(ctools_persistent_pool *pool);

/*
 * Shut down and free the persistent pool.
 * Waits for all pending work to complete, then terminates workers.
 */
void ctools_persistent_pool_destroy(ctools_persistent_pool *pool);

/*
 * Get a global persistent pool (singleton).
 * Lazily initialized with CTOOLS_IO_MAX_THREADS workers.
 * Returns NULL on initialization failure.
 */
ctools_persistent_pool *ctools_get_global_pool(void);

/*
 * Destroy the global persistent pool.
 * Call at plugin cleanup or program exit.
 */
void ctools_destroy_global_pool(void);

#endif /* CTOOLS_THREADS_H */
