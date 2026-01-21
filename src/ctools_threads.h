/*
 * ctools_threads.h - Thread pool utilities for ctools
 *
 * Provides a persistent thread pool for efficient parallel execution.
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
