#ifndef CTOOLS_ALLOC_H
#define CTOOLS_ALLOC_H
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

/* ============================================================================
   Cross-platform aligned memory allocation

   IMPORTANT: On Windows, _aligned_malloc() requires _aligned_free().
   On POSIX systems, posix_memalign() memory can be freed with regular free().
   Always use ctools_aligned_free() for memory allocated with ctools_aligned_alloc().
   ============================================================================ */

#if defined(_WIN32)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

/*
    Allocate memory aligned to specified boundary.
    Returns NULL on failure.
    MUST be freed with ctools_aligned_free().
*/
static inline void *ctools_aligned_alloc(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, alignment);
#else
    if (posix_memalign(&ptr, alignment, size) != 0)
    {
        return NULL;
    }
#endif
    return ptr;
}

/*
    Free memory allocated with ctools_aligned_alloc().
    Safe to call with NULL.
*/
static inline void ctools_aligned_free(void *ptr)
{
    if (ptr == NULL)
        return;
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/* ============================================================================
   Safe size arithmetic for memory allocation

   Prevents integer overflow when computing allocation sizes like N*K*sizeof(T).
   Returns 0 on success, -1 on overflow.
   ============================================================================ */

#include <stdint.h>

/*
    Multiply two size_t values safely, checking for overflow.
    Returns 0 on success (result stored in *out), -1 on overflow.
*/
static inline int ctools_safe_mul_size(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a)
    {
        *out = 0;
        return -1; /* overflow */
    }
    *out = a * b;
    return 0;
}

/*
    Compute a * b * c safely for allocation sizes.
    Typical use: ctools_safe_alloc_size(N, K, sizeof(double), &size)
    Returns 0 on success, -1 on overflow.
*/
static inline int ctools_safe_alloc_size(size_t a, size_t b, size_t c, size_t *out)
{
    size_t tmp;
    if (ctools_safe_mul_size(a, b, &tmp) != 0)
    {
        *out = 0;
        return -1;
    }
    if (ctools_safe_mul_size(tmp, c, out) != 0)
    {
        *out = 0;
        return -1;
    }
    return 0;
}

/*
    Allocate memory with overflow-checked size computation.
    Returns NULL on overflow or allocation failure.
    Usage: ptr = ctools_safe_malloc2(N, sizeof(double));
           ptr = ctools_safe_malloc3(N, K, sizeof(double));
*/
static inline void *ctools_safe_malloc2(size_t a, size_t b)
{
    size_t size;
    if (ctools_safe_mul_size(a, b, &size) != 0)
        return NULL;
    return malloc(size);
}

static inline void *ctools_safe_malloc3(size_t a, size_t b, size_t c)
{
    size_t size;
    if (ctools_safe_alloc_size(a, b, c, &size) != 0)
        return NULL;
    return malloc(size);
}

static inline void *ctools_safe_calloc2(size_t a, size_t b)
{
    size_t size;
    if (ctools_safe_mul_size(a, b, &size) != 0)
        return NULL;
    return calloc(1, size);
}

static inline void *ctools_safe_calloc3(size_t a, size_t b, size_t c)
{
    size_t size;
    if (ctools_safe_alloc_size(a, b, c, &size) != 0)
        return NULL;
    return calloc(1, size);
}

/*
    Overflow-safe aligned allocation for (count * element_size).
    Returns NULL on overflow or allocation failure.
    MUST be freed with ctools_aligned_free().
*/
static inline void *ctools_safe_aligned_alloc2(size_t alignment, size_t count, size_t element_size)
{
    size_t size;
    if (ctools_safe_mul_size(count, element_size, &size) != 0)
    {
        return NULL; /* Overflow */
    }
    return ctools_aligned_alloc(alignment, size);
}

/* ============================================================================
   Convenience macros for cache-line aligned allocation

   These macros provide a concise way to allocate memory aligned to the
   cache line boundary, which is the most common alignment requirement.
   ============================================================================ */

/*
    Allocate memory aligned to cache line boundary.
    Returns NULL on failure.
    MUST be freed with ctools_aligned_free().
*/
#define ctools_cacheline_alloc(size) ctools_aligned_alloc(CACHE_LINE_SIZE, (size))

/*
    Overflow-safe cache-line aligned allocation for (count * element_size).
    Returns NULL on overflow or allocation failure.
    MUST be freed with ctools_aligned_free().
*/
#define ctools_safe_cacheline_alloc2(count, elem_size)                                             \
    ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, (count), (elem_size))

#endif
