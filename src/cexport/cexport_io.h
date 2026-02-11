/*
 * cexport_io.h - Cross-platform I/O abstraction for parallel offset writes
 *
 * Provides a unified interface for high-performance file I/O that enables
 * parallel writes at specific offsets. This is critical for Windows where
 * stdio text mode has CRT locking and newline translation overhead.
 *
 * Backends:
 *   CEXPORT_IO_PWRITE - pwrite (POSIX) / WriteFile+OVERLAPPED (Windows)
 *   CEXPORT_IO_MMAP   - Memory-mapped writes for maximum throughput
 *
 * Performance optimizations:
 *   - Async mmap flush (MS_ASYNC) with optional final sync
 *   - Vectored I/O (pwritev/writev) for batched writes
 *   - Windows async batching with WaitForMultipleObjects
 *   - Pre-fault mmap pages (MADV_POPULATE_WRITE on Linux 5.14+)
 *   - Optional direct I/O (FILE_FLAG_NO_BUFFERING on Windows)
 *   - Deferred fsync for non-critical exports
 *
 * All writes are in binary mode with explicit line endings controlled
 * by the caller, ensuring correct offset math across platforms.
 */

#ifndef CEXPORT_IO_H
#define CEXPORT_IO_H

#include <stddef.h>
#include <stdint.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
/* ssize_t not defined on Windows */
#ifndef _SSIZE_T_DEFINED
#define _SSIZE_T_DEFINED
typedef ptrdiff_t ssize_t;
#endif
#endif

/* Maximum number of concurrent async writes for batching */
#define CEXPORT_IO_MAX_BATCH 64

/* I/O backend types */
typedef enum {
    CEXPORT_IO_PWRITE = 0,   /* pwrite/WriteFile with offsets (default) */
    CEXPORT_IO_MMAP   = 1    /* Memory-mapped write */
} cexport_io_backend_t;

/* I/O flags for fine-grained control */
typedef enum {
    CEXPORT_IO_FLAG_NONE       = 0,
    CEXPORT_IO_FLAG_DIRECT     = (1 << 0),  /* Bypass OS cache (large files) */
    CEXPORT_IO_FLAG_NOFSYNC    = (1 << 1),  /* Skip final fsync (faster, less durable) */
    CEXPORT_IO_FLAG_PREFAULT   = (1 << 2)   /* Pre-fault mmap pages */
} cexport_io_flags_t;

/* Cross-platform file handle */
typedef struct {
#ifdef _WIN32
    HANDLE handle;
    HANDLE mapping_handle;
    /* Async I/O batching */
    OVERLAPPED *batch_overlapped;
    HANDLE *batch_events;
    size_t batch_capacity;
#else
    int fd;
#endif
    size_t file_size;               /* Pre-allocated size */
    size_t actual_written;          /* Track actual bytes written */
    cexport_io_backend_t backend;
    cexport_io_flags_t flags;

    /* For mmap backend */
    char *mapped_data;

    char error_message[256];
} cexport_io_file;

/* Batch write descriptor for vectored I/O */
typedef struct {
    const void *buf;    /* Data buffer */
    size_t len;         /* Buffer length */
    size_t offset;      /* File offset */
} cexport_io_batch_entry;

/*
 * Initialize a cexport_io_file structure.
 * Call before any other operations.
 */
void cexport_io_init(cexport_io_file *file);

/*
 * Open file for writing in binary mode.
 * Creates or truncates the file.
 *
 * @param file     File structure to initialize
 * @param filename Path to output file
 * @param backend  I/O backend to use (PWRITE or MMAP)
 * @param flags    Optional flags (DIRECT, NOFSYNC, PREFAULT)
 * @return         0 on success, -1 on error (error_message set)
 */
int cexport_io_open(cexport_io_file *file, const char *filename,
                    cexport_io_backend_t backend, cexport_io_flags_t flags);

/*
 * Pre-allocate file to given size.
 * Required before parallel writes to establish valid offset range.
 * For mmap backend, also creates the memory mapping.
 *
 * @param file  Opened file
 * @param size  Total size to pre-allocate
 * @return      0 on success, -1 on error
 */
int cexport_io_presize(cexport_io_file *file, size_t size);

/*
 * Get direct pointer to mapped memory at offset (mmap backend only).
 * Allows zero-copy formatting directly into the mapped region.
 *
 * @param file    Opened and pre-sized file with MMAP backend
 * @param offset  Offset into mapped region
 * @return        Pointer to mapped memory, or NULL if not mmap backend
 */
char *cexport_io_get_mapped_ptr(cexport_io_file *file, size_t offset);

/*
 * Write data at specific offset (thread-safe for non-overlapping ranges).
 * Multiple threads can call this concurrently as long as their
 * [offset, offset+len) ranges do not overlap.
 *
 * @param file    Opened and pre-sized file
 * @param buf     Data to write
 * @param len     Number of bytes to write
 * @param offset  File offset to write at
 * @return        Bytes written on success, -1 on error
 */
ssize_t cexport_io_write_at(cexport_io_file *file, const void *buf,
                            size_t len, size_t offset);

/*
 * Batched write using vectored I/O (pwritev/writev or async batching).
 * More efficient than multiple cexport_io_write_at calls.
 *
 * @param file     Opened and pre-sized file
 * @param entries  Array of batch entries (buf, len, offset)
 * @param count    Number of entries (max CEXPORT_IO_MAX_BATCH)
 * @return         Total bytes written on success, -1 on error
 */
ssize_t cexport_io_write_batch(cexport_io_file *file,
                                const cexport_io_batch_entry *entries,
                                size_t count);

/*
 * Finalize and close file.
 * Truncates file to actual_size if smaller than pre-allocated size.
 * Flushes and unmaps memory for mmap backend.
 *
 * @param file        File to close
 * @param actual_size Final file size (for truncation)
 * @return            0 on success, -1 on error
 */
int cexport_io_close(cexport_io_file *file, size_t actual_size);

#endif /* CEXPORT_IO_H */
