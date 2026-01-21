/*
 * cexport_io.c - Cross-platform I/O implementation for parallel offset writes
 *
 * Performance optimizations implemented:
 *
 * 1. Async mmap flush - MS_ASYNC instead of MS_SYNC, with optional final sync
 * 2. Vectored I/O - pwritev() on POSIX, batched async WriteFile on Windows
 * 3. Windows async batching - WaitForMultipleObjects for concurrent writes
 * 4. Better madvise hints - MADV_SEQUENTIAL | MADV_WILLNEED, MADV_POPULATE_WRITE
 * 5. Optional direct I/O - FILE_FLAG_NO_BUFFERING on Windows for large files
 * 6. Deferred fsync - Skip final sync with CEXPORT_IO_FLAG_NOFSYNC
 * 7. Zero-copy mmap access - cexport_io_get_mapped_ptr() for direct formatting
 *
 * Windows Notes:
 * - FILE_FLAG_OVERLAPPED enables async I/O and offset-based writes
 * - FILE_FLAG_NO_BUFFERING bypasses cache but requires aligned writes
 * - SetFilePointerEx + SetEndOfFile pre-allocates without zeroing
 * - WaitForMultipleObjects batches completion waits for better throughput
 *
 * POSIX Notes:
 * - pwrite() is atomic and thread-safe for non-overlapping ranges
 * - pwritev() batches multiple buffers into single syscall
 * - mmap with MAP_SHARED + MS_ASYNC allows async writeback
 * - MADV_POPULATE_WRITE (Linux 5.14+) pre-faults pages to avoid stalls
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "cexport_io.h"

/* ============================================================================
   Common Functions
   ============================================================================ */

void cexport_io_init(cexport_io_file *file)
{
    if (file == NULL) return;
    memset(file, 0, sizeof(*file));
#ifdef _WIN32
    file->handle = INVALID_HANDLE_VALUE;
    file->mapping_handle = NULL;
    file->batch_overlapped = NULL;
    file->batch_events = NULL;
    file->batch_capacity = 0;
#else
    file->fd = -1;
#endif
    file->mapped_data = NULL;
    file->file_size = 0;
    file->actual_written = 0;
    file->backend = CEXPORT_IO_PWRITE;
    file->flags = CEXPORT_IO_FLAG_NONE;
    file->error_message[0] = '\0';
}

const char *cexport_io_line_ending(int use_crlf)
{
    return use_crlf ? "\r\n" : "\n";
}

int cexport_io_line_ending_len(int use_crlf)
{
    return use_crlf ? 2 : 1;
}

/* Common: get mapped pointer (works for both platforms) */
char *cexport_io_get_mapped_ptr(cexport_io_file *file, size_t offset)
{
    if (file == NULL || file->mapped_data == NULL) {
        return NULL;
    }
    if (offset >= file->file_size) {
        return NULL;
    }
    return file->mapped_data + offset;
}

/* ============================================================================
   Windows Implementation
   ============================================================================ */

#ifdef _WIN32

/*
 * Initialize async batch resources for Windows.
 * Pre-allocates OVERLAPPED structures and events for batched writes.
 */
static int win_init_batch(cexport_io_file *file, size_t capacity)
{
    if (capacity == 0 || capacity > CEXPORT_IO_MAX_BATCH) {
        capacity = CEXPORT_IO_MAX_BATCH;
    }

    file->batch_overlapped = (OVERLAPPED *)calloc(capacity, sizeof(OVERLAPPED));
    file->batch_events = (HANDLE *)calloc(capacity, sizeof(HANDLE));

    if (file->batch_overlapped == NULL || file->batch_events == NULL) {
        free(file->batch_overlapped);
        free(file->batch_events);
        file->batch_overlapped = NULL;
        file->batch_events = NULL;
        return -1;
    }

    /* Pre-create events for async I/O */
    for (size_t i = 0; i < capacity; i++) {
        file->batch_events[i] = CreateEvent(NULL, TRUE, FALSE, NULL);
        if (file->batch_events[i] == NULL) {
            /* Cleanup on failure */
            for (size_t j = 0; j < i; j++) {
                CloseHandle(file->batch_events[j]);
            }
            free(file->batch_overlapped);
            free(file->batch_events);
            file->batch_overlapped = NULL;
            file->batch_events = NULL;
            return -1;
        }
        file->batch_overlapped[i].hEvent = file->batch_events[i];
    }

    file->batch_capacity = capacity;
    return 0;
}

/*
 * Cleanup async batch resources.
 */
static void win_cleanup_batch(cexport_io_file *file)
{
    if (file->batch_events != NULL) {
        for (size_t i = 0; i < file->batch_capacity; i++) {
            if (file->batch_events[i] != NULL) {
                CloseHandle(file->batch_events[i]);
            }
        }
        free(file->batch_events);
        file->batch_events = NULL;
    }
    if (file->batch_overlapped != NULL) {
        free(file->batch_overlapped);
        file->batch_overlapped = NULL;
    }
    file->batch_capacity = 0;
}

int cexport_io_open(cexport_io_file *file, const char *filename,
                    cexport_io_backend_t backend, cexport_io_flags_t flags)
{
    if (file == NULL || filename == NULL) {
        return -1;
    }

    cexport_io_init(file);
    file->backend = backend;
    file->flags = flags;

    /*
     * Build file flags:
     * - FILE_FLAG_OVERLAPPED: enables async I/O and offset-based writes
     * - FILE_FLAG_SEQUENTIAL_SCAN: hints cache manager for sequential access
     * - FILE_FLAG_NO_BUFFERING: bypasses cache (optional, for large files)
     */
    DWORD dwFlags = FILE_ATTRIBUTE_NORMAL;

    if (backend == CEXPORT_IO_PWRITE) {
        dwFlags |= FILE_FLAG_OVERLAPPED;
    }

    dwFlags |= FILE_FLAG_SEQUENTIAL_SCAN;

    if (flags & CEXPORT_IO_FLAG_DIRECT) {
        /*
         * Direct I/O bypasses the file system cache.
         * Requires aligned writes (typically 4KB boundaries).
         * Best for very large files where cache pollution is a concern.
         */
        dwFlags |= FILE_FLAG_NO_BUFFERING | FILE_FLAG_WRITE_THROUGH;
    }

    file->handle = CreateFileA(
        filename,
        GENERIC_READ | GENERIC_WRITE,
        0,                              /* No sharing during write */
        NULL,                           /* Default security */
        CREATE_ALWAYS,                  /* Create or truncate */
        dwFlags,
        NULL                            /* No template */
    );

    if (file->handle == INVALID_HANDLE_VALUE) {
        DWORD err = GetLastError();
        snprintf(file->error_message, sizeof(file->error_message),
                 "Cannot create file '%s': Windows error %lu", filename, err);
        return -1;
    }

    /* Initialize async batch resources for pwrite backend */
    if (backend == CEXPORT_IO_PWRITE) {
        if (win_init_batch(file, CEXPORT_IO_MAX_BATCH) != 0) {
            CloseHandle(file->handle);
            file->handle = INVALID_HANDLE_VALUE;
            snprintf(file->error_message, sizeof(file->error_message),
                     "Failed to initialize async I/O batch");
            return -1;
        }
    }

    return 0;
}

int cexport_io_presize(cexport_io_file *file, size_t size)
{
    if (file == NULL || file->handle == INVALID_HANDLE_VALUE) {
        return -1;
    }

    if (size == 0) {
        file->file_size = 0;
        return 0;
    }

    /*
     * Pre-allocate file to specified size:
     * 1. SetFilePointerEx moves the file pointer
     * 2. SetEndOfFile extends the file to that position
     *
     * This avoids fragmentation on NTFS and ensures offsets are valid.
     */
    LARGE_INTEGER li;
    li.QuadPart = (LONGLONG)size;

    if (!SetFilePointerEx(file->handle, li, NULL, FILE_BEGIN)) {
        DWORD err = GetLastError();
        snprintf(file->error_message, sizeof(file->error_message),
                 "SetFilePointerEx failed: Windows error %lu", err);
        return -1;
    }

    if (!SetEndOfFile(file->handle)) {
        DWORD err = GetLastError();
        snprintf(file->error_message, sizeof(file->error_message),
                 "SetEndOfFile failed: Windows error %lu", err);
        return -1;
    }

    file->file_size = size;

    /* For mmap backend, create file mapping now */
    if (file->backend == CEXPORT_IO_MMAP) {
        DWORD size_high = (DWORD)(size >> 32);
        DWORD size_low = (DWORD)(size & 0xFFFFFFFF);

        file->mapping_handle = CreateFileMappingA(
            file->handle,
            NULL,           /* Default security */
            PAGE_READWRITE, /* Read/write access */
            size_high,
            size_low,
            NULL            /* No name */
        );

        if (file->mapping_handle == NULL) {
            DWORD err = GetLastError();
            snprintf(file->error_message, sizeof(file->error_message),
                     "CreateFileMapping failed: Windows error %lu", err);
            return -1;
        }

        file->mapped_data = (char *)MapViewOfFile(
            file->mapping_handle,
            FILE_MAP_WRITE,
            0, 0,           /* Offset: start of file */
            size            /* Map entire file */
        );

        if (file->mapped_data == NULL) {
            DWORD err = GetLastError();
            CloseHandle(file->mapping_handle);
            file->mapping_handle = NULL;
            snprintf(file->error_message, sizeof(file->error_message),
                     "MapViewOfFile failed: Windows error %lu", err);
            return -1;
        }

        /*
         * Pre-fault pages if requested.
         * Windows doesn't have MADV_POPULATE_WRITE, but we can use
         * PrefetchVirtualMemory (Windows 8+) or just touch pages.
         */
        if (file->flags & CEXPORT_IO_FLAG_PREFAULT) {
            /* Touch first byte of each 64KB region to fault in pages */
            volatile char touch;
            for (size_t off = 0; off < size; off += 65536) {
                touch = file->mapped_data[off];
                (void)touch;
            }
        }
    }

    return 0;
}

ssize_t cexport_io_write_at(cexport_io_file *file, const void *buf,
                            size_t len, size_t offset)
{
    if (file == NULL || buf == NULL) {
        return -1;
    }

    if (len == 0) {
        return 0;
    }

    /* Bounds check */
    if (offset + len > file->file_size) {
        snprintf(file->error_message, sizeof(file->error_message),
                 "Write at offset %zu + len %zu exceeds file size %zu",
                 offset, len, file->file_size);
        return -1;
    }

    if (file->backend == CEXPORT_IO_MMAP) {
        /* mmap backend: direct memory copy (zero syscall overhead) */
        memcpy(file->mapped_data + offset, buf, len);
        return (ssize_t)len;
    }

    /*
     * pwrite-style using OVERLAPPED with event for async completion.
     * Using event-based completion allows batching via WaitForMultipleObjects.
     */
    OVERLAPPED ov;
    memset(&ov, 0, sizeof(ov));
    ov.Offset = (DWORD)(offset & 0xFFFFFFFF);
    ov.OffsetHigh = (DWORD)(offset >> 32);
    /* Use NULL event for single synchronous write */

    DWORD written = 0;
    BOOL success = WriteFile(file->handle, buf, (DWORD)len, &written, &ov);

    if (!success) {
        DWORD err = GetLastError();
        if (err == ERROR_IO_PENDING) {
            /* Wait for completion */
            success = GetOverlappedResult(file->handle, &ov, &written, TRUE);
            if (!success) {
                err = GetLastError();
                snprintf(file->error_message, sizeof(file->error_message),
                         "GetOverlappedResult failed: Windows error %lu", err);
                return -1;
            }
        } else {
            snprintf(file->error_message, sizeof(file->error_message),
                     "WriteFile failed: Windows error %lu", err);
            return -1;
        }
    }

    return (ssize_t)written;
}

ssize_t cexport_io_write_batch(cexport_io_file *file,
                                const cexport_io_batch_entry *entries,
                                size_t count)
{
    if (file == NULL || entries == NULL || count == 0) {
        return -1;
    }

    if (file->backend == CEXPORT_IO_MMAP) {
        /* For mmap, just do sequential memcpy (already fast) */
        ssize_t total = 0;
        for (size_t i = 0; i < count; i++) {
            if (entries[i].offset + entries[i].len > file->file_size) {
                return -1;
            }
            memcpy(file->mapped_data + entries[i].offset,
                   entries[i].buf, entries[i].len);
            total += entries[i].len;
        }
        return total;
    }

    /* pwrite backend: use async batching */
    if (file->batch_overlapped == NULL || count > file->batch_capacity) {
        /* Fallback to sequential writes if batch not initialized */
        ssize_t total = 0;
        for (size_t i = 0; i < count; i++) {
            ssize_t written = cexport_io_write_at(file, entries[i].buf,
                                                   entries[i].len,
                                                   entries[i].offset);
            if (written < 0) return -1;
            total += written;
        }
        return total;
    }

    /*
     * Async batched writes:
     * 1. Issue all WriteFile calls with events
     * 2. WaitForMultipleObjects to wait for all completions
     * 3. Check results
     */
    for (size_t i = 0; i < count; i++) {
        ResetEvent(file->batch_events[i]);

        OVERLAPPED *ov = &file->batch_overlapped[i];
        ov->Offset = (DWORD)(entries[i].offset & 0xFFFFFFFF);
        ov->OffsetHigh = (DWORD)(entries[i].offset >> 32);
        ov->hEvent = file->batch_events[i];

        DWORD written = 0;
        BOOL success = WriteFile(file->handle, entries[i].buf,
                                  (DWORD)entries[i].len, &written, ov);

        if (!success && GetLastError() != ERROR_IO_PENDING) {
            snprintf(file->error_message, sizeof(file->error_message),
                     "WriteFile batch[%zu] failed: Windows error %lu",
                     i, GetLastError());
            return -1;
        }
    }

    /* Wait for all writes to complete */
    DWORD wait_result = WaitForMultipleObjects(
        (DWORD)count,
        file->batch_events,
        TRUE,       /* Wait for ALL */
        INFINITE
    );

    if (wait_result == WAIT_FAILED) {
        snprintf(file->error_message, sizeof(file->error_message),
                 "WaitForMultipleObjects failed: Windows error %lu",
                 GetLastError());
        return -1;
    }

    /* Verify all writes completed successfully */
    ssize_t total = 0;
    for (size_t i = 0; i < count; i++) {
        DWORD written = 0;
        if (!GetOverlappedResult(file->handle, &file->batch_overlapped[i],
                                  &written, FALSE)) {
            snprintf(file->error_message, sizeof(file->error_message),
                     "GetOverlappedResult batch[%zu] failed: Windows error %lu",
                     i, GetLastError());
            return -1;
        }
        if (written != entries[i].len) {
            snprintf(file->error_message, sizeof(file->error_message),
                     "Short write in batch[%zu]: expected %zu, got %lu",
                     i, entries[i].len, written);
            return -1;
        }
        total += written;
    }

    return total;
}

int cexport_io_close(cexport_io_file *file, size_t actual_size)
{
    int rc = 0;

    if (file == NULL) {
        return -1;
    }

    /* Unmap memory if using mmap backend */
    if (file->mapped_data != NULL) {
        /*
         * FlushViewOfFile with 0 means flush entire mapping.
         * This is async - doesn't wait for disk write.
         */
        if (!FlushViewOfFile(file->mapped_data, 0)) {
            rc = -1;
        }
        if (!UnmapViewOfFile(file->mapped_data)) {
            rc = -1;
        }
        file->mapped_data = NULL;
    }

    if (file->mapping_handle != NULL) {
        CloseHandle(file->mapping_handle);
        file->mapping_handle = NULL;
    }

    /* Cleanup async batch resources */
    win_cleanup_batch(file);

    /* Truncate to actual size if smaller than pre-allocated */
    if (file->handle != INVALID_HANDLE_VALUE) {
        if (actual_size < file->file_size) {
            LARGE_INTEGER li;
            li.QuadPart = (LONGLONG)actual_size;
            if (SetFilePointerEx(file->handle, li, NULL, FILE_BEGIN)) {
                SetEndOfFile(file->handle);
            }
        }

        /*
         * Flush file buffers before close.
         * Skip if NOFSYNC flag is set (faster, less durable).
         */
        if (!(file->flags & CEXPORT_IO_FLAG_NOFSYNC)) {
            FlushFileBuffers(file->handle);
        }

        if (!CloseHandle(file->handle)) {
            rc = -1;
        }
        file->handle = INVALID_HANDLE_VALUE;
    }

    file->file_size = 0;
    file->actual_written = 0;

    return rc;
}

#else /* POSIX */

/* ============================================================================
   POSIX Implementation (Linux, macOS, etc.)
   ============================================================================ */

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/uio.h>  /* For pwritev/writev */
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>

/* IOV_MAX fallback if not defined */
#ifndef IOV_MAX
#ifdef _XOPEN_IOV_MAX
#define IOV_MAX _XOPEN_IOV_MAX
#else
#define IOV_MAX 16  /* Conservative fallback */
#endif
#endif

/* Check for pwritev support */
#if defined(__linux__) || defined(__FreeBSD__) || defined(__APPLE__)
#define HAVE_PWRITEV 1
#endif

/* fdatasync is not available on macOS - use fsync instead */
#ifdef __APPLE__
#define fdatasync(fd) fsync(fd)
#endif

int cexport_io_open(cexport_io_file *file, const char *filename,
                    cexport_io_backend_t backend, cexport_io_flags_t flags)
{
    if (file == NULL || filename == NULL) {
        return -1;
    }

    cexport_io_init(file);
    file->backend = backend;
    file->flags = flags;

    /*
     * Open file for writing:
     * - O_WRONLY: write-only (O_RDWR needed for mmap)
     * - O_CREAT: create if doesn't exist
     * - O_TRUNC: truncate if exists
     * - 0644: standard file permissions
     */
    int open_flags = O_CREAT | O_TRUNC;
    if (backend == CEXPORT_IO_MMAP) {
        open_flags |= O_RDWR;  /* mmap needs read+write */
    } else {
        open_flags |= O_WRONLY;
    }

    /*
     * Direct I/O (O_DIRECT) bypasses the page cache.
     * Requires aligned buffers and offsets (typically 512 or 4096 bytes).
     * Not commonly used for CSV export due to alignment requirements.
     */
#ifdef O_DIRECT
    if (flags & CEXPORT_IO_FLAG_DIRECT) {
        open_flags |= O_DIRECT;
    }
#endif

    file->fd = open(filename, open_flags, 0644);
    if (file->fd < 0) {
        snprintf(file->error_message, sizeof(file->error_message),
                 "Cannot create file '%s': %s", filename, strerror(errno));
        return -1;
    }

    return 0;
}

int cexport_io_presize(cexport_io_file *file, size_t size)
{
    if (file == NULL || file->fd < 0) {
        return -1;
    }

    if (size == 0) {
        file->file_size = 0;
        return 0;
    }

    /*
     * ftruncate extends (or shrinks) the file to the specified size.
     * On most filesystems, the new space is zero-filled.
     * This ensures all offsets in [0, size) are valid for pwrite.
     */
    if (ftruncate(file->fd, (off_t)size) < 0) {
        snprintf(file->error_message, sizeof(file->error_message),
                 "ftruncate to %zu failed: %s", size, strerror(errno));
        return -1;
    }

    file->file_size = size;

    /* For mmap backend, create the mapping now */
    if (file->backend == CEXPORT_IO_MMAP) {
        file->mapped_data = (char *)mmap(
            NULL,                   /* Let OS choose address */
            size,
            PROT_READ | PROT_WRITE,
            MAP_SHARED,             /* Write through to file */
            file->fd,
            0                       /* Offset: start of file */
        );

        if (file->mapped_data == MAP_FAILED) {
            file->mapped_data = NULL;
            snprintf(file->error_message, sizeof(file->error_message),
                     "mmap of %zu bytes failed: %s", size, strerror(errno));
            return -1;
        }

        /*
         * Advise kernel about our access pattern.
         * Combining multiple hints for optimal performance.
         */
#ifdef MADV_SEQUENTIAL
        madvise(file->mapped_data, size, MADV_SEQUENTIAL);
#endif

#ifdef MADV_WILLNEED
        /* Tell kernel we'll need these pages soon */
        madvise(file->mapped_data, size, MADV_WILLNEED);
#endif

        /*
         * Pre-fault pages if requested.
         * MADV_POPULATE_WRITE (Linux 5.14+) faults in all pages for writing.
         * This avoids page fault latency during parallel formatting.
         */
        if (file->flags & CEXPORT_IO_FLAG_PREFAULT) {
#ifdef MADV_POPULATE_WRITE
            /* Best option: kernel pre-faults all pages */
            madvise(file->mapped_data, size, MADV_POPULATE_WRITE);
#else
            /* Fallback: manually touch first byte of each page */
            size_t page_size = 4096;
            volatile char touch;
            for (size_t off = 0; off < size; off += page_size) {
                touch = file->mapped_data[off];
                file->mapped_data[off] = touch;  /* Write to fault in for write */
            }
#endif
        }
    }

    return 0;
}

ssize_t cexport_io_write_at(cexport_io_file *file, const void *buf,
                            size_t len, size_t offset)
{
    if (file == NULL || buf == NULL) {
        return -1;
    }

    if (len == 0) {
        return 0;
    }

    /* Bounds check */
    if (offset + len > file->file_size) {
        snprintf(file->error_message, sizeof(file->error_message),
                 "Write at offset %zu + len %zu exceeds file size %zu",
                 offset, len, file->file_size);
        return -1;
    }

    if (file->backend == CEXPORT_IO_MMAP) {
        /* mmap backend: direct memory copy (no syscall) */
        memcpy(file->mapped_data + offset, buf, len);
        return (ssize_t)len;
    }

    /*
     * pwrite: atomic positioned write.
     * Single syscall for typical chunk sizes (2-5MB).
     * Modern kernels complete large writes atomically.
     */
    ssize_t written = pwrite(file->fd, buf, len, (off_t)offset);
    if (written < 0) {
        if (errno == EINTR) {
            /* Retry once on interrupt */
            written = pwrite(file->fd, buf, len, (off_t)offset);
        }
        if (written < 0) {
            snprintf(file->error_message, sizeof(file->error_message),
                     "pwrite at offset %zu failed: %s", offset, strerror(errno));
            return -1;
        }
    }

    return written;
}

ssize_t cexport_io_write_batch(cexport_io_file *file,
                                const cexport_io_batch_entry *entries,
                                size_t count)
{
    if (file == NULL || entries == NULL || count == 0) {
        return -1;
    }

    if (file->backend == CEXPORT_IO_MMAP) {
        /* For mmap, sequential memcpy is already optimal */
        ssize_t total = 0;
        for (size_t i = 0; i < count; i++) {
            if (entries[i].offset + entries[i].len > file->file_size) {
                return -1;
            }
            memcpy(file->mapped_data + entries[i].offset,
                   entries[i].buf, entries[i].len);
            total += entries[i].len;
        }
        return total;
    }

#ifdef HAVE_PWRITEV
    /*
     * pwritev: vectored positioned write.
     * Batches multiple buffers into a single syscall.
     *
     * NOTE: pwritev writes buffers contiguously from the given offset.
     * For non-contiguous writes, we still need individual pwrite calls.
     * But we can group contiguous ranges.
     */

    /* Check if all entries are contiguous */
    int contiguous = 1;
    for (size_t i = 1; i < count && contiguous; i++) {
        if (entries[i].offset != entries[i-1].offset + entries[i-1].len) {
            contiguous = 0;
        }
    }

    if (contiguous && count <= IOV_MAX) {
        /* All entries are contiguous - optimize with single syscall */
#ifdef __APPLE__
        /* macOS: pwritev only available on 11.0+, use combined buffer fallback */
        if (__builtin_available(macOS 11.0, *)) {
            struct iovec *iov = (struct iovec *)malloc(count * sizeof(struct iovec));
            if (iov != NULL) {
                for (size_t i = 0; i < count; i++) {
                    iov[i].iov_base = (void *)entries[i].buf;
                    iov[i].iov_len = entries[i].len;
                }
                ssize_t written = pwritev(file->fd, iov, (int)count,
                                           (off_t)entries[0].offset);
                free(iov);
                if (written < 0) {
                    snprintf(file->error_message, sizeof(file->error_message),
                             "pwritev failed: %s", strerror(errno));
                    return -1;
                }
                return written;
            }
        } else {
            /* Older macOS: combine buffers into single pwrite (avoids N syscalls) */
            size_t total_len = 0;
            for (size_t i = 0; i < count; i++) {
                total_len += entries[i].len;
            }
            char *combined = (char *)malloc(total_len);
            if (combined != NULL) {
                char *ptr = combined;
                for (size_t i = 0; i < count; i++) {
                    memcpy(ptr, entries[i].buf, entries[i].len);
                    ptr += entries[i].len;
                }
                ssize_t written = pwrite(file->fd, combined, total_len,
                                          (off_t)entries[0].offset);
                free(combined);
                if (written < 0) {
                    snprintf(file->error_message, sizeof(file->error_message),
                             "pwrite failed: %s", strerror(errno));
                    return -1;
                }
                return written;
            }
        }
#else
        /* Linux/FreeBSD: pwritev always available */
        struct iovec *iov = (struct iovec *)malloc(count * sizeof(struct iovec));
        if (iov != NULL) {
            for (size_t i = 0; i < count; i++) {
                iov[i].iov_base = (void *)entries[i].buf;
                iov[i].iov_len = entries[i].len;
            }
            ssize_t written = pwritev(file->fd, iov, (int)count,
                                       (off_t)entries[0].offset);
            free(iov);
            if (written < 0) {
                snprintf(file->error_message, sizeof(file->error_message),
                         "pwritev failed: %s", strerror(errno));
                return -1;
            }
            return written;
        }
#endif
    }
#endif

    /* Non-contiguous or allocation failed: individual writes */
    ssize_t total = 0;
    for (size_t i = 0; i < count; i++) {
        ssize_t written = cexport_io_write_at(file, entries[i].buf,
                                               entries[i].len,
                                               entries[i].offset);
        if (written < 0) return -1;
        total += written;
    }
    return total;
}

int cexport_io_close(cexport_io_file *file, size_t actual_size)
{
    int rc = 0;

    if (file == NULL) {
        return -1;
    }

    /* Unmap memory if using mmap backend */
    if (file->mapped_data != NULL) {
        /*
         * MS_ASYNC: schedule writeback but don't wait.
         * This is much faster than MS_SYNC for large files.
         * The actual disk write happens asynchronously.
         */
#ifdef MS_ASYNC
        if (msync(file->mapped_data, file->file_size, MS_ASYNC) < 0) {
            /* Non-fatal for async mode */
        }
#endif

        if (munmap(file->mapped_data, file->file_size) < 0) {
            rc = -1;
        }
        file->mapped_data = NULL;
    }

    /* Truncate to actual size if smaller than pre-allocated */
    if (file->fd >= 0) {
        if (actual_size < file->file_size) {
            if (ftruncate(file->fd, (off_t)actual_size) < 0) {
                /* Non-fatal: file may just be larger than needed */
            }
        }

        /*
         * Sync to disk before close.
         * Skip if NOFSYNC flag is set (faster, less durable).
         * Use fdatasync() instead of fsync() to skip metadata sync.
         * (On macOS, fdatasync is mapped to fsync above.)
         */
        if (!(file->flags & CEXPORT_IO_FLAG_NOFSYNC)) {
            fdatasync(file->fd);
        }

        if (close(file->fd) < 0) {
            rc = -1;
        }
        file->fd = -1;
    }

    file->file_size = 0;
    file->actual_written = 0;

    return rc;
}

#endif /* _WIN32 */
