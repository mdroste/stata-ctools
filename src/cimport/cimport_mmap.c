/*
 * cimport_mmap.c
 * Cross-platform memory-mapped file I/O for cimport
 */

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "cimport_mmap.h"

/* Need full context definition for file_data, file_size, error_message fields */
#include "cimport_context.h"

#ifdef CIMPORT_WINDOWS

int cimport_mmap_file(CImportContext *ctx, const char *filename)
{
    /* Open file */
    ctx->file_handle = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ,
                                    NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (ctx->file_handle == INVALID_HANDLE_VALUE) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot open file: error %lu", GetLastError());
        return -1;
    }

    /* Get file size */
    LARGE_INTEGER file_size;
    if (!GetFileSizeEx(ctx->file_handle, &file_size)) {
        CloseHandle(ctx->file_handle);
        ctx->file_handle = NULL;
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot get file size: error %lu", GetLastError());
        return -1;
    }

    ctx->file_size = (size_t)file_size.QuadPart;
    ctx->mmap_size = ctx->file_size;
    if (ctx->file_size == 0) {
        /* Empty file - set file_data to NULL and return success.
         * The caller should handle this case by creating an empty dataset. */
        CloseHandle(ctx->file_handle);
        ctx->file_handle = NULL;
        ctx->file_data = NULL;
        ctx->mmap_data = NULL;
        return 0;
    }

    /* Create file mapping */
    ctx->mapping_handle = CreateFileMappingA(ctx->file_handle, NULL, PAGE_READONLY,
                                              0, 0, NULL);
    if (ctx->mapping_handle == NULL) {
        CloseHandle(ctx->file_handle);
        ctx->file_handle = NULL;
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot create file mapping: error %lu", GetLastError());
        return -1;
    }

    /* Map view of file */
    ctx->file_data = (char *)MapViewOfFile(ctx->mapping_handle, FILE_MAP_READ,
                                            0, 0, ctx->file_size);
    if (ctx->file_data == NULL) {
        CloseHandle(ctx->mapping_handle);
        ctx->mapping_handle = NULL;
        CloseHandle(ctx->file_handle);
        ctx->file_handle = NULL;
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot map file: error %lu", GetLastError());
        return -1;
    }

    /* Store original mmap pointer for later unmapping
     * (file_data may be modified for BOM skip or encoding conversion) */
    ctx->mmap_data = ctx->file_data;

    return 0;
}

void cimport_munmap_file(CImportContext *ctx)
{
    /* Use mmap_data (original pointer) for unmapping, not file_data
     * which may have been modified for BOM skip or encoding conversion */
    if (ctx->mmap_data) {
        UnmapViewOfFile(ctx->mmap_data);
        ctx->mmap_data = NULL;
        ctx->file_data = NULL;
    }
    if (ctx->mapping_handle != NULL && ctx->mapping_handle != INVALID_HANDLE_VALUE) {
        CloseHandle(ctx->mapping_handle);
        ctx->mapping_handle = NULL;
    }
    if (ctx->file_handle != NULL && ctx->file_handle != INVALID_HANDLE_VALUE) {
        CloseHandle(ctx->file_handle);
        ctx->file_handle = NULL;
    }
}

#else /* Unix/POSIX */

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int cimport_mmap_file(CImportContext *ctx, const char *filename)
{
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot open file: %s", strerror(errno));
        return -1;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot stat file: %s", strerror(errno));
        return -1;
    }

    ctx->file_size = st.st_size;
    ctx->mmap_size = st.st_size;
    if (ctx->file_size == 0) {
        /* Empty file - set file_data to NULL and return success.
         * The caller should handle this case by creating an empty dataset. */
        close(fd);
        ctx->file_data = NULL;
        ctx->mmap_data = NULL;
        return 0;
    }

    ctx->file_data = mmap(NULL, ctx->file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);

    if (ctx->file_data == MAP_FAILED) {
        snprintf(ctx->error_message, sizeof(ctx->error_message),
                 "Cannot mmap file: %s", strerror(errno));
        ctx->file_data = NULL;
        ctx->mmap_data = NULL;
        return -1;
    }

    /* Store original mmap pointer for later unmapping
     * (file_data may be modified for BOM skip or encoding conversion) */
    ctx->mmap_data = ctx->file_data;

    /* Advise kernel about sequential access pattern */
#ifdef __linux__
    madvise(ctx->file_data, ctx->file_size, MADV_SEQUENTIAL | MADV_WILLNEED);
#elif defined(__APPLE__)
    madvise(ctx->file_data, ctx->file_size, MADV_SEQUENTIAL);
#endif

    return 0;
}

void cimport_munmap_file(CImportContext *ctx)
{
    /* Use mmap_data/mmap_size (original values) for unmapping, not file_data/file_size
     * which may have been modified for BOM skip or encoding conversion */
    if (ctx->mmap_data) {
        munmap(ctx->mmap_data, ctx->mmap_size);
        ctx->mmap_data = NULL;
        ctx->file_data = NULL;
    }
}

#endif /* CIMPORT_WINDOWS */
