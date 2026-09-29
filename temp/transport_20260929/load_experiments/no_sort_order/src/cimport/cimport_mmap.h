/*
 * cimport_mmap.h
 * Cross-platform memory-mapped file I/O for cimport
 *
 * Provides efficient read-only file access via memory mapping.
 * Supports both Windows (CreateFileMapping) and Unix (mmap) APIs.
 */

#ifndef CIMPORT_MMAP_H
#define CIMPORT_MMAP_H

#include <stdlib.h>
#include <stddef.h>

#ifdef _WIN32
#include <windows.h>
#define CIMPORT_WINDOWS 1
#endif

/* Opaque context for memory-mapped file (forward declaration) */
struct CImportContext;

/*
 * Memory-map a file for reading.
 *
 * @param ctx       Import context (receives file_data, file_size)
 * @param filename  Path to file
 * @return          0 on success, -1 on error (error_message set in ctx)
 */
int cimport_mmap_file(struct CImportContext *ctx, const char *filename);

/*
 * Unmap a previously mapped file and close handles.
 *
 * @param ctx  Import context with mapped file
 */
void cimport_munmap_file(struct CImportContext *ctx);

/*
 * Switch madvise hint for multi-threaded load phase.
 *
 * The initial mmap uses MADV_SEQUENTIAL for the single-pass parse.
 * Before multi-threaded cache build / SPI store (where multiple threads
 * read different columns from the same pages), switch to MADV_NORMAL
 * so the kernel doesn't aggressively discard pages behind one thread's
 * read position that another thread still needs.
 *
 * No-op on Windows.
 */
void cimport_madvise_normal(struct CImportContext *ctx);

#endif /* CIMPORT_MMAP_H */
