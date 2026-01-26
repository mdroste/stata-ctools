/*
 * cimport_xlsx_zip.h
 * ZIP archive interface for XLSX import
 *
 * Provides simple wrappers around miniz for XLSX file extraction.
 * XLSX files are ZIP archives containing XML files.
 */

#ifndef CIMPORT_XLSX_ZIP_H
#define CIMPORT_XLSX_ZIP_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* Forward declaration - actual type is miniz's mz_zip_archive */
typedef struct xlsx_zip_archive xlsx_zip_archive;

/* ============================================================================
 * ZIP Archive Operations
 * ============================================================================ */

/*
 * Open a ZIP archive from a file path.
 * Returns NULL on failure.
 */
xlsx_zip_archive *xlsx_zip_open_file(const char *filename);

/*
 * Open a ZIP archive from memory.
 * The memory must remain valid while the archive is open.
 * Returns NULL on failure.
 */
xlsx_zip_archive *xlsx_zip_open_memory(const void *data, size_t size);

/*
 * Close a ZIP archive and free resources.
 */
void xlsx_zip_close(xlsx_zip_archive *archive);

/*
 * Get the number of files in the archive.
 */
size_t xlsx_zip_get_num_files(xlsx_zip_archive *archive);

/*
 * Get the filename at a given index.
 * Returns the length written (excluding null terminator), or 0 on error.
 */
size_t xlsx_zip_get_filename(xlsx_zip_archive *archive, size_t file_index,
                             char *filename_buf, size_t buf_size);

/*
 * Find a file by name in the archive.
 * Returns file index, or (size_t)-1 if not found.
 */
size_t xlsx_zip_locate_file(xlsx_zip_archive *archive, const char *filename);

/*
 * Get the uncompressed size of a file.
 * Returns 0 on error.
 */
size_t xlsx_zip_get_file_size(xlsx_zip_archive *archive, size_t file_index);

/*
 * Extract a file to a newly allocated buffer.
 * Caller is responsible for freeing the returned buffer.
 * Sets *out_size to the extracted size.
 * Returns NULL on failure.
 */
void *xlsx_zip_extract_to_heap(xlsx_zip_archive *archive, size_t file_index,
                               size_t *out_size);

/*
 * Extract a file to a provided buffer.
 * Returns true on success, false on failure.
 */
bool xlsx_zip_extract_to_buffer(xlsx_zip_archive *archive, size_t file_index,
                                void *buffer, size_t buffer_size);

/*
 * Extract a file by name to a newly allocated buffer.
 * Convenience function combining locate + extract.
 * Returns NULL on failure.
 */
void *xlsx_zip_extract_file(xlsx_zip_archive *archive, const char *filename,
                            size_t *out_size);

/* ============================================================================
 * Common XLSX Archive Paths
 * ============================================================================ */

/* Core XLSX file paths within the archive */
#define XLSX_PATH_WORKBOOK      "xl/workbook.xml"
#define XLSX_PATH_SHARED_STRINGS "xl/sharedStrings.xml"
#define XLSX_PATH_STYLES        "xl/styles.xml"
#define XLSX_PATH_CONTENT_TYPES "[Content_Types].xml"
#define XLSX_PATH_RELS          "xl/_rels/workbook.xml.rels"

/* Worksheet path template (sheet1.xml, sheet2.xml, etc.) */
#define XLSX_WORKSHEET_FMT      "xl/worksheets/sheet%d.xml"

#endif /* CIMPORT_XLSX_ZIP_H */
