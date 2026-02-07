/*
 * cimport_xlsx_zip.c
 * ZIP archive interface for XLSX import
 *
 * Wrapper around miniz for extracting files from XLSX archives.
 */

#include "cimport_xlsx_zip.h"
#include "miniz/miniz.h"

#include <stdlib.h>
#include <string.h>

/* The xlsx_zip_archive is actually a miniz mz_zip_archive */
struct xlsx_zip_archive {
    mz_zip_archive zip;
    bool is_memory;  /* Track whether opened from memory */
};

xlsx_zip_archive *xlsx_zip_open_file(const char *filename)
{
    if (!filename) return NULL;

    xlsx_zip_archive *archive = (xlsx_zip_archive *)calloc(1, sizeof(xlsx_zip_archive));
    if (!archive) return NULL;

    mz_zip_zero_struct(&archive->zip);
    archive->is_memory = false;

    if (!mz_zip_reader_init_file(&archive->zip, filename, 0)) {
        free(archive);
        return NULL;
    }

    return archive;
}

xlsx_zip_archive *xlsx_zip_open_memory(const void *data, size_t size)
{
    if (!data || size == 0) return NULL;

    xlsx_zip_archive *archive = (xlsx_zip_archive *)calloc(1, sizeof(xlsx_zip_archive));
    if (!archive) return NULL;

    mz_zip_zero_struct(&archive->zip);
    archive->is_memory = true;

    if (!mz_zip_reader_init_mem(&archive->zip, data, size, 0)) {
        free(archive);
        return NULL;
    }

    return archive;
}

void xlsx_zip_close(xlsx_zip_archive *archive)
{
    if (!archive) return;

    mz_zip_reader_end(&archive->zip);
    free(archive);
}

size_t xlsx_zip_get_num_files(xlsx_zip_archive *archive)
{
    if (!archive) return 0;
    return mz_zip_reader_get_num_files(&archive->zip);
}

size_t xlsx_zip_get_filename(xlsx_zip_archive *archive, size_t file_index,
                             char *filename_buf, size_t buf_size)
{
    if (!archive || !filename_buf || buf_size == 0) return 0;

    mz_uint num_files = mz_zip_reader_get_num_files(&archive->zip);
    if (file_index >= num_files) return 0;

    mz_uint len = mz_zip_reader_get_filename(&archive->zip, (mz_uint)file_index,
                                              filename_buf, (mz_uint)buf_size);

    /* mz_zip_reader_get_filename returns length including null terminator */
    return (len > 0) ? (len - 1) : 0;
}

size_t xlsx_zip_locate_file(xlsx_zip_archive *archive, const char *filename)
{
    if (!archive || !filename) return (size_t)-1;

    int idx = mz_zip_reader_locate_file(&archive->zip, filename, NULL, 0);
    if (idx < 0) return (size_t)-1;

    return (size_t)idx;
}

size_t xlsx_zip_get_file_size(xlsx_zip_archive *archive, size_t file_index)
{
    if (!archive) return 0;

    mz_uint num_files = mz_zip_reader_get_num_files(&archive->zip);
    if (file_index >= num_files) return 0;

    mz_zip_archive_file_stat stat;
    if (!mz_zip_reader_file_stat(&archive->zip, (mz_uint)file_index, &stat)) {
        return 0;
    }

    return (size_t)stat.m_uncomp_size;
}

void *xlsx_zip_extract_to_heap(xlsx_zip_archive *archive, size_t file_index,
                               size_t *out_size)
{
    if (!archive) return NULL;

    mz_uint num_files = mz_zip_reader_get_num_files(&archive->zip);
    if (file_index >= num_files) return NULL;

    size_t size = 0;
    void *data = mz_zip_reader_extract_to_heap(&archive->zip, (mz_uint)file_index,
                                                &size, 0);

    if (out_size) *out_size = size;
    return data;
}

bool xlsx_zip_extract_to_buffer(xlsx_zip_archive *archive, size_t file_index,
                                void *buffer, size_t buffer_size)
{
    if (!archive || !buffer) return false;

    mz_uint num_files = mz_zip_reader_get_num_files(&archive->zip);
    if (file_index >= num_files) return false;

    return mz_zip_reader_extract_to_mem(&archive->zip, (mz_uint)file_index,
                                         buffer, buffer_size, 0) != 0;
}

void *xlsx_zip_extract_file(xlsx_zip_archive *archive, const char *filename,
                            size_t *out_size)
{
    if (!archive || !filename) return NULL;

    size_t idx = xlsx_zip_locate_file(archive, filename);
    if (idx == (size_t)-1) return NULL;

    return xlsx_zip_extract_to_heap(archive, idx, out_size);
}

/* ============================================================================
 * Streaming Extraction
 * ============================================================================ */

struct xlsx_zip_stream {
    mz_zip_reader_extract_iter_state *iter;
};

xlsx_zip_stream *xlsx_zip_stream_open(xlsx_zip_archive *archive, size_t file_index)
{
    if (!archive) return NULL;

    mz_uint num_files = mz_zip_reader_get_num_files(&archive->zip);
    if (file_index >= num_files) return NULL;

    xlsx_zip_stream *stream = (xlsx_zip_stream *)malloc(sizeof(xlsx_zip_stream));
    if (!stream) return NULL;

    stream->iter = mz_zip_reader_extract_iter_new(&archive->zip, (mz_uint)file_index, 0);
    if (!stream->iter) {
        free(stream);
        return NULL;
    }

    return stream;
}

size_t xlsx_zip_stream_read(xlsx_zip_stream *stream, void *buffer, size_t buffer_size)
{
    if (!stream || !stream->iter || !buffer || buffer_size == 0) return 0;
    return mz_zip_reader_extract_iter_read(stream->iter, buffer, buffer_size);
}

void xlsx_zip_stream_close(xlsx_zip_stream *stream)
{
    if (!stream) return;
    if (stream->iter) {
        mz_zip_reader_extract_iter_free(stream->iter);
    }
    free(stream);
}
