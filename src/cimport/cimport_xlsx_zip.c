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
 * Fast Extraction (libdeflate-accelerated)
 * ============================================================================ */

#include "libdeflate/libdeflate.h"

void *xlsx_zip_extract_fast(xlsx_zip_archive *archive, size_t file_index,
                            size_t *out_size)
{
    if (!archive) return NULL;

    mz_uint num_files = mz_zip_reader_get_num_files(&archive->zip);
    if (file_index >= num_files) return NULL;

    /* Get file stats */
    mz_zip_archive_file_stat stat;
    if (!mz_zip_reader_file_stat(&archive->zip, (mz_uint)file_index, &stat))
        return NULL;

    size_t uncomp_size = (size_t)stat.m_uncomp_size;
    size_t comp_size = (size_t)stat.m_comp_size;

    if (stat.m_method == 0) {
        /* STORE method — no compression, use standard extraction */
        return xlsx_zip_extract_to_heap(archive, file_index, out_size);
    }

    if (stat.m_method != MZ_DEFLATED) {
        /* Unknown method — fall back to miniz */
        return xlsx_zip_extract_to_heap(archive, file_index, out_size);
    }

    /* Extract raw compressed data using miniz */
    void *comp_buf = mz_zip_reader_extract_to_heap(&archive->zip,
                                                     (mz_uint)file_index,
                                                     &comp_size,
                                                     MZ_ZIP_FLAG_COMPRESSED_DATA);
    if (!comp_buf) {
        /* Fall back to standard miniz extraction */
        return xlsx_zip_extract_to_heap(archive, file_index, out_size);
    }

    /* Allocate output buffer */
    void *uncomp_buf = malloc(uncomp_size + 1);  /* +1 for potential NUL */
    if (!uncomp_buf) {
        free(comp_buf);
        return NULL;
    }

    /* Decompress with libdeflate */
    struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    if (!d) {
        free(comp_buf);
        free(uncomp_buf);
        return xlsx_zip_extract_to_heap(archive, file_index, out_size);
    }

    size_t actual_size = 0;
    enum libdeflate_result result = libdeflate_deflate_decompress(
        d, comp_buf, comp_size, uncomp_buf, uncomp_size, &actual_size);
    libdeflate_free_decompressor(d);
    free(comp_buf);

    if (result != LIBDEFLATE_SUCCESS) {
        free(uncomp_buf);
        /* Fall back to miniz on decompression failure */
        return xlsx_zip_extract_to_heap(archive, file_index, out_size);
    }

    if (out_size) *out_size = actual_size;
    return uncomp_buf;
}

void *xlsx_zip_extract_direct(xlsx_zip_archive *archive, size_t file_index,
                              size_t *out_size)
{
    if (!archive) return NULL;

    mz_zip_archive *pZip = &archive->zip;
    mz_uint num_files = mz_zip_reader_get_num_files(pZip);
    if (file_index >= num_files) return NULL;

    /* Get file metadata */
    mz_zip_archive_file_stat stat;
    if (!mz_zip_reader_file_stat(pZip, (mz_uint)file_index, &stat))
        return xlsx_zip_extract_fast(archive, file_index, out_size);

    size_t uncomp_size = (size_t)stat.m_uncomp_size;
    size_t comp_size = (size_t)stat.m_comp_size;

    /* Only handle DEFLATE — for STORE or unknown methods, use fast path */
    if (stat.m_method != MZ_DEFLATED)
        return xlsx_zip_extract_fast(archive, file_index, out_size);

    /* Read the 30-byte local file header to compute data offset.
     * ZIP local header: signature(4) + version(2) + flags(2) + method(2) +
     * modtime(2) + moddate(2) + crc32(4) + comp_size(4) + uncomp_size(4) +
     * filename_len(2) + extra_len(2) = 30 bytes */
    mz_uint8 local_hdr[30];
    if (pZip->m_pRead(pZip->m_pIO_opaque, stat.m_local_header_ofs,
                       local_hdr, 30) != 30)
        return xlsx_zip_extract_fast(archive, file_index, out_size);

    /* Validate local header signature: 0x04034b50 (little-endian) */
    if (local_hdr[0] != 0x50 || local_hdr[1] != 0x4B ||
        local_hdr[2] != 0x03 || local_hdr[3] != 0x04)
        return xlsx_zip_extract_fast(archive, file_index, out_size);

    mz_uint16 filename_len = (mz_uint16)(local_hdr[26] | (local_hdr[27] << 8));
    mz_uint16 extra_len = (mz_uint16)(local_hdr[28] | (local_hdr[29] << 8));
    mz_uint64 data_ofs = stat.m_local_header_ofs + 30 + filename_len + extra_len;

    /* Allocate and read compressed data directly from archive */
    void *comp_buf = malloc(comp_size);
    if (!comp_buf)
        return xlsx_zip_extract_fast(archive, file_index, out_size);

    if (pZip->m_pRead(pZip->m_pIO_opaque, data_ofs, comp_buf, comp_size)
        != comp_size) {
        free(comp_buf);
        return xlsx_zip_extract_fast(archive, file_index, out_size);
    }

    /* Allocate output buffer (+1 for NUL) */
    void *uncomp_buf = malloc(uncomp_size + 1);
    if (!uncomp_buf) {
        free(comp_buf);
        return NULL;
    }

    /* Decompress with libdeflate */
    struct libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    if (!d) {
        free(comp_buf);
        free(uncomp_buf);
        return xlsx_zip_extract_fast(archive, file_index, out_size);
    }

    size_t actual_size = 0;
    enum libdeflate_result result = libdeflate_deflate_decompress(
        d, comp_buf, comp_size, uncomp_buf, uncomp_size, &actual_size);
    libdeflate_free_decompressor(d);
    free(comp_buf);

    if (result != LIBDEFLATE_SUCCESS) {
        free(uncomp_buf);
        return xlsx_zip_extract_fast(archive, file_index, out_size);
    }

    if (out_size) *out_size = actual_size;
    return uncomp_buf;
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
