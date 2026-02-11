/*
 * cimport_encoding.h
 * Character encoding detection and conversion for cimport
 *
 * Supports encodings compatible with Stata's import delimited:
 * - UTF-8 (default)
 * - UTF-16LE, UTF-16BE
 * - ISO-8859-1 (Latin-1)
 * - Windows-1252 (CP1252)
 * - ISO-8859-15 (Latin-9)
 * - ASCII
 */

#ifndef CIMPORT_ENCODING_H
#define CIMPORT_ENCODING_H

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

/* Supported encodings */
typedef enum {
    CIMPORT_ENC_UNKNOWN = 0,
    CIMPORT_ENC_UTF8,
    CIMPORT_ENC_UTF8_BOM,
    CIMPORT_ENC_UTF16LE,
    CIMPORT_ENC_UTF16BE,
    CIMPORT_ENC_ISO_8859_1,     /* Latin-1 */
    CIMPORT_ENC_ISO_8859_15,    /* Latin-9 (with Euro sign) */
    CIMPORT_ENC_WINDOWS_1252,   /* CP1252 */
    CIMPORT_ENC_ASCII,
    CIMPORT_ENC_MACROMAN,       /* Mac OS Roman */
    CIMPORT_ENC_COUNT
} CImportEncoding;

/* Encoding detection result */
typedef struct {
    CImportEncoding encoding;
    int bom_length;             /* Length of BOM to skip (0 if none) */
    float confidence;           /* Detection confidence (0.0 - 1.0) */
} CImportEncodingDetection;

/*
 * Detect encoding from file data.
 * Uses BOM detection first, then heuristics if no BOM found.
 *
 * @param data      File data buffer
 * @param size      Size of data buffer
 * @return          Detection result with encoding and BOM length
 */
CImportEncodingDetection cimport_detect_encoding(const char *data, size_t size);

/*
 * Parse encoding name string to enum.
 * Accepts various formats: "utf-8", "UTF8", "utf8", "ISO-8859-1", "latin1", etc.
 *
 * @param name      Encoding name string
 * @return          Encoding enum value, or CIMPORT_ENC_UNKNOWN if not recognized
 */
CImportEncoding cimport_parse_encoding_name(const char *name);

/*
 * Get canonical name for encoding (for display/logging).
 *
 * @param enc       Encoding enum value
 * @return          Canonical name string (e.g., "UTF-8", "ISO-8859-1")
 */
const char *cimport_encoding_name(CImportEncoding enc);

/*
 * Check if encoding requires conversion to UTF-8.
 *
 * @param enc       Encoding enum value
 * @return          true if conversion needed, false if already UTF-8 compatible
 */
bool cimport_encoding_needs_conversion(CImportEncoding enc);

/*
 * Convert data from source encoding to UTF-8.
 * Allocates new buffer for converted data.
 *
 * @param data          Source data buffer
 * @param size          Size of source data
 * @param encoding      Source encoding
 * @param bom_skip      Bytes to skip at start (BOM)
 * @param out_data      Output: pointer to converted data (caller must free)
 * @param out_size      Output: size of converted data
 * @return              0 on success, -1 on error
 */
int cimport_convert_to_utf8(const char *data, size_t size, CImportEncoding encoding,
                            int bom_skip, char **out_data, size_t *out_size);

/*
 * Strip invalid UTF-8 sequences from data buffer.
 * Copies only valid UTF-8 bytes from src to dst.
 * dst must be at least src_size bytes.
 *
 * @param src       Source data buffer
 * @param src_size  Size of source data
 * @param dst       Output buffer (must be at least src_size bytes)
 * @return          Number of bytes written to dst
 */
size_t cimport_strip_invalid_utf8(const char *src, size_t src_size, char *dst);

/*
 * Check if data contains any invalid UTF-8 sequences.
 *
 * @param data      Data buffer
 * @param size      Size of data buffer
 * @return          true if invalid UTF-8 found, false if all valid
 */
bool cimport_has_invalid_utf8(const char *data, size_t size);

/*
 * SIMD-accelerated check if all bytes in data are ASCII (< 0x80).
 * Uses OR-reduction with SIMD vectors for high throughput.
 * Returns early on first non-ASCII byte found.
 *
 * @param data      Data buffer
 * @param size      Size of data buffer
 * @return          true if all bytes are ASCII, false otherwise
 */
bool cimport_is_ascii_simd(const char *data, size_t size);

#endif /* CIMPORT_ENCODING_H */
