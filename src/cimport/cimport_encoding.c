/*
 * cimport_encoding.c
 * Character encoding detection and conversion for cimport
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "cimport_encoding.h"

/* ============================================================================
 * Windows-1252 to Unicode mapping for bytes 0x80-0x9F
 * (these differ from ISO-8859-1)
 * ============================================================================ */

static const uint32_t cp1252_special[32] = {
    0x20AC, /* 0x80: Euro sign */
    0x0081, /* 0x81: undefined, use as-is */
    0x201A, /* 0x82: Single low-9 quotation mark */
    0x0192, /* 0x83: Latin small letter f with hook */
    0x201E, /* 0x84: Double low-9 quotation mark */
    0x2026, /* 0x85: Horizontal ellipsis */
    0x2020, /* 0x86: Dagger */
    0x2021, /* 0x87: Double dagger */
    0x02C6, /* 0x88: Modifier letter circumflex accent */
    0x2030, /* 0x89: Per mille sign */
    0x0160, /* 0x8A: Latin capital letter S with caron */
    0x2039, /* 0x8B: Single left-pointing angle quotation */
    0x0152, /* 0x8C: Latin capital ligature OE */
    0x008D, /* 0x8D: undefined */
    0x017D, /* 0x8E: Latin capital letter Z with caron */
    0x008F, /* 0x8F: undefined */
    0x0090, /* 0x90: undefined */
    0x2018, /* 0x91: Left single quotation mark */
    0x2019, /* 0x92: Right single quotation mark */
    0x201C, /* 0x93: Left double quotation mark */
    0x201D, /* 0x94: Right double quotation mark */
    0x2022, /* 0x95: Bullet */
    0x2013, /* 0x96: En dash */
    0x2014, /* 0x97: Em dash */
    0x02DC, /* 0x98: Small tilde */
    0x2122, /* 0x99: Trade mark sign */
    0x0161, /* 0x9A: Latin small letter s with caron */
    0x203A, /* 0x9B: Single right-pointing angle quotation */
    0x0153, /* 0x9C: Latin small ligature oe */
    0x009D, /* 0x9D: undefined */
    0x017E, /* 0x9E: Latin small letter z with caron */
    0x0178, /* 0x9F: Latin capital letter Y with diaeresis */
};

/* ISO-8859-15 (Latin-9) differences from ISO-8859-1 */
static const uint32_t latin9_special[8][2] = {
    {0xA4, 0x20AC}, /* Euro sign (replaces currency sign) */
    {0xA6, 0x0160}, /* S with caron */
    {0xA8, 0x0161}, /* s with caron */
    {0xB4, 0x017D}, /* Z with caron */
    {0xB8, 0x017E}, /* z with caron */
    {0xBC, 0x0152}, /* OE ligature */
    {0xBD, 0x0153}, /* oe ligature */
    {0xBE, 0x0178}, /* Y with diaeresis */
};

/* ============================================================================
 * Helper: Encode Unicode codepoint as UTF-8
 * ============================================================================ */

static int codepoint_to_utf8(uint32_t cp, char *out)
{
    if (cp < 0x80) {
        out[0] = (char)cp;
        return 1;
    } else if (cp < 0x800) {
        out[0] = (char)(0xC0 | (cp >> 6));
        out[1] = (char)(0x80 | (cp & 0x3F));
        return 2;
    } else if (cp < 0x10000) {
        out[0] = (char)(0xE0 | (cp >> 12));
        out[1] = (char)(0x80 | ((cp >> 6) & 0x3F));
        out[2] = (char)(0x80 | (cp & 0x3F));
        return 3;
    } else if (cp < 0x110000) {
        out[0] = (char)(0xF0 | (cp >> 18));
        out[1] = (char)(0x80 | ((cp >> 12) & 0x3F));
        out[2] = (char)(0x80 | ((cp >> 6) & 0x3F));
        out[3] = (char)(0x80 | (cp & 0x3F));
        return 4;
    }
    /* Invalid codepoint, use replacement character */
    out[0] = (char)0xEF;
    out[1] = (char)0xBF;
    out[2] = (char)0xBD;
    return 3;
}

/* ============================================================================
 * Encoding Detection
 * ============================================================================ */

CImportEncodingDetection cimport_detect_encoding(const char *data, size_t size)
{
    CImportEncodingDetection result = {CIMPORT_ENC_UTF8, 0, 0.5f};

    if (size < 2) {
        result.encoding = CIMPORT_ENC_ASCII;
        result.confidence = 0.3f;
        return result;
    }

    const unsigned char *bytes = (const unsigned char *)data;

    /* Check for BOM (Byte Order Mark) */
    if (size >= 3 && bytes[0] == 0xEF && bytes[1] == 0xBB && bytes[2] == 0xBF) {
        result.encoding = CIMPORT_ENC_UTF8_BOM;
        result.bom_length = 3;
        result.confidence = 1.0f;
        return result;
    }

    if (size >= 2 && bytes[0] == 0xFF && bytes[1] == 0xFE) {
        /* Could be UTF-16LE or UTF-32LE */
        if (size >= 4 && bytes[2] == 0x00 && bytes[3] == 0x00) {
            /* UTF-32LE - not commonly used, treat as UTF-16LE for now */
        }
        result.encoding = CIMPORT_ENC_UTF16LE;
        result.bom_length = 2;
        result.confidence = 1.0f;
        return result;
    }

    if (size >= 2 && bytes[0] == 0xFE && bytes[1] == 0xFF) {
        result.encoding = CIMPORT_ENC_UTF16BE;
        result.bom_length = 2;
        result.confidence = 1.0f;
        return result;
    }

    /* No BOM - use heuristics */

    /* Check for null bytes (indicates UTF-16 without BOM) */
    size_t null_count = 0;
    size_t null_even = 0;
    size_t null_odd = 0;
    size_t check_size = size > 4096 ? 4096 : size;

    for (size_t i = 0; i < check_size; i++) {
        if (bytes[i] == 0x00) {
            null_count++;
            if (i % 2 == 0) null_even++;
            else null_odd++;
        }
    }

    /* If many nulls in alternating positions, likely UTF-16 */
    if (null_count > check_size / 8) {
        if (null_odd > null_even * 2) {
            result.encoding = CIMPORT_ENC_UTF16LE;
            result.confidence = 0.7f;
            return result;
        } else if (null_even > null_odd * 2) {
            result.encoding = CIMPORT_ENC_UTF16BE;
            result.confidence = 0.7f;
            return result;
        }
    }

    /* Check UTF-8 validity and high-byte patterns */
    int utf8_valid = 1;
    int high_bytes = 0;
    int cp1252_special_chars = 0;
    int utf8_multibyte_sequences = 0;

    for (size_t i = 0; i < check_size && utf8_valid; i++) {
        unsigned char c = bytes[i];

        if (c >= 0x80) {
            high_bytes++;

            /* Check for Windows-1252 special characters (0x80-0x9F) */
            if (c >= 0x80 && c <= 0x9F) {
                cp1252_special_chars++;
            }

            /* Validate UTF-8 sequence */
            if ((c & 0xE0) == 0xC0) {
                /* 2-byte sequence */
                if (i + 1 < check_size && (bytes[i+1] & 0xC0) == 0x80) {
                    utf8_multibyte_sequences++;
                    i++;
                } else {
                    utf8_valid = 0;
                }
            } else if ((c & 0xF0) == 0xE0) {
                /* 3-byte sequence */
                if (i + 2 < check_size &&
                    (bytes[i+1] & 0xC0) == 0x80 &&
                    (bytes[i+2] & 0xC0) == 0x80) {
                    utf8_multibyte_sequences++;
                    i += 2;
                } else {
                    utf8_valid = 0;
                }
            } else if ((c & 0xF8) == 0xF0) {
                /* 4-byte sequence */
                if (i + 3 < check_size &&
                    (bytes[i+1] & 0xC0) == 0x80 &&
                    (bytes[i+2] & 0xC0) == 0x80 &&
                    (bytes[i+3] & 0xC0) == 0x80) {
                    utf8_multibyte_sequences++;
                    i += 3;
                } else {
                    utf8_valid = 0;
                }
            } else if ((c & 0xC0) == 0x80) {
                /* Unexpected continuation byte */
                utf8_valid = 0;
            } else {
                /* Invalid UTF-8 lead byte (0xFE, 0xFF) */
                utf8_valid = 0;
            }
        }
    }

    /* Determine encoding based on analysis */
    if (high_bytes == 0) {
        result.encoding = CIMPORT_ENC_ASCII;
        result.confidence = 0.95f;
        return result;
    }

    if (utf8_valid && utf8_multibyte_sequences > 0) {
        result.encoding = CIMPORT_ENC_UTF8;
        result.confidence = 0.9f;
        return result;
    }

    if (!utf8_valid) {
        /* Not valid UTF-8, likely Windows-1252 or ISO-8859-1 */
        if (cp1252_special_chars > 0) {
            result.encoding = CIMPORT_ENC_WINDOWS_1252;
            result.confidence = 0.8f;
        } else {
            result.encoding = CIMPORT_ENC_ISO_8859_1;
            result.confidence = 0.7f;
        }
        return result;
    }

    /* Default to UTF-8 */
    result.encoding = CIMPORT_ENC_UTF8;
    result.confidence = 0.6f;
    return result;
}

/* ============================================================================
 * Encoding Name Parsing
 * ============================================================================ */

CImportEncoding cimport_parse_encoding_name(const char *name)
{
    if (!name || !*name) return CIMPORT_ENC_UNKNOWN;

    /* Normalize: lowercase, remove hyphens/underscores */
    char normalized[64];
    int j = 0;
    for (int i = 0; name[i] && j < 63; i++) {
        char c = name[i];
        if (c != '-' && c != '_' && c != ' ') {
            normalized[j++] = (c >= 'A' && c <= 'Z') ? c + 32 : c;
        }
    }
    normalized[j] = '\0';

    /* Match against known encodings */
    if (strcmp(normalized, "utf8") == 0 ||
        strcmp(normalized, "utf8bom") == 0) {
        return CIMPORT_ENC_UTF8;
    }

    if (strcmp(normalized, "utf16le") == 0 ||
        strcmp(normalized, "utf16") == 0 ||
        strcmp(normalized, "unicode") == 0) {
        return CIMPORT_ENC_UTF16LE;
    }

    if (strcmp(normalized, "utf16be") == 0) {
        return CIMPORT_ENC_UTF16BE;
    }

    if (strcmp(normalized, "iso88591") == 0 ||
        strcmp(normalized, "latin1") == 0 ||
        strcmp(normalized, "l1") == 0) {
        return CIMPORT_ENC_ISO_8859_1;
    }

    if (strcmp(normalized, "iso885915") == 0 ||
        strcmp(normalized, "latin9") == 0 ||
        strcmp(normalized, "l9") == 0) {
        return CIMPORT_ENC_ISO_8859_15;
    }

    if (strcmp(normalized, "windows1252") == 0 ||
        strcmp(normalized, "cp1252") == 0 ||
        strcmp(normalized, "win1252") == 0) {
        return CIMPORT_ENC_WINDOWS_1252;
    }

    if (strcmp(normalized, "ascii") == 0 ||
        strcmp(normalized, "usascii") == 0) {
        return CIMPORT_ENC_ASCII;
    }

    if (strcmp(normalized, "macroman") == 0 ||
        strcmp(normalized, "macintosh") == 0 ||
        strcmp(normalized, "xmacroman") == 0) {
        return CIMPORT_ENC_MACROMAN;
    }

    return CIMPORT_ENC_UNKNOWN;
}

const char *cimport_encoding_name(CImportEncoding enc)
{
    switch (enc) {
        case CIMPORT_ENC_UTF8:       return "UTF-8";
        case CIMPORT_ENC_UTF8_BOM:   return "UTF-8 (BOM)";
        case CIMPORT_ENC_UTF16LE:    return "UTF-16LE";
        case CIMPORT_ENC_UTF16BE:    return "UTF-16BE";
        case CIMPORT_ENC_ISO_8859_1: return "ISO-8859-1";
        case CIMPORT_ENC_ISO_8859_15: return "ISO-8859-15";
        case CIMPORT_ENC_WINDOWS_1252: return "Windows-1252";
        case CIMPORT_ENC_ASCII:      return "ASCII";
        case CIMPORT_ENC_MACROMAN:   return "MacRoman";
        default:                     return "Unknown";
    }
}

bool cimport_encoding_needs_conversion(CImportEncoding enc)
{
    switch (enc) {
        case CIMPORT_ENC_UTF8:
        case CIMPORT_ENC_UTF8_BOM:
        case CIMPORT_ENC_ASCII:
            return false;
        default:
            return true;
    }
}

/* ============================================================================
 * Single-byte to UTF-8 Converters
 * ============================================================================ */

int cimport_latin1_to_utf8(unsigned char byte, char *out)
{
    if (byte < 0x80) {
        out[0] = (char)byte;
        return 1;
    }
    /* ISO-8859-1 bytes 0x80-0xFF map directly to Unicode codepoints 0x80-0xFF */
    out[0] = (char)(0xC0 | (byte >> 6));
    out[1] = (char)(0x80 | (byte & 0x3F));
    return 2;
}

int cimport_cp1252_to_utf8(unsigned char byte, char *out)
{
    if (byte < 0x80) {
        out[0] = (char)byte;
        return 1;
    }

    uint32_t codepoint;
    if (byte >= 0x80 && byte <= 0x9F) {
        /* Special Windows-1252 characters */
        codepoint = cp1252_special[byte - 0x80];
    } else {
        /* Same as ISO-8859-1 for 0xA0-0xFF */
        codepoint = byte;
    }

    return codepoint_to_utf8(codepoint, out);
}

int cimport_latin9_to_utf8(unsigned char byte, char *out)
{
    if (byte < 0x80) {
        out[0] = (char)byte;
        return 1;
    }

    /* Check for Latin-9 specific characters */
    for (int i = 0; i < 8; i++) {
        if (byte == latin9_special[i][0]) {
            return codepoint_to_utf8(latin9_special[i][1], out);
        }
    }

    /* Same as ISO-8859-1 for other characters */
    return cimport_latin1_to_utf8(byte, out);
}

/* Mac Roman to Unicode mapping (selected characters that differ from Latin-1) */
static const uint32_t macroman_to_unicode[128] = {
    0x00C4, 0x00C5, 0x00C7, 0x00C9, 0x00D1, 0x00D6, 0x00DC, 0x00E1,
    0x00E0, 0x00E2, 0x00E4, 0x00E3, 0x00E5, 0x00E7, 0x00E9, 0x00E8,
    0x00EA, 0x00EB, 0x00ED, 0x00EC, 0x00EE, 0x00EF, 0x00F1, 0x00F3,
    0x00F2, 0x00F4, 0x00F6, 0x00F5, 0x00FA, 0x00F9, 0x00FB, 0x00FC,
    0x2020, 0x00B0, 0x00A2, 0x00A3, 0x00A7, 0x2022, 0x00B6, 0x00DF,
    0x00AE, 0x00A9, 0x2122, 0x00B4, 0x00A8, 0x2260, 0x00C6, 0x00D8,
    0x221E, 0x00B1, 0x2264, 0x2265, 0x00A5, 0x00B5, 0x2202, 0x2211,
    0x220F, 0x03C0, 0x222B, 0x00AA, 0x00BA, 0x03A9, 0x00E6, 0x00F8,
    0x00BF, 0x00A1, 0x00AC, 0x221A, 0x0192, 0x2248, 0x2206, 0x00AB,
    0x00BB, 0x2026, 0x00A0, 0x00C0, 0x00C3, 0x00D5, 0x0152, 0x0153,
    0x2013, 0x2014, 0x201C, 0x201D, 0x2018, 0x2019, 0x00F7, 0x25CA,
    0x00FF, 0x0178, 0x2044, 0x20AC, 0x2039, 0x203A, 0xFB01, 0xFB02,
    0x2021, 0x00B7, 0x201A, 0x201E, 0x2030, 0x00C2, 0x00CA, 0x00C1,
    0x00CB, 0x00C8, 0x00CD, 0x00CE, 0x00CF, 0x00CC, 0x00D3, 0x00D4,
    0xF8FF, 0x00D2, 0x00DA, 0x00DB, 0x00D9, 0x0131, 0x02C6, 0x02DC,
    0x00AF, 0x02D8, 0x02D9, 0x02DA, 0x00B8, 0x02DD, 0x02DB, 0x02C7,
};

static int cimport_macroman_to_utf8(unsigned char byte, char *out)
{
    if (byte < 0x80) {
        out[0] = (char)byte;
        return 1;
    }
    return codepoint_to_utf8(macroman_to_unicode[byte - 0x80], out);
}

/* ============================================================================
 * UTF-16 to UTF-8 Conversion
 * ============================================================================ */

int cimport_utf16_to_utf8(const char *data, size_t size, bool is_little_endian,
                          char **out_data, size_t *out_size)
{
    const unsigned char *bytes = (const unsigned char *)data;

    /* Allocate output buffer (worst case: each UTF-16 unit -> 3 UTF-8 bytes) */
    /* Overflow check: (size/2) * 3 + 1 */
    size_t units = size / 2;
    if (units > (SIZE_MAX - 1) / 3) {
        return -1;  /* Overflow */
    }
    size_t max_out = units * 3 + 1;
    char *output = malloc(max_out);
    if (!output) return -1;

    size_t out_pos = 0;
    size_t i = 0;

    while (i + 1 < size) {
        uint32_t unit;
        if (is_little_endian) {
            unit = bytes[i] | ((uint32_t)bytes[i+1] << 8);
        } else {
            unit = ((uint32_t)bytes[i] << 8) | bytes[i+1];
        }
        i += 2;

        uint32_t codepoint;

        /* Check for surrogate pair */
        if (unit >= 0xD800 && unit <= 0xDBFF) {
            /* High surrogate - need low surrogate */
            if (i + 1 < size) {
                uint32_t low;
                if (is_little_endian) {
                    low = bytes[i] | ((uint32_t)bytes[i+1] << 8);
                } else {
                    low = ((uint32_t)bytes[i] << 8) | bytes[i+1];
                }

                if (low >= 0xDC00 && low <= 0xDFFF) {
                    /* Valid surrogate pair */
                    codepoint = 0x10000 + ((unit - 0xD800) << 10) + (low - 0xDC00);
                    i += 2;
                } else {
                    /* Invalid - use replacement character */
                    codepoint = 0xFFFD;
                }
            } else {
                codepoint = 0xFFFD;
            }
        } else if (unit >= 0xDC00 && unit <= 0xDFFF) {
            /* Unexpected low surrogate */
            codepoint = 0xFFFD;
        } else {
            codepoint = unit;
        }

        /* Encode to UTF-8 */
        out_pos += codepoint_to_utf8(codepoint, output + out_pos);
    }

    output[out_pos] = '\0';
    *out_data = output;
    *out_size = out_pos;
    return 0;
}

/* ============================================================================
 * Main Conversion Function
 * ============================================================================ */

int cimport_convert_to_utf8(const char *data, size_t size, CImportEncoding encoding,
                            int bom_skip, char **out_data, size_t *out_size)
{
    /* Skip BOM if present */
    data += bom_skip;
    size -= bom_skip;

    /* Handle UTF-16 separately */
    if (encoding == CIMPORT_ENC_UTF16LE) {
        return cimport_utf16_to_utf8(data, size, true, out_data, out_size);
    }
    if (encoding == CIMPORT_ENC_UTF16BE) {
        return cimport_utf16_to_utf8(data, size, false, out_data, out_size);
    }

    /* Single-byte encodings: allocate worst-case buffer */
    /* ISO-8859-1/Windows-1252: each byte can become up to 3 UTF-8 bytes */
    /* Overflow check: size * 3 + 1 */
    if (size > (SIZE_MAX - 1) / 3) {
        return -1;  /* Overflow */
    }
    size_t max_out = size * 3 + 1;
    char *output = malloc(max_out);
    if (!output) return -1;

    const unsigned char *src = (const unsigned char *)data;
    size_t out_pos = 0;

    for (size_t i = 0; i < size; i++) {
        unsigned char c = src[i];
        int bytes_written;

        switch (encoding) {
            case CIMPORT_ENC_ISO_8859_1:
                bytes_written = cimport_latin1_to_utf8(c, output + out_pos);
                break;

            case CIMPORT_ENC_ISO_8859_15:
                bytes_written = cimport_latin9_to_utf8(c, output + out_pos);
                break;

            case CIMPORT_ENC_WINDOWS_1252:
                bytes_written = cimport_cp1252_to_utf8(c, output + out_pos);
                break;

            case CIMPORT_ENC_MACROMAN:
                bytes_written = cimport_macroman_to_utf8(c, output + out_pos);
                break;

            default:
                /* ASCII or UTF-8 - copy as-is */
                output[out_pos] = (char)c;
                bytes_written = 1;
                break;
        }

        out_pos += bytes_written;
    }

    output[out_pos] = '\0';
    *out_data = output;
    *out_size = out_pos;
    return 0;
}
