/* BIFF8 in an OLE compound file. Files are published atomically only after the
   complete workbook and FAT/DIFAT chains have been written. */
#include "cexport_xls.h"
#include "cexport/cexport_parse.h"
#include "ctools_runtime.h"
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static void put16(unsigned char *p, unsigned v) {
    p[0] = (unsigned char)v;
    p[1] = (unsigned char)(v >> 8);
}
static void put32(unsigned char *p, uint32_t v) {
    for (int j = 0; j < 4; j++)
        p[j] = (unsigned char)(v >> (8 * j));
}
static void put64(unsigned char *p, uint64_t v) {
    for (int j = 0; j < 8; j++)
        p[j] = (unsigned char)(v >> (8 * j));
}
static int record(FILE *f, unsigned id, const void *data, unsigned size) {
    unsigned char h[4];
    if (size > 8224)
        return 109;
    put16(h, id);
    put16(h + 2, size);
    return fwrite(h, 1, 4, f) != 4 || (size && fwrite(data, 1, size, f) != size) ? 603 : 0;
}
static int bof(FILE *f, int sheet) {
    unsigned char b[16] = {0};
    put16(b, 0x600);
    put16(b + 2, sheet ? 0x10 : 5);
    put16(b + 4, 0xdbb);
    put16(b + 6, 1997);
    put32(b + 8, 0x41);
    put32(b + 12, 6);
    return record(f, 0x809, b, 16);
}
static int utf16(const char *s, unsigned char *out, size_t cap, size_t *count) {
    size_t n = 0;
    const unsigned char *p = (const unsigned char *)s;
    while (*p) {
        uint32_t cp;
        unsigned b = *p++;
        if (b < 0x80)
            cp = b;
        else {
            int k;
            if (b >= 0xc2 && b <= 0xdf) {
                cp = b & 31;
                k = 1;
            } else if (b >= 0xe0 && b <= 0xef) {
                cp = b & 15;
                k = 2;
            } else if (b >= 0xf0 && b <= 0xf4) {
                cp = b & 7;
                k = 3;
            } else
                return 198;
            for (int j = 0; j < k; j++) {
                if ((*p & 0xc0) != 0x80)
                    return 198;
                cp = (cp << 6) | (*p++ & 63);
            }
            if (cp > 0x10ffff || (cp >= 0xd800 && cp <= 0xdfff))
                return 198;
        }
        if (n + (cp >= 0x10000 ? 2 : 1) > cap)
            return 109;
        if (cp >= 0x10000) {
            cp -= 0x10000;
            put16(out + 2 * n++, 0xd800 + (cp >> 10));
            put16(out + 2 * n++, 0xdc00 + (cp & 1023));
        } else
            put16(out + 2 * n++, cp);
    }
    *count = n;
    return 0;
}
typedef struct {
    FILE *file;
    unsigned char data[8224];
    size_t length;
    unsigned count;
    long start;
    int first;
} xls_sst;
static int sst_flush(xls_sst *s) {
    int rc = record(s->file, s->first ? 0xfc : 0x3c, s->data, (unsigned)s->length);
    s->length = 0;
    s->first = 0;
    return rc;
}
static int sst_string(xls_sst *s, const char *text) {
    unsigned char *encoded = malloc(65534);
    if (!encoded)
        return 909;
    size_t n;
    int rc = utf16(text, encoded, 32767, &n);
    if (rc) {
        free(encoded);
        return rc;
    }
    if (sizeof(s->data) - s->length < 3 && (rc = sst_flush(s))) {
        free(encoded);
        return rc;
    }
    put16(s->data + s->length, (unsigned)n);
    s->data[s->length + 2] = 1;
    s->length += 3;
    s->count++;
    for (size_t i = 0; i < n; i++) {
        if (sizeof(s->data) - s->length < 2) {
            if ((rc = sst_flush(s))) {
                free(encoded);
                return rc;
            }
            s->data[0] = 1;
            s->length = 1;
        }
        memcpy(s->data + s->length, encoded + 2 * i, 2);
        s->length += 2;
    }
    free(encoded);
    return 0;
}
static int string_cell(FILE *f, unsigned row, unsigned col, unsigned index) {
    unsigned char b[10];
    put16(b, row);
    put16(b + 2, col);
    put16(b + 4, 0);
    put32(b + 6, index);
    return record(f, 0xfd, b, 10);
}
static int blank_cell(FILE *f, unsigned row, unsigned col) {
    unsigned char b[6] = {0};
    put16(b, row);
    put16(b + 2, col);
    return record(f, 0x201, b, 6);
}
static int string_value(int var, int obs, char **buffer, size_t *capacity) {
    size_t length = SF_var_is_strl(var) ? (size_t)SF_sdatalen(var, obs) : 2045;
    if (length > INT_MAX - 1 || (SF_var_is_strl(var) && SF_var_is_binary(var, obs)))
        return 109;
    if (length + 1 > *capacity) {
        char *tmp = realloc(*buffer, length + 1);
        if (!tmp)
            return 909;
        *buffer = tmp;
        *capacity = length + 1;
    }
    if (SF_var_is_strl(var)) {
        if (SF_strldata(var, obs, *buffer, (int)*capacity) != (int)length)
            return 609;
        (*buffer)[length] = 0;
        return 0;
    }
    return SF_sdata(var, obs, *buffer);
}
static void directory_entry(unsigned char *b, const char *name, unsigned type, unsigned child,
                            unsigned sector, uint64_t size) {
    memset(b, 0, 128);
    size_t n = strlen(name);
    for (size_t j = 0; j < n; j++)
        put16(b + 2 * j, (unsigned char)name[j]);
    put16(b + 64, (unsigned)(2 * (n + 1)));
    b[66] = (unsigned char)type;
    b[67] = 1;
    put32(b + 68, UINT32_MAX);
    put32(b + 72, UINT32_MAX);
    put32(b + 76, child);
    put32(b + 116, sector);
    put64(b + 120, size);
}
static int compound(FILE *out, FILE *biff, size_t bytes) {
    /* Avoid the mini stream by padding short workbooks to the 4096 byte cutoff. */
    size_t stream = bytes < 4096 ? 4096 : bytes, n = (stream + 511) / 512, nfat = 1, ndifat = 0;
    for (;;) {
        size_t next = (n + 1 + nfat + ndifat + 127) / 128,
               dif = next > 109 ? (next - 109 + 126) / 127 : 0;
        if (next == nfat && dif == ndifat)
            break;
        nfat = next;
        ndifat = dif;
    }
    size_t total = 1 + n + nfat + ndifat;
    if (total > UINT32_MAX)
        return 920;
    unsigned char b[512];
    memset(b, 0, 512);
    const unsigned char signature[8] = {0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1};
    memcpy(b, signature, 8);
    put16(b + 24, 0x3e);
    put16(b + 26, 3);
    put16(b + 28, 0xfffe);
    put16(b + 30, 9);
    put16(b + 32, 6);
    put32(b + 44, (uint32_t)nfat);
    put32(b + 48, 0);
    put32(b + 56, 4096);
    put32(b + 60, 0xfffffffe);
    put32(b + 68, ndifat ? (uint32_t)(1 + n + nfat) : 0xfffffffe);
    put32(b + 72, (uint32_t)ndifat);
    for (size_t i = 0; i < 109; i++)
        put32(b + 76 + 4 * i, i < nfat ? (uint32_t)(1 + n + i) : UINT32_MAX);
    if (fwrite(b, 1, 512, out) != 512)
        return 603;
    memset(b, 0, 512);
    directory_entry(b, "Root Entry", 5, 1, 0xfffffffe, 0);
    directory_entry(b + 128, "Workbook", 2, UINT32_MAX, 1, stream);
    if (fwrite(b, 1, 512, out) != 512)
        return 603;
    rewind(biff);
    size_t remaining = bytes;
    for (size_t i = 0; i < n; i++) {
        memset(b, 0, 512);
        size_t take = remaining < 512 ? remaining : 512;
        if (take && fread(b, 1, take, biff) != take)
            return 603;
        remaining -= take;
        if (fwrite(b, 1, 512, out) != 512)
            return 603;
    }
    for (size_t i = 0; i < nfat; i++) {
        for (size_t j = 0; j < 128; j++) {
            size_t sector = i * 128 + j;
            uint32_t next;
            if (sector >= total)
                next = UINT32_MAX;
            else if (!sector || sector == n)
                next = 0xfffffffe;
            else if (sector <= n)
                next = (uint32_t)sector + 1;
            else if (sector <= n + nfat)
                next = 0xfffffffd;
            else
                next = 0xfffffffc;
            put32(b + 4 * j, next);
        }
        if (fwrite(b, 1, 512, out) != 512)
            return 603;
    }
    for (size_t i = 0; i < ndifat; i++) {
        for (size_t j = 0; j < 127; j++) {
            size_t f = 109 + i * 127 + j;
            put32(b + 4 * j, f < nfat ? (uint32_t)(1 + n + f) : UINT32_MAX);
        }
        put32(b + 508, i + 1 < ndifat ? (uint32_t)(2 + n + nfat + i) : 0xfffffffe);
        if (fwrite(b, 1, 512, out) != 512)
            return 603;
    }
    return 0;
}
#include "cexport_xls_edit.inc"

ST_retcode cexport_xls_main(const char *args) {
    char *filename = NULL, *sheet = NULL, *opts = NULL, *missing = NULL;
    char **names = NULL;
    int *dates = NULL;
    int k = SF_nvars(), rc = 0, replace = 0, firstrow = 1, row0 = 0, col0 = 0;
    int action = 0, keep = 0, existing = 0, missing_numeric = 0;
    double missing_number = 0;
    xlsWorkBook *oldbook = NULL;
    char *text = NULL;
    size_t textcap = 0;
    unsigned stringindex = 0;
    size_t n = 0;
    FILE *biff = NULL, *out = NULL;
    cexport_output output;
    memset(&output, 0, sizeof(output));
    if ((rc = cexport_read_local("filename", &filename)))
        goto done;
    if ((rc = cexport_read_local("sheet", &sheet)))
        goto done;
    if (k < 1 || k > 256) {
        rc = 103;
        goto done;
    }
    if ((rc = cexport_read_local("missing", &missing)))
        goto done;
    char *missing_end;
    missing_number = strtod(missing, &missing_end);
    while (*missing_end && isspace((unsigned char)*missing_end))
        missing_end++;
    missing_numeric = missing_end != missing && !*missing_end && isfinite(missing_number);
    opts = strdup(args);
    names = calloc((size_t)k, sizeof(char *));
    dates = calloc((size_t)k, sizeof(int));
    if (!opts || !names || !dates) {
        rc = 909;
        goto done;
    }
    for (char *t = strtok(opts, " "); t; t = strtok(NULL, " ")) {
        if (!strcmp(t, "replace"))
            replace = 1;
        else if (!strcmp(t, "sheetmodify"))
            action = 1;
        else if (!strcmp(t, "sheetreplace"))
            action = 2;
        else if (!strcmp(t, "keepcellfmt"))
            keep = 1;
        else if (!strcmp(t, "nofirstrow"))
            firstrow = 0;
        else if (!strncmp(t, "cell=", 5)) {
            const char *p = t + 5;
            col0 = 0;
            while (*p >= 'A' && *p <= 'Z') {
                col0 = col0 * 26 + *p - 'A' + 1;
                p++;
            }
            col0--;
            row0 = atoi(p) - 1;
            if (row0 < 0 || col0 < 0) {
                rc = 198;
                goto done;
            }
        }
    }
    for (int i = SF_in1(); i <= SF_in2(); i++)
        if (SF_ifobs(i))
            n++;
    if (n + (size_t)row0 + (size_t)firstrow > 65536 || col0 + k > 256) {
        rc = 198;
        goto done;
    }
    for (int j = 0; j < k; j++) {
        int type;
        if ((rc = cexport_column_metadata((size_t)j, &names[j], &type, &dates[j])))
            goto done;
    }
    biff = tmpfile();
    if (!biff) {
        rc = 603;
        goto done;
    }
    if ((rc = bof(biff, 0)))
        goto done;
    unsigned char b[256] = {0};
    put16(b, 1200);
    if ((rc = record(biff, 0x42, b, 2)))
        goto done;
    memset(b, 0, 20);
    if ((rc = record(biff, 0xe0, b, 20)))
        goto done;
    put16(b + 2, 14);
    if ((rc = record(biff, 0xe0, b, 20)))
        goto done;
    put16(b + 2, 22);
    if ((rc = record(biff, 0xe0, b, 20)))
        goto done;
    xls_sst sst = {0};
    sst.file = biff;
    sst.first = 1;
    sst.length = 8;
    sst.start = ftell(biff) + 4;
    if (firstrow)
        for (int j = 0; j < k; j++)
            if ((rc = sst_string(&sst, names[j])))
                goto done;
    for (int i = SF_in1(); i <= SF_in2(); i++)
        if (SF_ifobs(i))
            for (int j = 0; j < k; j++)
                if (SF_var_is_string(j + 1)) {
                    if ((rc = string_value(j + 1, i, &text, &textcap)))
                        goto done;
                    if (*text && (rc = sst_string(&sst, text)))
                        goto done;
                } else if (*missing && !missing_numeric) {
                    double value;
                    if ((rc = (_stata_)->safevdata(j + 1, i, &value)))
                        goto done;
                    if (SF_is_missing(value) && (rc = sst_string(&sst, missing)))
                        goto done;
                }
    if ((rc = sst_flush(&sst)))
        goto done;
    long aftersst = ftell(biff);
    if (fseek(biff, sst.start, SEEK_SET)) {
        rc = 603;
        goto done;
    }
    put32(b, sst.count);
    put32(b + 4, sst.count);
    if (fwrite(b, 1, 8, biff) != 8 || fseek(biff, aftersst, SEEK_SET)) {
        rc = 603;
        goto done;
    }
    memset(b, 0, sizeof(b));
    size_t chars;
    rc = utf16(*sheet ? sheet : "Sheet1", b + 8, 31, &chars);
    if (rc)
        goto done;
    b[6] = (unsigned char)chars;
    b[7] = 1;
    long bound = ftell(biff) + 4;
    if ((rc = record(biff, 0x85, b, (unsigned)(8 + 2 * chars))) ||
        (rc = record(biff, 0x0a, NULL, 0)))
        goto done;
    long start = ftell(biff);
    if (start < 0 || start > UINT32_MAX) {
        rc = 603;
        goto done;
    }
    if (fseek(biff, bound, SEEK_SET)) {
        rc = 603;
        goto done;
    }
    put32(b, (uint32_t)start);
    if (fwrite(b, 1, 4, biff) != 4 || fseek(biff, start, SEEK_SET)) {
        rc = 603;
        goto done;
    }
    if ((rc = bof(biff, 1)))
        goto done;
    memset(b, 0, 14);
    put32(b, (uint32_t)row0);
    put32(b + 4, (uint32_t)(n + row0 + firstrow));
    put16(b + 8, (unsigned)col0);
    put16(b + 10, (unsigned)(col0 + k));
    if ((rc = record(biff, 0x200, b, 14)))
        goto done;
    if (firstrow)
        for (int j = 0; j < k; j++)
            if ((rc = string_cell(biff, (unsigned)row0, (unsigned)(col0 + j), stringindex++)))
                goto done;
    unsigned r = (unsigned)(row0 + firstrow);
    for (int i = SF_in1(); i <= SF_in2(); i++) {
        if (!SF_ifobs(i))
            continue;
        if ((i & 4095) == 0 && SF_poll()) {
            rc = 1;
            goto done;
        }
        for (int j = 0; j < k; j++) {
            if (SF_var_is_string(j + 1)) {
                if ((rc = string_value(j + 1, i, &text, &textcap)))
                    goto done;
                if (*text && (rc = string_cell(biff, r, (unsigned)(col0 + j), stringindex++)))
                    goto done;
                if (!*text && (rc = blank_cell(biff, r, (unsigned)(col0 + j))))
                    goto done;
            } else {
                double x;
                if ((rc = (_stata_)->safevdata(j + 1, i, &x)))
                    goto done;
                if (SF_is_missing(x)) {
                    if (missing_numeric)
                        x = missing_number;
                    else {
                        if ((rc = *missing
                                      ? string_cell(biff, r, (unsigned)(col0 + j), stringindex++)
                                      : blank_cell(biff, r, (unsigned)(col0 + j))))
                            goto done;
                        continue;
                    }
                }
                put16(b, r);
                put16(b + 2, (unsigned)(col0 + j));
                put16(b + 4, (unsigned)dates[j]);
                uint64_t bits;
                memcpy(&bits, &x, 8);
                put64(b + 6, bits);
                if ((rc = record(biff, 0x203, b, 14)))
                    goto done;
            }
        }
        r++;
    }
    if ((rc = record(biff, 0x0a, NULL, 0)))
        goto done;
    long bytes = ftell(biff);
    if (bytes < 0) {
        rc = 603;
        goto done;
    }
    if (!replace) {
        FILE *source = fopen(filename, "rb");
        if (source) {
            existing = 1;
            fclose(source);
        }
    }
    if (existing) {
        FILE *edited = NULL;
        rc = xe_edit(filename, sheet, action, keep, row0, col0, n + (size_t)firstrow, k, biff,
                     &edited, &oldbook);
        if (rc) {
            if (edited)
                fclose(edited);
            goto done;
        }
        fclose(biff);
        biff = edited;
        bytes = ftell(biff);
        if (bytes < 0) {
            rc = 603;
            goto done;
        }
    }
    if ((rc = cexport_output_prepare(&output, filename, replace || existing)))
        goto done;
    out = fopen(output.temporary, "wb");
    if (!out) {
        rc = 603;
        goto done;
    }
    if ((rc = existing ? xe_compound(out, filename, oldbook, biff, (size_t)bytes)
                       : compound(out, biff, (size_t)bytes)))
        goto done;
    if (fclose(out)) {
        out = NULL;
        rc = 603;
        goto done;
    }
    out = NULL;
    rc = cexport_output_commit(&output);
done:
    if (oldbook)
        xls_close_WB(oldbook);
    if (biff)
        fclose(biff);
    if (out)
        fclose(out);
    cexport_output_cleanup(&output);
    if (names)
        for (int j = 0; j < k; j++)
            free(names[j]);
    free(names);
    free(dates);
    free(opts);
    free(sheet);
    free(filename);
    free(text);
    free(missing);
    return rc;
}
