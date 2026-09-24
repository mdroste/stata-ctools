/*
 * cexport_parse.c
 * Argument parsing and setup utilities for cexport
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <fcntl.h>
#if defined(_WIN32)
#include <io.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "cexport_parse.h"
#include "cexport_context.h"

/* ========================================================================
   Context Management
   ======================================================================== */

void cexport_context_init(cexport_context *ctx)
{
    memset(ctx, 0, sizeof(*ctx));
    ctools_filtered_data_init(&ctx->filtered);

    /* Set defaults */
    ctx->delimiter = ',';
    ctx->write_header = true;
    ctx->quote_strings = false;
    ctx->quote_if_needed = true;
    ctx->verbose = false;

    /* Line ending defaults: LF (Unix-style) */
    ctx->use_crlf = false;
    ctx->line_ending[0] = '\n';
    ctx->line_ending[1] = '\0';
    ctx->line_ending[2] = '\0';
    ctx->line_ending_len = 1;

    /* I/O backend defaults */
    ctx->io_backend = CEXPORT_IO_PWRITE;
    ctx->io_flags = CEXPORT_IO_FLAG_NONE;
    ctx->use_parallel_io = true;

    /* Adaptive sizing defaults */
    ctx->actual_avg_row_size = 0;
    ctx->adaptive_chunk_size = CTOOLS_EXPORT_CHUNK_SIZE;
}

void cexport_context_cleanup(cexport_context *ctx)
{
    if (ctx->filename != NULL) {
        free(ctx->filename);
        ctx->filename = NULL;
    }

    if (ctx->fp != NULL) {
        fclose(ctx->fp);
        ctx->fp = NULL;
    }

    if (ctx->varnames != NULL) {
        for (size_t j = 0; j < ctx->nvars; j++) {
            if (ctx->varnames[j] != NULL) {
                free(ctx->varnames[j]);
            }
        }
        free(ctx->varnames);
        ctx->varnames = NULL;
    }

    if (ctx->vartypes != NULL) {
        free(ctx->vartypes);
        ctx->vartypes = NULL;
    }

    ctools_arena_free(&ctx->chunk_arena);
    ctools_filtered_data_free(&ctx->filtered);
}

/* ========================================================================
   Argument Parsing
   ======================================================================== */

int cexport_parse_args(cexport_context *ctx, const char *args)
{
    char *args_copy = strdup(args);
    if (args_copy == NULL) return -1;

    char *saveptr;
    char *token;
    int arg_idx = 0;
    if (cexport_read_local("filename", &ctx->filename) || !ctx->filename[0]) {
        free(args_copy);
        return 198;
    }

    token = strtok_r(args_copy, " \t", &saveptr);
    while (token != NULL) {
        if (arg_idx == 0) {
            /* Delimiter is the first option token. */
            if (strlen(token) == 1) {
                ctx->delimiter = token[0];
            } else if (strcmp(token, "tab") == 0) {
                ctx->delimiter = '\t';
            }
        } else {
            /* Options */
            if (strcmp(token, "replace") == 0) {
                ctx->replace = true;
            } else if (strcmp(token, "noheader") == 0) {
                ctx->write_header = false;
            } else if (strcmp(token, "quote") == 0) {
                ctx->quote_strings = true;
            } else if (strcmp(token, "noquoteif") == 0) {
                ctx->quote_if_needed = false;
            } else if (strcmp(token, "verbose") == 0) {
                ctx->verbose = true;
            } else if (strcmp(token, "crlf") == 0) {
                /* Windows-style line endings */
                ctx->use_crlf = true;
                ctx->line_ending[0] = '\r';
                ctx->line_ending[1] = '\n';
                ctx->line_ending[2] = '\0';
                ctx->line_ending_len = 2;
            } else if (strcmp(token, "mmap") == 0) {
                /* Use memory-mapped I/O backend */
                ctx->io_backend = CEXPORT_IO_MMAP;
            } else if (strcmp(token, "noparallel") == 0) {
                /* Disable parallel I/O (for debugging/comparison) */
                ctx->use_parallel_io = false;
            } else if (strcmp(token, "nofsync") == 0) {
                /* Skip final fsync for faster but less durable writes */
                ctx->io_flags |= CEXPORT_IO_FLAG_NOFSYNC;
            } else if (strcmp(token, "direct") == 0) {
                /* Direct I/O bypasses OS cache (for very large files) */
                ctx->io_flags |= CEXPORT_IO_FLAG_DIRECT;
            } else if (strcmp(token, "prefault") == 0) {
                /* Pre-fault mmap pages to avoid page fault latency */
                ctx->io_flags |= CEXPORT_IO_FLAG_PREFAULT;
            }
        }

        arg_idx++;
        token = strtok_r(NULL, " \t", &saveptr);
    }

    free(args_copy);

    if (ctx->filename == NULL) {
        SF_error("cexport: no output filename specified\n");
        return -1;
    }

    return 0;
}

/* User text is never tokenized. The ado supplies a byte length, allowing
 * truncation by SF_macro_use to be detected before any output is opened. */
ST_retcode cexport_read_local(const char *field, char **value)
{
    char name[96], length_text[32];
    int length;
    *value = NULL;
    snprintf(name, sizeof(name), "___cexport_%s_len", field);
    if (SF_macro_use(name, length_text, sizeof(length_text)) ||
        !ctools_safe_atoi(length_text, &length) || length < 0 || length > 1048576)
        return 198;
    char *buf = malloc((size_t)length + 1);
    if (!buf) return 920;
    memset(buf, 0, (size_t)length + 1);
    snprintf(name, sizeof(name), "___cexport_%s", field);
    if (SF_macro_use(name, buf, length + 1) || strlen(buf) != (size_t)length) {
        free(buf);
        return 198;
    }
    *value = buf;
    return 0;
}

ST_retcode cexport_column_metadata(size_t column, char **name, int *type, int *date)
{
    char macro[96], value[64];
    *name = NULL;
    snprintf(macro, sizeof(macro), "___cexport_name_%zu", column + 1);
    if (SF_macro_use(macro, value, sizeof(value)) || !value[0] || strlen(value) > 32)
        return 198;
    *name = strdup(value);
    if (!*name) return 920;
    snprintf(macro, sizeof(macro), "___cexport_type_%zu", column + 1);
    if (SF_macro_use(macro, value, sizeof(value)) || !ctools_safe_atoi(value, type) ||
        *type < 0 || *type > 5) return 198;
    snprintf(macro, sizeof(macro), "___cexport_date_%zu", column + 1);
    if (SF_macro_use(macro, value, sizeof(value)) || !ctools_safe_atoi(value, date) ||
        *date < 0 || *date > 2) return 198;
    return 0;
}

int cexport_load_varnames(cexport_context *ctx)
{
    ctx->nvars = SF_nvars();
    ctx->varnames = ctools_safe_calloc2(ctx->nvars, sizeof(char *));
    ctx->vartypes = ctools_safe_malloc2(ctx->nvars, sizeof(vartype_t));
    if (!ctx->varnames || !ctx->vartypes) return 920;
    for (size_t j = 0; j < ctx->nvars; j++) {
        int type, date;
        ST_retcode rc = cexport_column_metadata(j, &ctx->varnames[j], &type, &date);
        if (rc) return rc;
        ctx->vartypes[j] = (vartype_t)type;
    }
    return 0;
}

int cexport_load_vartypes(cexport_context *ctx)
{
    /* Names and types are loaded and validated together. */
    return ctx->vartypes ? 0 : 198;
}

/* Reserve a sibling file; publishing is a separate, atomic operation. */
ST_retcode cexport_output_prepare(cexport_output *out, const char *filename, bool replace)
{
    memset(out, 0, sizeof(*out));
    out->replace = replace;
    out->destination = strdup(filename);
    out->temporary = malloc(strlen(filename) + 20);
    if (!out->destination || !out->temporary) { cexport_output_cleanup(out); return 920; }
    sprintf(out->temporary, "%s.ctools.XXXXXX", filename);
#ifdef _WIN32
    int fd = -1;
    if (_mktemp_s(out->temporary, strlen(out->temporary) + 1) == 0)
        fd = _open(out->temporary, _O_CREAT | _O_EXCL | _O_BINARY | _O_WRONLY,
                   _S_IREAD | _S_IWRITE);
    if (fd >= 0) _close(fd);
#else
    int fd = mkstemp(out->temporary);
    if (fd >= 0) close(fd);
#endif
    if (fd < 0) { cexport_output_cleanup(out); return 603; }
    out->owns_temporary = true;
    return 0;
}

ST_retcode cexport_output_commit(cexport_output *out)
{
#ifdef _WIN32
    if (!MoveFileExA(out->temporary, out->destination,
            MOVEFILE_WRITE_THROUGH | (out->replace ? MOVEFILE_REPLACE_EXISTING : 0)))
        return GetLastError() == ERROR_ALREADY_EXISTS || GetLastError() == ERROR_FILE_EXISTS ? 602 : 603;
#else
    if (out->replace) {
        if (rename(out->temporary, out->destination)) return 603;
    } else {
        /* link() refuses an existing destination, including a racing writer. */
        if (link(out->temporary, out->destination)) return errno == EEXIST ? 602 : 603;
        unlink(out->temporary);
    }
#endif
    return 0;
}

void cexport_output_cleanup(cexport_output *out)
{
    if (out->owns_temporary) remove(out->temporary);
    free(out->temporary); free(out->destination);
    memset(out, 0, sizeof(*out));
}
