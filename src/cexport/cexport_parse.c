/*
 * cexport_parse.c
 * Argument parsing and setup utilities for cexport
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

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

    token = strtok_r(args_copy, " \t", &saveptr);
    while (token != NULL) {
        if (arg_idx == 0) {
            /* First arg: filename */
            ctx->filename = strdup(token);
            if (ctx->filename == NULL) {
                free(args_copy);
                return -1;
            }
        } else if (arg_idx == 1) {
            /* Second arg: delimiter */
            if (strlen(token) == 1) {
                ctx->delimiter = token[0];
            } else if (strcmp(token, "tab") == 0) {
                ctx->delimiter = '\t';
            }
        } else {
            /* Options */
            if (strcmp(token, "noheader") == 0) {
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

int cexport_load_varnames(cexport_context *ctx)
{
    size_t nvars = SF_nvars();
    ctx->varnames = (char **)ctools_safe_calloc2(nvars, sizeof(char *));
    if (ctx->varnames == NULL) return -1;

    ctx->nvars = nvars;

    /* Try to get variable names from macro */
    char varnames_buf[32768];
    size_t names_parsed = 0;
    if (SF_macro_use("CEXPORT_VARNAMES", varnames_buf, sizeof(varnames_buf)) == 0 &&
        strlen(varnames_buf) > 0) {
        /* Parse space-separated names */
        char *p = varnames_buf;
        for (size_t j = 0; j < nvars; j++) {
            /* Skip whitespace */
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\0') break;

            /* Find end of name */
            char *start = p;
            while (*p != ' ' && *p != '\t' && *p != '\0') p++;

            size_t len = p - start;
            ctx->varnames[j] = (char *)malloc(len + 1);
            if (ctx->varnames[j] == NULL) return -1;
            memcpy(ctx->varnames[j], start, len);
            ctx->varnames[j][len] = '\0';
            names_parsed++;
        }

        /* Fill remaining entries with fallback names if macro was short */
        for (size_t j = names_parsed; j < nvars; j++) {
            char namebuf[32];
            snprintf(namebuf, sizeof(namebuf), "v%zu", j + 1);
            ctx->varnames[j] = strdup(namebuf);
            if (ctx->varnames[j] == NULL) return -1;
        }
    } else {
        /* Fallback: use generic names */
        for (size_t j = 0; j < nvars; j++) {
            char namebuf[32];
            snprintf(namebuf, sizeof(namebuf), "v%zu", j + 1);
            ctx->varnames[j] = strdup(namebuf);
            if (ctx->varnames[j] == NULL) return -1;
        }
    }

    return 0;
}

int cexport_load_vartypes(cexport_context *ctx)
{
    size_t nvars = ctx->nvars;
    ctx->vartypes = (vartype_t *)ctools_safe_malloc2(nvars, sizeof(vartype_t));
    if (ctx->vartypes == NULL) return -1;

    /* Initialize to double (safest default - most precision) */
    for (size_t j = 0; j < nvars; j++) {
        ctx->vartypes[j] = VARTYPE_DOUBLE;
    }

    /* Try to get variable types from macro */
    char vartypes_buf[32768];
    if (SF_macro_use("CEXPORT_VARTYPES", vartypes_buf, sizeof(vartypes_buf)) == 0 &&
        strlen(vartypes_buf) > 0) {
        /* Parse space-separated type codes */
        char *p = vartypes_buf;
        for (size_t j = 0; j < nvars; j++) {
            /* Skip whitespace */
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\0') break;

            /* Parse type code */
            int type_code = (int)strtol(p, &p, 10);
            if (type_code >= 0 && type_code <= 5) {
                ctx->vartypes[j] = (vartype_t)type_code;
            }
        }
    }

    return 0;
}
