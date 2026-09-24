/*
 * cexport_parse.h
 * Argument parsing and setup utilities for cexport
 *
 * Contains functions for:
 * - Command argument parsing
 * - Variable name and type loading from Stata
 * - Context initialization and cleanup
 */

#ifndef CEXPORT_PARSE_H
#define CEXPORT_PARSE_H

#include "cexport_context.h"

typedef struct {
    char *destination;
    char *temporary;
    bool replace;
    bool owns_temporary;
} cexport_output;

ST_retcode cexport_read_local(const char *field, char **value);
ST_retcode cexport_column_metadata(size_t column, char **name, int *type, int *date);
ST_retcode cexport_output_prepare(cexport_output *out, const char *filename, bool replace);
ST_retcode cexport_output_commit(cexport_output *out);
void cexport_output_cleanup(cexport_output *out);

/* ========================================================================
   Context Management
   ======================================================================== */

/*
    Initialize a cexport context with default values.
    Must be called before any other operations on the context.

    @param ctx  Context to initialize
*/
void cexport_context_init(cexport_context *ctx);

/*
    Free all resources in context.
    Closes any open file handles and frees allocated memory.

    @param ctx  Context to cleanup
*/
void cexport_context_cleanup(cexport_context *ctx);

/* ========================================================================
   Argument Parsing
   ======================================================================== */

/*
    Parse command arguments.

    Format: delimiter [options...]; filename is a length-checked local
    Options: noheader, quote, noquoteif, verbose, crlf, mmap,
             noparallel, nofsync, direct, prefault

    @param ctx   Context to populate
    @param args  Space-separated argument string
    @return      0 on success, -1 on error
*/
int cexport_parse_args(cexport_context *ctx, const char *args);

/*
    Load variable names from Stata macro.
    The .ado supplies independent local metadata for every column.

    @param ctx  Context to populate (nvars and varnames)
    @return     0 on success, -1 on error
*/
int cexport_load_varnames(cexport_context *ctx);

/*
    Load variable types from Stata macro.
    The .ado supplies an independent local type code for each column:
    0=string, 1=byte, 2=int, 3=long, 4=float, 5=double

    @param ctx  Context to populate (vartypes array)
    @return     0 on success, -1 on error
*/
int cexport_load_vartypes(cexport_context *ctx);

#endif /* CEXPORT_PARSE_H */
