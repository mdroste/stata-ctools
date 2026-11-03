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

#include "formats/ctools_output.h"

ST_retcode cexport_read_local(const char *field, char **value);
ST_retcode cexport_column_metadata(size_t column, char **name, int *type, int *date);

/* Excel cell style codes that cexport.ado passes per column in the date
   field, assigned from display formats as export excel assigns them. The
   writers keep cell formats in this order: General, m/d/yy, a datetime with
   milliseconds, m/d/yy h:mm, h:mm:ss AM/PM, 0 for %tg, then 0, 0.00,
   #,##0, #,##0.00 and 0.00E+00 for other numbers. Codes 1-4 are dates
   (Excel serial values). Codes 1-5, the %t formats, keep their format under
   keepcellfmt and for a missing() number, as in export excel. */
#define CEXPORT_STYLE_COUNT 11
#define CEXPORT_STYLE_IS_DATE(code) ((code) >= 1 && (code) <= 4)
#define CEXPORT_STYLE_IS_TFORMAT(code) ((code) >= 1 && (code) <= 5)

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
    Load variable names and storage types from Stata macros.
    The .ado supplies independent local metadata for every column.

    @param ctx  Context to populate (nvars, varnames, vartypes)
    @return     0 on success, -1 on error
*/
int cexport_load_varnames(cexport_context *ctx);

#endif /* CEXPORT_PARSE_H */
