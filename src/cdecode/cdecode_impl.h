/*
 * cdecode_impl.h
 *
 * High-performance numeric-to-string decoding for Stata
 * Part of the ctools suite
 *
 * Counterpart to cencode - converts a labeled numeric variable
 * to a string variable containing the label text.
 */

#ifndef CDECODE_IMPL_H
#define CDECODE_IMPL_H

#include "stplugin.h"

/*
 * Scan a Stata label save file and return the max label length.
 * Used for pre-flight detection so the .ado can create the
 * destination string variable with the correct width.
 *
 * Stores _cdecode_maxlen scalar in Stata.
 *
 * @param args  "labelsfile=/path/to/file"
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cdecode_scan_main(const char *args);

/*
 * Main entry point for cdecode command.
 *
 * Called from ctools_plugin.c dispatcher.
 * Parses a Stata `label save` format file directly (no Mata needed).
 * Uses flat array lookup for dense label values, hash table for sparse.
 *
 * @param args  "maxlen=N labelsfile=/path/to/file"
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cdecode_main(const char *args);

#endif /* CDECODE_IMPL_H */
