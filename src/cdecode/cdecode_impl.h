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
 * Main entry point for cdecode command.
 *
 * Called from ctools_plugin.c dispatcher.
 *
 * Arguments format (space-separated):
 *   src_idx: 1-based index of source numeric variable
 *   dst_idx: 1-based index of destination string variable
 *   maxlen: maximum string length (0 = auto)
 *   labels=<encoded>: value-label pairs in format "value|label||value|label||..."
 *
 * @param args  Command arguments as space-separated string
 * @return      0 on success, Stata error code on failure
 */
ST_retcode cdecode_main(const char *args);

#endif /* CDECODE_IMPL_H */
