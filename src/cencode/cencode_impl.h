/*
 * cencode_impl.h
 *
 * High-performance string encoding for Stata
 * Part of the ctools suite
 *
 * Replaces "encode" with a C-accelerated implementation featuring:
 * - Parallel unique string collection using thread-local hash sets
 * - Lock-free parallel encoding phase
 * - Efficient string interning with arena allocator
 */

#ifndef CENCODE_IMPL_H
#define CENCODE_IMPL_H

#include "stplugin.h"

/*
 * Main entry point for cencode command.
 *
 * Arguments passed via args string:
 *   "var_idx gen_idx [label=name] [noextend]"
 *
 * Where:
 *   var_idx  - 1-based index of the source string variable
 *   gen_idx  - 1-based index of the destination numeric variable (pre-created by .ado)
 *
 * The plugin:
 *   1. Reads the source string variable
 *   2. Builds a mapping of unique strings to integers (1, 2, 3, ...)
 *   3. Writes the encoded values to the destination variable
 *   4. Stores the label mapping in macros for the .ado to read
 *
 * Returns:
 *   0 on success, Stata error code on failure
 */
ST_retcode cencode_main(const char *args);

#endif /* CENCODE_IMPL_H */
