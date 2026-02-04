#!/bin/bash
# Build standalone ASan test for civreghdfe

set -e

echo "Building civreghdfe ASan test harness..."

# Use clang with ASan
CC=clang

# Flags - force include mock header to define types before stplugin.h
CFLAGS="-fsanitize=address -fno-omit-frame-pointer -g -O1"
CFLAGS="$CFLAGS -include src/test_stplugin_mock.h"  # Force include mock first
CFLAGS="$CFLAGS -DSTPLUGIN_H"  # Then skip real stplugin.h via include guard
CFLAGS="$CFLAGS -DSD_FASTMODE"
CFLAGS="$CFLAGS -I/opt/homebrew/opt/libomp/include"

# OpenMP flags
OMPFLAGS="-Xpreprocessor -fopenmp"

# Include paths
INCLUDES="-Isrc -Isrc/civreghdfe -Isrc/creghdfe -Isrc/csort -Isrc/cmerge"

# Source files (excluding stplugin.c which has the real SPI)
SRCS="
    src/test_civreghdfe_asan.c
    src/ctools_types.c
    src/ctools_timer.c
    src/ctools_arena.c
    src/ctools_threads.c
    src/ctools_hash.c
    src/ctools_ols.c
    src/ctools_hdfe_utils.c
    src/ctools_data_io.c
    src/ctools_cleanup.c
    src/ctools_spi_checked.c
    src/civreghdfe/civreghdfe_impl.c
    src/civreghdfe/civreghdfe_estimate.c
    src/civreghdfe/civreghdfe_matrix.c
    src/civreghdfe/civreghdfe_vce.c
    src/civreghdfe/civreghdfe_tests.c
    src/creghdfe/creghdfe_hdfe.c
    src/creghdfe/creghdfe_solver.c
    src/creghdfe/creghdfe_utils.c
"

# Libraries
LIBS="-L/opt/homebrew/opt/libomp/lib -lomp -lm"

# Build
$CC $CFLAGS $OMPFLAGS $INCLUDES $SRCS $LIBS -o test_civreghdfe_asan

echo "Build complete: test_civreghdfe_asan"
echo ""
echo "Run with: ./test_civreghdfe_asan [num_iterations]"
echo "Default is 150 iterations"
