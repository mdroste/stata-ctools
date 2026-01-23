# ==============================================================================
#                         ctools Stata Plugin Build System
# ==============================================================================
#
# Builds high-performance Stata plugins for multiple platforms:
#   - macOS Apple Silicon: ctools_mac_arm.plugin    (with OpenMP)
#   - macOS Intel:         ctools_mac_x86.plugin    (with OpenMP)
#   - Windows (x64):       ctools_windows.plugin    (with OpenMP via llvm-mingw)
#   - Linux (x64):         ctools_linux.plugin      (with OpenMP)
#
# Usage:
#   make              - Build for current platform
#   make all          - Build all distribution plugins
#   make native       - Build optimized plugin for current architecture
#   make macos-arm    - Build macOS Apple Silicon plugin
#   make macos-intel  - Build macOS Intel plugin
#   make windows      - Build Windows x64 plugin
#   make linux        - Build Linux x64 plugin
#   make check        - Check build dependencies
#   make clean        - Remove compiled files
#   make help         - Show this help
#
# Requirements:
#   macOS:   Xcode Command Line Tools, Homebrew with libomp
#   Windows: Native MSVC/MinGW, or cross-compile with llvm-mingw
#   Linux:   GCC with OpenMP (libgomp)
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Directory Configuration
# ------------------------------------------------------------------------------
SRC_DIR = src
BUILD_DIR = build

# ------------------------------------------------------------------------------
# Platform Detection
# ------------------------------------------------------------------------------
ifeq ($(OS),Windows_NT)
    DETECTED_OS := Windows
    DETECTED_ARCH := x86_64
else
    UNAME_S := $(shell uname -s)
    UNAME_M := $(shell uname -m)
    ifeq ($(UNAME_S),Darwin)
        DETECTED_OS := macOS
        DETECTED_ARCH := $(UNAME_M)
    else ifeq ($(UNAME_S),Linux)
        DETECTED_OS := Linux
        DETECTED_ARCH := $(UNAME_M)
    else
        DETECTED_OS := Unknown
        DETECTED_ARCH := Unknown
    endif
endif

# ------------------------------------------------------------------------------
# Source Files (organized by subdirectory)
# ------------------------------------------------------------------------------
# Core plugin infrastructure
CORE_SRCS = $(SRC_DIR)/stplugin.c \
            $(SRC_DIR)/ctools_plugin.c

# Common utilities
COMMON_SRCS = $(SRC_DIR)/ctools_types.c \
              $(SRC_DIR)/ctools_timer.c \
              $(SRC_DIR)/ctools_error.c \
              $(SRC_DIR)/ctools_threads.c \
              $(SRC_DIR)/ctools_arena.c \
              $(SRC_DIR)/ctools_data_io.c \
              $(SRC_DIR)/ctools_hdfe_utils.c \
              $(SRC_DIR)/ctools_ols.c \
              $(SRC_DIR)/ctools_sort_radix_lsd.c \
              $(SRC_DIR)/ctools_sort_radix_msd.c \
              $(SRC_DIR)/ctools_sort_timsort.c \
              $(SRC_DIR)/ctools_sort_sample.c \
              $(SRC_DIR)/ctools_sort_counting.c \
              $(SRC_DIR)/ctools_sort_merge.c \
              $(SRC_DIR)/ctools_sort_ips4o.c

# csort module
CSORT_SRCS = $(SRC_DIR)/csort/csort_impl.c \
             $(SRC_DIR)/csort/csort_stream.c

# cimport module (uses shared ctools_arena from ctools_arena.c)
CIMPORT_SRCS = $(SRC_DIR)/cimport/cimport_impl.c \
               $(SRC_DIR)/cimport/cimport_mmap.c \
               $(SRC_DIR)/cimport/cimport_parse.c \
               $(SRC_DIR)/cimport/cimport_encoding.c

# cexport module
CEXPORT_SRCS = $(SRC_DIR)/cexport/cexport_impl.c \
               $(SRC_DIR)/cexport/cexport_io.c \
               $(SRC_DIR)/cexport/cexport_format.c \
               $(SRC_DIR)/cexport/cexport_parse.c

# cmerge module
CMERGE_SRCS = $(SRC_DIR)/cmerge/cmerge_impl.c \
              $(SRC_DIR)/cmerge/cmerge_memory.c \
              $(SRC_DIR)/cmerge/cmerge_radix_sort.c \
              $(SRC_DIR)/cmerge/cmerge_keys.c \
              $(SRC_DIR)/cmerge/cmerge_group_search.c \
              $(SRC_DIR)/cmerge/cmerge_join.c \
              $(SRC_DIR)/cmerge/cmerge_io.c


# creghdfe module
CREGHDFE_SRCS = $(SRC_DIR)/creghdfe/creghdfe_types.c \
                $(SRC_DIR)/creghdfe/creghdfe_utils.c \
                $(SRC_DIR)/creghdfe/creghdfe_hdfe.c \
                $(SRC_DIR)/creghdfe/creghdfe_solver.c \
                $(SRC_DIR)/creghdfe/creghdfe_ols.c \
                $(SRC_DIR)/creghdfe/creghdfe_vce.c \
                $(SRC_DIR)/creghdfe/creghdfe_regress.c \
                $(SRC_DIR)/creghdfe/creghdfe_impl.c

# cqreg module
CQREG_SRCS = $(SRC_DIR)/cqreg/cqreg_types.c \
             $(SRC_DIR)/cqreg/cqreg_linalg.c \
             $(SRC_DIR)/cqreg/cqreg_blas.c \
             $(SRC_DIR)/cqreg/cqreg_ipm.c \
             $(SRC_DIR)/cqreg/cqreg_fn.c \
             $(SRC_DIR)/cqreg/cqreg_irls.c \
             $(SRC_DIR)/cqreg/cqreg_sparsity.c \
             $(SRC_DIR)/cqreg/cqreg_vce.c \
             $(SRC_DIR)/cqreg/cqreg_hdfe.c \
             $(SRC_DIR)/cqreg/cqreg_regress.c \
             $(SRC_DIR)/cqreg/cqreg_impl.c

# cbinscatter module
CBINSCATTER_SRCS = $(SRC_DIR)/cbinscatter/cbinscatter_impl.c \
                   $(SRC_DIR)/cbinscatter/cbinscatter_bins.c \
                   $(SRC_DIR)/cbinscatter/cbinscatter_resid.c \
                   $(SRC_DIR)/cbinscatter/cbinscatter_fit.c

# civreghdfe module
CIVREGHDFE_SRCS = $(SRC_DIR)/civreghdfe/civreghdfe_matrix.c \
                  $(SRC_DIR)/civreghdfe/civreghdfe_estimate.c \
                  $(SRC_DIR)/civreghdfe/civreghdfe_vce.c \
                  $(SRC_DIR)/civreghdfe/civreghdfe_tests.c \
                  $(SRC_DIR)/civreghdfe/civreghdfe_impl.c

# cencode module
CENCODE_SRCS = $(SRC_DIR)/cencode/cencode_impl.c

# cwinsor module
CWINSOR_SRCS = $(SRC_DIR)/cwinsor/cwinsor_impl.c

# cdestring module
CDESTRING_SRCS = $(SRC_DIR)/cdestring/cdestring_impl.c

# cdecode module
CDECODE_SRCS = $(SRC_DIR)/cdecode/cdecode_impl.c

# All sources
SOURCES = $(CORE_SRCS) $(COMMON_SRCS) $(CSORT_SRCS) $(CIMPORT_SRCS) $(CEXPORT_SRCS) $(CMERGE_SRCS) $(CREGHDFE_SRCS) $(CQREG_SRCS) $(CBINSCATTER_SRCS) $(CIVREGHDFE_SRCS) $(CENCODE_SRCS) $(CWINSOR_SRCS) $(CDESTRING_SRCS) $(CDECODE_SRCS)

# Headers (for dependency tracking)
CORE_HEADERS = $(SRC_DIR)/stplugin.h \
               $(SRC_DIR)/ctools_types.h \
               $(SRC_DIR)/ctools_timer.h \
               $(SRC_DIR)/ctools_error.h \
               $(SRC_DIR)/ctools_threads.h \
               $(SRC_DIR)/ctools_arena.h \
               $(SRC_DIR)/ctools_hdfe_utils.h \
               $(SRC_DIR)/ctools_ols.h \
               $(SRC_DIR)/ctools_config.h \
               $(SRC_DIR)/ctools_spi.h

CSORT_HEADERS = $(SRC_DIR)/csort/csort_impl.h \
                $(SRC_DIR)/csort/csort_stream.h

CIMPORT_HEADERS = $(SRC_DIR)/cimport/cimport_impl.h \
                  $(SRC_DIR)/cimport/cimport_context.h \
                  $(SRC_DIR)/cimport/cimport_mmap.h \
                  $(SRC_DIR)/cimport/cimport_parse.h \
                  $(SRC_DIR)/cimport/cimport_encoding.h

CEXPORT_HEADERS = $(SRC_DIR)/cexport/cexport_impl.h \
                  $(SRC_DIR)/cexport/cexport_io.h \
                  $(SRC_DIR)/cexport/cexport_context.h \
                  $(SRC_DIR)/cexport/cexport_format.h \
                  $(SRC_DIR)/cexport/cexport_parse.h

CMERGE_HEADERS = $(SRC_DIR)/cmerge/cmerge_impl.h \
                 $(SRC_DIR)/cmerge/cmerge_memory.h \
                 $(SRC_DIR)/cmerge/cmerge_radix_sort.h \
                 $(SRC_DIR)/cmerge/cmerge_keys.h \
                 $(SRC_DIR)/cmerge/cmerge_group_search.h \
                 $(SRC_DIR)/cmerge/cmerge_join.h \
                 $(SRC_DIR)/cmerge/cmerge_io.h


CREGHDFE_HEADERS = $(SRC_DIR)/creghdfe/creghdfe_types.h \
                   $(SRC_DIR)/creghdfe/creghdfe_utils.h \
                   $(SRC_DIR)/creghdfe/creghdfe_hdfe.h \
                   $(SRC_DIR)/creghdfe/creghdfe_solver.h \
                   $(SRC_DIR)/creghdfe/creghdfe_ols.h \
                   $(SRC_DIR)/creghdfe/creghdfe_vce.h \
                   $(SRC_DIR)/creghdfe/creghdfe_regress.h \
                   $(SRC_DIR)/creghdfe/creghdfe_impl.h

CQREG_HEADERS = $(SRC_DIR)/cqreg/cqreg_types.h \
                $(SRC_DIR)/cqreg/cqreg_linalg.h \
                $(SRC_DIR)/cqreg/cqreg_blas.h \
                $(SRC_DIR)/cqreg/cqreg_ipm.h \
                $(SRC_DIR)/cqreg/cqreg_sparsity.h \
                $(SRC_DIR)/cqreg/cqreg_vce.h \
                $(SRC_DIR)/cqreg/cqreg_hdfe.h \
                $(SRC_DIR)/cqreg/cqreg_regress.h \
                $(SRC_DIR)/cqreg/cqreg_impl.h

CBINSCATTER_HEADERS = $(SRC_DIR)/cbinscatter/cbinscatter_types.h \
                      $(SRC_DIR)/cbinscatter/cbinscatter_bins.h \
                      $(SRC_DIR)/cbinscatter/cbinscatter_resid.h \
                      $(SRC_DIR)/cbinscatter/cbinscatter_fit.h \
                      $(SRC_DIR)/cbinscatter/cbinscatter_impl.h

CIVREGHDFE_HEADERS = $(SRC_DIR)/civreghdfe/civreghdfe_matrix.h \
                     $(SRC_DIR)/civreghdfe/civreghdfe_estimate.h \
                     $(SRC_DIR)/civreghdfe/civreghdfe_vce.h \
                     $(SRC_DIR)/civreghdfe/civreghdfe_tests.h \
                     $(SRC_DIR)/civreghdfe/civreghdfe_impl.h

CENCODE_HEADERS = $(SRC_DIR)/cencode/cencode_impl.h

CWINSOR_HEADERS = $(SRC_DIR)/cwinsor/cwinsor_impl.h

CDESTRING_HEADERS = $(SRC_DIR)/cdestring/cdestring_impl.h

CDECODE_HEADERS = $(SRC_DIR)/cdecode/cdecode_impl.h

# Include paths for subdirectories
INCLUDE_DIRS = -I$(SRC_DIR) -I$(SRC_DIR)/csort -I$(SRC_DIR)/cimport -I$(SRC_DIR)/cexport -I$(SRC_DIR)/cmerge -I$(SRC_DIR)/creshape -I$(SRC_DIR)/creghdfe -I$(SRC_DIR)/cqreg -I$(SRC_DIR)/cbinscatter -I$(SRC_DIR)/civreghdfe -I$(SRC_DIR)/cencode -I$(SRC_DIR)/cwinsor -I$(SRC_DIR)/cdestring -I$(SRC_DIR)/cdecode

# ------------------------------------------------------------------------------
# Output Configuration
# ------------------------------------------------------------------------------
PLUGIN_MAC_ARM    = $(BUILD_DIR)/ctools_mac_arm.plugin
PLUGIN_MAC_X86    = $(BUILD_DIR)/ctools_mac_x86.plugin
PLUGIN_WINDOWS    = $(BUILD_DIR)/ctools_windows.plugin
PLUGIN_LINUX      = $(BUILD_DIR)/ctools_linux.plugin
PLUGIN_NATIVE     = $(BUILD_DIR)/ctools.plugin

# ==============================================================================
#                              macOS Configuration
# ==============================================================================
CC_MAC = clang

# Find Homebrew libomp
LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null || echo "/opt/homebrew/opt/libomp")
LIBOMP_EXISTS := $(shell test -d $(LIBOMP_PREFIX) && echo yes || echo no)

# Base optimization flags for macOS
# SD_FASTMODE disables SPI bounds checking for max performance
MAC_BASE_FLAGS = -O3 -Wall -Wextra -fPIC -DSYSTEM=APPLEMAC -DSD_FASTMODE \
                 -ffast-math -funroll-loops -ftree-vectorize -flto \
                 -fno-strict-aliasing $(INCLUDE_DIRS)

# macOS Apple Silicon (arm64) with OpenMP and Accelerate
# Check for static libomp.a (preferred for bundling)
LIBOMP_STATIC_ARM := $(shell test -f $(LIBOMP_PREFIX)/lib/libomp.a && echo yes || echo no)

ifeq ($(LIBOMP_EXISTS),yes)
    CFLAGS_MAC_ARM = $(MAC_BASE_FLAGS) -arch arm64 \
                     -mmacosx-version-min=11.0 \
                     -mcpu=apple-m1 \
                     -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    ifeq ($(LIBOMP_STATIC_ARM),yes)
        # Static linking - bundle libomp into the plugin
        LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto $(LIBOMP_PREFIX)/lib/libomp.a \
                          -framework Accelerate
    else
        LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto -L$(LIBOMP_PREFIX)/lib -lomp \
                          -framework Accelerate
    endif
    MAC_ARM_HAS_OMP = yes
else
    CFLAGS_MAC_ARM = $(MAC_BASE_FLAGS) -arch arm64 \
                     -mmacosx-version-min=11.0 \
                     -mcpu=apple-m1
    LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto -framework Accelerate
    MAC_ARM_HAS_OMP = no
endif

# Check for x86_64 libomp (installed via Rosetta Homebrew)
LIBOMP_INTEL = /usr/local/opt/libomp
LIBOMP_INTEL_EXISTS := $(shell test -f $(LIBOMP_INTEL)/lib/libomp.dylib && \
                               file $(LIBOMP_INTEL)/lib/libomp.dylib 2>/dev/null | grep -q x86_64 && \
                               echo yes || echo no)
# Use lipo to check .a architecture (file command doesn't show arch for archives)
LIBOMP_STATIC_X86 := $(shell test -f $(LIBOMP_INTEL)/lib/libomp.a && \
                             lipo -info $(LIBOMP_INTEL)/lib/libomp.a 2>/dev/null | grep -q x86_64 && \
                             echo yes || echo no)

# macOS Intel (x86_64) with Accelerate
ifeq ($(LIBOMP_INTEL_EXISTS),yes)
    CFLAGS_MAC_X86 = $(MAC_BASE_FLAGS) -arch x86_64 \
                     -mmacosx-version-min=10.13 \
                     -march=x86-64 -mtune=haswell \
                     -Xpreprocessor -fopenmp -I$(LIBOMP_INTEL)/include
    ifeq ($(LIBOMP_STATIC_X86),yes)
        # Static linking - bundle libomp into the plugin
        LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto $(LIBOMP_INTEL)/lib/libomp.a \
                          -framework Accelerate
    else
        LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto -L$(LIBOMP_INTEL)/lib -lomp \
                          -framework Accelerate
    endif
    MAC_X86_HAS_OMP = yes
else
    CFLAGS_MAC_X86 = $(MAC_BASE_FLAGS) -arch x86_64 \
                     -mmacosx-version-min=10.13 \
                     -march=x86-64 -mtune=haswell
    LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto -framework Accelerate
    MAC_X86_HAS_OMP = no
endif

# ==============================================================================
#                             Windows Configuration
# ==============================================================================
# Try llvm-mingw first (has OpenMP), fall back to mingw-w64
LLVM_MINGW_PREFIX = /opt/llvm-mingw
LLVM_MINGW_EXISTS := $(shell test -x $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang && echo yes || echo no)

ifeq ($(DETECTED_OS),Windows)
    # Native Windows build - static link libgomp
    CC_WIN = gcc
    CFLAGS_WIN = -O3 -Wall -shared -DSYSTEM=STWIN32 -DSD_FASTMODE -fopenmp \
                 -ffast-math -funroll-loops -ftree-vectorize -fno-strict-aliasing $(INCLUDE_DIRS)
    LDFLAGS_WIN = -static-libgcc -Wl,-Bstatic -lgomp -Wl,-Bdynamic -lpthread
    WIN_HAS_OMP = yes
    WIN_OMP_STATIC = yes
else ifeq ($(LLVM_MINGW_EXISTS),yes)
    # Cross-compile with llvm-mingw (OpenMP requires DLL - bundle libomp.dll with plugin)
    CC_WIN = $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang
    CFLAGS_WIN = -O3 -Wall -shared -DSYSTEM=STWIN32 -DSD_FASTMODE -fopenmp \
                 -ffast-math -funroll-loops -ftree-vectorize -fno-strict-aliasing $(INCLUDE_DIRS)
    LDFLAGS_WIN = -fopenmp -lpthread
    WIN_HAS_OMP = yes
    WIN_OMP_STATIC = no
else
    # Cross-compile with mingw-w64 (no OpenMP)
    CC_WIN = x86_64-w64-mingw32-gcc
    MINGW_EXISTS := $(shell which x86_64-w64-mingw32-gcc 2>/dev/null && echo yes || echo no)
    CFLAGS_WIN = -O3 -Wall -Wno-unknown-pragmas -shared -DSYSTEM=STWIN32 -DSD_FASTMODE \
                 -ffast-math -funroll-loops -ftree-vectorize -fno-strict-aliasing $(INCLUDE_DIRS)
    LDFLAGS_WIN = -static-libgcc -lpthread
    WIN_HAS_OMP = no
    WIN_OMP_STATIC = no
endif

# ==============================================================================
#                              Linux Configuration
# ==============================================================================
ifeq ($(DETECTED_OS),Linux)
    CC_LINUX = gcc
else
    # Cross-compile attempt (usually won't work well)
    CC_LINUX = gcc
endif

LINUX_BASE_FLAGS = -O3 -Wall -Wextra -shared -fPIC -DSYSTEM=STUNIX -DSD_FASTMODE \
                   -ffast-math -funroll-loops -ftree-vectorize -flto -fno-strict-aliasing $(INCLUDE_DIRS)

# Check for OpenMP support on Linux
LINUX_HAS_OMP := $(shell which gcc >/dev/null 2>&1 && \
                         gcc -fopenmp -E - < /dev/null >/dev/null 2>&1 && \
                         echo yes || echo no)

# Note: Static linking of libgomp.a doesn't work on Linux for shared libraries
# because libgomp uses TLS relocations incompatible with PIC code.
# Linux plugins require libgomp.so at runtime.

# Check for OpenBLAS on Linux
LINUX_HAS_OPENBLAS := $(shell pkg-config --exists openblas 2>/dev/null && echo yes || \
                              (test -f /usr/include/cblas.h && test -f /usr/lib/libopenblas.so && echo yes) || \
                              echo no)

ifeq ($(LINUX_HAS_OMP),yes)
    ifeq ($(LINUX_HAS_OPENBLAS),yes)
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -fopenmp -DHAVE_OPENBLAS
        LDFLAGS_LINUX = -flto -fopenmp -lopenblas
    else
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -fopenmp
        LDFLAGS_LINUX = -flto -fopenmp
    endif
else
    ifeq ($(LINUX_HAS_OPENBLAS),yes)
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -DHAVE_OPENBLAS
        LDFLAGS_LINUX = -flto -lopenblas
    else
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS)
        LDFLAGS_LINUX = -flto
    endif
endif

# ==============================================================================
#                                  Targets
# ==============================================================================
.PHONY: all native macos-arm macos-intel windows linux clean help check rebuild

# Default: build for current platform
ifeq ($(DETECTED_OS),macOS)
    ifeq ($(DETECTED_ARCH),arm64)
default: macos-arm
    else
default: macos-intel
    endif
else ifeq ($(DETECTED_OS),Windows)
default: windows
else ifeq ($(DETECTED_OS),Linux)
default: linux
else
default: help
endif

# ------------------------------------------------------------------------------
# Help
# ------------------------------------------------------------------------------
help:
	@echo ""
	@echo "======================================================================"
	@echo "                   ctools Stata Plugin Build System                   "
	@echo "======================================================================"
	@echo ""
	@echo "  High-performance C plugins for Stata: csort, cmerge, cimport, cexport, creghdfe, cqreg"
	@echo ""
	@echo "  Usage:"
	@echo "    make              Build for current platform ($(DETECTED_OS) $(DETECTED_ARCH))"
	@echo "    make all          Build all distribution plugins"
	@echo "    make native       Build optimized plugin for current architecture"
	@echo "    make macos-arm    Build macOS Apple Silicon plugin"
	@echo "    make macos-intel  Build macOS Intel plugin"
	@echo "    make windows      Build Windows x64 plugin"
	@echo "    make linux        Build Linux x64 plugin"
	@echo "    make check        Check build dependencies"
	@echo "    make clean        Remove compiled files"
	@echo "    make rebuild      Clean and rebuild all"
	@echo ""
	@echo "  Detected: $(DETECTED_OS) $(DETECTED_ARCH)"
	@echo ""

# ------------------------------------------------------------------------------
# Dependency Check
# ------------------------------------------------------------------------------
check:
	@echo ""
	@echo "======================================================================"
	@echo "                     Build Environment Check                          "
	@echo "======================================================================"
	@echo ""
	@echo "  Detected OS:   $(DETECTED_OS)"
	@echo "  Architecture:  $(DETECTED_ARCH)"
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  macOS Apple Silicon (arm64)"
	@echo "----------------------------------------------------------------------"
	@printf "    clang:     " && (which clang >/dev/null 2>&1 && printf "OK\n" || printf "MISSING - install Xcode CLT\n")
ifeq ($(LIBOMP_EXISTS),yes)
	@echo "    libomp:    OK ($(LIBOMP_PREFIX))"
	@echo "    OpenMP:    Enabled"
else
	@echo "    libomp:    MISSING - run: brew install libomp"
	@echo "    OpenMP:    Disabled"
endif
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  macOS Intel (x86_64)"
	@echo "----------------------------------------------------------------------"
	@printf "    clang:     " && (which clang >/dev/null 2>&1 && printf "OK\n" || printf "MISSING\n")
ifeq ($(MAC_X86_HAS_OMP),yes)
	@echo "    libomp:    OK (/usr/local/opt/libomp)"
	@echo "    OpenMP:    Enabled"
else
	@echo "    libomp:    MISSING - run: arch -x86_64 /usr/local/bin/brew install libomp"
	@echo "    OpenMP:    Disabled"
endif
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  Windows x64 (cross-compile)"
	@echo "----------------------------------------------------------------------"
ifeq ($(LLVM_MINGW_EXISTS),yes)
	@echo "    llvm-mingw: OK ($(LLVM_MINGW_PREFIX))"
	@echo "    OpenMP:     Enabled"
else
	@printf "    mingw-w64:  " && (which x86_64-w64-mingw32-gcc >/dev/null 2>&1 && printf "OK (no OpenMP)\n" || printf "MISSING - brew install mingw-w64\n")
	@echo "    llvm-mingw: MISSING (optional, enables OpenMP)"
	@echo "                  Download: https://github.com/mstorsjo/llvm-mingw/releases"
	@echo "                  Extract to: /opt/llvm-mingw"
endif
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  Linux x64"
	@echo "----------------------------------------------------------------------"
ifeq ($(DETECTED_OS),Linux)
	@printf "    gcc:       " && (which gcc >/dev/null 2>&1 && printf "OK\n" || printf "MISSING\n")
ifeq ($(LINUX_HAS_OMP),yes)
	@echo "    OpenMP:    Enabled"
else
	@echo "    OpenMP:    MISSING - install libgomp-dev"
endif
else
	@echo "    (Cross-compile from $(DETECTED_OS) not fully supported)"
	@printf "    gcc:       " && (which gcc >/dev/null 2>&1 && printf "OK\n" || printf "MISSING\n")
endif
	@echo ""

# ------------------------------------------------------------------------------
# Build All Platforms
# ------------------------------------------------------------------------------
all: $(BUILD_DIR)
	@echo ""
	@echo "======================================================================"
	@echo "                   ctools Multi-Platform Build                        "
	@echo "  Building: csort, cmerge, cimport, cexport, creghdfe, cqreg          "
	@echo "======================================================================"
	@$(MAKE) --no-print-directory macos-arm
	@$(MAKE) --no-print-directory macos-intel
	@$(MAKE) --no-print-directory windows
	@$(MAKE) --no-print-directory linux
	@echo ""
	@echo "======================================================================"
	@echo "                        BUILD COMPLETE                                "
	@echo "======================================================================"
	@echo ""
	@echo "  Plugins created:"
	@ls -lh $(BUILD_DIR)/ctools_*.plugin 2>/dev/null | awk '{printf "    %-30s %s\n", $$9, $$5}'
	@echo ""

# Create build directory
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# ------------------------------------------------------------------------------
# macOS Apple Silicon (arm64)
# ------------------------------------------------------------------------------
macos-arm: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  [1/4] macOS Apple Silicon (arm64)"
	@echo "----------------------------------------------------------------------"
	@echo "    Compiler:  $(CC_MAC) -arch arm64"
ifeq ($(MAC_ARM_HAS_OMP),yes)
ifeq ($(LIBOMP_STATIC_ARM),yes)
	@echo "    OpenMP:    Enabled (static)"
else
	@echo "    OpenMP:    Enabled (dynamic)"
endif
	@echo "    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "    OpenMP:    Disabled (pthread only)"
	@echo "    Threads:   Parallel variable I/O"
endif
	@echo "    Features:  SD_FASTMODE, LTO, Apple M1 tuning"
ifeq ($(MAC_ARM_HAS_OMP),no)
	@echo "    Warning:   Install libomp for OpenMP: brew install libomp"
endif
	@$(CC_MAC) $(CFLAGS_MAC_ARM) -o $(PLUGIN_MAC_ARM) $(SOURCES) $(LDFLAGS_MAC_ARM)
	@echo "    Output:    $(PLUGIN_MAC_ARM)"
	@printf "    Size:      " && ls -lh $(PLUGIN_MAC_ARM) | awk '{print $$5}'
	@echo "----------------------------------------------------------------------"

# ------------------------------------------------------------------------------
# macOS Intel (x86_64)
# ------------------------------------------------------------------------------
macos-intel: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  [2/4] macOS Intel (x86_64)"
	@echo "----------------------------------------------------------------------"
	@echo "    Compiler:  $(CC_MAC) -arch x86_64"
ifeq ($(MAC_X86_HAS_OMP),yes)
ifeq ($(LIBOMP_STATIC_X86),yes)
	@echo "    OpenMP:    Enabled (static)"
else
	@echo "    OpenMP:    Enabled (dynamic)"
endif
	@echo "    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "    OpenMP:    Disabled (pthread only)"
	@echo "    Threads:   Parallel variable I/O"
endif
	@echo "    Features:  SD_FASTMODE, LTO, Haswell tuning"
	@$(CC_MAC) $(CFLAGS_MAC_X86) -o $(PLUGIN_MAC_X86) $(SOURCES) $(LDFLAGS_MAC_X86)
	@echo "    Output:    $(PLUGIN_MAC_X86)"
	@printf "    Size:      " && ls -lh $(PLUGIN_MAC_X86) | awk '{print $$5}'
	@echo "----------------------------------------------------------------------"

# ------------------------------------------------------------------------------
# Native build (for development)
# ------------------------------------------------------------------------------
native: $(BUILD_DIR)
ifeq ($(DETECTED_OS),macOS)
ifeq ($(DETECTED_ARCH),arm64)
	@echo "  Detected Apple Silicon - building arm64..."
	@$(MAKE) --no-print-directory macos-arm
	@cp $(PLUGIN_MAC_ARM) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
else
	@echo "  Detected Intel Mac - building x86_64..."
	@$(MAKE) --no-print-directory macos-intel
	@cp $(PLUGIN_MAC_X86) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
endif
else ifeq ($(DETECTED_OS),Linux)
	@echo "  Detected Linux - building native..."
	@$(MAKE) --no-print-directory linux
	@cp $(PLUGIN_LINUX) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
else ifeq ($(DETECTED_OS),Windows)
	@echo "  Detected Windows - building native..."
	@$(MAKE) --no-print-directory windows
	@cp $(PLUGIN_WINDOWS) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
endif

# ------------------------------------------------------------------------------
# Windows x64
# ------------------------------------------------------------------------------
windows: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  [3/4] Windows x64"
	@echo "----------------------------------------------------------------------"
	@echo "    Compiler:  $(CC_WIN)"
ifeq ($(WIN_HAS_OMP),yes)
	@echo "    OpenMP:    Enabled"
	@echo "    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "    OpenMP:    Disabled"
	@echo "    Threads:   Parallel variable I/O (pthread)"
endif
	@echo "    Features:  SD_FASTMODE, static linking"
ifeq ($(DETECTED_OS),macOS)
ifeq ($(LLVM_MINGW_EXISTS),no)
ifeq ($(MINGW_EXISTS),no)
	@echo "    ERROR: No Windows cross-compiler found"
	@echo "      Install mingw-w64:  brew install mingw-w64"
	@echo "      Or llvm-mingw for OpenMP support"
	@exit 1
endif
endif
endif
	@$(CC_WIN) $(CFLAGS_WIN) -o $(PLUGIN_WINDOWS) $(SOURCES) $(LDFLAGS_WIN) 2>&1 || \
		(echo "    Build failed. Check compiler installation." && exit 1)
	@echo "    Output:    $(PLUGIN_WINDOWS)"
	@printf "    Size:      " && ls -lh $(PLUGIN_WINDOWS) | awk '{print $$5}'
	@echo "----------------------------------------------------------------------"

# ------------------------------------------------------------------------------
# Linux x64
# ------------------------------------------------------------------------------
linux: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  [4/4] Linux x64"
	@echo "----------------------------------------------------------------------"
	@echo "    Compiler:  $(CC_LINUX)"
ifeq ($(DETECTED_OS),Linux)
ifeq ($(LINUX_HAS_OMP),yes)
	@echo "    OpenMP:    Enabled"
	@echo "    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "    OpenMP:    Disabled"
	@echo "    Threads:   Parallel variable I/O (pthread)"
endif
else
	@echo "    OpenMP:    Disabled (cross-compile)"
	@echo "    Note:      Cross-compile from $(DETECTED_OS) - build on Linux for best results"
endif
	@echo "    Features:  SD_FASTMODE, LTO"
ifeq ($(DETECTED_OS),Linux)
	@$(CC_LINUX) $(CFLAGS_LINUX) -o $(PLUGIN_LINUX) $(SOURCES) $(LDFLAGS_LINUX)
else
	@$(CC_LINUX) $(CFLAGS_LINUX) -o $(PLUGIN_LINUX) $(SOURCES) $(LDFLAGS_LINUX) 2>&1 || \
		(echo "    Cross-compile failed. Build on Linux instead." && exit 1)
endif
	@echo "    Output:    $(PLUGIN_LINUX)"
	@printf "    Size:      " && ls -lh $(PLUGIN_LINUX) | awk '{print $$5}'
	@echo "----------------------------------------------------------------------"

# ------------------------------------------------------------------------------
# Debug build with AddressSanitizer
# ------------------------------------------------------------------------------
ASAN_FLAGS = -fsanitize=address -fno-omit-frame-pointer -g -O1
TSAN_FLAGS = -fsanitize=thread -fno-omit-frame-pointer -g -O1

debug-asan: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  Debug build with AddressSanitizer"
	@echo "----------------------------------------------------------------------"
	$(CC_MAC) $(MAC_BASE_FLAGS) $(ASAN_FLAGS) -arch arm64 \
		-Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include \
		-o $(PLUGIN_MAC_ARM) $(SOURCES) \
		-bundle -arch arm64 -L$(LIBOMP_PREFIX)/lib -lomp -framework Accelerate
	@echo "    Output: $(PLUGIN_MAC_ARM) (with ASan)"

debug-tsan: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "----------------------------------------------------------------------"
	@echo "  Debug build with ThreadSanitizer"
	@echo "----------------------------------------------------------------------"
	$(CC_MAC) $(MAC_BASE_FLAGS) $(TSAN_FLAGS) -arch arm64 \
		-Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include \
		-o $(PLUGIN_MAC_ARM) $(SOURCES) \
		-bundle -arch arm64 -L$(LIBOMP_PREFIX)/lib -lomp -framework Accelerate
	@echo "    Output: $(PLUGIN_MAC_ARM) (with TSan)"

# ------------------------------------------------------------------------------
# Clean
# ------------------------------------------------------------------------------
clean:
	@echo ""
	@echo "  Cleaning build artifacts..."
ifeq ($(DETECTED_OS),Windows)
	@-del /Q $(BUILD_DIR)\ctools_*.plugin 2>NUL
	@-del /Q $(BUILD_DIR)\*.o 2>NUL
	@-del /Q $(SRC_DIR)\*.o 2>NUL
else
	@-rm -f $(BUILD_DIR)/ctools_*.plugin
	@-rm -f $(BUILD_DIR)/*.o
	@-rm -f $(SRC_DIR)/*.o
endif
	@echo "  Done."
	@echo ""

# ------------------------------------------------------------------------------
# Rebuild
# ------------------------------------------------------------------------------
rebuild: clean all
