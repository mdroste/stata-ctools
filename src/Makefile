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
# Core plugin infrastructure (in src/)
CORE_SRCS = stplugin.c \
            ctools_plugin.c

# Common utilities (in src/)
COMMON_SRCS = ctools_types.c \
              ctools_timer.c \
              ctools_error.c \
              ctools_threads.c \
              ctools_data_load.c \
              ctools_data_store.c \
              ctools_sort_radix_lsd.c

# csort module (in src/csort/)
CSORT_SRCS = csort/csort_impl.c

# cimport module (in src/cimport/)
CIMPORT_SRCS = cimport/cimport_impl.c

# cexport module (in src/cexport/)
CEXPORT_SRCS = cexport/cexport_impl.c

# cmerge module (in src/cmerge/)
CMERGE_SRCS = cmerge/cmerge_impl.c

# creghdfe module (in src/creghdfe/)
CREGHDFE_SRCS = creghdfe/creghdfe_types.c \
                creghdfe/creghdfe_utils.c \
                creghdfe/creghdfe_hdfe.c \
                creghdfe/creghdfe_solver.c \
                creghdfe/creghdfe_ols.c \
                creghdfe/creghdfe_vce.c \
                creghdfe/creghdfe_main.c \
                creghdfe/creghdfe_impl.c

# cqreg module (in src/cqreg/)
CQREG_SRCS = cqreg/cqreg_types.c \
             cqreg/cqreg_linalg.c \
             cqreg/cqreg_ipm.c \
             cqreg/cqreg_fn.c \
             cqreg/cqreg_irls.c \
             cqreg/cqreg_sparsity.c \
             cqreg/cqreg_vce.c \
             cqreg/cqreg_hdfe.c \
             cqreg/cqreg_main.c \
             cqreg/cqreg_impl.c

# All sources
SOURCES = $(CORE_SRCS) $(COMMON_SRCS) $(CSORT_SRCS) $(CIMPORT_SRCS) $(CEXPORT_SRCS) $(CMERGE_SRCS) $(CREGHDFE_SRCS) $(CQREG_SRCS)

# Headers (for dependency tracking)
CORE_HEADERS = stplugin.h \
               ctools_types.h \
               ctools_timer.h \
               ctools_error.h \
               ctools_threads.h \
               ctools_config.h

CSORT_HEADERS = csort/csort_impl.h

CIMPORT_HEADERS = cimport/cimport_impl.h

CEXPORT_HEADERS = cexport/cexport_impl.h

CMERGE_HEADERS = cmerge/cmerge_impl.h

CREGHDFE_HEADERS = creghdfe/creghdfe_types.h \
                   creghdfe/creghdfe_utils.h \
                   creghdfe/creghdfe_hdfe.h \
                   creghdfe/creghdfe_solver.h \
                   creghdfe/creghdfe_ols.h \
                   creghdfe/creghdfe_vce.h \
                   creghdfe/creghdfe_main.h \
                   creghdfe/creghdfe_impl.h

CQREG_HEADERS = cqreg/cqreg_types.h \
                cqreg/cqreg_linalg.h \
                cqreg/cqreg_ipm.h \
                cqreg/cqreg_sparsity.h \
                cqreg/cqreg_vce.h \
                cqreg/cqreg_hdfe.h \
                cqreg/cqreg_main.h \
                cqreg/cqreg_impl.h

HEADERS = $(CORE_HEADERS) $(CSORT_HEADERS) $(CIMPORT_HEADERS) $(CEXPORT_HEADERS) $(CMERGE_HEADERS) $(CREGHDFE_HEADERS) $(CQREG_HEADERS)

# Include paths for subdirectories
INCLUDE_DIRS = -I. -Icsort -Icimport -Icexport -Icmerge -Icreghdfe -Icqreg

# ------------------------------------------------------------------------------
# Output Configuration
# ------------------------------------------------------------------------------
BUILD_DIR = ../build

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

# macOS Apple Silicon (arm64) with OpenMP
ifeq ($(LIBOMP_EXISTS),yes)
    CFLAGS_MAC_ARM = $(MAC_BASE_FLAGS) -arch arm64 \
                     -mmacosx-version-min=11.0 \
                     -mcpu=apple-m1 \
                     -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto -L$(LIBOMP_PREFIX)/lib -lomp
    MAC_ARM_HAS_OMP = yes
else
    CFLAGS_MAC_ARM = $(MAC_BASE_FLAGS) -arch arm64 \
                     -mmacosx-version-min=11.0 \
                     -mcpu=apple-m1
    LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto
    MAC_ARM_HAS_OMP = no
endif

# Check for x86_64 libomp (installed via Rosetta Homebrew)
LIBOMP_INTEL_EXISTS := $(shell test -f /usr/local/opt/libomp/lib/libomp.dylib && \
                               file /usr/local/opt/libomp/lib/libomp.dylib 2>/dev/null | grep -q x86_64 && \
                               echo yes || echo no)

# macOS Intel (x86_64)
ifeq ($(LIBOMP_INTEL_EXISTS),yes)
    LIBOMP_INTEL = /usr/local/opt/libomp
    CFLAGS_MAC_X86 = $(MAC_BASE_FLAGS) -arch x86_64 \
                     -mmacosx-version-min=10.13 \
                     -march=x86-64 -mtune=haswell \
                     -Xpreprocessor -fopenmp -I$(LIBOMP_INTEL)/include
    LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto -L$(LIBOMP_INTEL)/lib -lomp
    MAC_X86_HAS_OMP = yes
else
    CFLAGS_MAC_X86 = $(MAC_BASE_FLAGS) -arch x86_64 \
                     -mmacosx-version-min=10.13 \
                     -march=x86-64 -mtune=haswell
    LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto
    MAC_X86_HAS_OMP = no
endif

# ==============================================================================
#                             Windows Configuration
# ==============================================================================
# Try llvm-mingw first (has OpenMP), fall back to mingw-w64
LLVM_MINGW_PREFIX = /opt/llvm-mingw
LLVM_MINGW_EXISTS := $(shell test -x $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang && echo yes || echo no)

ifeq ($(DETECTED_OS),Windows)
    # Native Windows build
    CC_WIN = gcc
    CFLAGS_WIN = -O3 -Wall -shared -DSYSTEM=STWIN32 -DSD_FASTMODE -fopenmp \
                 -ffast-math -funroll-loops -ftree-vectorize -fno-strict-aliasing $(INCLUDE_DIRS)
    LDFLAGS_WIN = -static -fopenmp
    WIN_HAS_OMP = yes
else ifeq ($(LLVM_MINGW_EXISTS),yes)
    # Cross-compile with llvm-mingw (has OpenMP support)
    CC_WIN = $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang
    CFLAGS_WIN = -O3 -Wall -shared -DSYSTEM=STWIN32 -DSD_FASTMODE -fopenmp \
                 -ffast-math -funroll-loops -ftree-vectorize -fno-strict-aliasing $(INCLUDE_DIRS)
    LDFLAGS_WIN = -fopenmp
    WIN_HAS_OMP = yes
else
    # Cross-compile with mingw-w64 (no OpenMP)
    CC_WIN = x86_64-w64-mingw32-gcc
    MINGW_EXISTS := $(shell which x86_64-w64-mingw32-gcc 2>/dev/null && echo yes || echo no)
    CFLAGS_WIN = -O3 -Wall -Wno-unknown-pragmas -shared -DSYSTEM=STWIN32 -DSD_FASTMODE \
                 -ffast-math -funroll-loops -ftree-vectorize -fno-strict-aliasing $(INCLUDE_DIRS)
    LDFLAGS_WIN = -static-libgcc
    WIN_HAS_OMP = no
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

ifeq ($(LINUX_HAS_OMP),yes)
    CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -fopenmp
    LDFLAGS_LINUX = -flto -fopenmp
else
    CFLAGS_LINUX = $(LINUX_BASE_FLAGS)
    LDFLAGS_LINUX = -flto
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
	@echo "╔══════════════════════════════════════════════════════════════════════════════╗"
	@echo "║                      ctools Stata Plugin Build System                        ║"
	@echo "╚══════════════════════════════════════════════════════════════════════════════╝"
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
	@echo "╔══════════════════════════════════════════════════════════════════════════════╗"
	@echo "║                        Build Environment Check                               ║"
	@echo "╚══════════════════════════════════════════════════════════════════════════════╝"
	@echo ""
	@echo "  Detected OS:   $(DETECTED_OS)"
	@echo "  Architecture:  $(DETECTED_ARCH)"
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  macOS Apple Silicon (arm64)                                                 │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@printf "│    clang:     " && (which clang >/dev/null 2>&1 && printf "✓ OK\n" || printf "✗ MISSING - install Xcode CLT\n")
ifeq ($(LIBOMP_EXISTS),yes)
	@echo "│    libomp:    ✓ OK ($(LIBOMP_PREFIX))"
	@echo "│    OpenMP:    ✓ Enabled"
else
	@echo "│    libomp:    ✗ MISSING - run: brew install libomp"
	@echo "│    OpenMP:    ✗ Disabled"
endif
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  macOS Intel (x86_64)                                                        │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@printf "│    clang:     " && (which clang >/dev/null 2>&1 && printf "✓ OK\n" || printf "✗ MISSING\n")
ifeq ($(MAC_X86_HAS_OMP),yes)
	@echo "│    libomp:    ✓ OK (/usr/local/opt/libomp)"
	@echo "│    OpenMP:    ✓ Enabled"
else
	@echo "│    libomp:    ✗ MISSING - run: arch -x86_64 /usr/local/bin/brew install libomp"
	@echo "│    OpenMP:    ✗ Disabled"
endif
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  Windows x64 (cross-compile)                                                 │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
ifeq ($(LLVM_MINGW_EXISTS),yes)
	@echo "│    llvm-mingw: ✓ OK ($(LLVM_MINGW_PREFIX))"
	@echo "│    OpenMP:     ✓ Enabled"
else
	@printf "│    mingw-w64:  " && (which x86_64-w64-mingw32-gcc >/dev/null 2>&1 && printf "✓ OK (no OpenMP)\n" || printf "✗ MISSING - brew install mingw-w64\n")
	@echo "│    llvm-mingw: ✗ MISSING (optional, enables OpenMP)"
	@echo "│                  Download: https://github.com/mstorsjo/llvm-mingw/releases"
	@echo "│                  Extract to: /opt/llvm-mingw"
endif
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  Linux x64                                                                   │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
ifeq ($(DETECTED_OS),Linux)
	@printf "│    gcc:       " && (which gcc >/dev/null 2>&1 && printf "✓ OK\n" || printf "✗ MISSING\n")
ifeq ($(LINUX_HAS_OMP),yes)
	@echo "│    OpenMP:    ✓ Enabled"
else
	@echo "│    OpenMP:    ✗ MISSING - install libgomp-dev"
endif
else
	@echo "│    (Cross-compile from $(DETECTED_OS) not fully supported)"
	@printf "│    gcc:       " && (which gcc >/dev/null 2>&1 && printf "✓ OK\n" || printf "✗ MISSING\n")
endif
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"
	@echo ""

# ------------------------------------------------------------------------------
# Build All Platforms
# ------------------------------------------------------------------------------
all: $(BUILD_DIR)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════════════════╗"
	@echo "║                      ctools Multi-Platform Build                             ║"
	@echo "╠══════════════════════════════════════════════════════════════════════════════╣"
	@echo "║  Building: csort · cmerge · cimport · cexport · creghdfe · cqreg             ║"
	@echo "╚══════════════════════════════════════════════════════════════════════════════╝"
	@$(MAKE) --no-print-directory macos-arm
	@$(MAKE) --no-print-directory macos-intel
	@$(MAKE) --no-print-directory windows
	@$(MAKE) --no-print-directory linux
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════════════════════╗"
	@echo "║                           BUILD COMPLETE                                     ║"
	@echo "╚══════════════════════════════════════════════════════════════════════════════╝"
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
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  [1/4] macOS Apple Silicon (arm64)                                           │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@echo "│    Compiler:  $(CC_MAC) -arch arm64"
ifeq ($(MAC_ARM_HAS_OMP),yes)
	@echo "│    OpenMP:    ✓ Enabled"
	@echo "│    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "│    OpenMP:    ✗ Disabled (pthread only)"
	@echo "│    Threads:   Parallel variable I/O"
endif
	@echo "│    Features:  SD_FASTMODE, LTO, Apple M1 tuning"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
ifeq ($(MAC_ARM_HAS_OMP),no)
	@echo "│    ⚠ Install libomp for OpenMP: brew install libomp"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
endif
	@$(CC_MAC) $(CFLAGS_MAC_ARM) -o $(PLUGIN_MAC_ARM) $(SOURCES) $(LDFLAGS_MAC_ARM)
	@echo "│    Output:    $(PLUGIN_MAC_ARM)"
	@printf "│    Size:      " && ls -lh $(PLUGIN_MAC_ARM) | awk '{print $$5}'
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"

# ------------------------------------------------------------------------------
# macOS Intel (x86_64)
# ------------------------------------------------------------------------------
macos-intel: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  [2/4] macOS Intel (x86_64)                                                  │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@echo "│    Compiler:  $(CC_MAC) -arch x86_64"
ifeq ($(MAC_X86_HAS_OMP),yes)
	@echo "│    OpenMP:    ✓ Enabled"
	@echo "│    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "│    OpenMP:    ✗ Disabled (pthread only)"
	@echo "│    Threads:   Parallel variable I/O"
endif
	@echo "│    Features:  SD_FASTMODE, LTO, Haswell tuning"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@$(CC_MAC) $(CFLAGS_MAC_X86) -o $(PLUGIN_MAC_X86) $(SOURCES) $(LDFLAGS_MAC_X86)
	@echo "│    Output:    $(PLUGIN_MAC_X86)"
	@printf "│    Size:      " && ls -lh $(PLUGIN_MAC_X86) | awk '{print $$5}'
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"

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
	@echo "  → Copied to $(PLUGIN_NATIVE) for local development"
else
	@echo "  Detected Intel Mac - building x86_64..."
	@$(MAKE) --no-print-directory macos-intel
	@cp $(PLUGIN_MAC_X86) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  → Copied to $(PLUGIN_NATIVE) for local development"
endif
else ifeq ($(DETECTED_OS),Linux)
	@echo "  Detected Linux - building native..."
	@$(MAKE) --no-print-directory linux
	@cp $(PLUGIN_LINUX) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  → Copied to $(PLUGIN_NATIVE) for local development"
else ifeq ($(DETECTED_OS),Windows)
	@echo "  Detected Windows - building native..."
	@$(MAKE) --no-print-directory windows
	@cp $(PLUGIN_WINDOWS) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  → Copied to $(PLUGIN_NATIVE) for local development"
endif

# ------------------------------------------------------------------------------
# Windows x64
# ------------------------------------------------------------------------------
windows: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  [3/4] Windows x64                                                           │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@echo "│    Compiler:  $(CC_WIN)"
ifeq ($(WIN_HAS_OMP),yes)
	@echo "│    OpenMP:    ✓ Enabled"
	@echo "│    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "│    OpenMP:    ✗ Disabled"
	@echo "│    Threads:   Parallel variable I/O (pthread)"
endif
	@echo "│    Features:  SD_FASTMODE, static linking"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
ifeq ($(DETECTED_OS),macOS)
ifeq ($(LLVM_MINGW_EXISTS),no)
ifeq ($(MINGW_EXISTS),no)
	@echo "│    ✗ ERROR: No Windows cross-compiler found"
	@echo "│      Install mingw-w64:  brew install mingw-w64"
	@echo "│      Or llvm-mingw for OpenMP support"
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"
	@exit 1
endif
endif
endif
	@$(CC_WIN) $(CFLAGS_WIN) -o $(PLUGIN_WINDOWS) $(SOURCES) $(LDFLAGS_WIN) 2>&1 || \
		(echo "│    ✗ Build failed. Check compiler installation." && \
		 echo "└──────────────────────────────────────────────────────────────────────────────┘" && \
		 exit 1)
	@echo "│    Output:    $(PLUGIN_WINDOWS)"
	@printf "│    Size:      " && ls -lh $(PLUGIN_WINDOWS) | awk '{print $$5}'
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"

# ------------------------------------------------------------------------------
# Linux x64
# ------------------------------------------------------------------------------
linux: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo "┌──────────────────────────────────────────────────────────────────────────────┐"
	@echo "│  [4/4] Linux x64                                                             │"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
	@echo "│    Compiler:  $(CC_LINUX)"
ifeq ($(DETECTED_OS),Linux)
ifeq ($(LINUX_HAS_OMP),yes)
	@echo "│    OpenMP:    ✓ Enabled"
	@echo "│    Threads:   Parallel variable I/O + parallel sorting"
else
	@echo "│    OpenMP:    ✗ Disabled"
	@echo "│    Threads:   Parallel variable I/O (pthread)"
endif
else
	@echo "│    OpenMP:    ✗ Disabled (cross-compile)"
	@echo "│    Note:      Cross-compile from $(DETECTED_OS) - build on Linux for best results"
endif
	@echo "│    Features:  SD_FASTMODE, LTO"
	@echo "├──────────────────────────────────────────────────────────────────────────────┤"
ifeq ($(DETECTED_OS),Linux)
	@$(CC_LINUX) $(CFLAGS_LINUX) -o $(PLUGIN_LINUX) $(SOURCES) $(LDFLAGS_LINUX)
else
	@$(CC_LINUX) $(CFLAGS_LINUX) -o $(PLUGIN_LINUX) $(SOURCES) $(LDFLAGS_LINUX) 2>&1 || \
		(echo "│    ⚠ Cross-compile failed. Build on Linux instead." && \
		 echo "└──────────────────────────────────────────────────────────────────────────────┘" && \
		 exit 1)
endif
	@echo "│    Output:    $(PLUGIN_LINUX)"
	@printf "│    Size:      " && ls -lh $(PLUGIN_LINUX) | awk '{print $$5}'
	@echo "└──────────────────────────────────────────────────────────────────────────────┘"

# ------------------------------------------------------------------------------
# Clean
# ------------------------------------------------------------------------------
clean:
	@echo ""
	@echo "  Cleaning build artifacts..."
ifeq ($(DETECTED_OS),Windows)
	@-del /Q $(BUILD_DIR)\ctools_*.plugin 2>NUL
	@-del /Q $(BUILD_DIR)\*.o 2>NUL
	@-del /Q *.o 2>NUL
else
	@-rm -f $(BUILD_DIR)/ctools_*.plugin
	@-rm -f $(BUILD_DIR)/*.o
	@-rm -f *.o
endif
	@echo "  Done."
	@echo ""

# ------------------------------------------------------------------------------
# Rebuild
# ------------------------------------------------------------------------------
rebuild: clean all
