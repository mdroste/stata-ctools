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
#   make                - Build for current platform
#   make all            - Build all distribution plugins
#   make native         - Build optimized plugin for current architecture
#   make macos-arm      - Build macOS Apple Silicon plugin
#   make macos-intel    - Build macOS Intel plugin
#   make windows        - Build Windows x64 plugin
#   make linux          - Build Linux x64 plugin
#   make windows-docker - Build Windows plugin via Docker
#   make linux-docker   - Build Linux plugin via Docker
#   make check          - Check build dependencies
#   make clean          - Remove compiled files
#   make help           - Show this help
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
# Source Files (auto-discovered from src/ directory tree)
# ------------------------------------------------------------------------------
# All .c files under src/ are compiled, all .h files tracked for dependencies,
# and all subdirectories are added as include paths. Adding new files or modules
# under src/ requires no changes to this Makefile.
#
# Third-party sources (libdeflate) are compiled separately with warnings
# suppressed to keep our own code warning-clean.
ALL_SOURCES  = $(shell find $(SRC_DIR) -name '*.c')
LIBDEFLATE_SRCS = $(shell find $(SRC_DIR)/cimport/libdeflate -name '*.c' 2>/dev/null)
SOURCES      = $(filter-out $(LIBDEFLATE_SRCS),$(ALL_SOURCES))
HEADERS      = $(shell find $(SRC_DIR) -name '*.h')
INCLUDE_DIRS = $(addprefix -I,$(sort $(dir $(ALL_SOURCES) $(HEADERS))))

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
        LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto -Wl,-S $(LIBOMP_PREFIX)/lib/libomp.a \
                          -framework Accelerate
    else
        LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto -Wl,-S -L$(LIBOMP_PREFIX)/lib -lomp \
                          -framework Accelerate
    endif
    MAC_ARM_HAS_OMP = yes
else
    CFLAGS_MAC_ARM = $(MAC_BASE_FLAGS) -arch arm64 \
                     -mmacosx-version-min=11.0 \
                     -mcpu=apple-m1
    LDFLAGS_MAC_ARM = -bundle -arch arm64 -flto -Wl,-S -framework Accelerate
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
                     -march=x86-64 -mtune=haswell -mavx2 \
                     -Xpreprocessor -fopenmp -I$(LIBOMP_INTEL)/include
    ifeq ($(LIBOMP_STATIC_X86),yes)
        # Static linking - bundle libomp into the plugin
        LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto -Wl,-S $(LIBOMP_INTEL)/lib/libomp.a \
                          -framework Accelerate
    else
        LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto -Wl,-S -L$(LIBOMP_INTEL)/lib -lomp \
                          -framework Accelerate
    endif
    MAC_X86_HAS_OMP = yes
else
    CFLAGS_MAC_X86 = $(MAC_BASE_FLAGS) -arch x86_64 \
                     -mmacosx-version-min=10.13 \
                     -march=x86-64 -mtune=haswell -mavx2
    LDFLAGS_MAC_X86 = -bundle -arch x86_64 -flto -Wl,-S -framework Accelerate
    MAC_X86_HAS_OMP = no
endif

# ==============================================================================
#                             Windows Configuration
# ==============================================================================
# Cross-compile: try llvm-mingw first (has OpenMP), fall back to mingw-w64
LLVM_MINGW_PREFIX ?= /opt/llvm-mingw

ifeq ($(DETECTED_OS),Windows)
    # Native Windows build (requires MSYS2/MinGW) - static link libgomp
    CC_WIN = gcc
    CFLAGS_WIN = -O3 -Wall -shared -DSYSTEM=STWIN32 -DSD_FASTMODE -fopenmp \
                 -ffast-math -funroll-loops -ftree-vectorize -flto \
                 -fdata-sections -ffunction-sections \
                 -fno-strict-aliasing -mavx2 $(INCLUDE_DIRS)
    LDFLAGS_WIN = -flto -Wl,--gc-sections -Wl,-S \
                  -static-libgcc -Wl,-Bstatic -lgomp -Wl,-Bdynamic -lpthread
    WIN_HAS_OMP = yes
    WIN_OMP_STATIC = yes
else
    LLVM_MINGW_EXISTS := $(shell test -x $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang && echo yes || echo no)
    MINGW_EXISTS := $(shell which x86_64-w64-mingw32-gcc 2>/dev/null && echo yes || echo no)

    ifeq ($(LLVM_MINGW_EXISTS),yes)
        # llvm-mingw: OpenMP support via bundled libomp.dll
        CC_WIN = $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang
        CFLAGS_WIN = -O3 -Wall -shared -DSYSTEM=STWIN32 -DSD_FASTMODE -fopenmp \
                     -ffast-math -funroll-loops -ftree-vectorize -flto \
                     -fdata-sections -ffunction-sections \
                     -fno-strict-aliasing -mavx2 $(INCLUDE_DIRS)
        LDFLAGS_WIN = -flto -Wl,--gc-sections -Wl,-S -fopenmp -lpthread
        WIN_HAS_OMP = yes
        WIN_OMP_STATIC = no
    else ifeq ($(MINGW_EXISTS),yes)
        # mingw-w64: no OpenMP
        CC_WIN = x86_64-w64-mingw32-gcc
        CFLAGS_WIN = -O3 -Wall -Wno-unknown-pragmas -shared -DSYSTEM=STWIN32 -DSD_FASTMODE \
                     -ffast-math -funroll-loops -ftree-vectorize -flto \
                     -fdata-sections -ffunction-sections \
                     -fno-strict-aliasing -mavx2 $(INCLUDE_DIRS)
        LDFLAGS_WIN = -flto -Wl,--gc-sections -Wl,-S -static-libgcc -lpthread
        WIN_HAS_OMP = no
        WIN_OMP_STATIC = no
    else
        # No cross-compiler available
        CC_WIN =
        CFLAGS_WIN =
        LDFLAGS_WIN =
        WIN_HAS_OMP = no
        WIN_OMP_STATIC = no
    endif
endif

# ==============================================================================
#                              Linux Configuration
# ==============================================================================
ifeq ($(DETECTED_OS),Linux)
    CC_LINUX = gcc
else
    # Check for a proper cross-compiler
    LINUX_CROSS := $(shell which x86_64-linux-gnu-gcc 2>/dev/null && echo yes || echo no)
    ifeq ($(LINUX_CROSS),yes)
        CC_LINUX = x86_64-linux-gnu-gcc
    else
        CC_LINUX =
    endif
endif

LINUX_BASE_FLAGS = -O3 -Wall -Wextra -shared -fPIC -DSYSTEM=STUNIX -DSD_FASTMODE \
                   -ffast-math -funroll-loops -ftree-vectorize -flto \
                   -fdata-sections -ffunction-sections \
                   -fno-strict-aliasing $(INCLUDE_DIRS)

# Add -mavx2 for x86_64 targets (native Linux or cross-compile)
ifeq ($(DETECTED_OS),Linux)
    LINUX_BASE_FLAGS += -mavx2
else ifneq ($(CC_LINUX),)
    LINUX_BASE_FLAGS += -mavx2
endif

# Check for OpenMP support (test the actual compiler, not just host gcc)
ifeq ($(DETECTED_OS),Linux)
    LINUX_HAS_OMP := $(shell gcc -fopenmp -E - < /dev/null >/dev/null 2>&1 && echo yes || echo no)
else ifneq ($(CC_LINUX),)
    LINUX_HAS_OMP := $(shell $(CC_LINUX) -fopenmp -E - < /dev/null >/dev/null 2>&1 && echo yes || echo no)
else
    LINUX_HAS_OMP = no
endif

# Note: Static linking of libgomp.a doesn't work on Linux for shared libraries
# because libgomp uses TLS relocations incompatible with PIC code.
# Linux plugins require libgomp.so at runtime.

# Check for OpenBLAS (native Linux only)
ifeq ($(DETECTED_OS),Linux)
    LINUX_HAS_OPENBLAS := $(shell pkg-config --exists openblas 2>/dev/null && echo yes || \
                                  (test -f /usr/include/cblas.h && test -f /usr/lib/libopenblas.so && echo yes) || \
                                  echo no)
else
    LINUX_HAS_OPENBLAS = no
endif

ifeq ($(LINUX_HAS_OMP),yes)
    ifeq ($(LINUX_HAS_OPENBLAS),yes)
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -fopenmp -DHAVE_OPENBLAS
        LDFLAGS_LINUX = -flto -Wl,--gc-sections -Wl,-S -fopenmp -lopenblas
    else
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -fopenmp
        LDFLAGS_LINUX = -flto -Wl,--gc-sections -Wl,-S -fopenmp
    endif
else
    ifeq ($(LINUX_HAS_OPENBLAS),yes)
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS) -DHAVE_OPENBLAS
        LDFLAGS_LINUX = -flto -Wl,--gc-sections -Wl,-S -lopenblas
    else
        CFLAGS_LINUX = $(LINUX_BASE_FLAGS)
        LDFLAGS_LINUX = -flto -Wl,--gc-sections -Wl,-S
    endif
endif

# ==============================================================================
#                                  Targets
# ==============================================================================
.PHONY: all native macos-arm macos-intel windows linux clean help check rebuild linux-docker windows-docker

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
	@echo "ctools plugin makefile"
	@echo ""
	@echo "  Builds the C plugin backend for ctools"
	@echo ""
	@echo "  Usage:"
	@echo "    make              Build for current platform ($(DETECTED_OS) $(DETECTED_ARCH))"
	@echo "    make all          Build all distribution plugins"
	@echo "    make native       Build optimized plugin for current architecture"
	@echo "    make macos-arm       Build macOS Apple Silicon plugin"
	@echo "    make macos-intel     Build macOS Intel plugin"
	@echo "    make windows         Build Windows x64 plugin"
	@echo "    make linux           Build Linux x64 plugin"
	@echo "    make windows-docker  Build Windows plugin via Docker"
	@echo "    make linux-docker    Build Linux plugin via Docker"
	@echo "    make check           Check build dependencies"
	@echo "    make clean           Remove compiled files"
	@echo "    make rebuild         Clean and rebuild all"
	@echo ""
	@echo "  Detected: $(DETECTED_OS) $(DETECTED_ARCH)"
	@echo ""

# ------------------------------------------------------------------------------
# Dependency Check
# ------------------------------------------------------------------------------
check:
	@echo ""
	@echo "ctools build environment check"
	@echo ""
	@echo "  Detected OS:   $(DETECTED_OS)"
	@echo "  Architecture:  $(DETECTED_ARCH)"
	@echo ""
	@echo "  macOS Apple Silicon (arm64)"
	@printf "    clang:     " && (which clang >/dev/null 2>&1 && printf "OK\n" || printf "MISSING - install Xcode CLT\n")
ifeq ($(LIBOMP_EXISTS),yes)
	@echo "    libomp:    OK ($(LIBOMP_PREFIX))"
	@echo "    OpenMP:    Enabled"
else
	@echo "    libomp:    MISSING - run: brew install libomp"
	@echo "    OpenMP:    Disabled"
endif
	@echo ""
	@echo "  macOS Intel (x86_64)"
	@printf "    clang:     " && (which clang >/dev/null 2>&1 && printf "OK\n" || printf "MISSING\n")
ifeq ($(MAC_X86_HAS_OMP),yes)
	@echo "    libomp:    OK (/usr/local/opt/libomp)"
	@echo "    OpenMP:    Enabled"
else
	@echo "    libomp:    MISSING - run: arch -x86_64 /usr/local/bin/brew install libomp"
	@echo "    OpenMP:    Disabled"
endif
	@echo ""
	@echo "  Windows x64"
ifeq ($(DETECTED_OS),Windows)
	@printf "    gcc:        " && (which gcc >/dev/null 2>&1 && printf "OK (MSYS2/MinGW)\n" || printf "MISSING - install MSYS2\n")
	@echo "    OpenMP:     Enabled"
else
	@printf "    llvm-mingw: " && (test -x $(LLVM_MINGW_PREFIX)/bin/x86_64-w64-mingw32-clang && printf "OK ($(LLVM_MINGW_PREFIX)) - OpenMP enabled\n" || printf "MISSING (optional, enables OpenMP)\n")
	@printf "    mingw-w64:  " && (which x86_64-w64-mingw32-gcc >/dev/null 2>&1 && printf "OK (no OpenMP)\n" || printf "MISSING\n")
ifneq ($(CC_WIN),)
	@echo "    Status:     Ready ($(CC_WIN))"
else
	@echo "    Status:     No cross-compiler found"
	@echo "                  Use: make windows-docker, brew install mingw-w64,"
	@echo "                  or install llvm-mingw from https://github.com/mstorsjo/llvm-mingw/releases"
endif
endif
	@echo ""
	@echo "  Linux x64"
ifeq ($(DETECTED_OS),Linux)
	@printf "    gcc:       " && (which gcc >/dev/null 2>&1 && printf "OK\n" || printf "MISSING\n")
ifeq ($(LINUX_HAS_OMP),yes)
	@echo "    OpenMP:    Enabled"
else
	@echo "    OpenMP:    MISSING - install libgomp-dev"
endif
else ifneq ($(CC_LINUX),)
	@echo "    Cross:     $(CC_LINUX)"
else
	@echo "    Cross:     No cross-compiler found (no x86_64-linux-gnu-gcc)"
	@echo "               Use: make linux-docker"
endif
	@echo ""
	@echo "  Docker (cross-compilation)"
	@printf "    docker:    " && (which docker >/dev/null 2>&1 && printf "OK\n" || printf "MISSING - https://docs.docker.com/get-docker/\n")
	@echo "    Targets:   make linux-docker, make windows-docker"
	@echo ""

# ------------------------------------------------------------------------------
# Build All Platforms
# ------------------------------------------------------------------------------
all: $(BUILD_DIR)
	@echo ""
	@echo "Building ctools plugins..."
	@$(MAKE) --no-print-directory macos-arm
	@$(MAKE) --no-print-directory macos-intel
	@$(MAKE) --no-print-directory windows
	@$(MAKE) --no-print-directory linux
	@echo ""
	@echo "Build complete."
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
	@echo " Building MacOS plugin (arm64)"
	@echo "    Compiler:  $(CC_MAC) -arch arm64"
ifeq ($(MAC_ARM_HAS_OMP),yes)
ifeq ($(LIBOMP_STATIC_ARM),yes)
	@echo "    OpenMP:    Enabled (static)"
else
	@echo "    OpenMP:    Enabled (dynamic)"
endif
else
	@echo "    OpenMP:    Disabled (pthread only)"
endif
ifeq ($(MAC_ARM_HAS_OMP),no)
	@echo "    Warning:   Install libomp for OpenMP: brew install libomp"
endif
	@for f in $(LIBDEFLATE_SRCS); do \
		d=$$(basename $$(dirname $$f)); b=$$(basename $$f .c); \
		$(CC_MAC) $(CFLAGS_MAC_ARM) -w -c $$f -o $(BUILD_DIR)/ld_$${d}_$${b}_arm.o; \
	done
	@$(CC_MAC) $(CFLAGS_MAC_ARM) -o $(PLUGIN_MAC_ARM) $(SOURCES) $(BUILD_DIR)/ld_*_arm.o $(LDFLAGS_MAC_ARM)
	@rm -f $(BUILD_DIR)/ld_*_arm.o
	@echo "    Output:    $(PLUGIN_MAC_ARM)"
	@printf "    Size:      " && ls -lh $(PLUGIN_MAC_ARM) | awk '{print $$5}'

# ------------------------------------------------------------------------------
# macOS Intel (x86_64)
# ------------------------------------------------------------------------------
macos-intel: $(BUILD_DIR) $(SOURCES) $(HEADERS)
	@echo ""
	@echo " Building MacOS plugin (x86_64)"
	@echo "    Compiler:  $(CC_MAC) -arch x86_64"
ifeq ($(MAC_X86_HAS_OMP),yes)
ifeq ($(LIBOMP_STATIC_X86),yes)
	@echo "    OpenMP:    Enabled (static)"
else
	@echo "    OpenMP:    Enabled (dynamic)"
endif
else
	@echo "    OpenMP:    Disabled (pthread only)"
endif
	@for f in $(LIBDEFLATE_SRCS); do \
		d=$$(basename $$(dirname $$f)); b=$$(basename $$f .c); \
		$(CC_MAC) $(CFLAGS_MAC_X86) -w -c $$f -o $(BUILD_DIR)/ld_$${d}_$${b}_x86.o; \
	done
	@$(CC_MAC) $(CFLAGS_MAC_X86) -o $(PLUGIN_MAC_X86) $(SOURCES) $(BUILD_DIR)/ld_*_x86.o $(LDFLAGS_MAC_X86)
	@rm -f $(BUILD_DIR)/ld_*_x86.o
	@echo "    Output:    $(PLUGIN_MAC_X86)"
	@printf "    Size:      " && ls -lh $(PLUGIN_MAC_X86) | awk '{print $$5}'

# ------------------------------------------------------------------------------
# Native build (for development)
# ------------------------------------------------------------------------------
native: $(BUILD_DIR)
ifeq ($(DETECTED_OS),macOS)
ifeq ($(DETECTED_ARCH),arm64)
	@echo "  Detected Apple Silicon - building macOS arm64 plugin..."
	@$(MAKE) --no-print-directory macos-arm
	@cp $(PLUGIN_MAC_ARM) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
else
	@echo "  Detected Intel Mac - building macOS x86_64 plugin..."
	@$(MAKE) --no-print-directory macos-intel
	@cp $(PLUGIN_MAC_X86) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
endif
else ifeq ($(DETECTED_OS),Linux)
	@echo "  Detected Linux - building Linux plugin..."
	@$(MAKE) --no-print-directory linux
	@cp $(PLUGIN_LINUX) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
else ifeq ($(DETECTED_OS),Windows)
	@echo "  Detected Windows - building Windows plugin..."
	@$(MAKE) --no-print-directory windows
	@cp $(PLUGIN_WINDOWS) $(PLUGIN_NATIVE)
	@echo ""
	@echo "  -> Copied to $(PLUGIN_NATIVE) for local development"
endif

# ------------------------------------------------------------------------------
# Windows x64
# ------------------------------------------------------------------------------
windows: $(BUILD_DIR) $(SOURCES) $(HEADERS)
ifeq ($(CC_WIN),)
	@echo ""
	@echo "  ERROR: No Windows cross-compiler found."
	@echo ""
	@echo "  Options:"
	@echo "    make windows-docker              Build via Docker (recommended)"
	@echo "    brew install mingw-w64           Install cross-compiler (macOS)"
	@echo "    apt install mingw-w64            Install cross-compiler (Linux)"
	@echo ""
	@exit 1
else
	@echo ""
	@echo " Building Windows plugin (x86_64)"
	@echo "    Compiler:  $(CC_WIN)"
ifeq ($(WIN_HAS_OMP),yes)
	@echo "    OpenMP:    Enabled"
else
	@echo "    OpenMP:    Disabled"
endif
	@for f in $(LIBDEFLATE_SRCS); do \
		d=$$(basename $$(dirname $$f)); b=$$(basename $$f .c); \
		$(CC_WIN) $(CFLAGS_WIN) -w -c $$f -o $(BUILD_DIR)/ld_$${d}_$${b}_win.o; \
	done
	@$(CC_WIN) $(CFLAGS_WIN) -o $(PLUGIN_WINDOWS) $(SOURCES) $(BUILD_DIR)/ld_*_win.o $(LDFLAGS_WIN)
	@rm -f $(BUILD_DIR)/ld_*_win.o
	@echo "    Output:    $(PLUGIN_WINDOWS)"
	@printf "    Size:      " && ls -lh $(PLUGIN_WINDOWS) | awk '{print $$5}'
endif

# ------------------------------------------------------------------------------
# Linux x64
# ------------------------------------------------------------------------------
linux: $(BUILD_DIR) $(SOURCES) $(HEADERS)
ifeq ($(CC_LINUX),)
	@echo ""
	@echo "  ERROR: No Linux cross-compiler found."
	@echo ""
	@echo "  Options:"
	@echo "    make linux-docker                Build via Docker (recommended)"
	@echo "                                     Or build natively on a Linux machine"
	@echo ""
	@exit 1
else
	@echo ""
	@echo " Building Linux plugin (x86_64)"
ifneq ($(DETECTED_OS),Linux)
	@echo "    Compiler:  $(CC_LINUX) (cross-compile)"
else
	@echo "    Compiler:  $(CC_LINUX)"
endif
ifeq ($(LINUX_HAS_OMP),yes)
	@echo "    OpenMP:    Enabled"
else
	@echo "    OpenMP:    Disabled"
endif
	@for f in $(LIBDEFLATE_SRCS); do \
		d=$$(basename $$(dirname $$f)); b=$$(basename $$f .c); \
		$(CC_LINUX) $(CFLAGS_LINUX) -w -c $$f -o $(BUILD_DIR)/ld_$${d}_$${b}_linux.o; \
	done
	@$(CC_LINUX) $(CFLAGS_LINUX) -o $(PLUGIN_LINUX) $(SOURCES) $(BUILD_DIR)/ld_*_linux.o $(LDFLAGS_LINUX)
	@rm -f $(BUILD_DIR)/ld_*_linux.o
	@echo "    Output:    $(PLUGIN_LINUX)"
	@printf "    Size:      " && ls -lh $(PLUGIN_LINUX) | awk '{print $$5}'
endif

# ------------------------------------------------------------------------------
# Docker Cross-Compilation
# ------------------------------------------------------------------------------
DOCKER = docker
DOCKER_LINUX_IMAGE = ctools-build-linux
DOCKER_WINDOWS_IMAGE = ctools-build-windows

linux-docker: $(BUILD_DIR)
	@command -v $(DOCKER) >/dev/null 2>&1 || \
		{ echo ""; echo "  ERROR: Docker is not installed."; echo "    Install: https://docs.docker.com/get-docker/"; echo ""; exit 1; }
	@echo ""
	@echo " Building Linux plugin via Docker"
	@echo "    Preparing container..."
	@printf 'FROM ubuntu:22.04\nRUN apt-get update -qq && apt-get install -y -qq gcc make libgomp-dev libopenblas-dev pkg-config > /dev/null 2>&1\n' | \
		$(DOCKER) build --platform linux/amd64 -q -t $(DOCKER_LINUX_IMAGE) - > /dev/null
	@echo "    Compiling (x86_64)..."
	@$(DOCKER) run --rm --platform linux/amd64 -v "$$(pwd)":/src -w /src $(DOCKER_LINUX_IMAGE) make linux
	@echo ""

windows-docker: $(BUILD_DIR)
	@command -v $(DOCKER) >/dev/null 2>&1 || \
		{ echo ""; echo "  ERROR: Docker is not installed."; echo "    Install: https://docs.docker.com/get-docker/"; echo ""; exit 1; }
	@echo ""
	@echo " Building Windows plugin via Docker"
	@echo "    Preparing container..."
	@printf 'FROM ubuntu:22.04\nRUN apt-get update -qq && apt-get install -y -qq gcc-mingw-w64-x86-64 make > /dev/null 2>&1\n' | \
		$(DOCKER) build -q -t $(DOCKER_WINDOWS_IMAGE) - > /dev/null
	@echo "    Compiling (cross-compile)..."
	@$(DOCKER) run --rm -v "$$(pwd)":/src -w /src $(DOCKER_WINDOWS_IMAGE) make windows
	@echo ""

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
