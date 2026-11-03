# Installation and building

## Compatibility and Requirements

ctools is compatible with Stata 14.1+ (plugin interface version 3.0). The [platform requirements](development/PLATFORMS.md) list supported processors, operating systems, and required libraries: macOS 11 or later on Apple Silicon or Intel; Windows 10 or later on x86-64; x86-64 Linux with glibc 2.17 or later (RHEL/CentOS 7 and newer); the OpenMP runtime is built into the plugin.

**Note:** ctools does not support datasets exceeding 2^31−1 (2,147,483,647) observations. This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of Stata's plugin API for C and can only be addressed with an internal Stata update. ctools will gracefully exit with an error if your dataset exceeds this limit.


## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

### Updating an installation

Restart Stata after replacing a loaded plugin. The interpolation command is
`cinterpolate`; `cinterpolate, cubic` implements SSC `cipolate`-style interpolation.
When upgrading from a ctools installation that supplied `cipolate.ado` and
`cipolate.sthlp`, locate those files with `which cipolate`, check their headers,
and remove only the stale ctools copies. See the [migration instructions](USAGE.md#interpolation-cinterpolate).

### Manual Installation

1. Download `ctools-<version>.zip` from [GitHub Releases](https://github.com/mdroste/stata-ctools/releases) and extract it.
2. In Stata, run `sysdir` to find your PERSONAL directory, or `adopath` to choose another directory already on the search path.
3. Copy the extracted files directly into that directory. The ado files, help files, and all four plugins belong together; there is no `build/` wrapper in the ZIP. Restart Stata if replacing an installed plugin.

Alternatively, run `net install ctools, from("/absolute/path/to/extracted/archive") replace`.
From a source checkout, verify the extracted package with
`python3 validation/run.py --checks package /path/to/extracted/archive --integrity`.

GitHub Releases are created explicitly and receive `ctools-<version>.zip` after
the build checks pass. The version comes from `config/release.json`. Pushes to
`main` run CI and produce artifacts. Updating the installable plugins on `main`
requires a manual **Build Plugins** run on `main` with **publish_plugins** selected.
Actions runs retain the ZIP in the `ctools-<version>` artifact, including
candidate builds from development branches and pull requests. Pull-request
Windows plugins are unsigned, and Mac plugins lack Developer ID signing and
notarization. Downloading that
artifact adds an outer Actions ZIP; extract it to obtain the install ZIP.

A source checkout's `build/` is not necessarily an installable distribution: it
may contain only the locally compiled binary. After building from source, use
`make package BUILD_REVISION=<revision> PACKAGE_DIR=dist/ctools-<revision>` to
produce a complete package for the host platform. Install from that staged
directory. Its manifest explicitly identifies the platform and `BUILD_INFO.json`
records the release version, date, source revision, and file checksums. Full release
packages contain all four platform plugins, including the signed Windows plugin
and Developer ID signed, notarized Mac plugins. Mac plugins use Apple's online
notarization tickets; first use may require internet access. See
[macOS signing](development/MACOS_SIGNING.md) for verification and setup details.


## Troubleshooting

### The plugin does not match the installed ado files

Every ctools command first checks that the installed plugin and ado files use
compatible interfaces. A plugin from another release stops the command, before any data or
file changes, with an error that names both and asks you to reinstall:

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

Then restart Stata: a plugin stays loaded until Stata exits, so a reinstall in a
running session does not replace it. After a manual installation, copy every file
from one release ZIP. Old plugins report no protocol.

### The plugin cannot be loaded

The error lists the plugin ctools selected, Stata's `c(os)` and `c(machine_type)`,
the path `findfile` found (or that the file is missing), and then the message from
Stata's plugin loader; the command exits with the loader's return code.

- **File not found:** reinstall ctools. A manual installation must include all four
  `ctools_*.plugin` files in the same directory as the ado files.
- **Linux:** the plugin needs x86-64 and glibc 2.17 or later (`ldd --version`); its
  OpenMP runtime is linked in, so no `libgomp` is required. A message such as
  ``version `GLIBC_2.xx' not found`` names the missing piece; `ldd -r` on the plugin
  path from the error lists every unresolved library and symbol.
- **macOS:** release plugins are signed with a Developer ID and notarized, and the
  first load may need internet access to check the notarization ticket. A plugin
  built locally, or one Gatekeeper has quarantined, can fail with a code-signature
  or "not valid mach-o file" message. `codesign --verify --verbose` on the plugin
  path checks the signature; reinstall from the release rather than overriding
  Gatekeeper for a file of unknown origin.
- **Windows:** the release plugin carries an Authenticode signature (see the file's
  Properties, Digital Signatures). SmartScreen does not screen plugins that Stata
  loads, but endpoint-security software can block a DLL; ask your administrator to
  allow the signed plugin.

## Building from Source Files

Release settings are defined in `config/release.json`. Run `python3 scripts/sync_release.py` after changing it; CI checks the generated copies. `ctools, version` reports version information for the installed ado files and plugin.


You probably do not need to compile the ctools plugin yourself; GitHub automatically builds plugins for Windows, Mac, and Linux. If you want or need to build plugins from source, follow the [platform build instructions](development/PLATFORMS.md), including the deployment-compatible macOS OpenMP runtime.

```bash
make              # Build for current platform
make all          # Build all platform plugins (requires every toolchain)
make package      # Build and stage an installable host-platform package
make check        # Check build dependencies
make clean        # Remove compiled files
```

If you have a workstation/server CPU with lots of cache, you might want to try playing around with the settings in [src/ctools_config.h](../src/ctools_config.h).

See [docs/development/README.md](development/README.md) for additional information on ctools' architecture and core logic.



## Release validation

Publication requires successful builds for Mac ARM, Mac Intel, Windows, and Linux,
the native CI checks, platform dependency checks, and a verified, timestamped
Windows signature, and verified Developer ID signatures and notarization for both
Mac plugins. See [Windows signing](development/WINDOWS_SIGNING.md) and
[macOS signing](development/MACOS_SIGNING.md) for credential and GitHub configuration.
GitHub Actions does not require self-hosted runners or
run licensed Stata tests. Stata correctness suites remain available for local testing.

Exports normally write a temporary file and rename it after writing succeeds. Exports containing several files retain the originals until all new files are in place and restore them if an error occurs. A system crash can leave only some files replaced. If restoration fails, ctools retains the backups and reports their locations. On filesystems that cannot rename or link files without replacing an existing destination, ctools creates a new destination and copies into it. That filename is visible during copying and is removed if copying fails.
