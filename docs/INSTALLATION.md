# Installation and building

## Compatibility and Requirements

ctools is compatible with Stata 14.1+ (plugin interface version 3.0). The [platform contract](development/PLATFORMS.md) specifies CPU/OS baselines and runtime dependencies, including Linux libgomp.

**Note:** ctools does not support datasets exceeding 2^31 (~2.147 billion) observations. This is a [known limitation](https://github.com/mcaceresb/stata-gtools/issues/43) of Stata's plugin API for C and can only be addressed with an internal Stata update. ctools will gracefully exit with an error if your dataset exceeds this limit.


## Installation

### From GitHub (recommended)

```stata
net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace
```

### Manual Installation

1. Download `ctools-<version>.zip` from [GitHub Releases](https://github.com/mdroste/stata-ctools/releases) and extract it.
2. In Stata, run `sysdir` to find your PERSONAL directory, or `adopath` to choose another directory already on the search path.
3. Copy the extracted files directly into that directory. The ado files, help files, and all four plugins belong together; there is no `build/` wrapper in the ZIP. Restart Stata if replacing an installed plugin.

Alternatively, run `net install ctools, from("/absolute/path/to/extracted/archive") replace`.
From a source checkout, verify the extracted package with
`python3 -m validation.checks.check_package /path/to/extracted/archive --integrity`.

Each main build that passes all publication gates creates a GitHub prerelease
tagged `v<version>-build.<run-number>`, containing `ctools-<version>.zip`.
The version comes from `release.json`; the build number distinguishes builds of
the same version. Manually created GitHub Releases also receive this ZIP.
Actions runs retain the ZIP in the `ctools-<version>` artifact, including
candidate builds that have not passed licensed validation. Downloading that
artifact adds an outer Actions ZIP; extract it to obtain the install ZIP.

A source checkout's `build/` is not necessarily an installable distribution: it
may contain only the locally compiled binary. After building from source, use
`make package BUILD_REVISION=<revision> PACKAGE_DIR=dist/ctools-<revision>` to
produce a complete package for the host platform. Install from that staged
directory. Its manifest explicitly identifies the platform and `BUILD_INFO.json`
records the revision and file checksums. Full release packages contain all four
platform plugins; publication requires the licensed correctness gate.


## Building from Source Files

Release identity is defined in `release.json`. Run `python3 scripts/sync_release.py` after changing it; CI checks the generated copies. `ctools, version` reports the installed ado and plugin identities.


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



## Licensed release validation

Publication requires a successful Stata correctness run on each fresh candidate
plugin: Mac ARM, Mac Intel, Windows, and Linux. Configure repository variable
`CTOOLS_STATA_RUNNERS` as a JSON object mapping `mac_arm`, `mac_x86`, `windows`,
and `linux` to trusted self-hosted runner label arrays. Each runner needs Python
3, the documented reference commands and `CTOOLS_STATA_EXECUTABLE` set to the
absolute licensed Stata executable path. The runner's OS/architecture and loaded
plugin revision are checked against the candidate. Missing configuration permits
candidate artifact builds but blocks public plugin commits and release assets.
Untrusted pull requests never run on these licensed runners. Windows uses Stata's
[noninteractive `/e` batch mode](https://www.stata.com/support/faqs/windows/batch-mode/).

Exports normally publish with atomic rename. Multi-file SHP and XPORT exports
retain originals until all outputs are published and roll back ordinary failures;
this is not a crash-atomic filesystem transaction. Backups are retained with a
reported path if rollback itself fails. On POSIX filesystems lacking both exclusive
rename and hard links, the no-overwrite fallback copies into an exclusively created
private destination and removes it on failure. Its name is visible while copying.
