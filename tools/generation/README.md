# Generated import tables

These scripts create the character encoding and locale tables used by `cimport`.
The generated files are stored in `src/cimport/generated/` and included in the
repository, so building ctools requires neither Java nor a network connection.
The handwritten import code is in `src/cimport/`.

From the repository root, check that the generated files are up to date:

```sh
python3 tools/generation/generate_cimport_charsets.py --check
python3 tools/generation/generate_cimport_locales.py --check
```

Omit `--check` to regenerate the files.

## Character encoding tables

The inputs are copies of ICU 67.1's `release-67-1` source files. Their SHA-256
checksums are recorded in `src/vendor/icu_charset/SHA256SUMS`; the Unicode
license is included with the copies. The script checks both input files before
reading their tables. It reads Java source as text.

## Locale tables

The locale data were extracted from OpenJDK 17.0.14 (2025-01-21 LTS).
`inputs/locales.json` records the Java version and checksums of the extraction
program and its UTF-8 output. To check the saved data against that Java version:

```sh
python3 tools/generation/generate_cimport_locales.py --java /path/to/jdk17/bin/java --check
```

This check requires the recorded Java version. Changes to the input data need a
separate compatibility review because they can change how ctools imports files.
Git preserves the input files' line endings so their checksums match on every
platform.
