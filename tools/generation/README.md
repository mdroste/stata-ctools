# Generated import tables

Ordinary builds consume checked-in files under `src/cimport/generated/`. They
need neither Java nor a network connection. Handwritten import implementation
fragments remain in `src/cimport/`.

Run the deterministic, read-only drift checks:

```sh
python3 tools/generation/generate_cimport_charsets.py --check
python3 tools/generation/generate_cimport_locales.py --check
```

Omit `--check` to regenerate outputs. ICU input sources are the exact 67.1
`release-67-1` files, pinned by `src/vendor/icu_charset/SHA256SUMS`, with the
Unicode license alongside the input copies. The generator verifies both source
checksums before parsing them; it does not execute Java.

Locale inputs are the extracted facts from OpenJDK 17.0.14 (2025-01-21 LTS).
`inputs/locales.json` records the exact runtime identification and SHA-256 of
both the Java extractor source and its UTF-8 facts. To independently verify
those facts with that runtime, pass `--java /path/to/jdk17/bin/java --check`.
A newer JDK is deliberately rejected. Updating inputs is a separately reviewed
compatibility change; do not refresh them to make a failing check pass.
Git attributes preserve these pinned input bytes without platform newline conversion.
