# Character encoding detection from ICU

`src/cimport/cimport_charset.inc` implements Unicode ICU 67.1's default
character encoding detectors in C. The generated tables in
`src/cimport/generated/cimport_charset_tables.inc` contain byte maps,
three-character language patterns and common multibyte characters. Building
and using ctools requires neither Java nor the ICU libraries.

The implementation follows ICU's
[Java source files](https://github.com/unicode-org/icu/tree/release-67-1/icu4j/main/classes/core/src/com/ibm/icu/text):
`CharsetDetector.java`, `CharsetRecog_UTF8.java`, `CharsetRecog_Unicode.java`,
`CharsetRecog_2022.java`, `CharsetRecog_mbcs.java` and `CharsetRecog_sbcs.java`.
The Unicode license is included in `LICENSE` and the distribution's
`ctools-codecs-LICENSE.txt`.

The C implementation uses ICU's 8,000-byte input sample, zero padding for UTF-16
detection, confidence calculations, detector order and rules for ties. It omits
the EBCDIC detectors that ICU disables by default and the optional HTML tag
filter. When ISO-8859-8-I scores highest, it uses UTF-8, matching Stata's
behavior when Java rejects that decoder name.

From the repository root, regenerate the tables using the saved source files:

```sh
python3 tools/generation/generate_cimport_charsets.py
```

To read another copy of the same ICU sources, pass its directory as an argument.
The script checks the files against `SHA256SUMS`. See
[generated import tables](../../../tools/generation/README.md) for input versions
and verification commands.
