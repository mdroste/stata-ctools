# C adaptation of ICU charset recognition

`src/cimport/cimport_charset.inc` adapts Unicode ICU 67.1's default charset
recognizers into portable C. `cimport_charset_tables.inc` contains its public
byte maps, language trigrams and common multibyte character tables. Java and
ICU libraries are not build or runtime dependencies.

Source: https://github.com/unicode-org/icu/tree/release-67-1/icu4j/main/classes/core/src/com/ibm/icu/text
Files: `CharsetDetector.java`, `CharsetRecog_UTF8.java`,
`CharsetRecog_Unicode.java`, `CharsetRecog_2022.java`, `CharsetRecog_mbcs.java`,
and `CharsetRecog_sbcs.java`. The Unicode license is included in `LICENSE` and
the distribution's `ctools-codecs-LICENSE.txt`.

The C implementation retains the 8,000-byte InputStream sample, zero padding
for UTF-16 recognition, confidence calculations, recognizer order and tie
behavior. Default-disabled EBCDIC recognizers and optional HTML tag filtering
are omitted. A winning ISO-8859-8-I match falls back to UTF-8, following the
installed native command's behavior when Java rejects this decoder name.

Regenerate the frequency tables from the public source directory:

```sh
python3 validation/generate_cimport_charsets.py /path/to/icu/text
```
