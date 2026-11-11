# Writing .dta files from the plugin

Status: the writer is implemented (`src/formats/ctools_dta.c`). The opt-in
`tempdta` options of cimport (`cimport fulldta` for the whole dataset,
`cimport strldta` for the strL columns) and of cmerge (merges with strL
variables) use it. csort keeps native `sort` for strL data, which a temporary
file cannot beat (see [Rollout](#rollout)). All measurements below come from StataNow 19.5 MP
(4-core license) on a Ryzen 9 9950X (16 cores, 32 threads, 94 GB RAM) under
Windows 11, with files in `%TEMP%` (the folder `tempfile` uses), on
2026-11-08. Plugin timings use a zig-built test plugin without OpenMP.
Linux and macOS have not been measured.

## Summary

The plugin would write a complete format 118 (or 119) `.dta` to a
`tempfile`, and the ado file would run `use, clear` on it. This gives the
plugin everything the SPI cannot do: strL values of any length (binary
included), all storage types, new variables, variable and value labels,
formats, characteristics, notes, the data label and a sort marker, all in one
step.

It is not faster across the board. Stata's `use` reads numeric data at about
3 GB/s. The existing SPI store writes into existing variables 3.5 to 6 times
faster than that. The DTA path wins when:

- the command produces a new dataset or new variables, where today the ado
  first creates every variable, with `generate` or Mata's `st_addvar`
  (measured for `generate`: 10 doubles on 20M rows take 1.13 s, more than
  `use` of the whole dataset at 0.52 s);
- the output contains strLs, binary strings or metadata the SPI cannot write.

For in-place updates of existing numeric or `str#` columns (`cwinsor`,
`cinterpolate`, the regressions' generated variables), the SPI should stay.

The writer's design follows from three measured facts:

1. Stata reads strLs in step with the data, in row-major order (details
   below). GSOs must appear in exactly the order of the cells that define
   them. Parallel work must therefore be split by row block, not by column.
2. On NTFS, one thread appending in order is the fastest way to write a
   file: 4.4 GB/s, where parallel positioned writes to the same file reach
   only 2.3 to 2.9 GB/s. Workers should build blocks in parallel, and one
   writer should append them in order.
3. Both the write and `use` scale with file size. Narrow storage types
   (byte instead of double is 8 times fewer bytes) matter more than any
   kernel tuning.

## Why: what the SPI cannot do today

| Limitation | Current workaround |
|---|---|
| `SF_sstore` on a strL returns 0 and stores nothing | cimport writes strL columns through a blob file read in Mata (`_cimport_read_blobs`) or through a strL-only `.dta` (`cimport strldta`); cio stores strL values one column at a time with `st_sstore` |
| No binary strings or values over 2,045 bytes | loaders refuse them (error 109); `_ctools_strlkey` substitutes numeric keys |
| `csort` and `cmerge` cannot write strLs | both fall back to native `sort` and `merge` whenever a strL is present |
| Every value is a double, and the plugin cannot create variables | the ado runs `set obs` and creates every output variable (`generate` or `st_addvar`) before the plugin fills them |
| No value labels, characteristics or data label | separate plugin calls per label (`cio label`), a Mata-built empty 118 file for the data label (`_cio_label_stub`) |

`cimport strldta` (`src/cimport/cimport_impl.c`, `cimport_strldta`) already
proves the approach: it writes strL columns as a `.dta`, and the cimport ado
loads it with `use`. The writer proposed here generalizes it.

## Measurements

### Moving data into Stata (20M observations)

| Operation | File size | Time | Per cell |
|---|---:|---:|---:|
| `use`, 10 doubles | 1.6 GB | 0.52 s | 2.6 ns |
| `save`, 10 doubles | 1.6 GB | 0.74 s | |
| SPI store, 10 doubles (`csort` store phase) | | 0.149 s | 0.75 ns |
| SPI load, 10 doubles (`csort` load phase) | | 0.163 s | |
| `generate double z = .`, 10 times | | 1.13 s | 5.7 ns |
| `use`, 10 bytes | 200 MB | 0.27 s | 1.3 ns |
| SPI store, 10 bytes | | 0.16 s | 0.8 ns |
| `use`, 5 `str16` | 1.6 GB | 0.89 s | 8.9 ns |
| SPI store, 5 `str16` | | 0.147 s | 1.5 ns |

### strLs (2M observations of distinct 40-byte text)

| Operation | File size | Time |
|---|---:|---:|
| `use`, empty memory beforehand | 138 MB | 0.29 s (145 ns per value) |
| `use`, 2M strLs already in memory | 138 MB | 0.51 s |
| `clear` of 2M strLs | | 0.22 s |
| `use`, 1,000 distinct values shared by reference | 16 MB | 0.29 s |
| Mata `st_sstore` into an existing strL | | 0.46 s |
| `use`, the same values as `str40` | 80 MB | 0.025 s |

Two conclusions follow. Sharing repeated strL values reduces file size and
write time, but not load time. And the ado should `clear` before the plugin
runs whenever the result replaces the data, so that freeing the old strLs
does not add to the cost of `use`.

### Writing 1.6 GB to the Windows file cache (Win32 calls via Python ctypes, 8 MB writes, best of 3)

| Pattern | Throughput |
|---|---:|
| One thread, sequential append | 4.4 GB/s |
| `SetEndOfFile` presize, then in-order positioned writes | 3.5 GB/s |
| Presize, then writes in reverse order | 2.05 GB/s |
| Presize, 4 threads, each writing a contiguous quarter | 2.3 GB/s |
| Presize, 4 threads, interleaved chunks | 2.7 GB/s |
| No presize, 4 threads, interleaved chunks | 2.9 GB/s |
| Presize, 8 threads, interleaved chunks | 2.4 GB/s |
| 4 threads, each writing its own file (not usable here) | 10.4 GB/s |

Presizing itself is instant (about 1 ms). Writes that land beyond the current
end of valid data are what cost time. Chunk size also matters: 1 MB and 8 MB
writes ran at 3.6 to 3.9 GB/s, 32 MB writes at 1.9 GB/s.
`FILE_ATTRIBUTE_TEMPORARY` did not change throughput; its benefit is keeping
a short-lived file out of the lazy writer.

`cimport strldta` currently presizes the file and writes blocks from several
threads at their offsets, which is the slower pattern above.

### The writer itself (fresh files, best of 3)

| Data | File | Writer | Of which in `WriteFile` | Stata's `save` |
|---|---:|---:|---:|---:|
| 20M × 10 doubles | 1.6 GB | 0.42 s | 0.41 s | 0.74 s |
| 20M × 5 `str16` | 1.6 GB | 0.42 s | 0.41 s | 0.83 s |
| 20M × 10 bytes | 200 MB | 0.077 s | 0.065 s | 0.27 s |
| 2M distinct 40-byte strLs | 138 MB | 0.040 s | 0.035 s | 0.83 s |

The writer is limited by the single appender's copy into the file cache. It
reaches about 3.8 GB/s for large files, against 4.7 GB/s for a bare append
from Python with nothing else running. With one thread and no overlap,
1.6 GB of doubles took 0.54 to 0.70 s, and `str16` 1.6 s.

### End to end: a new 20M × 10 dataset

| Route | Steps | Time |
|---|---|---:|
| SPI (today) | 10 × `generate` + SPI store | about 1.28 s |
| DTA, doubles | write 1.6 GB (0.42 s) + `use` (0.52 s) | about 0.94 s |
| DTA, values that fit in a byte | write 200 MB (0.08 s) + `use` (0.27 s) | about 0.35 s |
| In-place update of 10 existing doubles: SPI | store | 0.15 s |
| In-place update of 10 existing doubles: DTA | write + `use` | about 0.94 s (do not use) |

## Format facts the writer relies on

Everything here was checked against files Stata 19 wrote itself, and against
hand-built files that Stata either loaded or rejected. Stata 19's `save`
writes format 118 when there are at most 32,767 variables. Stata 19 MP also
loaded a hand-built 119 file with three variables.

### Layout (format 118, little-endian "LSF")

```
<stata_dta><header><release>118</release><byteorder>LSF</byteorder>
  <K>u16</K><N>u64</N><label>u16 n, n bytes</label><timestamp>u8 n (0 or 17), n bytes</timestamp></header>
<map>14 × u64 file offsets</map>
<variable_types>K × u16</variable_types>
<varnames>K × 129 bytes</varnames>
<sortlist>(K+1) × u16, 1-based variable numbers, 0-terminated</sortlist>
<formats>K × 57 bytes</formats>
<value_label_names>K × 129 bytes</value_label_names>
<variable_labels>K × 321 bytes</variable_labels>
<characteristics>{ <ch> u32 n, varname[129], charname[129], text with terminating 0 </ch> }</characteristics>
<data>N records of W bytes, fields packed in variable order, no alignment</data>
<strls>{ "GSO" u32 v, u64 o, u8 t, u32 n, n bytes }</strls>
<value_labels>{ <lbl> u32 n, name[129], 3 pad, table </lbl> }</value_labels>
</stata_dta>
```

- Format 119 differs in three places: K is u32, sortlist entries are u32,
  and a strL cell packs v in 3 bytes and o in 5 (118 uses 2 and 6).
- In 118 the header's N field starts at byte 79 (81 in 119). Its position is
  fixed, so a streaming writer can patch N and the map at the end.
- Type codes: `str1`..`str2045` = 1..2045, strL = 32768, double = 65526,
  float = 65527, long = 65528, int = 65529, byte = 65530.
- Default formats as Stata writes them: byte and int `%8.0g`, long `%12.0g`,
  float `%9.0g`, double `%10.0g`, `str#` `%#s` with a minimum of `%9s`,
  strL `%9s`.
- `str#` cells are NUL-padded. A value of full width has no terminator.
- Notes are characteristics: `var[note0]` holds the count and
  `var[note1]`… hold the text.
- The value label table is `u32 count, u32 textlen, u32 offsets[count],
  i32 values[count], text`. Stata writes entries sorted by value; unsorted
  tables were not tested.
- The timestamp may be empty. Stata writes garbage after the sortlist's
  terminating 0; the writer should write zeros.

### Missing values (bit patterns confirmed against Stata's output)

| Type | Valid range | `.` | `.a` … `.z` |
|---|---|---|---|
| byte | −127..100 | 101 | 102..127 |
| int | −32,767..32,740 | 32,741 | 32,742..32,767 |
| long | −2,147,483,647..2,147,483,620 | 2,147,483,621 | +1 each |
| float | ±1.70141173319e38 | 0x7F000000 | +0x800 each |
| double | ±8.98846567431157e307 | 0x7FE0000000000000 | +2^40 each |

Stata's double range is symmetric: `c(mindouble)` is −8.98846567431157e307,
not −DBL_MAX. A file can hold −1.797e308, and `use` loads it, but the value
lies outside what Stata itself produces.

The writer reduces each missing value to an index m (0 for `.`, 1..26 for
`.a`..`.z`), computed as `(bits − 0x7FE0000000000000) >> 40` for doubles.
NaN, ±Inf, values of 2^1023 or more in magnitude, and patterns at or above
2^1023 that are not canonical all map to `.`. Then it encodes m in the
target type.

### The strL rule (the constraint that shapes the design)

Each strL cell holds (v,o). Stata reads the strls section in step with the
cells of the data, in row-major order (observation by observation, and
variable by variable within each observation). Each cell is one of:

1. `(0,0)`, an empty string; or
2. its own position: then the next GSO in the strls section must carry
   exactly that key; or
3. an earlier position in row-major order, belonging to a strL variable whose
   GSO is already defined: this cross-link shares one stored value.

What Stata accepted and rejected:

| Probe | Result |
|---|---|
| A cell in variable 2 pointing to (1,1) in the same observation (cross-variable sharing) | loaded |
| A cell pointing to an earlier observation of a later variable, (1,2) → (2,1) | loaded |
| Cells of later observations pointing to (1,1) (what Stata itself writes for shared values) | loaded |
| GSOs in reverse order, or ordered variable by variable | rejected: "obs/var contained unexpected values" |
| A pointer to a later observation, or to a later variable of the same observation | rejected: "crosslink a strl to another strl that had not yet been defined" |
| A pointer to a numeric variable's position | rejected: "crosslink a strl to a non-strl variable" |
| A GSO whose own cell does not reference it, or an extra unreferenced GSO | rejected |
| A self-referencing cell with no GSO | rejected: "marker for the start of a strL was not found" |

GSO types: 130 is text, stored with its terminating 0 counted in n; 129 is
binary, stored exactly. Both a binary value with an embedded NUL and an
empty text GSO (n = 1) load correctly.

### Sort marker

Stata trusts the sortlist without checking it: a file claiming
`sorted by x` with unsorted data loads as sorted. The writer may set it only
when the plugin has actually established the order, such as csort's output.
A wrong sortlist silently breaks later `by`, `merge` and `tsset` results.

## Writer design

### Module and API

`src/formats/ctools_dta.h` declares the writer, and `ctools_dta.c`
implements it (owning prefix `ctools_dta_`). The writer never calls the SPI.
Callers describe columns, and the writer reads them in place; it never
copies a whole column.

```c
typedef enum { CTOOLS_DTA_BYTE = 1, CTOOLS_DTA_INT, CTOOLS_DTA_LONG, CTOOLS_DTA_FLOAT,
               CTOOLS_DTA_DOUBLE, CTOOLS_DTA_STR, CTOOLS_DTA_STRL } ctools_dta_type;

typedef enum {
    CTOOLS_DTA_SRC_F64 = 1,  /* const double *, Stata's encoding (numeric types) */
    CTOOLS_DTA_SRC_CSTR,     /* const char *const *, NUL-terminated; NULL reads as "" */
    CTOOLS_DTA_SRC_FLAT,     /* fixed-width slots: data + row * stride */
    CTOOLS_DTA_SRC_SPAN,     /* const ctools_dta_span *: binary-safe, any length */
    CTOOLS_DTA_SRC_PRODUCER  /* filled a row group at a time by spec->producer */
} ctools_dta_source;

typedef struct { const char *ptr; size_t len; int binary; } ctools_dta_span;

typedef struct {
    const char *name;
    ctools_dta_type type;
    int width;                 /* str#: 1..2045, or 0 for the longest value */
    const char *format, *label, *value_label;  /* NULL format: the type's default */
    ctools_dta_source source;
    const void *data;
    size_t stride;             /* FLAT */
    const int64_t *rows;       /* optional gather; a negative entry gives . or "" */
    unsigned flags;            /* CTOOLS_DTA_COL_TEXT: no NUL bytes, skip the scan;
                                  CTOOLS_DTA_COL_SHARE: repeated source rows share a strL */
} ctools_dta_column;

int ctools_dta_write(const char *path, const ctools_dta_spec *spec, ctools_dta_result *result);
```

`ctools_dta_spec` holds the row count, the columns, the format release (0
chooses 118, or 119 past 32,767 variables), the data label, the timestamp,
the sortlist, characteristics (notes included), value labels, an optional
producer, a thread count, a temporary-file flag and the producer keep
budget. `ctools_dta_result` returns the file size, the number of strL
values, the total time, the time spent in file writes and a message on
failure.

The `rows` gather is per column, so a merge can read master and using
columns through different maps. It accepts `int64_t` maps with −1 for "no
row", the same convention as `ctools_stream_var_permuted`.

A producer supplies rows in groups: `fill(ctx, group, phase, out)` sets one
array per producer column (`double *` or `ctools_dta_span *`), and
`release` frees it. fill runs on worker threads, concurrently for different
groups. Since every data record comes before every strL in the file, a group
is needed twice. The writer keeps groups filled in the data pass, up to
`keep_bytes` (default 1 GiB) by the producer's `bytes` estimates, and fills
the rest again for the strL pass. cimport's producer decodes its strL fields
from the mapped CSV exactly as the old writer did.

Not implemented yet: a `NARROWEST` type (the smallest exact storage type),
strL sharing by content, a `stata_data` adapter, and a streaming mode where
N is not known in advance. N would be patched like the map is.

### Pipeline

```
plan  ──►  header + metadata  ──►  data units  ──►  strL units  ──►  value labels  ──►  patch the map
(when       (one buffer,           (built by the pool a wave at a time, appended in order
 needed)     appended)              by one thread as each unit finishes)
```

1. **Plan**, only for in-memory columns that need it: the longest value of
   each `str#` column of width 0, and the strL bytes of each range of 65,536
   rows. This sizes the strL units.
2. **Header and metadata**: one buffer, with the map zeroed. These sections
   are small (about 640 bytes per variable).
3. **Data units**: R rows with R × W about 1 MiB, or the producer's groups.
   A worker encodes column by column into the unit, writing each column's
   cells with stride W at its field offset. There is one loop per storage
   type over constant-size `memcpy`. Integer types refuse fractions and
   values out of range with r(459), and float refuses values past the float
   range. strL cells are computed, not looked up: `(0,0)` when empty, else
   the cell's own (v,o).
4. **strL units**: the same rows, in order, about 4 MiB of copied bytes each
   (or the producer's groups). For each row, then each strL variable, a
   nonempty value emits a 20-byte header and its payload: text (type 130)
   with its terminating NUL, or binary (type 129) when the span says so or
   the bytes hold a NUL. A payload of 64 KiB or more is written straight
   from source memory, not copied. Concatenating the units in order
   satisfies the strL rule.
5. **Value labels and the closing tag**, then one positioned rewrite of the
   112-byte map, inside the part of the file already written.

### Threads and I/O

- **Ordered pipeline.** The pool (`ctools_get_global_pool`) builds a wave
  of units into one of two sets of reusable task buffers. A persistent
  appending thread takes the waves in order and writes each unit as soon as
  its worker signals it is done. It overlaps the rest of that wave and the
  next one. A set is reused only after its previous wave is written.
- **Few encoders for in-memory columns.** Encoding is several times faster
  than the append, and extra encoders only take memory bandwidth from the
  appender's copy into the file cache. Waves of in-memory columns hold 4
  units (`DTA_IO_WAVE_UNITS`). With 16 encoders, 1.6 GB of doubles took
  0.54 s; with this limit, 0.42 to 0.46 s, the time `WriteFile` alone took.
  Producer passes (cimport's parsing) use every thread.
- **Windows:** `CreateFileW` with `FILE_FLAG_SEQUENTIAL_SCAN` (plus
  `FILE_ATTRIBUTE_TEMPORARY` when `spec->temporary` is set), no presize,
  `WriteFile` in pieces of at most 4 MiB, no flush. Overwriting an existing
  cached file instead of writing a new one cost about 20% (3.8 against
  4.7 GB/s), another reason to write each result to a fresh `tempfile`.
- **POSIX:** `open(O_WRONLY|O_CREAT|O_TRUNC|O_CLOEXEC)`, `write` in the same
  pieces, `pwrite` for the map, no `fsync`. Unmeasured.
- Worker and appending threads never call the SPI. On any failure the file
  is closed and removed: r(198) for an invalid spec, r(459) for a value that
  does not fit its type, r(603) or r(693) for the file, r(920) for memory.
  The message says which variable and observation.

### Sharing strL values

- **By source row** (implemented, column flag `CTOOLS_DTA_COL_SHARE`): when
  `rows` repeats a source row, as in m:1 merges or bootstrap expansion, the
  first output observation using a row defines the value, and later ones
  cross-link to it. `first[j] = min{i : rows[i] = j}` takes one pass before
  the plan. This is what Stata's own files hold: native merge's saved m:1
  result (2M rows) had 2.4M stored values and 1.42M cross-links. Here it
  matters for load time too. Without sharing, cmerge's file held 3.82M
  stored values and `use` took 0.98 s from empty memory, against 0.63 s for
  Stata's file. (A synthetic file whose 2M cells shared 1,000 values loaded
  no faster than 2M distinct ones, so the gain depends on the data.)
- **By content** (not implemented): a concurrent hash on the bytes that
  stores the minimum row-major position seen (CAS-min), in one pass before
  the data units. Since positions are totally ordered, the result does not
  depend on thread timing, and every cross-link points backwards as the
  strL rule requires.

### Validation the writer does itself

Stata trusts the file (it loaded a false sort marker without complaint), so
before writing anything the writer refuses:

- invalid or duplicate names: 1 to 32 characters (128 bytes) of letters,
  digits and `_`, not starting with a digit, not a reserved word or `str#`;
- a format that does not start with `%`, is longer than 56 bytes, or does
  not match the type (`…s` for strings only);
- labels too long for their fields;
- a sortlist with a repeated or out-of-range variable;
- K > 32,767 in format 118, or more observations than the strL cell can
  number (2^48 in 118, 2^40 in 119);
- invalid characteristics or value labels (duplicate values or names).

While writing, it refuses a `str#` value longer than its width and a strL
value over 2,000,000,000 bytes (r(459)).

### Testing

`validation/stata/dta.do` (45 checks) drives the internal plugin command
`dta` (`src/formats/ctools_dta_plugin.c`):

- `dta selftest` writes every storage type with all 27 missing values and
  the boundary values, full-width `str1` and `str2045`, and strLs: empty,
  1, 2,045, 2,046 and 100,000 bytes, binary with and without NULs, and
  Unicode. It also writes a value label given out of order, notes,
  characteristics, a data label and a sort marker. Both formats, with 1, 2
  and the default number of threads, at 1, 3,000 and 50,000 rows.
- `dta save` round trips random datasets (0 to 40,000 rows, binary strLs
  from `fileread()`) and compares them with `cf` and `describe`.
- `dta producer` covers producer groups beside an in-memory column, kept and
  filled again, and checks that every fill is released.
- `dta error` covers the refusals, and that no file is left behind.
- `dta wide` writes 32,768 variables as format 119.
- A cimport test compares a 70,000-row CSV with two strL columns against
  `import delimited` on each route: default (Mata), `tempdta` (`fulldta`,
  since the strL text outweighs the other columns), and `tempdta` with
  `global CTOOLS_cimport_fault "fulldta"`. That validation hook in
  cimport.ado makes `fulldta` decline, so the import takes `strldta`. Each
  route must also leave `c(changed)` at 1.

## Ado-side integration

```stata
tempfile __ctools_out
* the path travels in a local, as cimport already does (spaces, Unicode)
plugin call ctools_plugin `vars', "<cmd> ... dta"
quietly use `"`__ctools_out'"', clear
capture erase `"`__ctools_out'"'
```

- Run `clear` (or drop what will be replaced) before the plugin call when
  the result replaces the data, so `use` does not also pay to free the old
  dataset.
- Restore `c(filename)` and `c(filedate)` after `use`, as cimport already
  does with `S_FN` and `S_FNDATE`.
- `use` resets `c(changed)` to 0. Each command must match its native
  counterpart's `c(changed)` (see the cimport note below).
- Value labels, characteristics and the data label belong in the file, so
  they arrive with the data, rather than in separate calls.
- Format 118 loads in Stata 14 and later, matching the ado's `version 14.1`.
  Use 119 only when K > 32,767, which needs Stata 15 MP or later. Check
  `c(maxvar)` before writing.
- Import routes through a temporary `.dta` are opt-in (decided 2026-11-08):
  every import would write its whole result to disk, and writes wear SSDs.
  `FILE_ATTRIBUTE_TEMPORARY` asks Windows to keep the file in cache, but
  that is a hint, not a guarantee. cimport's `tempdta` option enables them.
  It tries `cimport fulldta` (every column through the writer, filled by the
  store workers), then `cimport strldta`, then the Mata transfer, which is
  also the default. Both temporary-file routes run only when the strL text
  is at least as large as the other columns' data.
- Under the ado's `version 14.1`, neither characteristics nor labels set
  `c(changed)` after `use`. cimport adds and drops a variable to match
  native (about 1 ms).

## Where to use it

| Command | Today | With the writer |
|---|---|---|
| `cimport` | by default `set obs`, every variable created in Mata (`_cimport_addvars`), SPI fill, strLs through Mata | done, as option `tempdta`: the whole dataset in one file and one `use` (`fulldta`), or the strL columns (`strldta`) |
| `cio` import (SAS, SPSS, XPORT, Excel) | variables created one at a time, `st_sstore` strLs, per-label plugin calls, a Mata-built stub for the data label | data, strLs, value labels, formats, variable labels and the data label in one file |
| `csort` with any strL variable | native `sort` | keep native `sort`: it moves 8-byte strL references, and `use` of the result alone costs more (see Rollout) |
| `cmerge` with kept strLs (merge refuses strL keys) | native `merge` | done, as option `tempdta`: plugin merge, output written with sharing by source row |
| commands that expand or rebuild the dataset (`cbsample`, `crangejoin`) | per command | candidates once the writer exists |
| in-place numeric or `str#` updates (`cwinsor`, `cinterpolate`, predictions) | SPI | keep the SPI: 3.5 to 6 times faster |

## Rollout

1. **Writer and tests.** Done: see [Testing](#testing). There are no C unit
   tests with the SPI mocked yet. Sharing, `NARROWEST` and Unicode labels at
   their byte limits are not covered, and the Python probe tools are not in
   the repository.
2. **Switch `cimport strldta` to the writer.** Done. On a 2M-row CSV with
   two strL columns (about 210 MB of strLs), the step takes 0.10 to 0.11 s,
   against 0.085 to 0.09 s for the old writer, which presized the file and
   wrote from every thread. The `use` that follows takes 0.5 to 0.8 s either
   way, and both produce the same `datasignature`. The single appender sets
   the pace here.
3. **Write the whole cimport result**, as an option. Done (`tempdta`, see
   [Ado-side integration](#ado-side-integration)), by the cimport session.
   On a 500 MB CSV of 75,784 rows × 131 columns (28 strL), it reported
   1.41 s for the old code, 1.14 s by default and 0.97 s with `tempdta`. The
   writer took about 0.24 s for about 640 MB, and `use`, about 0.5 s, is the
   floor.
4. **strL support for `csort` and `cmerge`.**
   - **csort: not worth it.** Native `sort` is as fast with a strL present
     as without one (2.17 s and 2.56 s for 10M rows with 6 doubles), because
     Stata moves only an 8-byte reference per strL. `use` of the sorted
     result alone took 2.35 s, so a temporary file cannot beat the fallback.
   - **cmerge: done, as option `tempdta`.** Phase 1 reads the using strL
     columns whole into the plugin's using cache. Execute then writes the
     merged dataset to a `.dta` instead of storing through the SPI:
     - master columns are gathered through the output's master rows, and
       new using columns through its using rows, both shared by source row;
     - keys, shared variables (update and replace rules) and `_merge` are
       assembled first.

     The ado keeps every schema step it already had (recasts, the template
     append for labels, notes and characteristics), then uses the file and
     puts back formats, labels, value labels, characteristics and the data
     label from a Mata snapshot. Before the `use`, rollback works as before.
     After it, the merged data stand.
   - **Measured** (m:1, master with key, 5 doubles and a 40-byte strL; using
     with 3 doubles and a strL, a fifth of the rows):

     | Master rows | native `merge` | `cmerge ..., tempdta` |
     |---|---:|---:|
     | 2M | 2.0 to 2.4 s | 2.2 to 2.7 s |
     | 10M | 14.2 to 17.3 s | 7.6 to 9.3 s |

     Native merge slows sharply with strLs (2.75 s without them at 10M). At
     10M rows the `tempdta` time is mostly fixed costs:
     - reading the master strLs through the SPI, one value at a time (1.8
       to 3.6 s);
     - freeing the old master and `use` (about 5 s);
     - the ado's template append and recasts (1.5 s), which a `.dta` route
       does not need. Avoiding them means re-implementing what `append` does
       for notes and characteristics, so it is left for later.
   - **SPI strL reads stay serial.** Reading rows of one strL variable from
     several threads corrupted Stata's heap (exit 0xC0000374), so values
     are read one at a time. Each value is read into a 256-byte slot first,
     and only longer values are measured and read again. That cut reading
     2M values from 1.57 s (three SPI calls each) to 0.42 s. The shared
     loaders read strL variables in parallel only one variable per thread.
5. **`cio` import** with labels and metadata in the file.
6. **Optional:** measure on Linux and macOS, `NARROWEST` typing for commands
   without a native type rule to match, content-based sharing, and a test of
   parallel mapped-file writes. A single-threaded mapped fill measured no
   faster than appending; parallel fills were not tried.

## Reading .dta files in the plugin (later)

The same module can parse 118, 119, and Stata 18's 120 and 121 (alias
variables). That would give the plugin full strL and binary access and the
real storage types. But the round trip starts with `save`, which took 0.74 s
for the 1.6 GB dataset, where the SPI load took 0.16 s. So the read path
pays only for strLs, binary data or metadata, the same trade-off as writing.

## Open questions

- Does `use` validate variable names, formats or value label order? The
  writer validates them itself either way.
- A cleaner way to leave `c(changed)` = 1 after `use` where native commands
  do (cimport adds and drops a variable).
- Is reading different strL variables from different threads safe? The
  shared loaders do it, but rows of one variable read concurrently corrupted
  Stata's heap.
- Throughput of ordered append, and of Stata's `use`, on Linux (where
  `c(tmpdir)` may be tmpfs) and on macOS (APFS).
- `use` of our files in Stata 14 to 18. Only Stata 19 was tested.
