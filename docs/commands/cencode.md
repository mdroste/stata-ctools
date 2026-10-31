# cencode

High-performance C-accelerated string encoding for Stata datasets.

## Overview

`cencode` is a drop-in replacement for Stata's `encode` command. It converts string variables to numeric variables with value labels, using an optimized streaming algorithm.

## Syntax

```stata
cencode varlist [if] [in], {generate(newvarlist) | replace} [options]
```

## Options

| Option | Description |
|--------|-------------|
| `generate(newvarlist)` | One distinct new output name per source variable; required unless `replace` is specified |
| `replace` | Replace each source with its encoded numeric variable; cannot be combined with `generate()` |
| `label(name)` | Name for the value label (default: same as newvar) |
| `noextend` | Do not extend existing value label |
| `verbose` | Display timing breakdown |
| `threads(#)` | Maximum workers; zero uses the ctools thread setting (automatic cap: 16) |

## Examples

```stata
* Basic usage - encode a string variable
sysuse auto, clear
cencode make, generate(make_code)

* View the encoded values and labels
tabulate make_code
label list make_code

* Encode with a custom label name
cencode make, generate(make_num) label(car_makes)

* Encode with verbose timing output
cencode make, generate(make_verbose) verbose

* Encode only foreign cars
cencode make if foreign == 1, generate(foreign_make)
```

## Algorithm

`cencode` uses a two-pass streaming algorithm for fixed-width string variables:

### Pass 1: Read and Code
- Each thread reads its share of the observations directly from Stata
- A thread-local hash table gives each distinct string a provisional code; only those strings and a 4-byte code per observation are kept
- The threads' distinct strings are merged and sorted bytewise; with an existing value label, labeled strings keep their (smallest) value and new ones are numbered after its largest value, as `encode` does

### Pass 2: Write Codes
- Each thread writes the final codes for its observations to the output variable
- Without `if`/`in`, the output variable is created without initial values and every observation is written

`encode` handles strL inputs. It also handles `noextend` with an unlabeled string, existing labels with extended missing values or text longer than 2,045 bytes, and codes beyond the `long` range.

### Speedup Tricks
1. **Streaming two-pass design**: Only K unique strings and a 4-byte code per observation are stored in memory, not the strings themselves
2. **Thread-local FNV-1a hash tables**: O(1) average lookup per observation without locks; tables are merged once at the end
3. **Arena allocator**: Unique strings are bulk-allocated in a contiguous arena with O(1) bump-pointer allocation and a single bulk free, eliminating per-string `malloc`/`free` overhead
4. **Alphabetical ordering**: After collecting uniques, codes are assigned in sorted order (1, 2, 3, ...) matching Stata's `encode` convention

## Performance

See the [published profiler results](../benchmarks/readme_profile/results.md)
and [run metadata](../benchmarks/readme_profile/manifest.json) for measured
encoding times. Reproduce them with the [benchmark runner](../performance/README.md).
The streaming C path keeps a four-byte code per observation and storage for the
unique strings, so memory use depends on both the number of rows and the number
and length of distinct values.

## Timing Breakdown

When using the `verbose` option, you'll see:

```
cencode timing breakdown:
-------------------------------------------------------
  C plugin internals:
    Argument parsing:         0.0000 sec
    Read and code strings:    0.0343 sec
    Merge unique values:      0.0000 sec
    Sort unique:              0.0000 sec
    Encode values:            0.0376 sec
    Store labels:             0.0000 sec
  -----------------------------------------------------
    C plugin total:           0.0719 sec
  -----------------------------------------------------
  Stata overhead:
    Pre-plugin setup:         0.0070 sec
    Plugin call overhead:     0.0001 sec
    Post-plugin labels:       0.0040 sec
-------------------------------------------------------
    Wall clock total:         0.0830 sec
-------------------------------------------------------

  Unique values encoded:    100
```

## Stored Results

`cencode` stores the following in `r()`:

| Scalar | Description |
|--------|-------------|
| `r(N_unique)` | Number of unique strings encoded for the last source variable |
| `r(N_vars)` | Number of variables encoded |

## Limitations

- Maximum of 65,536 unique values (Stata value label limit)
- Empty strings and observations excluded by `if`/`in` become numeric missing values, including with `replace`

## Comparison with Native `encode`

| Feature | Stata `encode` | `cencode` |
|---------|----------------|-----------|
| Multiple source variables | No | Yes |
| Replace source variables | No | Yes |
| String coding | Native command | Hash tables plus bytewise sorting of unique strings |

## See Also

- [ctools Overview](../../README.md)
- Stata's `encode` command: `help encode`
- Stata's `decode` command for reverse operation

## Validation and failure behavior

The complete source and output specification is validated even with zero observations. Generated targets must be distinct new names. A failed multi-variable operation restores source variables, outputs, and value labels; successful replace commits all requested variables together.
