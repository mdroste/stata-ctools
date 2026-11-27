# cencode

Convert strings to labeled numeric variables.

## Overview

`cencode` converts string variables to numeric variables with value labels. It accepts multiple variables and can replace the originals.

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
| `threads(#)` | Maximum threads; zero uses the ctools thread setting (automatic cap: 16) |

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

`cencode` uses a two-pass streaming algorithm:

### Pass 1: Read and Code
- Each thread reads its share of the observations directly from Stata; `strL` values are read on Stata's own thread, the only one the plugin interface allows, chunk after chunk
- A thread-local hash table gives each distinct string a provisional code; only those strings and a 4-byte code per observation are kept
- The threads' distinct strings are merged and sorted bytewise; with an existing value label, labeled strings keep their (smallest) value and new ones are numbered after its largest value, as `encode` does

### Pass 2: Write Codes
- Each thread writes the final codes for its observations to the output variable
- The output variable is created without initial values and every observation is written, missing outside `if` and `in`. The outputs are created before the `if` marker, so the marker is the last variable and Stata drops it without a pass over the data

### Value labels
- A label of up to a few thousand mappings is installed through Mata's `st_vlmodify()`. That function copies the whole label for each mapping it adds, so its time grows with the square of the label's size (40,000 mappings of 45 bytes take about 6 seconds). A larger label comes back from the plugin as a 0-observation dataset that `append` loads in milliseconds; a stable sort then restores the sort order, which `append` clears. The plugin chooses by an estimate of both costs
- An existing label too long for a macro reaches the plugin in a file

`encode` handles `strL` values longer than 2,045 bytes or holding binary data. It also handles `noextend` with an unlabeled string, existing labels with extended missing values or text longer than 2,045 bytes, and codes beyond the `long` range.

### Implementation details
1. **Streaming two-pass design**: The C implementation stores each distinct string once and a four-byte code for each observation
2. **Thread-local hash tables**: open addressing over power-of-two capacities, with a word-at-a-time hash; O(1) average lookup per observation without locks; tables are merged once at the end
3. **Arena allocator**: Unique strings are copied into a chain of large blocks with O(1) bump-pointer allocation and a bulk free, eliminating per-string `malloc`/`free` overhead
4. **Alphabetical ordering**: Codes are assigned in alphabetical order (1, 2, 3, ...) after collecting the distinct strings, matching Stata's `encode` convention

## Performance

See the [published profiler results](../benchmarks/readme_profile/results.md)
and [run metadata](../benchmarks/readme_profile/manifest.json) for measured
encoding times. Reproduce them with the [benchmark runner](../performance/README.md).
The C implementation keeps a four-byte code per observation and storage for the
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

`encode` leaves `r()` as it was. `cencode` keeps the results already in `r()` too, also when it fails, and adds the following:

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

The command checks input variables and output names even when there are no observations. Output names must be distinct and unused. With `replace`, it updates the requested variables only after all conversions succeed. A failed operation restores variables and value labels.
