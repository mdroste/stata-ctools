# csql

Run SQL on the data in memory and on Parquet, CSV, JSON and Stata files,
through [DuckDB](https://duckdb.org), without Python, R or a database server.

`csql` is a new ctools command, not a replacement for a Stata command.

## Overview

```
                      +-----------------------------+
 data in memory ----> |                             | ----> list (rows())
 'big.dta'      ----> |   DuckDB, inside Stata:     | ----> data in memory (clear)
 '*.parquet'    ----> |   filters, joins, GROUP BY, | ----> .dta / .parquet /
 '*.csv', JSON  ----> |   spills to disk past its   |       .csv / .json (saving())
 database file  ----> |   memory limit              |
                      +-----------------------------+
```

DuckDB is an analytical SQL engine that runs inside the Stata process. It
reads only the columns (and, in Parquet files, the row groups) a query needs,
uses every core, and moves joins, aggregations and sorts that outgrow memory to
disk. Files much larger than memory can be filtered, joined and aggregated;
only the result is loaded into Stata.

## Syntax

```stata
csql [using filename], query(sql) | file(filename) [clear saving(filename [, replace])
    rows(#) nolabel nodatetime strl(#) database(filename) readonly memory(limit)
    tempdir(path) threads(#) library(filename) verbose]

csql, install [replace from(location) library(filename)]
csql, version [library(filename)]
```

The rows of the last statement are listed (`rows()`, 20 by default), loaded
with `clear`, or saved with `saving()`. See `help csql` for every option.

## Tables

| In SQL | Reads |
| --- | --- |
| `stata` | the data in memory |
| `data` | the `using` file; the data in memory without one |
| `'file.dta'`, `read_dta('file.dta')` | a Stata dataset, formats 113–119, read in place |
| `'file.parquet'`, `read_parquet('dir/*.parquet')` | Parquet files; patterns read many |
| `'file.csv'`, `read_csv(...)` | delimited text |
| `'file.json'`, `read_json(...)` | JSON |

Anything else DuckDB reads works too: tables of a `database()` file,
`read_xlsx()`, files on the web or in S3 (DuckDB installs the extensions it
needs the first time).

## Examples

```stata
csql, install                       // once: downloads DuckDB 2.0.1 (about 20 MB)

* Aggregate a Parquet file larger than memory
csql using "transactions.parquet", query( ///
    "SELECT firm_id, year, SUM(sales) AS sales ///
     FROM data GROUP BY firm_id, year") clear

* Join Parquet files
csql, query("SELECT a.*, b.industry ///
    FROM read_parquet('sales.parquet') a ///
    JOIN read_parquet('firms.parquet') b ///
    ON a.firm_id = b.firm_id") clear

* Join the data in memory with a Stata file, without loading the file
csql, query("SELECT s.*, f.industry FROM stata s LEFT JOIN 'firms.dta' f USING (firm_id)") clear

* Subset a big .dta file
csql using "panel.dta", query("SELECT id, year, wage FROM data WHERE year >= 2010") clear

* Convert a folder of CSV files to one Parquet file, larger than memory if need be
csql, query("SELECT * FROM read_csv('raw/*.csv')") saving("all.parquet", replace)

* Keep tables between calls
csql, database("work.duckdb") query("CREATE OR REPLACE TABLE firms AS FROM 'firms.dta'")
csql, database("work.duckdb") query("SELECT industry, count(*) FROM firms GROUP BY ALL")
```

Names that are SQL keywords (`foreign`, `order`, ...) must be quoted, and the
SQL then goes in compound quotes: ``query(`"SELECT "foreign" FROM stata"')``.
Stata expands `$` and backquoted macros before `csql` sees the SQL; `file()`
reads SQL that is not expanded.

## Types, labels and formats

The data in memory and `.dta` files read as `TINYINT`, `SMALLINT`, `INTEGER`,
`FLOAT`, `DOUBLE` and `VARCHAR`; `%td` and `%tc` variables as `DATE` and
`TIMESTAMP` (unless `nodatetime`). Missing values (`.`, `.a`–`.z`) are `NULL`;
strings are never `NULL` (`""` is the empty string). SQL compares names
without regard to case, so a variable whose name differs from an earlier one
only in case (`X` after `x`) reads as `X_2`.

Results get the narrowest Stata type that holds their values, as `compress`
would: integers as `byte`/`int`/`long` by range, else `double`; `DATE` as
`long` `%td`; `TIMESTAMP` as `double` `%tc`; `VARCHAR` as `str#` or `strL`;
`BLOB` as binary `strL`; `ENUM` as integer codes with a value label. Lists,
structs, maps, intervals and UUIDs become text. A result column with the name
of a variable the query read takes its variable label, display format and,
for integers, its value label (`nolabel` turns this off). Invalid Stata names
are changed and keep the original as the variable label.

Integers beyond 2^53 lose precision as `double`; NaN and infinities become
missing. `csql` reports how many values each affects.

## Speed

Measured on 20 million rows (Windows, 16 threads, `firm_id`, `year`, `sales`,
a `str10` category):

| Task | Stata | csql |
| --- | ---: | ---: |
| Sum by `firm_id year` in memory | `collapse` 2.5 s | 0.5 s |
| The same from a Parquet file | — | 0.7 s |
| Group-by over a `.dta` file, not loaded | — | 0.1–0.2 s |
| Load the rows of one year from a `.dta` file | `use if` 0.54 s | 0.22–0.26 s |
| Load all rows from Parquet | (`use` of `.dta` 0.27 s) | 1.1 s |

Loading a whole file into Stata is bound by Stata's `use`, which `csql` calls
on the `.dta` file it writes; `csql` is for the work before that.

## How it works

`csql` loads the official DuckDB library through its C API. It registers the
data in memory and the `.dta` reader as table functions that read only the
requested columns, decode in parallel, and keep row order. The rows of the last
statement are measured once (integer ranges, string lengths) and written by
ctools' parallel `.dta` writer, which Stata then loads with `use`; Parquet,
CSV and JSON results are written by DuckDB's `COPY` as they are produced.
Break interrupts a running query. See the [developer
guide](../development/README.md#sql-through-duckdb-csql).

## Installing DuckDB

The library is not in the ctools package (it is 60–100 MB per platform).
`csql, install` downloads the pinned DuckDB release for the platform from
[github.com/duckdb/duckdb/releases](https://github.com/duckdb/duckdb/releases/tag/v2.0.1),
checks its SHA-256 checksum, and saves the library in PLUS. Without internet
access, download `duckdb-shared-libs-<platform>.tar.gz` elsewhere and use
`csql, install from(folder)`. Another DuckDB (1.1 or later) works through
`library()` or `global CSQL_DUCKDB`.

DuckDB's Linux library needs glibc 2.25 and GCC 6's C++ runtime or newer (RHEL
8, Ubuntu 18.04 and later; not CentOS 7); its macOS library needs macOS 11.

## Citation

Raasveldt, M., and H. Mühleisen. 2019. DuckDB: An embeddable analytical
database. *Proceedings of the 2019 International Conference on Management of
Data (SIGMOD '19)*, 1981–1984. https://doi.org/10.1145/3299869.3320212
