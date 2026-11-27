{smcl}
{* *! version 1.4.0 20261127}{...}
{viewerjumpto "Syntax" "csql##syntax"}{...}
{viewerjumpto "Description" "csql##description"}{...}
{viewerjumpto "Options" "csql##options"}{...}
{viewerjumpto "Remarks" "csql##remarks"}{...}
{viewerjumpto "Installing DuckDB" "csql##install"}{...}
{viewerjumpto "Examples" "csql##examples"}{...}
{viewerjumpto "Stored results" "csql##results"}{...}
{viewerjumpto "Acknowledgments" "csql##acknowledgments"}{...}
{title:Title}

{phang}
{bf:csql} {hline 2} SQL queries on Stata data and on Parquet, CSV, JSON and .dta files, through DuckDB


{marker syntax}{...}
{title:Syntax}

{phang2}
Run a query, and list its result, load it, or save it

{p 12 16 2}
{cmd:csql} [{cmd:using} {it:{help filename}}]{cmd:,} {opt q:uery(sql)} [{it:options}]

{p 12 16 2}
{cmd:csql} [{cmd:using} {it:{help filename}}]{cmd:,} {opt f:ile(filename)} [{it:options}]

{phang2}
Install DuckDB, or show the version in use

{p 12 16 2}
{cmd:csql,} {opt inst:all} [{opt replace} {opt from(location)} {opt lib:rary(filename)}]

{p 12 16 2}
{cmd:csql,} {opt ver:sion} [{opt lib:rary(filename)}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Query}
{p2coldent:* {opt q:uery(sql)}}the SQL to run: one statement or several separated by semicolons{p_end}
{p2coldent:* {opt f:ile(filename)}}read the SQL from {it:filename} instead{p_end}

{syntab:Result}
{synopt:{opt clear}}replace the data in memory with the rows of the last statement{p_end}
{synopt:{cmd:saving(}{it:filename} [{cmd:, replace}]{cmd:)}}save those rows to a .dta, .parquet, .csv or .json file{p_end}
{synopt:{opt row:s(#)}}rows to list when the result is neither loaded nor saved; default is {cmd:rows(20)}{p_end}
{synopt:{opt nolab:el}}do not carry variable and value labels into the result{p_end}
{synopt:{opt nodate:time}}read {cmd:%td} and {cmd:%tc} variables as numbers, not as SQL dates and timestamps{p_end}
{synopt:{opt strl(#)}}store strings longer than {it:#} bytes as strL; default is {cmd:strl(2045)}{p_end}

{syntab:Engine}
{synopt:{opt data:base(filename)}}open a DuckDB database file; default is a temporary database in memory{p_end}
{synopt:{opt read:only}}open {opt database()} read-only{p_end}
{synopt:{opt mem:ory(limit)}}DuckDB's memory limit, such as {cmd:memory(16GB)}; default is 80% of RAM{p_end}
{synopt:{opt temp:dir(path)}}where DuckDB spills data that do not fit in memory{p_end}
{synopt:{opt thr:eads(#)}}number of threads; default is the ctools setting ({help ctools##threads:ctools, threads()}){p_end}
{synopt:{opt lib:rary(filename)}}the DuckDB library to use{p_end}
{synopt:{opt ver:bose}}report timings{p_end}
{synoptline}
{p 4 6 2}* {opt query()} or {opt file()} is required.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:csql} runs SQL in {browse "https://duckdb.org":DuckDB}, an analytical
database engine that works inside the Stata process. It needs no server,
Python or R: DuckDB is one library file, which {cmd:csql, install} downloads
once (see {help csql##install:Installing DuckDB}).

{pstd}
A query can read the data in memory, Stata datasets on disk, and Parquet, CSV
and JSON files, singly or many at a time, and join them. DuckDB reads only the
columns and row groups a query needs, runs on all cores, and moves
aggregations, joins and sorts that do not fit in memory to disk. Files much
larger than memory can therefore be filtered, aggregated and joined, and only
the result, which must fit, is loaded into Stata.

{pstd}
The rows of the last statement are listed in the Results window (the first
{opt rows()} of them), loaded into memory with {opt clear}, or saved with
{opt saving()}. A saved Parquet, CSV or JSON file is written by DuckDB as the
rows are produced, so it may be larger than memory.

{pstd}
SQL names a query's tables in its {cmd:FROM} clause:

{p2colset 9 36 38 2}{...}
{p2col:{cmd:stata}}the data in memory{p_end}
{p2col:{cmd:data}}the {cmd:using} file; the data in memory when there is no {cmd:using} file{p_end}
{p2col:{cmd:'file.dta'}, {cmd:read_dta('file.dta')}}a Stata dataset (formats 113 to 119, Stata 8 and later){p_end}
{p2col:{cmd:'file.parquet'}, {cmd:read_parquet('*.parquet')}}Parquet files; a pattern reads all matching files{p_end}
{p2col:{cmd:'file.csv'}, {cmd:read_csv('file.csv')}}delimited text, with types detected{p_end}
{p2col:{cmd:'file.json'}, {cmd:read_json('file.json')}}JSON, one record per line or an array{p_end}
{p2colreset}{...}

{pstd}
and anything else DuckDB reads, such as the tables of {opt database()}, Excel
files ({cmd:read_xlsx()}) or files on the web and in cloud storage (DuckDB
downloads the extensions it needs). The {browse "https://duckdb.org/docs/stable/sql/introduction":DuckDB documentation}
describes the SQL.


{marker options}{...}
{title:Options}

{dlgtab:Query}

{phang}
{opt query(sql)} gives the SQL. Several statements are separated by
semicolons; they run in order, and the rows of the last one are the result.
Write the SQL within compound double quotes, {cmd:query(`"}...{cmd:"')}, when it
holds double quotes, as quoted names such as {cmd:"foreign"} do. Stata expands
macros in the SQL before {cmd:csql} sees it: write {cmd:\$} for a dollar sign.
In a do-file, {cmd:///} continues the SQL on the next line.

{phang}
{opt file(filename)} reads the SQL from a text file (UTF-8), which may hold
any number of statements and comments and is not subject to macro expansion.

{dlgtab:Result}

{phang}
{opt clear} loads the rows of the last statement into memory, replacing the
data there, as {cmd:use, clear} does. Without {opt clear} and
{opt saving()}, {cmd:csql} lists the rows and leaves the data in memory as
they were.

{phang}
{cmd:saving(}{it:filename} [{cmd:, replace}]{cmd:)} saves the rows to
{it:filename}. The extension sets the format: {cmd:.dta} (or none),
{cmd:.parquet}, {cmd:.csv} or {cmd:.json}. {opt replace} allows an existing
file to be overwritten. With {opt clear} as well, the saved .dta file is then
loaded; {opt clear} cannot be combined with the other formats.

{phang}
{opt rows(#)} sets how many rows to list when the result is neither loaded
nor saved. {cmd:csql} still counts all of them.

{phang}
{opt nolabel} stops {cmd:csql} from carrying variable and value labels from
the data the query read into the result (see
{help csql##labels:Labels and formats}).

{phang}
{opt nodatetime} makes variables with a {cmd:%td} or {cmd:%tc} format read as
plain numbers (days or milliseconds since 1960) instead of SQL {cmd:DATE} and
{cmd:TIMESTAMP} values.

{phang}
{opt strl(#)} stores string columns with values longer than {it:#} bytes as
strL. Columns up to 2,045 bytes otherwise become {cmd:str}{it:#}, the
smallest that holds their longest value.

{dlgtab:Engine}

{phang}
{opt database(filename)} opens (or creates) a DuckDB database file, so that
tables created by one {cmd:csql} call remain for the next. By default each call
works in a new database in memory that disappears when it ends.

{phang}
{opt readonly} opens {opt database()} for reading only, so that other
programs can read it at the same time.

{phang}
{opt memory(limit)} sets how much memory DuckDB may use, for example
{cmd:memory(8GB)} or {cmd:memory(500MB)}. Past it, DuckDB moves data to
{opt tempdir()}. The data in memory in Stata do not count toward the limit.

{phang}
{opt tempdir(path)} sets where DuckDB writes data that do not fit in memory.
The default is Stata's temporary folder ({cmd:c(tmpdir)}), or the folder of
{opt database()}.

{phang}
{opt threads(#)} sets the number of threads for this call. The default is
the ctools setting; see {help ctools##threads:ctools, threads()}.

{phang}
{opt library(filename)} names the DuckDB library to load, for example a
newer release than the one {cmd:csql, install} installs. The global macro
{cmd:CSQL_DUCKDB} does the same for the rest of the session. With
{cmd:csql, install}, {opt library()} is where the library is installed.

{phang}
{opt verbose} reports how long each stage took.

{phang}
{opt install}, {opt replace} and {opt from(location)} are described under
{help csql##install:Installing DuckDB}. {opt version} reports the DuckDB
version and the library file in use.


{marker remarks}{...}
{title:Remarks}

{pstd}
Remarks are presented under the following headings:

{phang2}{help csql##tables:The data in memory and Stata datasets}{p_end}
{phang2}{help csql##types:Types}{p_end}
{phang2}{help csql##labels:Labels and formats}{p_end}
{phang2}{help csql##large:Data larger than memory}{p_end}
{phang2}{help csql##sql:Writing SQL in Stata}{p_end}

{marker tables}{...}
{title:The data in memory and Stata datasets}

{pstd}
{cmd:stata} has one column per variable, in dataset order, and one row per
observation. {cmd:csql} reads only the variables a query uses and does not
copy the dataset. Numeric missing values ({cmd:.} and {cmd:.a} to {cmd:.z})
read as SQL {cmd:NULL}; strings are never {cmd:NULL}, and {cmd:""} reads as the
empty string, so {cmd:WHERE name = ''} selects observations where {cmd:name}
is empty. {cmd:if} and {cmd:in} are not allowed; use {cmd:WHERE} and
{cmd:LIMIT}.

{pstd}
A Stata dataset on disk is read in place, without loading it, and may be
larger than memory: {cmd:SELECT ... FROM 'big.dta'} reads only the columns
the query uses. {cmd:read_dta('file.dta', decode := true)} returns the text of
value labels instead of the values of labeled variables (unlabeled values
become their numbers as text).

{pstd}
A selection or filter without {cmd:ORDER BY} keeps the order of the rows it
reads, as {cmd:keep if} does. Joins and aggregations return rows in an order
of DuckDB's choosing; add {cmd:ORDER BY} when order matters. Queries that do
not need order run faster when they say so, with
{cmd:SET preserve_insertion_order = false;} before them.

{marker types}{...}
{title:Types}

{pstd}
Stata's storage types read as {cmd:TINYINT} (byte), {cmd:SMALLINT} (int),
{cmd:INTEGER} (long), {cmd:FLOAT} (float), {cmd:DOUBLE} (double) and
{cmd:VARCHAR} (str{it:#} and strL). Variables with a {cmd:%td} format read as
{cmd:DATE} and with a {cmd:%tc} format as {cmd:TIMESTAMP}, unless
{opt nodatetime}; other {cmd:%t} formats ({cmd:%tm}, {cmd:%tq}, ...) read as
numbers.

{pstd}
Results are stored in the narrowest Stata type that holds every value, as
{helpb compress} would: integers as {cmd:byte}, {cmd:int} or {cmd:long} by
their range, else {cmd:double}; {cmd:BOOLEAN} as {cmd:byte} (0 or 1);
{cmd:FLOAT} as {cmd:float}; {cmd:DOUBLE} and {cmd:DECIMAL} as {cmd:double};
{cmd:DATE} as {cmd:long} with format {cmd:%td}; {cmd:TIMESTAMP} (any precision;
with time zone, in UTC) as {cmd:double} with format {cmd:%tc}; {cmd:TIME} as
{cmd:double} with format {cmd:%tcHH:MM:SS}; {cmd:VARCHAR} as {cmd:str}{it:#} or
strL; {cmd:BLOB} as a binary strL; and {cmd:ENUM} as integer codes 1, 2, ...
with a value label that holds its values. Other types (lists, structs, maps,
intervals, UUIDs) are converted to text. {cmd:NULL} becomes missing ({cmd:.})
or {cmd:""}.

{pstd}
Stata cannot hold every SQL value. Integers beyond 2^53 lose precision in a
{cmd:double}, and NaN, infinities, doubles beyond {cmd:c(maxdouble)} and
{cmd:FLOAT} values beyond {cmd:c(maxfloat)} become missing; {cmd:csql} reports how many values were affected. Column names that
are not valid Stata names are changed (for example {cmd:count_star()} becomes
{cmd:count_star__}), and the variable label then holds the original name.

{marker labels}{...}
{title:Labels and formats}

{pstd}
A result column with the name of a variable the query read, from the data in
memory or a Stata dataset, takes that variable's variable label and display
format and, when it holds integers, its value label. When several sources
have the variable, the first one the query reads provides them. Labels do not
follow renaming ({cmd:price AS p}), and an aggregate named like a source
variable ({cmd:SUM(sales) AS sales}) keeps that variable's label. Specify
{opt nolabel} to carry no labels.

{pstd}
An {cmd:ENUM} column's value label has the column's name. When two value
labels would share a name (an {cmd:ENUM} named like an inherited label, or
labels of the same name but different contents from two sources), the
inherited label, or the one from the first column that uses it, keeps the
name and the other is named {it:name}{cmd:_1} (then {it:name}{cmd:_2}, ...),
a name no other label has. Same-named labels with the same contents are
written once.

{marker large}{...}
{title:Data larger than memory}

{pstd}
DuckDB streams through files and keeps only what a query needs in memory.
Filters, column selection and aggregations over Parquet files read only the
parts of the files they need. Hash tables for joins and aggregations, and
sorts, move to {opt tempdir()} when they outgrow {opt memory()}. Save large
results as Parquet with {opt saving()}; they need not fit in memory. Results
loaded with {opt clear} or saved as .dta must fit, as must the data in
memory.

{pstd}
A folder of Parquet files written by year, say, is read as one table with
{cmd:read_parquet('sales/*/*.parquet', hive_partitioning := true)}, and a
filter on the partition column skips the other folders.

{marker sql}{...}
{title:Writing SQL in Stata}

{pstd}
Names that are SQL keywords, such as {cmd:foreign}, {cmd:year} in some
contexts, or {cmd:order}, must be quoted, {cmd:"foreign"}, and the SQL then
goes within compound quotes:

{phang2}{cmd:. csql, query(`"SELECT "foreign", avg(price) FROM stata GROUP BY "foreign""')}{p_end}

{pstd}
Strings in SQL take single quotes ({cmd:'cat3'}), as do file names. Long
queries read best from a file with {opt file()}.


{marker install}{...}
{title:Installing DuckDB}

{pstd}
{cmd:csql} uses the official DuckDB library, which is not part of the ctools
package because of its size. {cmd:csql, install} downloads DuckDB
{cmd:2.0.1} for your computer from
{browse "https://github.com/duckdb/duckdb/releases/tag/v2.0.1":github.com/duckdb/duckdb}
(about 20 MB), checks it against the SHA-256 checksum pinned in {cmd:csql.ado},
and saves the library in your PLUS folder ({cmd:c(sysdir_plus)}, subfolder
{cmd:d}). It is needed once per computer. {opt replace} installs it again.

{pstd}
On a computer without internet access, download
{cmd:duckdb-shared-libs-windows-amd64.tar.gz},
{cmd:duckdb-shared-libs-osx-arm64.tar.gz} (Apple Silicon),
{cmd:duckdb-shared-libs-osx-amd64.tar.gz} (Intel Mac) or
{cmd:duckdb-shared-libs-linux-amd64.tar.gz} from that page elsewhere, copy it
over, and type {cmd:csql, install from(}{it:folder or file}{cmd:)}. A library file
({cmd:duckdb.dll}, {cmd:libduckdb.dylib} or {cmd:libduckdb.so}) given to
{opt from()} is copied as it is, without a checksum.

{pstd}
Another DuckDB release (1.1 or later) can be used with {opt library()} or
{cmd:global CSQL_DUCKDB}; {cmd:csql} is tested with 2.0.1. DuckDB's Linux
library needs glibc 2.25 or later and the C++ runtime of GCC 6 or later
(RHEL 8, Ubuntu 18.04 and newer; not CentOS 7, which ctools otherwise
supports). Its macOS library needs macOS 11 or later.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. csql, install}{p_end}
{phang2}{cmd:. sysuse auto}{p_end}

{pstd}List the most expensive cars{p_end}
{phang2}{cmd:. csql, query("SELECT make, price, mpg FROM stata ORDER BY price DESC LIMIT 5")}{p_end}

{pstd}Means by group, loaded into memory; {cmd:foreign} keeps its labels{p_end}
{phang2}{cmd:. csql, query(`"SELECT "foreign", count(*) AS n, avg(price) AS price FROM stata GROUP BY ALL"') clear}{p_end}

{pstd}Aggregate a Parquet file larger than memory{p_end}
{phang2}{cmd:. csql using "transactions.parquet", query("SELECT firm_id, year, SUM(sales) AS sales FROM data GROUP BY firm_id, year") clear}{p_end}

{pstd}Join two Parquet files{p_end}
{phang2}{cmd:. csql, query("SELECT a.*, b.industry FROM read_parquet('sales.parquet') a JOIN read_parquet('firms.parquet') b ON a.firm_id = b.firm_id") clear}{p_end}

{pstd}Join the data in memory with a Stata dataset on disk, without loading it{p_end}
{phang2}{cmd:. csql, query("SELECT s.*, f.industry FROM stata s LEFT JOIN 'firms.dta' f USING (firm_id)") clear}{p_end}

{pstd}Subset a large .dta file, as {cmd:use if} does{p_end}
{phang2}{cmd:. csql using "big.dta", query("SELECT id, year, income FROM data WHERE year >= 2010") clear}{p_end}

{pstd}Convert a folder of CSV files to one Parquet file{p_end}
{phang2}{cmd:. csql, query("SELECT * FROM read_csv('raw/*.csv')") saving("all.parquet", replace)}{p_end}

{pstd}Keep tables in a database file between calls{p_end}
{phang2}{cmd:. csql, database("work.duckdb") query("CREATE OR REPLACE TABLE firms AS SELECT * FROM 'firms.dta'")}{p_end}
{phang2}{cmd:. csql, database("work.duckdb") query("SELECT industry, count(*) FROM firms GROUP BY ALL")}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:csql} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}rows of the last statement (listed, loaded or saved){p_end}
{synopt:{cmd:r(k)}}columns of the last statement{p_end}
{synopt:{cmd:r(changed)}}rows changed, when the last statement changes a table{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(version)}}the DuckDB version{p_end}
{synopt:{cmd:r(library)}}the DuckDB library file{p_end}
{synopt:{cmd:r(saving)}}the file saved by {opt saving()}{p_end}
{synopt:{cmd:r(renamed)}}variables whose names were changed to valid Stata names{p_end}
{p2colreset}{...}


{* ctools-attribution:begin (generated by scripts/sync_attribution.py; edit config/attribution.json)}{...}
{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}
{cmd:csql} is a ctools command, not a reimplementation of a Stata command. Its
queries run in {browse "https://duckdb.org":DuckDB}, developed by Mark
Raasveldt, Hannes Mühleisen, and the DuckDB contributors, and released by the
DuckDB Foundation under the MIT license.

{pstd}
{cmd:csql} is tested with this version of DuckDB:

{p2colset 9 27 29 2}{...}
{p2col:Software}DuckDB{p_end}
{p2col:Authors}Mark Raasveldt, Hannes Mühleisen, and the DuckDB
contributors{p_end}
{p2col:Version}2.0.1{p_end}
{p2col:Source}{browse "https://github.com/duckdb/duckdb/releases/tag/v2.0.1":GitHub (duckdb/duckdb)}{p_end}
{p2col:License}MIT, Stichting DuckDB Foundation{p_end}
{p2colreset}{...}

{pstd}
ctools ships DuckDB's C API header ({cmd:duckdb.h}). The library itself is not
part of ctools: {cmd:csql, install} downloads it from DuckDB's own releases,
and {cmd:library()} loads another DuckDB release (1.1 or later).
{help ctools_licenses:help ctools_licenses} reproduces the DuckDB license.

{pstd}
If you use {cmd:csql} in published work, please cite DuckDB:

{phang2}Raasveldt, M., and H. Mühleisen. 2019. DuckDB: An embeddable
analytical database. In
{it:Proceedings of the 2019 International Conference on Management of Data (SIGMOD '19)},
1981-1984. New York: ACM.
{browse "https://doi.org/10.1145/3299869.3320212":https://doi.org/10.1145/3299869.3320212}.{p_end}

{pstd}
Please report bugs and questions about {cmd:csql} at
{browse "https://github.com/mdroste/stata-ctools/issues":github.com/mdroste/stata-ctools/issues}.
{* ctools-attribution:end}{...}


{title:Also see}

{pstd}{help use}, {help collapse}, {help joinby}, {help cimport}, {help cexport}, {help ctools}
