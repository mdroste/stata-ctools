* Debug cimport failures

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"

* Add build directory to adopath
adopath ++ "build"

* ==== TEST 1: Space delimiter ====
di _n "=== SPACE DELIMITER TEST ==="
file open fh using "temp/space_delim.csv", write replace
file write fh "id name value" _n
file write fh "1 Alpha 100" _n
file write fh "2 Beta 200" _n
file write fh "3 Gamma 300" _n
file close fh

di "File contents:"
type "temp/space_delim.csv"

di _n "Stata import delimited:"
import delimited using "temp/space_delim.csv", delimiters(" ") clear
desc, short
list

di _n "cimport:"
cimport delimited using "temp/space_delim.csv", delimiters(" ") clear verbose
desc, short
list

* ==== TEST 2: All numeric headers ====
di _n "=== ALL NUMERIC HEADERS TEST ==="
file open fh using "temp/numeric_headers.csv", write replace
file write fh "1,2,3,4,5" _n
file write fh "10,20,30,40,50" _n
file write fh "11,21,31,41,51" _n
file close fh

di "File contents:"
type "temp/numeric_headers.csv"

di _n "Stata import delimited:"
import delimited using "temp/numeric_headers.csv", clear
desc, short
list

di _n "cimport:"
cimport delimited using "temp/numeric_headers.csv", clear verbose
desc, short
list

* ==== TEST 3: Log file format ====
di _n "=== LOG FILE FORMAT TEST ==="
file open fh using "temp/log_format.csv", write replace
file write fh "timestamp,level,message" _n
file write fh "2023-01-15 10:30:45.123,INFO,Application started" _n
file write fh "2023-01-15 10:30:46.456,DEBUG,Loading configuration" _n
file write fh `"2023-01-15 10:30:47.789,WARN,"Deprecated API used""' _n
file write fh `"2023-01-15 10:30:48.012,ERROR,"Connection failed: timeout""' _n
file close fh

di "File contents:"
type "temp/log_format.csv"

di _n "Stata import delimited:"
import delimited using "temp/log_format.csv", clear
desc, short
list

di _n "cimport:"
cimport delimited using "temp/log_format.csv", clear verbose
desc, short
list

* ==== TEST 4: Config data ====
di _n "=== CONFIG DATA TEST ==="
file open fh using "temp/config_data.csv", write replace
file write fh "key,value,type" _n
file write fh "max_connections,100,int" _n
file write fh "timeout_ms,30000,int" _n
file write fh "enable_cache,true,bool" _n
file write fh "api_endpoint,https://api.example.com/v2,string" _n
file write fh "rate_limit,1.5,float" _n
file close fh

di "File contents:"
type "temp/config_data.csv"

di _n "Stata import delimited:"
import delimited using "temp/config_data.csv", clear
desc, short
list

di _n "cimport:"
cimport delimited using "temp/config_data.csv", clear verbose
desc, short
list

* ==== TEST 5: Delimiter at start/end ====
di _n "=== DELIMITER AT START/END TEST ==="
file open fh using "temp/delim_edges.csv", write replace
file write fh ",a,b," _n
file write fh ",1,2," _n
file write fh ",3,4," _n
file close fh

di "File contents:"
type "temp/delim_edges.csv"

di _n "Stata import delimited:"
import delimited using "temp/delim_edges.csv", clear
desc, short
list

di _n "cimport:"
cimport delimited using "temp/delim_edges.csv", clear verbose
desc, short
list
