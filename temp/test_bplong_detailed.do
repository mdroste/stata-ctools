adopath ++ build

sysuse bplong, clear
export delimited using temp/test_bplong.csv, replace

* Now import with cimport and check raw values
cimport delimited using temp/test_bplong.csv, clear
describe
list agegrp in 1/5

* Compare types
import delimited using temp/test_bplong.csv, clear
describe
list agegrp in 1/5
