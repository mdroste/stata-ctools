version 17
clear all
set more off
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
foreach f in good_sst truncated_sst cdata comment escaped {
    foreach cmd in import cimport {
        capture noisily `cmd' excel using "/tmp/ctools-io-audit/`f'.xlsx", clear
        local rc=_rc
        di "RESULT `cmd' `f' RC=`rc' N=" _N " K=" c(k)
        capture noisily list, noobs
    }
}
foreach f in good truncated {
    capture noisily cimport shp using "/tmp/ctools-io-audit/`f'.shp", clear
    local rc=_rc
    di "RESULT cimport `f'.shp RC=`rc' N=" _N
    capture noisily list _ID _X _Y, noobs
}
clear
set obs 4
gen str30 s="a"+char(1)+"b" in 1
replace s="a"+char(13)+"b" in 2
replace s="a_x000D_b" in 3
replace s="a_x0041_b" in 4
export excel using "/tmp/ctools-io-audit/native.xlsx", replace
cexport excel using "/tmp/ctools-io-audit/custom.xlsx", replace
foreach f in native custom {
    foreach cmd in import cimport {
        capture noisily `cmd' excel using "/tmp/ctools-io-audit/`f'.xlsx", clear
        local rc=_rc
        di "RESULT `cmd' `f' RC=`rc'"
        capture gen str100 hex=ustrtohex(A)
        capture noisily list, noobs
    }
}
exit, clear
