* args: directory containing the plugin and ado files, output log path
args plugin_dir log_path
clear all
set more off
adopath ++ "`plugin_dir'"
log using "`log_path'", text replace
* Duplicates that do not match the opposite side; empty-side bypasses.
foreach side in master using {
    foreach empty in 0 1 {
        foreach strings in 0 1 {
            clear
            set obs 3
            gen long k = cond(_n == 1,1,2)
            gen str12 payload = "duplicate"
            if `strings' tostring k, replace
            tempfile duplicate other
            save `duplicate'
            clear
            set obs 1
            gen long k = 99
            gen str12 other = "unique"
            if `strings' tostring k, replace
            if `empty' keep if 0
            save `other'
            if "`side'" == "master" {
                use `duplicate', clear
                capture cmerge 1:1 k using `other'
            }
            else {
                use `other', clear
                capture cmerge 1:1 k using `duplicate'
            }
            assert _rc == 459
        }
    }
}
* Borrowed key strings, using-only keys, duplicates and shared payload.
clear
input str4 k double x str10 s
"b" 2 "bee"
"" 0 "empty"
"a" 1 "aye"
"b" 3 "bee2"
end
tempfile m u ref
save `m'
clear
input str4 k double z str10 t
"d" 4 "dee"
"a" 1 "aye"
"b" 2 "bee"
"" 0 "empty"
end
save `u'
use `m', clear
merge m:1 k using `u'
sort k x
save `ref'
use `m', clear
cmerge m:1 k using `u'
sort k x
cf _all using `ref'
* Sorted inputs with inserted using-only rows, shared columns, and duplicate
* master groups. Compare every value against merge; include both key types.
foreach strings in 0 1 {
    foreach sorted_master in 0 1 {
        foreach shared in 0 1 {
            clear
            set obs 3
            gen long k = cond(_n == 2, 1, 3)
            gen double x = cond(_n == 1, 30, cond(_n == 2, 10, .))
            gen str12 s = cond(_n == 1, "three", cond(_n == 2, "one", ""))
            if `strings' tostring k, replace
            if `sorted_master' sort k
            tempfile master using reference
            save `master'
            clear
            set obs 2
            gen long k = _n + 1
            gen double x = cond(_n == 1, 20, 33)
            gen str12 s = cond(_n == 1, "two", "using3")
            gen double z = 100 * k
            if !`shared' drop x s
            if `strings' tostring k, replace
            sort k
            save `using'
            foreach opts in "" "update" "update replace" {
                use `master', clear
                local genopt ""
                if "`opts'" != "" local genopt "nogenerate"
                merge m:1 k using `using', `opts' `genopt'
                sort k x s
                save `reference', replace
                use `master', clear
                cmerge m:1 k using `using', `opts' `genopt'
                sort k x s
                cf _all using `reference'
            }
        }
    }
}
di "BIGDATA_EDGE_PASS"
log close
exit, clear
