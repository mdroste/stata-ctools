/*******************************************************************************
    benchmark_ips4o.do
    Performance comparison of IPS4o vs other sort algorithms
*******************************************************************************/

clear all
set more off

* Setup
adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

* Results storage
tempname results
postfile `results' str20 algorithm obs double time using "benchmark_ips4o_results.dta", replace

* Test sizes
local sizes 10000 50000 100000 500000 1000000 2000000

di as txt _n "======================================================================"
di as txt "IPS4o Performance Benchmark"
di as txt "======================================================================"

foreach n of local sizes {
    di as txt _n "Testing with `n' observations..."
    di as txt "--------------------------------------------------"

    * Generate random data
    clear
    set seed 12345
    set obs `n'
    gen double value = rnormal()

    * Test each algorithm
    foreach alg in lsd ips4o sample {
        preserve

        timer clear 1
        timer on 1
        quietly csort value, algorithm(`alg')
        timer off 1

        quietly timer list 1
        local t = r(t1)

        di as txt "  `alg':" _col(15) %8.3f `t' " seconds"

        post `results' ("`alg'") (`n') (`t')

        restore
    }

    * Also test native Stata sort for comparison
    preserve

    timer clear 1
    timer on 1
    quietly sort value
    timer off 1

    quietly timer list 1
    local t = r(t1)

    di as txt "  stata:" _col(15) %8.3f `t' " seconds"

    post `results' ("stata") (`n') (`t')

    restore
}

postclose `results'

* Display summary table
di as txt _n "======================================================================"
di as txt "Summary: Time (seconds) by algorithm and dataset size"
di as txt "======================================================================"

use "benchmark_ips4o_results.dta", clear
reshape wide time, i(obs) j(algorithm) string

* Calculate speedup ratios
gen double ips4o_vs_stata = timestata / timeips4o
gen double ips4o_vs_lsd = timelsd / timeips4o

format time* %8.3f
format *_vs_* %8.2f

list obs timelsd timeips4o timesample timestata ips4o_vs_stata ips4o_vs_lsd, noobs table

di as txt _n "IPS4o speedup vs Stata: Higher = faster"
di as txt "IPS4o speedup vs LSD: Higher = faster"

* Cleanup
erase "benchmark_ips4o_results.dta"

di as txt _n "======================================================================"
di as txt "Benchmark complete"
di as txt "======================================================================"
