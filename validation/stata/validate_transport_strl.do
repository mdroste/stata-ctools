/* Real-SPI strL regression. Supply the standalone benchmark_transport plugin.
   The driver must invoke this through the stata shell alias and finish with exit, clear. */
args transport_plugin
clear all
set more off
program transport_io, plugin using("`transport_plugin'")
forvalues threads = 1(7)8 {
    foreach width in 8 244 2045 {
        clear
        quietly set obs 25000
        quietly gen strL text1 = "x"
        forvalues power = 1/11 {
            quietly replace text1 = text1 + text1
        }
        quietly replace text1 = substr(string(_n,"%06.0f") + text1,1,`width')
        quietly gen double number = cond(mod(_n,19)==0,.z,_n/7)
        quietly gen strL text2 = "café_" + string(_n)
        quietly replace text1 = "" if mod(_n,17)==0
        quietly gen str12 fixed = string(_n)
        quietly datasignature
        local signature "`r(datasignature)'"
        foreach vars in "text1" "text1 text2" "number text1 fixed text2" "text1 number text2 fixed" {
            foreach rep in 1 2 3 {
                plugin call transport_io `vars', "strl" "`threads'" "read" "1"
                plugin call transport_io `vars' if mod(_n,3)!=0 in 3/24999, "strl" "`threads'" "filteredread" "1"
            }
        }
        quietly datasignature
        assert "`r(datasignature)'" == "`signature'"
    }
}
/* Cross production scheduling thresholds with one to three columns. */
clear
quietly set obs 250003
quietly gen strL text = "café_" + string(_n)
quietly gen strL text2 = "second_" + string(_n)
quietly gen double number = cond(mod(_n,19)==0,.z,_n/7)
foreach threads in 1 8 12 {
    foreach vars in "text" "text text2" "text text2 number" {
        local __ctools_strw "64,64,0"
        plugin call transport_io `vars', "large_strl" "`threads'" "read" "1"
        plugin call transport_io `vars' if mod(_n,3)!=0, "large_strl" "`threads'" "filteredread" "1"
        local __ctools_strw ""
        plugin call transport_io `vars', "large_strl" "`threads'" "read" "1"
    }
}
/* Bounds and write rejection must remain explicit errors. */
clear
quietly set obs 3
quietly gen strL text = "x"
forvalues power = 1/11 {
    quietly replace text = text + text
}
quietly replace text = substr(text,1,2046)
capture plugin call transport_io text, "oversize" "8" "read" "0"
assert _rc == 5
quietly replace text = "short"
capture plugin call transport_io text, "write" "8" "identity" "0"
assert _rc == 5
display "TRANSPORT_STRL_COMPLETE"
