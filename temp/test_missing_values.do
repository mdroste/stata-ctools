// Test how Stata treats extended missing values
clear
di "Regular missing value: " .
di "Extended missing .a: " .a
di "Extended missing .b: " .b
di "Extended missing .z: " .z

di ""
di "Testing if . < .a: " (. < .a)
di "Testing if .a < .b: " (.a < .b)
di "Testing if .a == .a: " (.a == .a)
di "Testing if .a == .b: " (.a == .b)

di ""
di "Numeric representations:"
di "." _col(10) "= " %21.14e .
di ".a" _col(10) "= " %21.14e .a
di ".b" _col(10) "= " %21.14e .b
di ".c" _col(10) "= " %21.14e .c
di ".z" _col(10) "= " %21.14e .z

di ""
di "Missing value tests for different extended missing values:"
di "missing(.) = " missing(.)
di "missing(.a) = " missing(.a)
di "missing(.b) = " missing(.b)

