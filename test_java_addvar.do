* Test Java variable creation vs Mata _st_addvar
clear all
set more off

local classpath "/Users/Mike/Documents/GitHub/stata-ctools"

* Create test dataset with many observations
set obs 50000000

* Test 1: Mata _st_addvar (8 variables)
di as text _n "Test 1: Mata _st_addvar (8 vars)..."
timer clear 1
timer on 1
mata: _st_addvar(("double","double","double","long","long","str20","str20","byte"), ("v1","v2","v3","v4","v5","v6","v7","v8"))
timer off 1
timer list 1
drop v1-v8

* Test 2: Java sequential (8 variables, one at a time)
di as text _n "Test 2: Java AddVar sequential (8 vars)..."
timer clear 2
timer on 2
javacall AddVar3 create, classpath(`classpath') args(double v1)
javacall AddVar3 create, classpath(`classpath') args(double v2)
javacall AddVar3 create, classpath(`classpath') args(double v3)
javacall AddVar3 create, classpath(`classpath') args(long v4)
javacall AddVar3 create, classpath(`classpath') args(long v5)
javacall AddVar3 create, classpath(`classpath') args(str20 v6)
javacall AddVar3 create, classpath(`classpath') args(str20 v7)
javacall AddVar3 create, classpath(`classpath') args(byte v8)
timer off 2
timer list 2
drop v1-v8

* Test 3: Java batch (single javacall)
di as text _n "Test 3: Java AddVarBatch (1 call for 8 vars)..."
timer clear 3
timer on 3
javacall AddVarBatch create, classpath(`classpath') args(double+double+double+long+long+str20+str20+byte v1+v2+v3+v4+v5+v6+v7+v8)
timer off 3
timer list 3
drop v1-v8

* Test 4: Sequential gen for comparison
di as text _n "Test 4: Sequential gen (8 vars)..."
timer clear 4
timer on 4
qui gen double v1 = .
qui gen double v2 = .
qui gen double v3 = .
qui gen long v4 = .
qui gen long v5 = .
qui gen str20 v6 = ""
qui gen str20 v7 = ""
qui gen byte v8 = .
timer off 4
timer list 4

di as text _n "Done!"
