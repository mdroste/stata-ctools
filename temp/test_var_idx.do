capture do "validation/validate_setup.do"
* Test what index the plugin actually sees
sysuse census, clear
di "state is at position 1 in allvars"
di "state2 is at position 2 in allvars"

* Plugin call passes only 2 vars: state2 and test
* So inside plugin: state2 = var 1, test = var 2
quietly gen long test = .
capture noisily plugin call ctools_plugin state2 test, "cencode 1 2 label=test"
di _rc

* Compare: using wrong indices (2 and 13 or similar)
drop test
quietly gen long test = .
capture noisily plugin call ctools_plugin state2 test, "cencode 2 13 label=test"
di _rc
