clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/screen_candidate.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/candidate"
foreach categories in 10 1000 10000 {
    clear
    set obs 200000
    gen str12 x = "value_" + string(mod(_n-1,`categories'), "%05.0f")
    cencode x, gen(y) verbose threads(8)
    encode x, gen(z)
    assert y==z
    decode y, gen(ytext)
    assert ytext==x
}
clear
set obs 1000000
set seed 81523
gen double x=rnormal()
gen long g=mod(_n,5000)
cwinsor x, by(g) suffix(_w) verbose threads(8)
local rc=0
log close
exit, clear
