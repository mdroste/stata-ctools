gettoken __case 0 : 0
if "`__case'" == "" local __case chdfe
if "`__case'" == "chdfe" {
do "validation/support/validate_setup.do"
quietly {
clear
set seed 74629
set obs 3000
gen long order=_n
gen int a=ceil(runiform()*50)
gen int b=ceil(runiform()*30)
gen double t=rnormal()
gen double u=rnormal()
gen double y=2*a-b+t+rnormal()
gen double x=a+2*b+rnormal()
gen double w=1+mod(_n,4)
gen str8 s="a"+string(a)
replace y=. if mod(_n,53)==0
replace a=. if mod(_n,59)==0
replace a=999 if _n==7
foreach absorb in "a" "a b" "a#b" "a##c.t b" "a#c.t b" "a##c.(t u) b" {
    foreach wt in "" "[aw=w]" "[fw=w]" "[pw=w]" {
        * The hdfe 3 default CG can stagnate above its requested tolerance.
        * Plain projections agree with explicit dummy-variable OLS here.
        hdfe y x if order>5 `wt', absorb(`absorb') gen(ref_) sample(refsample) tol(1e-12) accel(none)
        local df=e(df_a)
        local nh=e(N_hdfe)
        chdfe y x if order>5 `wt', absorb(`absorb') gen(got_) sample(gotsample) tol(1e-12)
        assert e(df_a)==`df' & e(N_hdfe)==`nh'
        assert refsample==gotsample
        assert_var_equal got_y ref_y $DEFAULT_SIGFIGS "hdfe `absorb' `wt' y"
        assert_var_equal got_x ref_x $DEFAULT_SIGFIGS "hdfe `absorb' `wt' x"
        assert order==_n
        drop got_* ref_* refsample gotsample
    }
}
test_pass "weighted projections, interactions, slopes, singleton samples, and degrees of freedom"
foreach dof in "all" "none" "firstpair" "pairwise" "clusters" "continuous" {
    * hdfe 3.2.9 parses firstpair but its Mata routine expects firstpairs.
    * With two effects, pairwise is the same adjustment and is a usable oracle.
    local refdof = cond("`dof'"=="firstpair","pairwise","`dof'")
    hdfe y x, absorb(a b) gen(ref_) clustervars(a) dof(`refdof') keepsingletons tol(1e-12) accel(none)
    local df=e(df_a)
    local m=e(mobility)
    chdfe y x, absorb(a b) gen(got_) clustervars(a) dof(`dof') keepsingletons tol(1e-12)
    assert e(df_a)==`df' & e(mobility)==`m'
    assert_var_equal got_y ref_y $DEFAULT_SIGFIGS "cluster nesting and dof(`dof')"
    drop got_* ref_*
}
test_pass "cluster nesting and degrees-of-freedom options"
chdfe y x, absorb(s b) gen(ref_) tol(1e-12) keepsingletons
foreach trans in sym kac cim {
    foreach accel in none cg sd aitken {
        if "`trans'"=="kac" & "`accel'"=="cg" continue
        chdfe y x, absorb(s b) gen(got_) tol(1e-12) keepsingletons ///
            transform(`trans') acceleration(`accel') threads(4) pool(1) verbose
        assert_var_equal got_y ref_y $DEFAULT_SIGFIGS "transform `trans' acceleration `accel'"
        drop got_*
    }
}
drop ref_*
test_pass "solver choices, strings, pool size and threads"
chdfe y x, absorb(a b) gen(ref_) sample(refsample) groupvar(mobility) tol(1e-12)
* Clear retains only requested variables, so compare through saved data below.
preserve
keep if refsample
keep order ref_y ref_x
tempfile expected
save `expected'
restore
preserve
chdfe y x, absorb(a b) clear keepvars(order) keepids tol(1e-12)
confirm variable __ID1__ __ID2__
merge 1:1 order using `expected', assert(match) nogen
assert_var_equal y ref_y $DEFAULT_SIGFIGS "clear sample and values"
restore
drop ref_* refsample mobility
* Failure must leave all user variables and observations unchanged.
capture chdfe y x, absorb(a b) clear maxit(1) tol(1e-16)
assert _rc==430
assert order==_n
confirm variable y x a b t u w s
gen double xy = 1
capture chdfe y, absorb(a) gen(x)
assert _rc==110
assert xy==1
drop xy
capture chdfe y if 0, absorb(a) gen(z_)
assert _rc==2000
capture confirm variable z_y
assert _rc==111
test_pass "clear output, saved IDs, empty sample and rollback"
* Factor variables and time-series operators have the reference output names.
replace y=1 if missing(y)
tsset order
hdfe y i.b L.x, absorb(a) gen(ref_) keepsingletons tol(1e-12)
chdfe y i.b L.x, absorb(a) gen(got_) keepsingletons tol(1e-12)
unab refvars : ref_*
foreach v of local refvars {
    local dest : subinstr local v "ref_" "got_"
    assert_var_equal `v' `dest' $DEFAULT_SIGFIGS "factor and time-series expansion"
}
test_pass "factor and time-series varlist"
drop ref_* got_*
chdfe y x, absorb(a b) gen(base_) tol(1e-12)
foreach diagnostic in "verbose(1)" "v(0)" "verbosity(1)" "timeit" {
    chdfe y x, absorb(a b) gen(test_) tol(1e-12) `diagnostic' maxit(.)
    assert test_y==base_y & test_x==base_x
    drop test_*
}
chdfe, version
assert "`e(version)'"!=""
test_pass "diagnostics, unlimited iterations, and version"
clear
input double(a b y expected)
20 4 1 2
20 5 2 2
10 1 3 1
10 2 4 1
21 5 5 2
21 4 6 2
11 2 7 1
11 1 8 1
end
hdfe y, absorb(a b) gen(ref_) groupvar(refgroup) accel(none)
local df=e(df_a)
chdfe y, absorb(a b) gen(got_) groupvar(gotgroup)
assert e(df_a)==`df'
assert gotgroup==expected & refgroup==expected
assert_var_equal got_y ref_y $DEFAULT_SIGFIGS "disconnected mobility groups"
chdfe y, absorb(a) gen(one_) groupvar(no_pair)
capture confirm variable no_pair
assert _rc==111
test_pass "mobility group numbering and absent eligible pair"
clear
set seed 68481
set obs 1000
gen int a=mod(_n,11)
gen int b=mod(_n,13)
gen double t=rnormal()
gen double y=rnormal()+a+b
gen double w=1+runiform()
chdfe y [aw=w], absorb(a##c.t b) gen(base_) tol(1e-12)
foreach scale in 1e-200 1e200 {
    gen double z=y*`scale'
    gen double scaled_t=t*`scale'
    gen double scaled_w=w*`scale'
    chdfe z [aw=scaled_w], absorb(a##c.scaled_t b) gen(out_) tol(1e-12)
    replace out_z=out_z/`scale'
    assert_var_equal out_z base_y $DEFAULT_SIGFIGS "outcome, slope and weight scaling"
    drop z scaled_t scaled_w out_z
}
test_pass "extreme finite scales"
* Clear keeps source variables used by expanded terms and simple weights.
gen long id=_n
tsset id
chdfe y L.t i.b [aw=w], absorb(a) clear keepids
confirm variable id t b w __ID1__
confirm variable __L__t
assert _N==999
test_pass "clear factor/time sources and weight retention"

}
noisily print_summary "chdfe"
global CTOOLS_CASE_COMPLETE "chdfe"
}
