*===============================================================================
* FILE: speed benchmark
*===============================================================================


timer clear
clear all

sysuse auto
csort price

local N 25000000
local G1 10000
local G2 10000000

local N 250000
local G1 10000
local G2 10000


*-------------------------------------------------------------------------------
* Create test datasets
*-------------------------------------------------------------------------------

* Main 25M obs test dataset
clear
set obs `N'
forval i=1/4 {
	gen float`i' = rnormal()
	gen int`i' = ceil(runiform()*`G1')
	gen bigint`i' = ceil(runiform()*`G2')
}
gen string1 = char(runiformint(65,90)) +  string(runiformint(0,9)) + char(runiformint(65,90)) + char(runiformint(65,90)) + string(runiformint(0,9)) + char(runiformint(65,90))
gen string2 = char(runiformint(65,90)) +  string(runiformint(0,9)) + char(runiformint(65,90)) + char(runiformint(65,90)) + string(runiformint(0,9)) + char(runiformint(65,90))
gen y = 1 + float1 + 2*float2 + 3*float3 + sin(int1) + cos(int2) + rnormal()
save temp, replace

* Auxiliary dataset for merging on big ints
clear
set obs `G2'
gen bigint1 = _n
gen z1 = rnormal()
save temp2, replace

* Auxiliary dataset for merging on small ints
clear
set obs `G1'
gen int1 = _n
gen z2 = rnormal()
save temp3, replace


*-------------------------------------------------------------------------------
* Sort benchmarks
*-------------------------------------------------------------------------------

use temp, clear

* Sort on float
timer on 1
sort float1
timer off 1
timer on 2
csort float2
timer off 2

* Sort on small int1
timer on 3
sort int1
timer off 3
timer on 4
csort int2
timer off 4

* Sort on big int1
timer on 5
sort bigint1
timer off 5
timer on 6
csort bigint2
timer off 6

* Sort on string
timer on 7
sort string1
timer off 7
timer on 8
csort string2
timer off 8

*-------------------------------------------------------------------------------
* Binscatter benchmarks
*-------------------------------------------------------------------------------

* Simple binscatter
timer on 9
binscatter y float1
timer off 9
timer on 10
cbinscatter y float1
timer off 10

* Absorbing binscatter
timer on 11
binscatter y float1, controls(float2 float3) absorb(int1)
timer off 11
timer on 12
cbinscatter y float1, controls(float2 float3) absorb(int1)
timer off 12


*-------------------------------------------------------------------------------
* reghdfe benchmarks
*-------------------------------------------------------------------------------

* basic reghdfe
timer on 13
reghdfe y float1 float2 float3, absorb(int1 int2)
timer off 13
timer on 14
creghdfe y float1 float2 float3, absorb(int1 int2)
timer off 14

* basic ivreghdfe
timer on 15
ivreghdfe y (float1=float4) float2 float3, absorb(int1 int2)
timer off 15
timer on 16
civreghdfe y (float1=float4) float2 float3, absorb(int1 int2)
timer off 16

*-------------------------------------------------------------------------------
* merge benchmarks
*-------------------------------------------------------------------------------

* m:1 big int merge
use temp, clear
timer on 17
merge m:1 bigint1 using temp2
timer off 17
use temp, clear
timer on 18
cmerge m:1 bigint1 using temp2
timer off 18

timer list

