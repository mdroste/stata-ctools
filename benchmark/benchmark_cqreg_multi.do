* Benchmark cqreg vs qreg - Multiple test cases
* Goal: Match standard errors across different specifications

clear all
set more off
adopath + "../build"

sysuse auto, clear

di _newline "{hline 78}"
di "TEST 1: price mpg"
di "{hline 78}"
quietly qreg price mpg
matrix qreg_V1 = e(V)
local qreg_se1_1 = sqrt(qreg_V1[1,1])
local qreg_se1_2 = sqrt(qreg_V1[2,2])

quietly cqreg price mpg
matrix cqreg_V1 = e(V)
local cqreg_se1_1 = sqrt(cqreg_V1[1,1])
local cqreg_se1_2 = sqrt(cqreg_V1[2,2])

di "mpg SE:   qreg=" %10.4f `qreg_se1_1' "  cqreg=" %10.4f `cqreg_se1_1' "  diff=" %10.4f `qreg_se1_1'-`cqreg_se1_1' "  ratio=" %8.4f `cqreg_se1_1'/`qreg_se1_1'
di "_cons SE: qreg=" %10.4f `qreg_se1_2' "  cqreg=" %10.4f `cqreg_se1_2' "  diff=" %10.4f `qreg_se1_2'-`cqreg_se1_2' "  ratio=" %8.4f `cqreg_se1_2'/`qreg_se1_2'

di _newline "{hline 78}"
di "TEST 2: price weight"
di "{hline 78}"
quietly qreg price weight
matrix qreg_V2 = e(V)
local qreg_se2_1 = sqrt(qreg_V2[1,1])
local qreg_se2_2 = sqrt(qreg_V2[2,2])

quietly cqreg price weight
matrix cqreg_V2 = e(V)
local cqreg_se2_1 = sqrt(cqreg_V2[1,1])
local cqreg_se2_2 = sqrt(cqreg_V2[2,2])

di "weight SE: qreg=" %10.4f `qreg_se2_1' "  cqreg=" %10.4f `cqreg_se2_1' "  diff=" %10.4f `qreg_se2_1'-`cqreg_se2_1' "  ratio=" %8.4f `cqreg_se2_1'/`qreg_se2_1'
di "_cons SE:  qreg=" %10.4f `qreg_se2_2' "  cqreg=" %10.4f `cqreg_se2_2' "  diff=" %10.4f `qreg_se2_2'-`cqreg_se2_2' "  ratio=" %8.4f `cqreg_se2_2'/`qreg_se2_2'

di _newline "{hline 78}"
di "TEST 3: price mpg weight"
di "{hline 78}"
quietly qreg price mpg weight
matrix qreg_V3 = e(V)
local qreg_se3_1 = sqrt(qreg_V3[1,1])
local qreg_se3_2 = sqrt(qreg_V3[2,2])
local qreg_se3_3 = sqrt(qreg_V3[3,3])

quietly cqreg price mpg weight
matrix cqreg_V3 = e(V)
local cqreg_se3_1 = sqrt(cqreg_V3[1,1])
local cqreg_se3_2 = sqrt(cqreg_V3[2,2])
local cqreg_se3_3 = sqrt(cqreg_V3[3,3])

di "mpg SE:    qreg=" %10.4f `qreg_se3_1' "  cqreg=" %10.4f `cqreg_se3_1' "  diff=" %10.4f `qreg_se3_1'-`cqreg_se3_1' "  ratio=" %8.4f `cqreg_se3_1'/`qreg_se3_1'
di "weight SE: qreg=" %10.4f `qreg_se3_2' "  cqreg=" %10.4f `cqreg_se3_2' "  diff=" %10.4f `qreg_se3_2'-`cqreg_se3_2' "  ratio=" %8.4f `cqreg_se3_2'/`qreg_se3_2'
di "_cons SE:  qreg=" %10.4f `qreg_se3_3' "  cqreg=" %10.4f `cqreg_se3_3' "  diff=" %10.4f `qreg_se3_3'-`cqreg_se3_3' "  ratio=" %8.4f `cqreg_se3_3'/`qreg_se3_3'

di _newline "{hline 78}"
di "TEST 4: mpg price (dependent variable swap)"
di "{hline 78}"
quietly qreg mpg price
matrix qreg_V4 = e(V)
local qreg_se4_1 = sqrt(qreg_V4[1,1])
local qreg_se4_2 = sqrt(qreg_V4[2,2])

quietly cqreg mpg price
matrix cqreg_V4 = e(V)
local cqreg_se4_1 = sqrt(cqreg_V4[1,1])
local cqreg_se4_2 = sqrt(cqreg_V4[2,2])

di "price SE: qreg=" %10.4f `qreg_se4_1' "  cqreg=" %10.4f `cqreg_se4_1' "  diff=" %10.4f `qreg_se4_1'-`cqreg_se4_1' "  ratio=" %8.4f `cqreg_se4_1'/`qreg_se4_1'
di "_cons SE: qreg=" %10.4f `qreg_se4_2' "  cqreg=" %10.4f `cqreg_se4_2' "  diff=" %10.4f `qreg_se4_2'-`cqreg_se4_2' "  ratio=" %8.4f `cqreg_se4_2'/`qreg_se4_2'

di _newline "{hline 78}"
di "SUMMARY"
di "{hline 78}"
