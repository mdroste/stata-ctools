/*******************************************************************************
 * debug_ivstats.do
 *
 * Compare underidentification and weak ID test statistics between
 * ivreghdfe and civreghdfe
 ******************************************************************************/

clear all
set more off

* Add build directory to adopath
adopath ++ "build"

* Load auto dataset
sysuse auto, clear

di as text ""
di as text "=========================================="
di as text "Test 1: Single endogenous, single instrument"
di as text "=========================================="

* Run ivreghdfe
ivreghdfe price (mpg = weight), absorb(foreign)
di as text ""
di as text "ivreghdfe results:"
di as text "  N = " e(N)
di as text "  idstat (underid) = " e(idstat)
di as text "  iddf = " e(iddf)
di as text "  idp = " e(idp)
di as text "  widstat (Cragg-Donald F) = " e(widstat)
di as text "  First-stage F = " e(rkf)

local iv_idstat = e(idstat)
local iv_widstat = e(widstat)
local iv_rkf = e(rkf)

* Run civreghdfe
civreghdfe price (mpg = weight), absorb(foreign)
di as text ""
di as text "civreghdfe results:"
di as text "  N = " e(N)
di as text "  idstat (underid) = " e(idstat)
di as text "  iddf = " e(iddf)
di as text "  idp = " e(idp)
di as text "  cd_f (Cragg-Donald F) = " e(cd_f)
di as text "  First-stage F = " e(F_first1)

local civ_idstat = e(idstat)
local civ_cd_f = e(cd_f)
local civ_F1 = e(F_first1)

di as text ""
di as text "COMPARISON:"
di as text "  underid stat: ivreghdfe=" `iv_idstat' " civreghdfe=" `civ_idstat' " diff=" (`iv_idstat' - `civ_idstat')
di as text "  weak ID F:    ivreghdfe=" `iv_widstat' " civreghdfe=" `civ_cd_f' " diff=" (`iv_widstat' - `civ_cd_f')

di as text ""
di as text "=========================================="
di as text "Test 2: Single endogenous, two instruments"
di as text "=========================================="

sysuse auto, clear
ivreghdfe price (mpg = weight length), absorb(foreign)
di as text ""
di as text "ivreghdfe results:"
di as text "  idstat = " e(idstat)
di as text "  widstat = " e(widstat)
di as text "  rkf = " e(rkf)

local iv_idstat = e(idstat)
local iv_widstat = e(widstat)

civreghdfe price (mpg = weight length), absorb(foreign)
di as text ""
di as text "civreghdfe results:"
di as text "  idstat = " e(idstat)
di as text "  cd_f = " e(cd_f)
di as text "  F_first1 = " e(F_first1)

local civ_idstat = e(idstat)
local civ_cd_f = e(cd_f)

di as text ""
di as text "COMPARISON:"
di as text "  underid stat: ivreghdfe=" `iv_idstat' " civreghdfe=" `civ_idstat' " diff=" (`iv_idstat' - `civ_idstat')
di as text "  weak ID F:    ivreghdfe=" `iv_widstat' " civreghdfe=" `civ_cd_f' " diff=" (`iv_widstat' - `civ_cd_f')

di as text ""
di as text "=========================================="
di as text "Test 3: With robust VCE"
di as text "=========================================="

sysuse auto, clear
ivreghdfe price (mpg = weight), absorb(foreign) vce(robust)
di as text ""
di as text "ivreghdfe results (robust):"
di as text "  idstat = " e(idstat)
di as text "  widstat = " e(widstat)

local iv_idstat = e(idstat)
local iv_widstat = e(widstat)

civreghdfe price (mpg = weight), absorb(foreign) vce(robust)
di as text ""
di as text "civreghdfe results (robust):"
di as text "  idstat = " e(idstat)
di as text "  cd_f = " e(cd_f)
di as text "  widstat = " e(widstat)
di as text "  kp_f = " e(kp_f)

local civ_idstat = e(idstat)
local civ_cd_f = e(cd_f)
local civ_widstat = e(widstat)

di as text ""
di as text "COMPARISON:"
di as text "  underid stat: ivreghdfe=" `iv_idstat' " civreghdfe=" `civ_idstat' " diff=" (`iv_idstat' - `civ_idstat')
di as text "  weak ID F (KP): ivreghdfe=" `iv_widstat' " civreghdfe=" `civ_widstat' " diff=" (`iv_widstat' - `civ_widstat')
di as text "  weak ID F (CD): ivreghdfe CD=" 106.9201 " civreghdfe=" `civ_cd_f'

di as text ""
di as text "=========================================="
di as text "Test 4: Two endogenous variables"
di as text "=========================================="

sysuse auto, clear
ivreghdfe price (mpg weight = length turn displacement), absorb(foreign)
di as text ""
di as text "ivreghdfe results (2 endog):"
di as text "  idstat = " e(idstat)
di as text "  widstat = " e(widstat)

local iv_idstat = e(idstat)
local iv_widstat = e(widstat)

civreghdfe price (mpg weight = length turn displacement), absorb(foreign)
di as text ""
di as text "civreghdfe results (2 endog):"
di as text "  idstat = " e(idstat)
di as text "  cd_f = " e(cd_f)

local civ_idstat = e(idstat)
local civ_cd_f = e(cd_f)

di as text ""
di as text "COMPARISON:"
di as text "  underid stat: ivreghdfe=" `iv_idstat' " civreghdfe=" `civ_idstat' " diff=" (`iv_idstat' - `civ_idstat')
di as text "  weak ID F:    ivreghdfe=" `iv_widstat' " civreghdfe=" `civ_cd_f' " diff=" (`iv_widstat' - `civ_cd_f')

di as text ""
di as text "=========================================="
di as text "All e() scalars from ivreghdfe"
di as text "=========================================="

sysuse auto, clear
ivreghdfe price (mpg = weight), absorb(foreign)
ereturn list
