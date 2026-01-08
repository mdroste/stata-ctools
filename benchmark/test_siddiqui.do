clear all
set more off
sysuse auto, clear

local q = 0.5
local h = 0.23141524

* Fit quantile regression at q-h, q, and q+h
local q_lo = `q' - `h'
local q_hi = `q' + `h'

display "q = " `q'
display "h = " `h'
display "q_lo = " `q_lo'
display "q_hi = " `q_hi'

* Get mean of X
quietly sum mpg
local x_bar = r(mean)
display "X_bar (mean mpg) = " `x_bar'

* Fit at q_lo
quietly qreg price mpg, quantile(`q_lo')
local b0_lo = _b[_cons]
local b1_lo = _b[mpg]
local Qhat_lo = `b0_lo' + `b1_lo' * `x_bar'
display ""
display "=== Quantile regression at q-h = " `q_lo' " ==="
display "b[_cons] = " `b0_lo'
display "b[mpg] = " `b1_lo'
display "Q_hat(q-h | X_bar) = " `Qhat_lo'

* Fit at q_hi
quietly qreg price mpg, quantile(`q_hi')
local b0_hi = _b[_cons]
local b1_hi = _b[mpg]
local Qhat_hi = `b0_hi' + `b1_hi' * `x_bar'
display ""
display "=== Quantile regression at q+h = " `q_hi' " ==="
display "b[_cons] = " `b0_hi'
display "b[mpg] = " `b1_hi'
display "Q_hat(q+h | X_bar) = " `Qhat_hi'

* Compute Siddiqui sparsity
local s_siddiqui = (`Qhat_hi' - `Qhat_lo') / (2 * `h')
display ""
display "=== Siddiqui (fitted) sparsity estimate ==="
display "Q_hat(q+h) - Q_hat(q-h) = " (`Qhat_hi' - `Qhat_lo')
display "Sparsity = [Q(q+h) - Q(q-h)] / (2h) = " `s_siddiqui'

* Compare to qreg fitted method
display ""
display "=== Stata qreg IID fitted ==="
quietly qreg price mpg
display "e(sparsity) = " e(sparsity)
display ""
display "Ratio: Siddiqui / qreg = " `s_siddiqui' / e(sparsity)
