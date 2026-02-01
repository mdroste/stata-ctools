* Test replace option and varlist support for cencode and cdecode
* These are ctools-specific extensions

clear all
set more off

* Load ctools
adopath ++ build

* =============================================================================
* Test cencode replace option (single variable)
* =============================================================================

di _n "Testing cencode replace option (single variable)..."

* Create test data
clear
set obs 100
gen str20 category = ""
replace category = "Apple" if _n <= 25
replace category = "Banana" if _n > 25 & _n <= 50
replace category = "Cherry" if _n > 50 & _n <= 75
replace category = "Date" if _n > 75

* Verify it's a string variable
confirm string variable category

* Use replace option
cencode category, replace

* Verify it's now numeric with value labels
confirm numeric variable category
local lblname : value label category
assert "`lblname'" == "category"

* Check values are correct
tab category
assert category[1] == 1  // Apple
assert category[26] == 2 // Banana
assert category[51] == 3 // Cherry
assert category[76] == 4 // Date

di as text "cencode replace option (single): PASSED"

* =============================================================================
* Test cdecode replace option (single variable)
* =============================================================================

di _n "Testing cdecode replace option (single variable)..."

* Now decode it back
cdecode category, replace

* Verify it's now a string
confirm string variable category

* Check values are correct
assert category[1] == "Apple"
assert category[26] == "Banana"
assert category[51] == "Cherry"
assert category[76] == "Date"

di as text "cdecode replace option (single): PASSED"

* =============================================================================
* Test cencode varlist support with generate()
* =============================================================================

di _n "Testing cencode varlist support with generate()..."

* Create test data with multiple string variables
clear
set obs 50
gen str10 fruit = ""
replace fruit = "Apple" if _n <= 25
replace fruit = "Banana" if _n > 25

gen str10 color = ""
replace color = "Red" if _n <= 25
replace color = "Yellow" if _n > 25

gen str10 size = ""
replace size = "Small" if _n <= 25
replace size = "Large" if _n > 25

* Encode all three at once
cencode fruit color size, generate(fruit_code color_code size_code)

* Verify all are numeric with correct labels
confirm numeric variable fruit_code
confirm numeric variable color_code
confirm numeric variable size_code

assert fruit_code[1] == 1  // Apple
assert fruit_code[26] == 2 // Banana
assert color_code[1] == 1  // Red
assert color_code[26] == 2 // Yellow
assert size_code[1] == 2   // Small (alphabetical: Large=1, Small=2)
assert size_code[26] == 1  // Large

di as text "cencode varlist with generate(): PASSED"

* =============================================================================
* Test cdecode varlist support with generate()
* =============================================================================

di _n "Testing cdecode varlist support with generate()..."

* Decode all three at once
cdecode fruit_code color_code size_code, generate(fruit_str color_str size_str)

* Verify all are strings with correct values
confirm string variable fruit_str
confirm string variable color_str
confirm string variable size_str

assert fruit_str[1] == "Apple"
assert fruit_str[26] == "Banana"
assert color_str[1] == "Red"
assert color_str[26] == "Yellow"
assert size_str[1] == "Small"
assert size_str[26] == "Large"

di as text "cdecode varlist with generate(): PASSED"

* =============================================================================
* Test cencode varlist support with replace
* =============================================================================

di _n "Testing cencode varlist support with replace..."

* Use the original string variables (fruit, color, size still exist)
confirm string variable fruit
confirm string variable color
confirm string variable size

* Replace all three at once
cencode fruit color size, replace

* Verify all are now numeric
confirm numeric variable fruit
confirm numeric variable color
confirm numeric variable size

assert fruit[1] == 1  // Apple
assert fruit[26] == 2 // Banana

di as text "cencode varlist with replace: PASSED"

* =============================================================================
* Test cdecode varlist support with replace
* =============================================================================

di _n "Testing cdecode varlist support with replace..."

* Replace back to strings
cdecode fruit color size, replace

* Verify all are now strings
confirm string variable fruit
confirm string variable color
confirm string variable size

assert fruit[1] == "Apple"
assert fruit[26] == "Banana"

di as text "cdecode varlist with replace: PASSED"

* =============================================================================
* Test error handling
* =============================================================================

di _n "Testing error handling..."

* Create fresh test data
clear
set obs 10
gen str10 var1 = "A"
gen str10 var2 = "B"

* Try to use both generate and replace (should fail)
cap cencode var1, generate(var1_code) replace
assert _rc == 198

* Try with no option (should fail)
cap cencode var1
assert _rc == 198

* Try with mismatched varlist and generate counts (should fail)
cap cencode var1 var2, generate(only_one)
assert _rc == 198

di as text "Error handling: PASSED"

* =============================================================================
* All tests passed
* =============================================================================

di _n as text "{hline 50}"
di as text "All replace option and varlist tests PASSED"
di as text "{hline 50}"
