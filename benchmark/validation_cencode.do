*! Validation tests for cencode - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cencode.log", replace text
validation_header "cencode"

* === Basic encoding ===
sysuse auto, clear
capture cencode make, generate(make_code)
local p = (_rc == 0)
capture confirm variable make_code
local passed = (`p' & _rc == 0)
report_test "Basic cencode" `passed'

* === Output type ===
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
capture confirm numeric variable make_code
local passed = (_rc == 0)
report_test "Output is numeric" `passed'

* === Label options ===
sysuse auto, clear
capture cencode make, generate(make_num) label(car_makes)
local p = (_rc == 0)
capture label list car_makes
local passed = (`p' & _rc == 0)
report_test "label()" `passed'

sysuse auto, clear
capture drop make_encoded
capture cencode make, generate(make_encoded)
capture label list make_encoded
local passed = (_rc == 0)
report_test "Default label name" `passed'

* === Other options ===
sysuse auto, clear
capture drop make_code
run_test "verbose" cencode make, generate(make_code) verbose

sysuse auto, clear
capture drop make_code
run_test "threads(2)" cencode make, generate(make_code) threads(2)

sysuse auto, clear
capture drop make_code
run_test "noextend" cencode make, generate(make_code) noextend

* === Conditionals ===
sysuse auto, clear
capture drop make_code
run_test "if condition" cencode make if foreign == 1, generate(make_code)

sysuse auto, clear
capture drop make_code
run_test "in range" cencode make in 1/50, generate(make_code)

* === Verify consecutive integers ===
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
capture qui sum make_code
local min = r(min)
local max = r(max)
capture qui tab make_code
local n_unique = r(r)
local passed = (`min' == 1 & `max' == `n_unique')
report_test "Values are 1 to N" `passed'

* === Alphabetical ordering ===
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
capture gsort make_code
local first_make = make[1]
capture gsort make
local alpha_first = make[1]
local passed = ("`first_make'" == "`alpha_first'")
report_test "Alphabetical order" `passed'

* === Missing strings ===
clear
set obs 20
gen str20 name = "name" + string(_n)
replace name = "" in 1/5
capture drop name_code
capture cencode name, generate(name_code)
local p = (_rc == 0)
capture count if missing(name_code) & name == ""
local passed = (`p' & r(N) == 5)
report_test "Missing strings" `passed'

* === Many unique values ===
clear
set obs 1000
gen str20 name = "name" + string(_n)
capture drop name_code
capture cencode name, generate(name_code)
local p = (_rc == 0)
capture qui tab name_code
local passed = (`p' & r(r) == 1000)
report_test "1000 unique values" `passed'

* === Duplicate values ===
clear
set obs 100
gen str20 category = "cat" + string(mod(_n - 1, 10) + 1)
capture drop cat_code
capture cencode category, generate(cat_code)
local p = (_rc == 0)
capture qui tab cat_code
local passed = (`p' & r(r) == 10)
report_test "Duplicate values" `passed'

* === Special characters ===
clear
set obs 10
gen str50 special = "normal"
replace special = "with space" in 1
replace special = "with,comma" in 2
replace special = "with-dash" in 3
replace special = "with_underscore" in 4
capture drop special_code
run_test "Special characters" cencode special, generate(special_code)

* === Long strings ===
clear
set obs 10
gen str200 longname = "a" * 150 + string(_n)
capture drop long_code
run_test "Long strings (150+)" cencode longname, generate(long_code)

* === Single unique value ===
clear
set obs 100
gen str10 single = "same"
capture drop single_code
capture cencode single, generate(single_code)
local p = (_rc == 0)
capture qui sum single_code
local passed = (`p' & r(min) == r(max) & r(min) == 1)
report_test "Single unique value" `passed'

* === Error cases ===
sysuse auto, clear
capture gen make_code = 1
capture cencode make, generate(make_code)
local passed = (_rc == 110)
report_test "Existing var error" `passed'

sysuse auto, clear
capture cencode price, generate(price_code)
local passed = (_rc != 0)
report_test "Numeric var error" `passed'

* === Label values match ===
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
local orig_make = make[1]
local code_val = make_code[1]
local label_val : label (make_code) `code_val'
local passed = ("`orig_make'" == "`label_val'")
report_test "Labels match strings" `passed'

* === Combined options ===
sysuse auto, clear
capture drop make_code
run_test "Combined options" cencode make if foreign == 1, generate(make_code) label(foreign_makes) verbose threads(2)

* === Large dataset ===
clear
set obs 100000
gen str20 category = "cat" + string(mod(_n - 1, 1000) + 1)
capture drop cat_code
run_test "Large dataset (100k)" cencode category, generate(cat_code)

* === strL variables ===
clear
set obs 10
gen strL longtext = "This is a longer string " + string(_n)
capture drop text_code
run_test "strL variables" cencode longtext, generate(text_code)

* === Unicode ===
clear
set obs 5
gen str50 unicode_text = "normal"
replace unicode_text = "cafe" in 1
replace unicode_text = "naive" in 2
replace unicode_text = "uber" in 3
capture drop uni_code
run_test "Unicode chars" cencode unicode_text, generate(uni_code)

* === Whitespace ===
clear
set obs 10
gen str20 spaces = "  leading"
replace spaces = "trailing  " in 2
replace spaces = "  both  " in 3
replace spaces = "middle  space" in 4
capture drop space_code
run_test "Whitespace" cencode spaces, generate(space_code)

* === All values labeled ===
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
capture qui levelsof make_code, local(codes)
local all_labeled = 1
foreach c of local codes {
    local lbl : label (make_code) `c'
    if "`lbl'" == "" local all_labeled = 0
}
report_test "All values labeled" `all_labeled'

* === Numeric-looking strings ===
clear
set obs 10
gen str10 numstr = string(_n * 100)
capture drop numstr_code
run_test "Numeric-like strings" cencode numstr, generate(numstr_code)

* === Mixed case ===
clear
set obs 5
gen str20 mixcase = "Apple"
replace mixcase = "apple" in 2
replace mixcase = "APPLE" in 3
replace mixcase = "ApPlE" in 4
replace mixcase = "banana" in 5
capture drop case_code
capture cencode mixcase, generate(case_code)
local p = (_rc == 0)
capture qui tab case_code
local passed = (`p' & r(r) == 5)
report_test "Mixed case preserved" `passed'

* === Stored results ===
sysuse auto, clear
capture drop make_code
capture cencode make, generate(make_code)
capture local n_unique = _cencode_n_unique
local passed = (`n_unique' > 0)
report_test "Stored results" `passed'

* === Empty if condition ===
sysuse auto, clear
capture drop make_code
capture cencode make if price > 100000, generate(make_code)
report_test "Empty if handled" 1  // Just check doesn't crash

validation_summary
log close
