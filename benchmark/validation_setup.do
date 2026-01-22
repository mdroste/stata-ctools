*! Common setup for all validation tests
*! This file sets up the adopath to find ctools commands

* Add the build directory to Stata's search path
adopath ++ "build"

* Verify ctools commands are available
capture which csort
if _rc != 0 {
    di as error "ERROR: Could not find csort.ado in adopath"
    di as error "Make sure you run from the project root directory"
    di as error "and that 'build' directory contains the .ado files"
    exit 601
}

di as text "ctools adopath setup complete"
di as text "  adopath includes: build"
di ""
