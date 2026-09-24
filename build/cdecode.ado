*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

program define cdecode, rclass
    version 14.1
    preserve
    capture noisily _cdecode_impl `0'
    local rc = _rc
    if `rc' {
        restore
        exit `rc'
    }
    return add
    restore, not
end

program define _cdecode_impl, rclass
    version 14.1
    syntax varlist [if] [in], [Generate(string) replace MAXLength(integer 0) Verbose THReads(integer 0)]

    if ("`generate'" == "") == ("`replace'" == "") {
        di as error "cdecode: specify either generate() or replace"
        exit 198
    }
    if `maxlength' < 0 | `threads' < 0 {
        di as error "cdecode: maxlength() and threads() must be nonnegative"
        exit 198
    }
    local nvars : word count `varlist'
    if "`generate'" != "" {
        local nout : word count `generate'
        local unique : list uniq generate
        local nunique : word count `unique'
        if `nout' != `nvars' | `nunique' != `nout' {
            di as error "cdecode: generate() requires one distinct name per input variable"
            exit 198
        }
        confirm new variable `generate'
    }
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc exit 108
        local label : value label `v'
        if "`label'" == "" {
            di as error "cdecode: `v' has no value label attached"
            exit 182
        }
    }

    * Keep labeled missing values in the sample. Native decode retrieves literal
    * label text directly, including strL labels, without parsing label-save code.
    marksample touse, novarlist
    local maxopt ""
    if `maxlength' > 0 local maxopt "maxlength(`maxlength')"
    local staged ""
    foreach v of local varlist {
        tempvar decoded
        quietly decode `v' if `touse', generate(`decoded') `maxopt'
        local staged `staged' `decoded'
    }

    local i = 0
    foreach v of local varlist {
        local ++i
        local decoded : word `i' of `staged'
        if "`replace'" != "" {
            order `decoded', before(`v')
            drop `v'
            rename `decoded' `v'
        }
        else {
            local target : word `i' of `generate'
            rename `decoded' `target'
        }
    }
    if "`verbose'" != "" {
        di as text "cdecode: decoded `nvars' variable(s) with Stata's value-label engine"
        di as text "threads() is retained for syntax compatibility; decoding is managed by Stata"
    }
    return scalar N_vars = `nvars'
end
