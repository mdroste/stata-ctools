*! Private helper for cmerge; maintained in package/ado.
program define _cmerge_label_mergevar
    args mergevar
    capture label list _merge
    if _rc {
        label define _merge 1 "Master only (1)" 2 "Using only (2)" ///
            3 "Matched (3)" 4 "Missing updated (4)" 5 "Nonmissing conflict (5)"
    }
    label values `mergevar' _merge
    label variable `mergevar' "Matching result from merge"
end
