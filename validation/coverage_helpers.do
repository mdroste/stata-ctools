/* Differential option tests. A case passes only after both commands succeed
   and their observable results agree; errors never count as successful parity. */
capture program drop coverage_fixture
program define coverage_fixture
    args variant
    clear
    set seed `=71431+`variant''
    set obs 600
    gen long id = _n
    gen byte g = mod(_n,20)
    gen byte h = mod(floor((_n-1)/20),15)
    gen byte group = mod(_n,3)
    gen double x = rnormal()
    gen double z = rnormal()
    gen double z2 = rnormal()
    gen double endog = .8*z+.5*z2+rnormal()
    gen double y = 1+.7*x+.4*endog+.03*g+rnormal()
    gen double count_y = rpoisson(exp(.3+.2*x+.1*z+.01*g))
    gen double w = .5+runiform()
    gen byte fw = 1+mod(_n,3)
    gen double ps = .05+.9*runiform()
    gen byte treat = runiform()<ps
    gen double time = mod(_n,31)
    gen double lo = time-2
    gen double hi = time+2
    gen str12 s = cond(mod(_n,5)==0,"",cond(mod(_n,2),"café","東京"))
    gen str30 tokens = cond(mod(_n,7)==0,""," 1, 2%,3x ")
    gen str30 number = string(y,"%21.12g")
    label define groups 0 "zero" 1 "one" 2 "two" 99 "unused"
    label values group groups
    if `variant' == 2 {
        replace x = exp(x/2)
        replace y = 1e6*y
        replace g = mod(_n,7)
        replace number = string(y,"%21.12g")
        replace tokens = "東京::café::" if mod(_n,7)==0
    }
    if `variant' == 3 {
        replace y = 1e-6*y
        replace y = .a if mod(_n,37)==0
        replace x = . if mod(_n,53)==0
        replace s = " quote " + char(34) + " x " if mod(_n,11)==0
        replace number = string(y,"%21.12g")
        replace group = .b if mod(_n,41)==0
    }
end

mata:
struct coverage_snapshot {
    string rowvector names, labels, formats
    string matrix value_labels, notes
    string matrix strings
    real matrix numbers
    real rowvector isnumeric
}
struct coverage_snapshot scalar coverage_snapshot_data() {
    struct coverage_snapshot scalar s
    real scalar j,k
    real colvector label_values
    string colvector label_text
    string scalar label_name
    string colvector note_names
    string rowvector note_vars
    s.names=J(1,st_nvar(),""); s.labels=s.names; s.formats=s.names
    s.value_labels=J(0,3,"")
    s.notes=J(0,3,"")
    note_vars="_dta"
    for(j=1;j<=st_nvar();j++) note_vars=note_vars,st_varname(j)
    for(j=1;j<=cols(note_vars);j++) {
        note_names=sort(st_dir("char",note_vars[j],"note*"),1)
        for(k=1;k<=rows(note_names);k++) s.notes=s.notes\(note_vars[j],note_names[k],st_global(note_vars[j]+"["+note_names[k]+"]"))
    }
    s.isnumeric=J(1,st_nvar(),0)
    s.numbers=J(st_nobs(),st_nvar(),.)
    s.strings=J(st_nobs(),st_nvar(),"")
    for(j=1;j<=st_nvar();j++) {
        s.names[j]=st_varname(j); s.labels[j]=st_varlabel(j); s.formats[j]=st_varformat(j)
        label_name=st_varvaluelabel(j)
        if(label_name!="") {
            st_vlload(label_name,label_values,label_text)
            for(k=1;k<=rows(label_values);k++) s.value_labels=s.value_labels\(s.names[j],strofreal(label_values[k],"%21x"),label_text[k])
        }
        s.isnumeric[j]=st_isnumvar(j)
        if(s.isnumeric[j]) s.numbers[,j]=st_data(.,j)
        else s.strings[,j]=st_sdata(.,j)
    }
    return(s)
}
real scalar coverage_isnumeric_equal(real matrix a, real matrix b) {
    real matrix scale, ok
    if(rows(a)!=rows(b) | cols(a)!=cols(b)) return(0)
    if(!rows(a) | !cols(a)) return(1)
    scale=rowmax((abs(vec(a)),abs(vec(b))))
    ok=(vec(a):==vec(b)) :| ((vec(a):<.):&(vec(b):<.):&
        (abs(vec(a)-vec(b)):<=1e-7:*scale))
    return(all(ok))
}
/* Predictions are sums and differences (xb = sum x*b, residuals = y - xbd), so
   their rounding and solver error scales with the column, not with an element
   that cancels to near zero: compare at 1e-7 of the column's largest value.
   Missingness must still match exactly. */
real scalar coverage_prediction_equal(real colvector a, real colvector b) {
    real colvector present
    real scalar scale
    if(rows(a)!=rows(b)) return(0)
    if(!rows(a)) return(1)
    present=(a:<.):&(b:<.)
    if(any((a:<.):!=(b:<.))) return(0)
    if(any((!present):&(a:!=b))) return(0)
    if(!any(present)) return(1)
    scale=max(abs(select(a,present))\abs(select(b,present)))
    return(all((!present) :| (abs(a-b):<=1e-7*scale)))
}
void coverage_compare_data(struct coverage_snapshot scalar ref) {
    struct coverage_snapshot scalar got
    real scalar same, j
    got=coverage_snapshot_data()
    same=0
    if(cols(got.names)==cols(ref.names) & rows(got.numbers)==rows(ref.numbers)) {
        same=all(got.names:==ref.names) & all(got.isnumeric:==ref.isnumeric)
        if(same) same=coverage_isnumeric_equal(ref.numbers,got.numbers) &
            (rows(ref.strings)==0 | all(ref.strings:==got.strings)) & all(ref.labels:==got.labels) & all(ref.formats:==got.formats)
        if(same) {
            same=(rows(ref.value_labels)==rows(got.value_labels))
            if(same & rows(ref.value_labels)>0) same=all(ref.value_labels:==got.value_labels)
            if(same) same=(rows(ref.notes)==rows(got.notes))
            if(same & rows(ref.notes)>0) same=all(ref.notes:==got.notes)
        }
    }
    if(!same & cols(got.names)==cols(ref.names) & rows(got.numbers)==rows(ref.numbers)) {
        printf("DATA_DETAIL values=%g strings=%g names=%g labels=%g\n", coverage_isnumeric_equal(ref.numbers,got.numbers), (rows(ref.strings)==0 | all(ref.strings:==got.strings)), all(ref.names:==got.names), all(ref.labels:==got.labels))
        for(j=1;j<=cols(ref.names);j++) {
            if(ref.formats[j]!=got.formats[j]) printf("DATA_FORMAT %s ref=%s got=%s\n",ref.names[j],ref.formats[j],got.formats[j])
            if(ref.labels[j]!=got.labels[j]) printf("DATA_LABEL %s ref=%s got=%s\n",ref.names[j],ref.labels[j],got.labels[j])
        }
        if(rows(ref.notes)!=rows(got.notes)) printf("DATA_NOTES ref=%g got=%g\n",rows(ref.notes),rows(got.notes))
        if(rows(ref.value_labels)!=rows(got.value_labels)) printf("DATA_VALUE_LABELS ref=%g got=%g\n",rows(ref.value_labels),rows(got.value_labels))
    }
    st_local("same",strofreal(same))
}
end

capture program drop coverage_data
program define coverage_data
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string) [SORTby(varlist) REFAFTER(string)]
    preserve
    tempfile original
    quietly save `original'
    capture quietly `reference'
    local nr = _rc
    if !`nr' & `"`refafter'"' != "" {
        capture quietly `refafter'
        local nr = _rc
    }
    if !`nr' {
        if "`sortby'" != "" sort `sortby', stable
        mata: coverage_ref=coverage_snapshot_data()
    }
    quietly use `original', clear
    capture quietly `candidate'
    local cr = _rc
    local same = 0
    if !`cr' & !`nr' {
        if "`sortby'" != "" sort `sortby', stable
        mata: coverage_compare_data(coverage_ref)
    }
    restore
    if `nr' | `cr' | !`same' {
        coverage_failure "`name'" "reference rc=`nr'; candidate rc=`cr'; data equal=`same'"
    }
    else test_pass "`name'"
end

capture program drop coverage_estimate
program define coverage_estimate
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string) ///
        [PREDict(string) SCAlars(string) REQUIRED(string)]
    preserve
    tempname rb rV cb cV
    tempfile original
    quietly save `original'
    local before = $TESTS_FAILED
    capture quietly `reference'
    local nr = _rc
    if `nr' {
        restore
        coverage_failure "`name'" "reference rc=`nr'"
        exit
    }
    matrix `rb'=e(b)
    matrix `rV'=e(V)
    local names : colnames `rb'
    local vrows : rownames `rV'
    local vcols : colnames `rV'
    local reference_N = e(N)
    tempvar esample
    gen byte `esample'=e(sample)
    mata: coverage_sample=st_data(.,st_local("esample"))
    local scalar_items `scalars' `required'
    local scalar_index=0
    foreach item of local scalar_items {
        local ++scalar_index
        local mapped : subinstr local item ":" " ", all
        local refkey : word 1 of `mapped'
        local ref_scalar_`scalar_index'=e(`refkey')
    }
    foreach item of local predict {
        tempvar p
        capture quietly predict double `p', `item'
        local prc_`item' = _rc
        if !`prc_`item'' mata: coverage_pred_`item'=st_data(.,st_local("p"))
    }
    quietly use `original', clear
    capture quietly `candidate'
    local cr = _rc
    if `cr' {
        restore
        coverage_failure "`name'" "candidate rc=`cr'"
        exit
    }
    matrix `cb'=e(b)
    matrix `cV'=e(V)
    local cnames : colnames `cb'
    local cvrows : rownames `cV'
    local cvcols : colnames `cV'
    if "`vrows'" != "`cvrows'" | "`vcols'" != "`cvcols'" coverage_failure "`name' V names" "covariance names/order differ"
    if "`names'" != "`cnames'" coverage_failure "`name' names" "coefficient names/order differ"
    if `reference_N' != e(N) coverage_failure "`name' N" "e(N) differs"
    gen byte `esample'=e(sample)
    mata: st_local("same",strofreal(all(coverage_sample:==st_data(.,st_local("esample")))))
    if !`same' coverage_failure "`name' sample" "e(sample) differs"
    assert_matrix_equal `rb' `cb' $DEFAULT_SIGFIGS "`name' b"
    assert_matrix_equal `rV' `cV' $DEFAULT_SIGFIGS "`name' V"
    local scalar_index=0
    foreach item of local scalar_items {
        local ++scalar_index
        local mapped : subinstr local item ":" " ", all
        local refkey : word 1 of `mapped'
        local candkey : word 2 of `mapped'
        if "`candkey'"=="" local candkey `refkey'
        local expected=`ref_scalar_`scalar_index''
        local actual=e(`candkey')
        local must_exist : list item in required
        if missing(`expected') | missing(`actual') {
            if `must_exist' | `expected'!=`actual' coverage_failure "`name' `item'" "required statistic missing or missingness differs"
        }
        else {
            sigfigs `expected' `actual'
            if r(sigfigs)<$DEFAULT_SIGFIGS coverage_failure "`name' `item'" "less than 7 significant digits"
        }
    }
    foreach item of local predict {
        tempvar p
        capture quietly predict double `p', `item'
        local crc = _rc
        if `prc_`item'' | `crc' coverage_failure "`name' predict `item'" "prediction failed"
        else {
            mata: st_local("same",strofreal(coverage_prediction_equal(coverage_pred_`item',st_data(.,st_local("p")))))
            if !`same' coverage_failure "`name' predict `item'" "prediction values/missingness differ"
        }
    }
    restore
    if $TESTS_FAILED == `before' test_pass "`name' results and predictions"
end

capture program drop coverage_failure
program define coverage_failure
    args name reason
    noi di as error "OPTION_FAILURE: `name': `reason'"
    test_fail "`name'" "`reason'"
end

capture program drop coverage_matching
program define coverage_matching
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string)
    preserve
    tempfile original
    quietly save `original'
    capture quietly `reference'
    local nr=_rc
    if `nr' {
        restore
        coverage_failure "`name'" "reference rc=`nr'"
        exit
    }
    local att=r(att)
    local se=r(seatt)
    sort id
    replace _weight=0 if missing(_weight)
    mata: coverage_matches=st_data(.,tokens("_pscore _weight _support _treated"))
    quietly use `original', clear
    capture quietly `candidate'
    local cr=_rc
    if `cr' {
        restore
        coverage_failure "`name'" "candidate rc=`cr'"
        exit
    }
    local actual=r(att)
    assert_psmatch_inference `=r(att_se)' `=r(att_t)' `se' `att' $DEFAULT_SIGFIGS "`name'"
    sigfigs `att' `actual'
    if r(sigfigs)<$DEFAULT_SIGFIGS coverage_failure "`name'" "ATT differs"
    sort id
    replace _weight=0 if missing(_weight)
    mata: st_local("same",strofreal(coverage_isnumeric_equal(coverage_matches,st_data(.,tokens("_pscore _weight _support _treated")))))
    if !`same' coverage_failure "`name'" "scores, weights, or support differ"
    else test_pass "`name' scores weights support"
    restore
end

capture program drop coverage_sampling
program define coverage_sampling
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string)
    local stratified = strpos(`"`reference'"',"strata(group)") | strpos(`"`reference'"',"by(group)")
    preserve
    tempfile original
    quietly save `original'
    capture quietly `reference'
    local nr=_rc
    if !`nr' {
        if !`stratified' replace group=0
        * Compare exact draws per group, and the eligibility of retained rows.
        if _N {
            contract group, freq(draws)
            sort group
            mata: coverage_draws=st_data(.,tokens("group draws"))
        }
        else mata: coverage_draws=J(0,2,.)
    }
    quietly use `original', clear
    capture quietly `candidate'
    local cr=_rc
    local same=0
    if !`cr' {
        if strpos(`"`candidate'"',"cluster(g)") & _N {
            tempfile sampled
            quietly save `sampled'
            contract g id, freq(copies)
            bysort g: gen members=_N
            bysort g: gen equal_copies=(copies==copies[1])
            capture assert members==3 & equal_copies
            local cluster_rc=_rc
            quietly use `sampled', clear
            if `cluster_rc' local cr=9
        }
        if strpos(`"`candidate'"',"csample ")==1 {
            capture isid id
            if _rc & _N local cr=9
        }
        if !`stratified' replace group=0
        capture assert inrange(id,1,600) & id==floor(id)
        if _rc local cr=9
        if _N {
            contract group, freq(draws)
            sort group
            if !`nr' mata: st_local("same",strofreal(coverage_isnumeric_equal(coverage_draws,st_data(.,tokens("group draws")))))
        }
        else if !`nr' mata: st_local("same",strofreal(rows(coverage_draws)==0))
    }
    restore
    if `nr' | `cr' | !`same' coverage_failure "`name'" "reference rc=`nr'; candidate rc=`cr'; group counts equal=`same'"
    else test_pass "`name' sampling counts"
end

capture program drop coverage_destring_reference
program define coverage_destring_reference
    syntax varname [if] [in], GENerate(name) [*]
    marksample selected, novarlist
    tempvar original
    clonevar `original'=`varlist'
    quietly replace `varlist'="" if !`selected'
    destring `varlist', gen(`generate') `options'
    quietly replace `varlist'=`original'
end

/* RAWfields: the native file with every quoted field unquoted, the raw fields
   that cexport's documented noquoteif extension writes (native export has no
   such option; validate_option_contracts.do checks the raw-field contract). */
capture program drop coverage_export_readback
program define coverage_export_readback
    syntax, CTOOLS(integer) KIND(string) [OPTS(string asis) RAWfields]
    tempfile output
    local extension=cond("`kind'"=="excel","xlsx","csv")
    local path "`output'.`extension'"
    local command=cond(`ctools',"cexport","export")
    `command' `kind' id y s group using "`path'", replace `opts'
    if "`rawfields'"!="" mata: coverage_csv_rawfields(st_local("path"))
    if "`kind'"=="excel" import excel using "`path'", clear
    else import delimited using "`path'", clear varnames(nonames) asdouble
    erase "`path'"
end

mata:
void coverage_csv_rawfields(string scalar path) {
    real scalar fh, p, quoted
    string scalar text, out, q
    fh=fopen(path,"r"); text=fread(fh,100000000); fclose(fh)
    q=char(34); out=""; quoted=0
    // Outside a field's quotes every quote opens one; inside, a doubled quote
    // is a literal quote and a single quote closes the field.
    while((p=strpos(text,q))>0) {
        out=out+substr(text,1,p-1)
        text=substr(text,p+1,.)
        if(quoted & substr(text,1,1)==q) {
            out=out+q
            text=substr(text,2,.)
        }
        else quoted=!quoted
    }
    out=out+text
    unlink(path)
    fh=fopen(path,"w"); fwrite(fh,out); fclose(fh)
}
end

capture program drop coverage_bins
program define coverage_bins
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string)
    preserve
    tempfile original
    quietly save `original'
    coverage_prepare_bins_reference
    local precise_reference : subinstr local reference "binscatter " "coverage_binscatter ", all
    capture quietly `precise_reference'
    local nr=_rc
    if !`nr' {
        mata: raw_bins=st_matrix("e(coverage_bindata)"); coverage_bins_ref=J(0,3,.)
        mata: for(j=1;j<=cols(raw_bins)/2;j++) coverage_bins_ref=coverage_bins_ref\(J(rows(raw_bins),1,j),raw_bins[,(2*j-1,2*j)])
        mata: coverage_bins_ref=select(coverage_bins_ref,coverage_bins_ref[,2]:<.)
        mata: coverage_bins_ref=sort(coverage_bins_ref,(1,2))
    }
    quietly use `original', clear
    capture quietly `candidate'
    local cr=_rc
    local same=0
    if !`cr' & !`nr' {
        mata: coverage_bins_got=st_matrix("e(bindata)")
        mata: coverage_bins_got=select(coverage_bins_got,coverage_bins_got[,2]:<.)
        mata: coverage_bins_got=sort(coverage_bins_got[,(1,3,4)],(1,2))
        * binscatter stores controls()/absorb() residuals with predict's default
        * float type, so its bin coordinates carry rounding of 6e-8 of the
        * column's scale: groups compare exactly, coordinates column-wise
        mata: st_local("same",strofreal(rows(coverage_bins_ref)==rows(coverage_bins_got)))
        if `same' mata: st_local("same",strofreal(all(coverage_bins_ref[,1]:==coverage_bins_got[,1]) & coverage_prediction_equal(coverage_bins_ref[,2],coverage_bins_got[,2]) & coverage_prediction_equal(coverage_bins_ref[,3],coverage_bins_got[,3])))
    }
    restore
    if `nr' | `cr' | !`same' coverage_failure "`name'" "reference rc=`nr'; candidate rc=`cr'; bins equal=`same'"
    else test_pass "`name' bin coordinates"
end

capture program drop coverage_weighted_bootstrap
program define coverage_weighted_bootstrap
    syntax, CTOOLS(integer) [*]
    gen double draw_weight=.
    local cmd=cond(`ctools',"cbsample","bsample")
    `cmd', weight(draw_weight) `options'
    assert draw_weight>=0 & draw_weight<. & draw_weight==floor(draw_weight)
    drop if draw_weight==0
    expand draw_weight
end

capture program drop coverage_graph
program define coverage_graph
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string) [TEXT(string)]
    tempfile native_svg candidate_svg
    capture quietly {
        `reference'
        graph export "`native_svg'.svg", as(svg) replace
    }
    local nr=_rc
    capture quietly {
        `candidate'
        graph export "`candidate_svg'.svg", as(svg) replace
    }
    local cr=_rc
    local same=0
    if !`nr' & !`cr' {
        mata: fh=fopen(st_local("native_svg")+".svg","r"); coverage_svg_ref=fread(fh,10000000); fclose(fh)
        mata: fh=fopen(st_local("candidate_svg")+".svg","r"); coverage_svg_got=fread(fh,10000000); fclose(fh)
        local same=1
        foreach tag in circle polygon {
            mata: st_local("ns",strofreal((strlen(coverage_svg_ref)-strlen(subinstr(coverage_svg_ref,"<"+st_local("tag"),"")))/(strlen(st_local("tag"))+1)))
            mata: st_local("cs",strofreal((strlen(coverage_svg_got)-strlen(subinstr(coverage_svg_got,"<"+st_local("tag"),"")))/(strlen(st_local("tag"))+1)))
            if `ns'!=`cs' local same=0
        }
        foreach token of local text {
            mata: st_local("present",strofreal(strpos(coverage_svg_ref,st_local("token"))>0 & strpos(coverage_svg_got,st_local("token"))>0))
            if !`present' local same=0
        }
    }
    capture erase "`native_svg'.svg"
    capture erase "`candidate_svg'.svg"
    if `nr' | `cr' | !`same' coverage_failure "`name'" "reference rc=`nr'; candidate rc=`cr'; rendered properties equal=`same'"
    else test_pass "`name' rendered properties"
end

/* Restore the named first-stage/reduced-form estimate, never the main IV fit.
   Clear previous estimates so a stale native result cannot mask a missing one. */
capture program drop coverage_saved_iv
program define coverage_saved_iv, eclass
    syntax, CTOOLS(integer) STAGE(string)
    estimates clear
    local command=cond(`ctools',"civreghdfe","ivreghdfe")
    if "`stage'"=="first" {
        `command' y x (endog=z z2), absorb(g) savefirst savefprefix(coverage_)
        estimates restore coverage_endog
    }
    else {
        `command' y x (endog=z z2), absorb(g) saverf saverfprefix(coverage_)
        estimates restore coverage_y
    }
end

capture program drop coverage_iv_firstref
program define coverage_iv_firstref, eclass
    ivreghdfe `0'
    tempname first
    matrix `first'=e(first)
    ereturn scalar F_first1=`first'[rownumb(`first',"F"),1]
    ereturn scalar partial_r2_1=`first'[rownumb(`first',"pr2"),1]
end

/* Reference estimators lose accuracy when y is far from unit scale:
   ivreg2's CUE optimize() tests convergence with mreldif(), which is absolute
   for coefficients far below one, and reghdfe's saved residuals are off by
   1e-7 of their range for y ~ 1e-6 (an explicit-dummy regress agrees with
   creghdfe to 1e-13). Both estimators are equivariant to the units of y:
   estimate with y in units of its standard deviation, then rescale b, V, the
   scale-dependent scalars and the saved residuals. */
capture program drop coverage_unitscale
program define coverage_unitscale, eclass
    gettoken command 0 : 0
    gettoken depvar : 0
    quietly summarize `depvar'
    local s = r(sd)
    tempvar original
    quietly clonevar `original' = `depvar'
    quietly replace `depvar' = `depvar'/`s'
    capture noisily `command' `0'
    local rc = _rc
    quietly replace `depvar' = `original'
    if `rc' exit `rc'
    tempname b V
    matrix `b' = e(b)*`s'
    matrix `V' = e(V)*(`s'^2)
    ereturn repost b=`b' V=`V'
    foreach stat in rss mss tss tss_within {
        if !missing(e(`stat')) ereturn scalar `stat' = e(`stat')*`s'^2
    }
    if !missing(e(rmse)) ereturn scalar rmse = e(rmse)*`s'
    foreach stat in ll ll_0 {
        if !missing(e(`stat')) ereturn scalar `stat' = e(`stat') - e(N)*ln(`s')
    }
    if "`e(resid)'" != "" quietly replace `e(resid)' = `e(resid)'*`s'
end

/* ivreghdfe 1.1.4 prints its main table with mata _coef_table(), which reads
   only the plus and level locals of DispMain: the eform() and factor display
   options it parses and passes in as dispopt are dropped (its commented-out
   ereturn display call, and ivreg2's DispMain, apply them). The reference is
   that intended display: the reference estimates, then ereturn display. */
capture program drop coverage_iv_replay
program define coverage_iv_replay
    syntax, ESTimate(string asis) [DISPlay(string asis)]
    quietly `estimate'
    ereturn display, `display'
end

/* Inspect option-specific output. Numerical estimands are checked separately;
   this catches options that are accepted but ignored by the display layer. */
capture program drop coverage_display
program define coverage_display
    syntax, REFerence(string asis) CANDidate(string asis) NAME(string) ///
        [HAS(string) ABSENT(string) TABLE]
    tempfile reference_log candidate_log
    tempname reference_handle candidate_handle
    log using "`reference_log'.log", text replace name(`reference_handle')
    capture noisily `reference'
    local nr=_rc
    log close `reference_handle'
    log using "`candidate_log'.log", text replace name(`candidate_handle')
    capture noisily `candidate'
    local cr=_rc
    log close `candidate_handle'
    mata: fh=fopen(st_local("reference_log")+".log","r"); ref_display=fread(fh,10000000); fclose(fh)
    mata: fh=fopen(st_local("candidate_log")+".log","r"); got_display=fread(fh,10000000); fclose(fh)
    mata: st_local("logs_present",strofreal(length(ref_display)==1 & length(got_display)==1))
    if !`logs_present' {
        capture erase "`reference_log'.log"
        capture erase "`candidate_log'.log"
        coverage_failure "`name'" "empty display capture; invoke this helper outside quietly blocks"
        exit
    }
    local same=1
    if "`has'"!="" {
        mata: st_local("present",strofreal(regexm(strlower(ref_display),strlower(st_local("has"))) & regexm(strlower(got_display),strlower(st_local("has")))))
        if !`present' local same=0
    }
    if "`absent'"!="" {
        mata: st_local("missing",strofreal(!regexm(strlower(ref_display),strlower(st_local("absent"))) & !regexm(strlower(got_display),strlower(st_local("absent")))))
        if !`missing' local same=0
    }
    if "`table'"!="" {
        mata: reference_table=coverage_display_table(ref_display); candidate_table=coverage_display_table(got_display)
        mata: st_local("table_same",strofreal(strlen(reference_table)>0 & reference_table==candidate_table))
        if !`table_same' {
            local same=0
            noi di "REFERENCE_TABLE:"
            mata: printf("%s\n",reference_table)
            noi di "CANDIDATE_TABLE:"
            mata: printf("%s\n",candidate_table)
        }
    }
    capture erase "`reference_log'.log"
    capture erase "`candidate_log'.log"
    if `nr' | `cr' | !`same' coverage_failure "`name'" "reference rc=`nr'; candidate rc=`cr'; displayed properties equal=`same'"
    else test_pass "`name' displayed properties"
end
mata:
string scalar coverage_display_table(string scalar text) {
    string rowvector lines
    string scalar result,line
    real scalar active,i
    lines=tokens(text,char(10)); active=0; result=""
    for(i=1;i<=cols(lines);i++) {
        line=strtrim(stritrim(lines[i]))
        if(strpos(line,"Std. err.") | strpos(line,"Std. Err.")) active=1
        if(active) {
            if(regexm(line,"^[-]+$")) break
            if(strpos(line,"|")) result=result+line+char(10)
        }
    }
    return(result)
}
end

capture program drop coverage_export_bytes
program define coverage_export_bytes
    syntax, NAME(string) [OPTS(string) CRLF]
    tempfile refcsv gotcsv
    capture quietly export delimited id y s group using "`refcsv'", replace `opts'
    local nr=_rc
    capture quietly cexport delimited id y s group using "`gotcsv'", replace `opts' `crlf'
    local cr=_rc
    local same=0
    if !`nr' & !`cr' {
        mata: fh=fopen(st_local("refcsv"),"r"); ref_bytes=fread(fh,10000000); fclose(fh)
        mata: fh=fopen(st_local("gotcsv"),"r"); got_bytes=fread(fh,10000000); fclose(fh)
        if "`crlf'"!="" {
            mata: ref_bytes=subinstr(ref_bytes,char(13)+char(10),char(10)); ref_bytes=subinstr(ref_bytes,char(10),char(13)+char(10))
        }
        mata: st_local("same",strofreal(ref_bytes==got_bytes))
    }
    if `nr' | `cr' | !`same' coverage_failure "`name'" "reference rc=`nr'; candidate rc=`cr'; CSV bytes equal=`same'"
    else test_pass "`name' CSV bytes"
end

/* Native encode/decode have no replace option: implement that extension by
   generating, dropping, and renaming, preserving original variable order. */
capture program drop coverage_codec_replace
program define coverage_codec_replace
    syntax varname, COMMAND(string)
    unab original_order : _all
    tempvar encoded
    if "`command'"=="encode" encode `varlist', gen(`encoded') label(`varlist')
    else decode `varlist', gen(`encoded')
    drop `varlist'
    rename `encoded' `varlist'
    order `original_order'
end

capture program drop coverage_merge_preserved
program define coverage_merge_preserved
    syntax using/
    tempvar original_order
    gen long `original_order'=_n
    merge 1:1 id using "`using'"
    sort `original_order', stable
    drop `original_order'
end

/* binscatter's CSV uses file write's short numeric format. Expose its original
   in-memory matrix in a namespaced test copy, without changing any computation.
   Exact anchors make reference-source changes fail instead of silently adapting. */
capture program drop coverage_prepare_bins_reference
program define coverage_prepare_bins_reference
    capture program list coverage_binscatter
    if !_rc exit
    quietly findfile binscatter.ado
    local source_path `"`r(fn)'"'
    tempfile reference_copy
    tempname src dst
    file open `src' using `"`source_path'"', read text
    file open `dst' using "`reference_copy'", write text replace
    file read `src' line
    local declarations=0
    local exposures=0
    while r(eof)==0 {
        if strpos(`"`macval(line)'"',"program define binscatter,") local ++declarations
        foreach identifier in binscatter means_in_boundaries fastxtile characterize_unique_vals_sorted {
            local replacement "coverage_`identifier'"
            if "`identifier'"=="characterize_unique_vals_sorted" local replacement "coverage_unique_vals"
            local line : subinstr local line "`identifier'" "`replacement'", all
        }
        file write `dst' `"`macval(line)'"' _n
        if strpos(`"`macval(line)'"',"ereturn scalar N =") {
            local expose="ereturn matrix coverage_bindata = "+char(96)+"y1_scatterpts"+char(39)
            file write `dst' `"`macval(expose)'"' _n
            local ++exposures
        }
        file read `src' line
    }
    file close `src'
    file close `dst'
    assert `declarations'==1 & `exposures'==1
    quietly run "`reference_copy'"
end
