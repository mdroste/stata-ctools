*! Binary statistical writers; serialization is provided by C.
program define _cio_export, rclass
    version 14.1
    gettoken kind 0 : 0
    local options "REPLACE"
    if "`kind'"=="spss" local options "`options' NOVALLabels"
    if "`kind'"=="sasxport5" local options "`options' REName VALLabfile(string)"
    if "`kind'"=="sasxport8" local options "`options' VALLabfile"
    if "`kind'"=="dbase" local options "`options' VERsion(string) DATAFmt ORIGDBFDATE"
    capture syntax [varlist] [if] [in] using/ [, `options']
    if _rc {
        local 0 `"using `0'"'
        syntax [if] [in] using/ [, `options']
    }
    if "`varlist'" == "" unab varlist : _all
    local ext ".sas7bdat"
    if "`kind'" == "dbase" local ext ".dbf"
    if "`kind'" == "spss" local ext ".sav"
    if "`kind'" == "sasxport5" local ext ".xpt"
    if "`kind'" == "sasxport8" local ext ".v8xpt"
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`base'"',".") local using `"`using'`ext'"'
    if "`replace'" == "" confirm new file `"`using'"'
    * Native's only explicit legacy version is version(3); IV is the default.
    if "`version'" != "" & "`version'" != "3" exit 198
    if "`version'" == "3" local version "III"
    local __cio_rename = ("`rename'" != "")
    local __cio_version "`version'"
    local __cio_origdbfdate = ("`origdbfdate'" != "")
    foreach field in year month day {
        local __cio_dbf_`field' : char _dta[dbf_`field']
    }
    local __cio_datafmt = ("`datafmt'" != "")
    local __cio_xpfoutput ""
    if "`kind'"=="sasxport5" {
        if "`vallabfile'"=="" local vallabfile "xpf"
        foreach full in xpf sascode both none {
            if "`vallabfile'" != "" & substr("`full'",1,strlen("`vallabfile'")) == "`vallabfile'" local vallabfile "`full'"
        }
        if !inlist("`vallabfile'","xpf","sascode","both","none") exit 198
        if inlist("`vallabfile'","xpf","both") {
            mata: st_local("__cio_xpfoutput",pathjoin(pathgetparent(st_local("using")),"formats.xpf"))
        }
    }
    local __cio_sascode ""
    if ("`kind'"=="sasxport8" & "`vallabfile'"!="") | ("`kind'"=="sasxport5" & inlist("`vallabfile'","sascode","both")) {
        mata: st_local("__cio_sascode",pathrmsuffix(st_local("using"))+".sas")
    }
    mata: st_local("__cio_cwd",pwd())
    mata: st_local("__cio_transportpath",pathisabs(st_local("using"))?st_local("using"):pathjoin(pwd(),st_local("using")))
    local __cio_filename `"`using'"'
    local __cio_replace = ("`replace'" != "")
    local __cio_datalabel : data label
    mata: st_local("__cio_table",strupper(pathrmsuffix(pathbasename(st_local("using")))))
    local j=0
    foreach v of local varlist {
        local ++j
        local __cio_name_`j' "`v'"
        local __cio_type_`j' : type `v'
        local __cio_format_`j' : format `v'
        local __cio_label_`j' : variable label `v'
        local __cio_labelset_`j' ""
        local __cio_lcount_`j' 0
        if "`novallabels'" == "" {
            local __cio_labelset_`j' : value label `v'
            if "`__cio_labelset_`j''" != "" {
                mata: _cio_label_metadata(`j',st_local("__cio_labelset_`j'"))
            }
        }
    }
    local __cio_nvar `j'
    local referenced ""
    forvalues i=1/`j' {
        local referenced `referenced' `__cio_labelset_`i''
    }
    quietly label dir
    local allsets `r(names)'
    if inlist("`kind'","sasxport5","sasxport8") {
        foreach set of local allsets {
            if !`: list set in referenced' {
                local ++j
                local __cio_labelset_`j' "`set'"
                mata: _cio_label_metadata(`j',st_local("set"))
            }
        }
    }
    local __cio_nmeta `j' 
    local exportvars `varlist'
    if "`kind'"=="dbase" & "`datafmt'"!="" {
        local j=0
        local formatted ""
        foreach v of local varlist {
            local ++j
            local fmt : format `v'
            if !strpos("`: type `v''","str") & substr("`fmt'",1,2)!="%t" & "`fmt'"!="%d" {
                if !regexm("`fmt'","^%-?[0-9]+[.][0-9]+[fge]$") exit 120
                tempvar dbffmt`j'
                quietly gen str2045 `dbffmt`j'' = strtrim(string(`v',"`fmt'"))
                local formatted `formatted' `dbffmt`j''
                local __cio_fmtvar_`j' = `__cio_nvar' + `: word count `formatted''
            }
        }
        local exportvars `exportvars' `formatted'
    }
    marksample touse, novarlist
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    plugin call ctools_plugin `exportvars' if `touse', "cio export `kind'"
    return local filename `"`using'"'
    quietly count if `touse'
    return scalar N = r(N)
    return scalar k = `__cio_nvar'
    di as text `"file `using' saved"'
end
mata:
void _cio_label_metadata(real scalar j,string scalar set) {
    real colvector values
    string colvector text
    real scalar k
    st_vlload(set,values,text)
    st_local(sprintf("__cio_lcount_%g",j),strofreal(rows(values)))
    for (k=1;k<=rows(values);k++) {
        st_local(sprintf("__cio_lv_%g_%g",j,k),strofreal(values[k],"%24.17e"))
        st_local(sprintf("__cio_lt_%g_%g",j,k),text[k])
    }
}
end
