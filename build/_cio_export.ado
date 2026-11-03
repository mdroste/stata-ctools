*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
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
    * dBase limits, checked as export dbase does: the variable count of the
    * selection and of the dataset, names, string widths, then the
    * observation width of the selection (c(width) after keeping it).
    if "`kind'"=="dbase" {
        if `: word count `varlist'' > 255 | c(k) > 255 {
            di as error "too many variables specified"
            di as error "    Number of variables must be less than 255."
            exit 103
        }
        local width 0
        foreach v of local varlist {
            if strlen("`v'") > 10 {
                di as error "`v': variable name must be less than 10 characters"
                exit 693
            }
            local type : type `v'
            if "`type'" == "strL" {
                di as error "variable `v' cannot be a strL"
                exit 693
            }
            if substr("`type'",1,3) == "str" & real(substr("`type'",4,.)) > 254 {
                di as error "variable `v' cannot exceed 254 characters in length"
                exit 108
            }
            local width = `width' + cond(substr("`type'",1,3)=="str", real(substr("`type'",4,.)), cond("`type'"=="byte",1,cond("`type'"=="int",2,cond("`type'"=="double",8,4))))
        }
        if `width' > 4000 {
            di as error "invalid observation width"
            di as error "    Observation width must be less than 4000."
            exit 902
        }
    }
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
            * formats.xpf beside the file, using only Stata 14.1 Mata functions.
            mata: st_local("__cio_xpfoutput",substr(st_local("using"),1,strlen(st_local("using"))-strlen(pathbasename(st_local("using"))))+"formats.xpf")
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
    * Native names a version 5 member after the file name up to its first
    * period, unchanged otherwise (the plugin keeps eight bytes).
    if "`kind'"=="sasxport5" {
        mata: st_local("__cio_table",strupper(pathbasename(st_local("using"))))
        mata: st_local("__cio_table",substr(st_local("__cio_table"),1,strpos(st_local("__cio_table")+".",".")-1))
    }
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
    * Version 5 transport writes, and checks the names of, only the value
    * labels that variables use, as export sasxport5 does.
    if "`kind'"=="sasxport8" {
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
    * SPSS and SAS version 8 transport timestamps exclude leap seconds, and
    * both store monthly, quarterly, weekly and half-yearly dates as the date
    * of the period's first day. Keep source metadata and convert only the
    * values passed to the writer; version 5 transport retains Stata units.
    if inlist("`kind'","spss","sasxport8") {
        local exportvars ""
        foreach v of local varlist {
            local fmt : format `v'
            local unit ""
            if regexm("`fmt'","^%-?t([Cmqwh])") local unit = regexs(1)
            if "`unit'" != "" {
                local convert = cond("`unit'"=="C","cofC",cond("`unit'"=="m","dofm",cond("`unit'"=="q","dofq",cond("`unit'"=="w","dofw","dofh"))))
                tempvar ordinary_time
                quietly gen double `ordinary_time' = `convert'(`v') if !missing(`v')
                quietly replace `ordinary_time' = `v' if missing(`v')
                local exportvars `exportvars' `ordinary_time'
            }
            else local exportvars `exportvars' `v'
        }
    }
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
    * Stata answers an if one observation at a time. Without one the plugin
    * takes every observation in range; an if is marked and passed last, as a
    * sample marker the plugin reads in parallel. The plugin reports the count.
    local touse ""
    local __cio_select "skipif"
    if `"`if'"' != "" {
        marksample touse, novarlist
        local __cio_select "touse"
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture noisily plugin call ctools_plugin `exportvars' `touse' `in', "cio export `kind' `__cio_select'"
    * A write-protected destination (plugin r(608)) gives native's code and text
    if _rc == 608 {
        if inlist("`kind'","spss","sasxport8") di as error `"file `using' already exists"'
        else if "`kind'"=="dbase" di as error `"file `using' could not be opened"'
        else di as error `"file `using' is read-only; cannot be modified or erased"'
        exit cond(inlist("`kind'","spss","sasxport8","dbase"),602,608)
    }
    if _rc exit _rc
    return local filename `"`using'"'
    return scalar N = `__cio_N'
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
