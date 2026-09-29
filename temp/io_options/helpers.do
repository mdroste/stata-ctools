/* Exact values (including binary strL), variable and dataset metadata, and
   label definitions. Native commands are reference fixtures, never fallback. */
clear all
set more off
set linesize 255
adopath ++ "build"
args pluginpath
if `"`pluginpath'"' != "" adopath ++ `"`pluginpath'"'
global IO_FORMATS_PASSED 0
global IO_FORMATS_FAILED 0

mata:
struct cio_snapshot {
    string matrix values, metadata, labels, chars
    string scalar datalabel, sortedby
}
struct cio_snapshot scalar cio_snapshot_data() {
    struct cio_snapshot scalar s
    real scalar j, k
    string colvector sets, names, text
    real colvector numbers
    s.values=J(st_nobs(),st_nvar(),""); s.metadata=J(st_nvar(),5,"")
    s.labels=J(0,3,""); s.chars=J(0,2,"")
    s.datalabel=st_local("datalabel"); s.sortedby=st_local("sortedby")
    for(j=1;j<=st_nvar();j++) {
        if(st_isnumvar(j)) s.values[,j]=strofreal(st_data(.,j),"%21x")
        else s.values[,j]=st_sdata(.,j)
        s.metadata[j,]=(st_varname(j),st_vartype(j),st_varformat(j),st_varlabel(j),st_varvaluelabel(j))
    }
    sets=sort(tokens(st_local("labelsets"))',1)
    for(j=1;j<=rows(sets);j++) {
        st_vlload(sets[j],numbers,text)
        for(k=1;k<=rows(numbers);k++)s.labels=s.labels\(sets[j],strofreal(numbers[k],"%21x"),text[k])
    }
    names=sort(st_dir("char","_dta","*"),1)
    for(j=1;j<=rows(names);j++)s.chars=s.chars\(names[j],st_global("_dta["+names[j]+"]"))
    return(s)
}
void cio_long_formats(struct cio_snapshot scalar s) {
    real scalar j
    for(j=1;j<=rows(s.metadata);j++)if(s.metadata[j,2]=="strL")s.metadata[j,3]="%9s"
}
real scalar cio_same_matrix(string matrix a,string matrix b) {
    if(rows(a)!=rows(b) || cols(a)!=cols(b))return(0)
    if(!rows(a) || !cols(a))return(1)
    return(all(a:==b))
}
void cio_compare(struct cio_snapshot scalar ref) {
    struct cio_snapshot scalar got
    real scalar same, j
    got=cio_snapshot_data()
    same=cio_same_matrix(ref.values,got.values) && cio_same_matrix(ref.metadata,got.metadata) && cio_same_matrix(ref.labels,got.labels) && cio_same_matrix(ref.chars,got.chars) && ref.datalabel==got.datalabel && ref.sortedby==got.sortedby
    if(!same) {
        printf("DETAIL: values=%g metadata=%g labels=%g chars=%g datalabel=%g sorted=%g\n",cio_same_matrix(ref.values,got.values),cio_same_matrix(ref.metadata,got.metadata),cio_same_matrix(ref.labels,got.labels),cio_same_matrix(ref.chars,got.chars),ref.datalabel==got.datalabel,ref.sortedby==got.sortedby)
        if(!cio_same_matrix(ref.values,got.values)) {
            for(j=1;j<=rows(ref.values);j++)printf("REF values: %s\n",invtokens(ref.values[j,]))
            for(j=1;j<=rows(got.values);j++)printf("GOT values: %s\n",invtokens(got.values[j,]))
        }
        if(!cio_same_matrix(ref.metadata,got.metadata)) {
            for(j=1;j<=rows(ref.metadata);j++)printf("REF metadata: %s\n",invtokens(ref.metadata[j,]))
            for(j=1;j<=rows(got.metadata);j++)printf("GOT metadata: %s\n",invtokens(got.metadata[j,]))
        }
        if(!cio_same_matrix(ref.labels,got.labels)) {
            for(j=1;j<=rows(ref.labels);j++)printf("REF labels: %s\n",invtokens(ref.labels[j,]))
            for(j=1;j<=rows(got.labels);j++)printf("GOT labels: %s\n",invtokens(got.labels[j,]))
        }
        if(!cio_same_matrix(ref.chars,got.chars)) {
            for(j=1;j<=rows(ref.chars);j++)printf("REF chars: %s\n",invtokens(ref.chars[j,]))
            for(j=1;j<=rows(got.chars);j++)printf("GOT chars: %s\n",invtokens(got.chars[j,]))
        }
    }
    st_local("same",strofreal(same))
}
end

program define cio_record
    args name rc
    if `rc' {
        di as error "FAIL: `name' (rc=`rc')"
        global IO_FORMATS_FAILED = $IO_FORMATS_FAILED+1
    }
    else {
        di "PASS: `name'"
        global IO_FORMATS_PASSED = $IO_FORMATS_PASSED+1
    }
end

program define cio_workbook_test
    syntax using/, NAME(string)
    local filepath `"`using'"'
    local using ""
    quietly import excel "`filepath'", describe
    local n=r(N_worksheet)
    forvalues j=1/`n' {
        local sheet`j' `"`r(worksheet_`j')'"'
        local range`j' "`r(range_`j')'"
    }
    local datalabel : data label
    local sortedby : sortedby
    quietly label dir
    local labelsets `r(names)'
    mata: cio_ref=cio_snapshot_data()
    capture noisily cimport excel "`filepath'", describe
    local rc=_rc
    if !`rc' {
        if r(N_worksheet)!=`n' local rc=9
        forvalues j=1/`n' {
            if `"`r(worksheet_`j')'"'!=`"`sheet`j''"' | "`r(range_`j')'"!="`range`j''" local rc=9
        }
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_compare(cio_ref)
        if !`same' local rc=9
    }
    cio_record "`name'" `rc'
end

program define cio_import_test
    syntax using/, KIND(string) NAME(string) [OPTS(string asis) SELECT(string asis) LONGSTRING]
    local filepath `"`using'"'
    local using ""
    local clause using
    if "`kind'"=="sasxport5" local clause ""
    clear
    capture noisily import `kind' `select' `clause' "`filepath'", clear `opts' 
    local rc=_rc
    if `rc' {
        di "REFERENCE FAILED: `name' (rc=`rc')"
        global IO_FORMATS_FAILED = $IO_FORMATS_FAILED+1
        exit
    }
    local datalabel : data label
    local sortedby : sortedby
    quietly label dir
    local labelsets `r(names)'
    mata: cio_ref=cio_snapshot_data()
    if "`longstring'"!="" {
        * The SPI cannot assign Stata's invalid native formats for long strings.
        mata: cio_long_formats(cio_ref)
        di "KNOWN GAP: `name': native invalid long-string display format; C uses %9s"
    }
    capture noisily cimport `kind' `select' using `"`filepath'"', clear `opts'
    local rc=_rc
    if !`rc' {
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_compare(cio_ref)
        if !`same' local rc=9
    }
    cio_record "`name'" `rc'
end

program define cio_export_test
    syntax, KIND(string) EXT(string) NAME(string) [OPTS(string asis) IOPTS(string asis) SELECT(string asis)]
    if "`select'"=="" local select "*"
    if "`kind'"=="shp" local select ""
    local clause using
    if "`kind'"=="sasxport5" local clause ""
    tempfile source
    quietly save `source' 
    capture noisily export `kind' `select' using "temp/io_allformats/native/output.`ext'", replace `opts'
    local rc=_rc
    if `rc' {
        di "REFERENCE FAILED: `name' (rc=`rc')"
        global IO_FORMATS_FAILED = $IO_FORMATS_FAILED+1
        exit
    }
    capture noisily cexport `kind' `select' using "temp/io_allformats/result/output.`ext'", replace `opts'
    local rc=_rc
    if !`rc' {
        quietly import `kind' `clause' "temp/io_allformats/native/output.`ext'", clear `iopts'
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_ref=cio_snapshot_data()
        quietly import `kind' `clause' "temp/io_allformats/result/output.`ext'", clear `iopts'
        local datalabel : data label
        local sortedby : sortedby
        quietly label dir
        local labelsets `r(names)'
        mata: cio_compare(cio_ref)
        if !`same' local rc=9
    }
    cio_record "`name'" `rc'
    use `source', clear
end
