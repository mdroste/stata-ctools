*! Binary statistical readers; parsing and values are provided by C.
program define _cio_import, rclass
    version 14.1
    gettoken kind 0 : 0, parse(" ,")
    local options CLEAR
    if inlist("`kind'","sas","spss","sasxport8","dbase","xls","xlsx") local options "`options' CASE(string)"
    if "`kind'"=="sas" local options "`options' ENCoding(string) BCAT(string)"
    if "`kind'"=="sasxport8" local options "`options' VLABFile(string)"
    if "`kind'"=="spss" local options "`options' ENCoding(string) ZSAV"
    if "`kind'"=="sasxport5" local options "`options' Describe NOVALlabels Member(string)"
    if inlist("`kind'","xls","xlsx") local options "`options' SHEET(string) CELLRange(string) FIRSTrow ALLString ALLSTRINGFMT(string) DEScribe"
    * Parse selectors as text: their variables belong to the incoming dataset.
    syntax [anything(everything equalok)] [, `options']
    if "`kind'"=="sasxport5" & "`describe'"!="" & ("`clear'"!="" | "`novallabels'"!="") exit 198
    local rest `"`anything'"'
    local prefix ""
    local found 0
    while `"`rest'"'!="" {
        gettoken token rest : rest, bind quotes
        if `"`token'"'=="using" {
            local found 1
            continue, break
        }
        local prefix `"`prefix' `token'"'
    }
    if `found' {
        gettoken using extra : rest
        local selection `"`prefix' `extra'"'
    }
    else {
        gettoken using selection : anything
    }
    if `"`using'"'=="" exit 198
    local rest = strtrim(`"`selection'"')
    local wanted ""
    local selectors ""
    while `"`rest'"'!="" {
        gettoken token rest : rest, bind quotes
        if inlist(`"`token'"',"if","in") {
            local selectors `"`token' `rest'"'
            continue, break
        }
        local wanted `"`wanted' `token'"'
    }
    local wanted = strtrim(`"`wanted'"')
    if !inlist("`kind'","sas","spss","xls","xlsx") & `"`wanted'"'!="" exit 198
    if !inlist("`kind'","sas","spss") & `"`selectors'"'!="" exit 198
    local ext ".sas7bdat"
    if "`kind'" == "xls" local ext ".xls"
    if "`kind'" == "xlsx" local ext ".xlsx"
    if "`kind'" == "shp" local ext ".shp"
    if "`kind'" == "dbase" local ext ".dbf"
    if "`kind'" == "spss" local ext ".sav"
    if "`kind'" == "spss" & "`zsav'" != "" local ext ".zsav"
    if "`kind'" == "sasxport5" local ext ".xpt"
    if "`kind'" == "sasxport8" local ext ".v8xpt"
    mata: st_local("base",pathbasename(st_local("using")))
    if !strpos(`"`base'"',".") local using `"`using'`ext'"'
    confirm file `"`using'"'
    if "`case'" == "" local case = cond("`kind'"=="sasxport5","lower","preserve")
    if "`case'"!="" {
        if substr("lower",1,strlen("`case'"))=="`case'" local case "lower"
        if substr("upper",1,strlen("`case'"))=="`case'" local case "upper"
        if strlen("`case'")>=3 & substr("preserve",1,strlen("`case'"))=="`case'" local case "preserve"
    }
    if inlist("`case'","low","u","pre") {
        if "`case'" == "low" local case "lower"
        if "`case'" == "u" local case "upper"
        if "`case'" == "pre" local case "preserve"
    }
    if !inlist("`case'","preserve","lower","upper") exit 198
    if "`clear'" == "" & "`describe'" == "" {
        if "`kind'"=="sasxport5" {
            if c(changed) exit 4
        }
        else if (_N | c(k)) exit 4
    }
    if `"`bcat'"' != "" {
        mata: st_local("base",pathbasename(st_local("bcat")))
        if !strpos(`"`base'"',".") local bcat `"`bcat'.sas7bcat"'
        confirm file `"`bcat'"'
    }
    if "`kind'"=="sasxport8" & `"`vlabfile'"'!="" {
        mata: st_local("vlabbase",pathbasename(st_local("vlabfile")))
        if !strpos(`"`vlabbase'"',".") local vlabfile `"`vlabfile'.v8xpt"'
        confirm file `"`vlabfile'"'
    }
    * Native statistical/DBF/SHP readers clear before parsing the file, even
    * when parsing fails. Excel validates its container before replacing data.
    if "`describe'"=="" & !inlist("`kind'","xls","xlsx") quietly clear
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local __cio_member "`member'"
    local __cio_sheet `"`sheet'"'
    local __cio_cellrange "`cellrange'"
    local __cio_firstrow = ("`firstrow'" != "")
    local __cio_allstring = ("`allstring'" != "")
    local __cio_allstringfmt `"`allstringfmt'"'
    local __cio_filename `"`using'"'
    local __cio_encoding `"`encoding'"'
    local __cio_bcat `"`bcat'"'
    local __cio_xpf ""
    local __cio_rawstorage 0
    if "`kind'" == "sasxport5" & "`novallabels'" == "" & "`describe'" == "" {
        mata: st_local("xpfdir",pathgetparent(st_local("using")))
        foreach candidate in formats.xpf FORMATS.xpf {
            capture confirm file `"`xpfdir'/`candidate'"'
            if !_rc local __cio_xpf `"`xpfdir'/`candidate'"'
        }
    }
    if "`describe'"!="" & inlist("`kind'","xls","xlsx") {
        capture noisily plugin call ctools_plugin, "cio workbook `kind'"
        local rc=_rc
        if `rc' {
            macro drop CTOOLS_CIO_*
            exit `rc'
        }
        local sheets $CTOOLS_CIO_nworksheets
        return scalar N_worksheet=`sheets'
        forvalues j=1/`sheets' {
            mata: st_local("worksheet",st_global("CTOOLS_CIO_worksheet_"+st_local("j")))
            mata: st_local("range",st_global("CTOOLS_CIO_range_"+st_local("j")))
            return local worksheet_`j' `"`worksheet'"'
            return local range_`j' "`range'"
            di `"`worksheet'  `range'"'
        }
        macro drop CTOOLS_CIO_*
        exit
    }
    if "`kind'"=="sasxport8" & `"`vlabfile'"'!="" {
        local __cio_filename `"`vlabfile'"'
        local __cio_rawstorage 1
        capture noisily plugin call ctools_plugin, "cio scan sasxport8"
        local rc=_rc
        if `rc' {
            capture plugin call ctools_plugin, "cio clear"
            macro drop CTOOLS_CIO_*
            exit `rc'
        }
        * The native command loads its companion before reading the main file.
        * Its second read fails before compression or case conversion, leaving
        * the companion's declared storage and original names in memory.
        if $CTOOLS_CIO_nvar>0 {
            local nobs $CTOOLS_CIO_nobs
            local nvar $CTOOLS_CIO_nvar
            local nlabels $CTOOLS_CIO_nlabels
            capture noisily _cio_import_load sasxport8 `nobs' `nvar' `nlabels', case(preserve)
            local rc=_rc
            capture plugin call ctools_plugin, "cio clear"
            macro drop CTOOLS_CIO_*
            if `rc' exit `rc'
            return scalar N=_N
            return scalar k=c(k)
            return scalar width=c(width)
            return scalar changed=c(changed)
            local datalabel : data label
            return local datalabel `"`macval(datalabel)'"'
            di as error "no; dataset in memory has changed since last saved"
            exit 4
        }
        local __cio_rawstorage 0
        local __cio_filename `"`using'"'
    }
    capture noisily plugin call ctools_plugin, "cio scan `kind'"
    local rc = _rc
    if `rc' {
        macro drop CTOOLS_CIO_*
        exit `rc'
    }
    local nobs = $CTOOLS_CIO_nobs
    local nvar = $CTOOLS_CIO_nvar
    local nlabels = $CTOOLS_CIO_nlabels
    if "`describe'" != "" {
        local nmembers $CTOOLS_CIO_nmembers
        local members = strlower("$CTOOLS_CIO_members")
        local descriptions "`member'"
        if "`kind'"=="sasxport5" & "`member'"=="" local descriptions "`members'"
        foreach item of local descriptions {
            if "`kind'"=="sasxport5" {
                local __cio_member "`item'"
                plugin call ctools_plugin, "cio scan `kind'"
                local nobs $CTOOLS_CIO_nobs
                local nvar $CTOOLS_CIO_nvar
                local size $CTOOLS_CIO_size
            }
            di "`item': `nvar' variables, `nobs' observations"
            forvalues j=1/`nvar' {
                plugin call ctools_plugin, "cio column `j'"
                di "${CTOOLS_CIO_name}  ${CTOOLS_CIO_type}  ${CTOOLS_CIO_label}"
            }
        }
        plugin call ctools_plugin, "cio clear"
        macro drop CTOOLS_CIO_*
        return scalar N = `nobs'
        return scalar k = `nvar'
        if "`kind'"=="sasxport5" {
            return scalar size = `size'
            if "`member'"=="" {
                return scalar n_members = `nmembers'
                return local members "`members'"
            }
        }
        exit
    }
    * Preserve also protects data from allocation or SPI failures after scan.
    preserve
    local loadcase "`case'"
    if inlist("`kind'","sas","spss") | (inlist("`kind'","xls","xlsx") & "`firstrow'"=="") local loadcase preserve
    capture noisily _cio_import_load `kind' `nobs' `nvar' `nlabels', case(`loadcase') wanted(`"`wanted'"') allstringfmt(`"`allstringfmt'"') `novallabels' `allstring'
    local rc = _rc
    capture plugin call ctools_plugin, "cio clear"
    if `rc' {
        macro drop CTOOLS_CIO_*
        restore
        exit `rc'
    }
    macro drop CTOOLS_CIO_*
    if `"`selectors'"' != "" quietly keep `selectors'
    if `"`wanted'"' != "" & !inlist("`kind'","xls","xlsx") keep `wanted'
    if inlist("`kind'","sas","spss") {
        if "`case'" != "preserve" rename _all, `case'
        quietly compress
    }
    restore, not
    return scalar N = _N
    return scalar k = c(k)
    return scalar width = c(width)
    return scalar changed = 1
    local datalabel : data label
    return local datalabel `"`macval(datalabel)'"'
    return local filename `"`using'"'
end

program define _cio_import_load
    version 14.1
    gettoken kind 0 : 0, parse(" ,")
    gettoken nobs 0 : 0, parse(" ,")
    gettoken nvar 0 : 0, parse(" ,")
    gettoken nlabels 0 : 0, parse(" ,")
    syntax, CASE(string) [WANTED(string) NOVALLabels ALLSTRINGFMT(string) ALLString]
    * Plugin registrations belong to the invoking ado program.
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    clear
    quietly set obs `=max(1,`nobs')'
    if `nvar' {
    forvalues j=1/`nvar' {
        plugin call ctools_plugin, "cio column `j'"
        if inlist("`kind'","xls","xlsx") {
            mata: st_local("vname",_cio_excel_name(st_global("CTOOLS_CIO_name"),st_global("CTOOLS_CIO_excel_column")))
        }
        else mata: st_local("vname",ustrtoname(st_global("CTOOLS_CIO_name")))
        if "`case'" == "lower" & !inlist("`kind'","xls","xlsx") local vname = ustrlower("`vname'")
        if "`case'" == "upper" & !inlist("`kind'","xls","xlsx") local vname = ustrupper("`vname'")
        if inlist("`kind'","xls","xlsx") {
            if "`case'" == "lower" local vname = strlower("`vname'")
            if "`case'" == "upper" local vname = strupper("`vname'")
        }
        capture confirm new variable `vname'
        if _rc {
            if inlist("`kind'","xls","xlsx") local vname "${CTOOLS_CIO_excel_column}"
            else local vname "v`j'"
        }
        if $CTOOLS_CIO_external local external `external' `j'
        if $CTOOLS_CIO_numerictext local textnumbers `textnumbers' `j'
        local sourcecol`j' "${CTOOLS_CIO_excel_column}"
        local sourcevar`j' "`vname'"
        local type ${CTOOLS_CIO_type}
        mata: (void) st_addvar(st_local("type"),st_local("vname"))
        mata: st_varlabel(st_local("vname"),st_global("CTOOLS_CIO_label"))
        mata: st_varformat(st_local("vname"),st_global("CTOOLS_CIO_format"))
        if ("`novallabels'" == "" | "`kind'"=="sasxport5") & inlist("`kind'","spss","sasxport5") & "${CTOOLS_CIO_labelset}" != "" {
            local set ${CTOOLS_CIO_labelset}
            * Native SPSS labels use their dictionary ordinal, beginning at 0.
            local pos : list posof "`set'" in sets
            if !`pos' {
                local sets `sets' `set'
                local pos : word count `sets'
            }
            local labelname "labels`=`pos'-1'"
            if "`kind'" == "sasxport5" local labelname = ustrlower(ustrtoname("`set'"))
            mata: st_varvaluelabel(st_local("vname"),st_local("labelname"))
        }
    }
    }
    if `nobs' {
    foreach hdridx of local external {
        forvalues i=1/`nobs' {
            local offset 0
            mata: __ctools_cio_blob = ""
            while 1 {
                plugin call ctools_plugin, "cio blob `hdridx' `i' `offset'"
                mata: __ctools_cio_blob = __ctools_cio_blob + _cio_hex_bytes(st_global("CTOOLS_CIO_blobhex"))
                local offset = $CTOOLS_CIO_blobnext
                if !`offset' continue, break
            }
            mata: st_sstore(`i',`hdridx',__ctools_cio_blob)
        }
        mata: mata drop __ctools_cio_blob
    }
    }
    if "`kind'" == "shp" {
        foreach field in type xmin ymin xmax ymax zmin zmax mmin mmax num_records {
            mata: st_global("_dta[shp_"+st_local("field")+"]",st_global("CTOOLS_CIO_shp_"+st_local("field")))
        }
        char _dta[dta_shp_ver] "1"
        char _dta[shp_code] "9994"
        char _dta[shp_version] "1000"
    }
    if "`kind'"=="dbase" {
        foreach field in version year month day {
            mata: st_global("_dta[dbf_"+st_local("field")+"]",st_global("CTOOLS_CIO_dbf_"+st_local("field")))
        }
    }
    if `nobs' & `nvar' plugin call ctools_plugin *, "cio load"
    if `nobs' & `"`allstringfmt'"'!="" {
        foreach j of local textnumbers {
            tempvar rawnumber
            quietly gen double `rawnumber' = .
            plugin call ctools_plugin `rawnumber', "cio textnumbers `j'"
            local target `sourcevar`j''
            quietly recast strL `target'
            mata: _cio_format_text(st_local("target"),st_local("rawnumber"),st_local("allstringfmt"))
            quietly drop `rawnumber'
            quietly compress `target'
            mata: st_varformat(st_local("target"),"%"+strofreal(st_vartype(st_local("target"))=="strL" ? 9 : max((9,max(strlen(st_sdata(.,st_local("target")))))))+"s")
        }
    }
    if !`nobs' quietly drop in 1
    if "`kind'" == "shp" {
        capture confirm variable shape_order
        if !_rc sort _ID
    }
    if "`novallabels'" == "" & `nlabels' {
        forvalues i=1/`nlabels' {
            plugin call ctools_plugin, "cio label `i'"
            local set ${CTOOLS_CIO_labelset}
            local pos : list posof "`set'" in sets
            if "`kind'" == "spss" local labelname "labels`=`pos'-1'"
            else mata: st_local("labelname",ustrlower(ustrtoname(st_global("CTOOLS_CIO_labelset"))))
            mata: st_vlmodify(st_local("labelname"),strtoreal(st_global("CTOOLS_CIO_labelvalue")),st_global("CTOOLS_CIO_labeltext"))
        }
    }
    if inlist("`kind'","xls","xlsx") & `"`wanted'"'!="" {
        local rest `"`wanted'"'
        local desired ""
        local selected ""
        local position 0
        while strtrim(`"`rest'"')!="" {
            gettoken name rest : rest, parse(" =")
            if `: list name in desired' exit 198
            local ++position
            gettoken equal tail : rest, parse(" =")
            local column ""
            if "`equal'"=="=" {
                gettoken column rest : tail
                local column=upper("`column'")
            }
            else {
                mata: st_local("column",_cio_excel_column(`position'))
            }
            local index 0
            forvalues j=1/`nvar' {
                if "`sourcecol`j''"=="`column'" local index `j'
            }
            tempvar picked`position'
            if !`index' {
                if !`nvar' exit 198
                mata: st_local("columnindex",strofreal(_cio_excel_index(st_local("column"))))
                mata: st_local("columnlimit",strofreal(_cio_excel_index(st_local("sourcecol`nvar'"))+1))
                if `columnindex'<1 | `columnindex'>`columnlimit' exit 198
                if "`allstring'"!="" quietly gen str1 `picked`position''=""
                else {
                    quietly gen byte `picked`position''=.
                    format `picked`position'' %10.0g
                }
            }
            else clonevar `picked`position''=`sourcevar`index''
            local selected `selected' `picked`position''
            local desired `desired' `name'
        }
        keep `selected'
        local j=0
        foreach v of local selected {
            local ++j
            local name : word `j' of `desired'
            rename `v' `name'
        }
    }
    mata: st_local("__data_label",st_global("CTOOLS_CIO_datalabel"))
    label data `"`macval(__data_label)'"'
end

mata:
real scalar _cio_excel_index(string scalar column) {
    real scalar j, n, code
    n=0
    for(j=1;j<=strlen(column);j++) {
        code=ascii(substr(column,j,1))-64
        if(code<1 || code>26)return(.)
        n=n*26+code
        if(n>16384)return(.)
    }
    return(n)
}
string scalar _cio_excel_column(real scalar index) {
    string scalar column
    column=""
    while(index>0) {
        column=char(65+mod(index-1,26))+column
        index=floor((index-1)/26)
    }
    return(column)
}
string scalar _cio_hex_bytes(string scalar hex) {
    real scalar i, n
    real rowvector bytes
    string scalar digits
    digits="0123456789abcdef"; n=strlen(hex)/2; bytes=J(1,n,0)
    for(i=1;i<=n;i++) bytes[i]=(strpos(digits,substr(hex,2*i-1,1))-1)*16+strpos(digits,substr(hex,2*i,1))-1
    return(char(bytes))
}
end

mata:
void _cio_format_text(string scalar target, string scalar source, string scalar fmt) {
    real colvector values, selected
    values=st_data(.,source)
    selected=selectindex(values:<.)
    if (rows(selected)) st_sstore(selected,target,strofreal(values[selected],fmt))
}
end

mata:
string scalar _cio_excel_name(string scalar source, string scalar fallback) {
    string scalar name
    name=ustrregexra(source,"[^\p{L}\p{N}_]","")
    name=usubstr(name,1,32)
    if (name=="" || ustrregexm(name,"^[0-9]")) return(fallback)
    if (!st_isvarname(name)) return(fallback)
    return(name)
}
end
