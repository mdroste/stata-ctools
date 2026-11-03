*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
*! Binary statistical readers; parsing and values are provided by C.
program define _cio_import, rclass
    version 14.1
    gettoken kind 0 : 0, parse(" ,")
    local options CLEAR SNAPshot FAST
    if inlist("`kind'","sas","spss","sasxport8","dbase","xls","xlsx") local options "`options' CASE(string)"
    if "`kind'"=="sas" local options "`options' ENCoding(string) BCAT(string)"
    if "`kind'"=="sasxport8" local options "`options' VLABFile(string)"
    if "`kind'"=="spss" local options "`options' ENCoding(string) ZSAV"
    if "`kind'"=="sasxport5" local options "`options' Describe NOVALlabels Member(string)"
    if inlist("`kind'","xls","xlsx") local options "`options' SHEET(string) CELLRange(string) FIRSTrow ALLString ALLSTRINGFMT(string) DEScribe"
    * Parse selectors as text: their variables belong to the incoming dataset.
    syntax [anything(everything equalok)] [, `options']
    if "`kind'"=="sasxport5" & "`describe'"!="" & ("`clear'"!="" | "`novallabels'"!="") exit 198
    if "`snapshot'" != "" & "`fast'" != "" {
        di as error "options snapshot and fast may not be combined"
        exit 198
    }
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
    * Explicit snapshot protects the caller before any destructive operation.
    * preserve is program-scoped: early scan/load errors and Break restore it.
    local snap 0
    local snapframe ""
    if "`snapshot'" != "" & "`describe'" == "" {
        preserve
        local snap 1
    }
    * The plugin is checked before the data are cleared.
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    * Native statistical/DBF/SHP readers clear before parsing the file, even
    * when parsing fails. Excel validates its container before replacing data.
    if "`describe'"=="" & !inlist("`kind'","xls","xlsx") quietly clear
    local __cio_member "`member'"
    local __cio_sheet `"`sheet'"'
    local __cio_cellrange "`cellrange'"
    local __cio_firstrow = ("`firstrow'" != "")
    local __cio_maxvar = c(maxvar)
    local __cio_maxmemory = c(max_memory)
    local __cio_allstring = ("`allstring'" != "")
    local __cio_allstringfmt `"`allstringfmt'"'
    local __cio_filename `"`using'"'
    local __cio_encoding `"`encoding'"'
    local __cio_bcat `"`bcat'"'
    local __cio_xpf ""
    local __cio_rawstorage 0
    if "`kind'" == "sasxport5" & "`novallabels'" == "" & "`describe'" == "" {
        * The file's directory with its trailing separator ("" in the working
        * directory), from pathbasename(), which Stata 14.1 has.
        mata: st_local("xpfdir",substr(st_local("using"),1,strlen(st_local("using"))-strlen(pathbasename(st_local("using")))))
        foreach candidate in formats.xpf FORMATS.xpf {
            capture confirm file `"`xpfdir'`candidate'"'
            if !_rc local __cio_xpf `"`xpfdir'`candidate'"'
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
        * Sheet names come from the file: they reach r() and the screen
        * without passing through macro expansion.
        forvalues j=1/`sheets' {
            mata: st_local("worksheet",st_global("CTOOLS_CIO_worksheet_"+st_local("j")))
            mata: st_local("range",st_global("CTOOLS_CIO_range_"+st_local("j")))
            return local worksheet_`j' : copy local worksheet
            return local range_`j' : copy local range
            mata: printf("%s  %s\n", st_local("worksheet"), st_local("range"))
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
            return local datalabel : copy local datalabel
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
        * Member names, variable names and labels come from the file: they are
        * handled in Mata and reach r() and the screen without macro expansion.
        mata: st_local("members", strlower(st_global("CTOOLS_CIO_members")))
        mata: st_local("descriptions", st_local("member"))
        if "`kind'"=="sasxport5" & "`member'"=="" local descriptions : copy local members
        mata: st_local("ndescriptions", strofreal(cols(tokens(st_local("descriptions")))))
        * An error or Break in the listing still releases the scanned file
        capture noisily {
            forvalues m=1/`ndescriptions' {
                mata: st_local("item", tokens(st_local("descriptions"))[`m'])
                if "`kind'"=="sasxport5" {
                    local __cio_member : copy local item
                    plugin call ctools_plugin, "cio scan `kind'"
                    local nobs $CTOOLS_CIO_nobs
                    local nvar $CTOOLS_CIO_nvar
                    local size $CTOOLS_CIO_size
                }
                mata: printf("%s: %s variables, %s observations\n", st_local("item"), st_local("nvar"), st_local("nobs"))
                forvalues j=1/`nvar' {
                    plugin call ctools_plugin, "cio column `j'"
                    mata: printf("%s  %s  %s\n", st_global("CTOOLS_CIO_name"), st_global("CTOOLS_CIO_type"), st_global("CTOOLS_CIO_label"))
                }
            }
        }
        local rc = _rc
        capture plugin call ctools_plugin, "cio clear"
        macro drop CTOOLS_CIO_*
        if `rc' exit `rc'
        return scalar N = `nobs'
        return scalar k = `nvar'
        if "`kind'"=="sasxport5" {
            return scalar size = `size'
            if "`member'"=="" {
                return scalar n_members = `nmembers'
                return local members : copy local members
            }
        }
        exit
    }
    * A failure after the scan (allocation, SPI, Break) leaves the data in
    * memory as they were. Only an Excel workbook reaches this point with data
    * in memory (the other readers cleared them above, as native ones do);
    * from Stata 16 they are copied into a frame of this program's own, in
    * memory (the copy fails before anything changes if memory is short), and
    * restored from it with frame copy (see cimport.ado for why the load is
    * not built in a frame that is then renamed). Option snapshot preserves
    * the data instead (on disk above c(max_preservemem)), as do loads before
    * Stata 16 and loads into an empty dataset. Option fast keeps no copy and
    * undoes nothing. preserve and the frame belong to this program, so these
    * steps are inline.
    if "`fast'" == "" & !`snap' {
        if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
            tempname snapframe
            capture noisily frame copy `c(frame)' `snapframe'
        }
        else {
            capture noisily preserve
            local snap 1
        }
        * A failed copy still releases the scanned file
        if _rc {
            local rc = _rc
            capture plugin call ctools_plugin, "cio clear"
            macro drop CTOOLS_CIO_*
            exit `rc'
        }
    }
    local loadcase "`case'"
    if inlist("`kind'","sas","spss") | (inlist("`kind'","xls","xlsx") & "`firstrow'"=="") local loadcase preserve
    capture noisily _cio_import_load `kind' `nobs' `nvar' `nlabels', case(`loadcase') wanted(`"`wanted'"') allstringfmt(`"`allstringfmt'"') `novallabels' `allstring' filename(`"`using'"')
    local rc = _rc
    capture plugin call ctools_plugin, "cio clear"
    macro drop CTOOLS_CIO_*
    if !`rc' {
        capture noisily {
            if `"`selectors'"' != "" quietly keep `selectors'
            if `"`wanted'"' != "" & !inlist("`kind'","xls","xlsx") keep `wanted'
            if inlist("`kind'","sas","spss") {
                if "`case'" != "preserve" rename _all, `case'
                quietly compress
            }
            * Validation hook: fail as Break would once the data are loaded
            if "$CTOOLS_cimport_fault" == "load" error 1
        }
        local rc = _rc
    }
    if `rc' {
        nobreak {
            if "`snapframe'" != "" frame copy `snapframe' `c(frame)', replace
            else if `snap' restore
        }
        exit `rc'
    }
    if "`snapframe'" != "" frame drop `snapframe'
    else if `snap' restore, not
    return scalar N = _N
    return scalar k = c(k)
    return scalar width = c(width)
    return scalar changed = 1
    local datalabel : data label
    return local datalabel : copy local datalabel
    return local filename `"`using'"'
end

program define _cio_import_load
    version 14.1
    gettoken kind 0 : 0, parse(" ,")
    gettoken nobs 0 : 0, parse(" ,")
    gettoken nvar 0 : 0, parse(" ,")
    gettoken nlabels 0 : 0, parse(" ,")
    syntax, CASE(string) [WANTED(string) NOVALLabels ALLSTRINGFMT(string) ALLString FILEname(string)]
    * Plugin registrations belong to the invoking ado program.
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    clear
    * The file's data label. label data cannot quote one holding an
    * unbalanced compound quote (it fails or stores other text); that label
    * comes from an empty dataset written to hold it. Later steps keep it.
    mata: st_local("__data_label",st_global("CTOOLS_CIO_datalabel"))
    capture label data `"`macval(__data_label)'"'
    local __label_rc = _rc
    local __label_got : data label
    mata: st_local("__label_rc",strofreal(`__label_rc' | st_local("__label_got")!=usubstr(st_local("__data_label"),1,80)))
    if `__label_rc' {
        tempfile __label_stub
        mata: _cio_label_stub(st_local("__label_stub"),st_local("__data_label"))
        quietly use `"`__label_stub'"', clear
    }
    quietly set obs `=max(1,`nobs')'
    mata: __ctools_cio_sets = J(1,0,"")
    if `nvar' {
    if inlist("`kind'","xls","xlsx") {
        * Resolve the entire header before creating variables: a later column's
        * fallback can require moving an earlier header to its own column letter.
        forvalues j=1/`nvar' {
            plugin call ctools_plugin, "cio column `j'"
            mata: st_local("__cio_candidate_`j'",_cio_excel_name(st_global("CTOOLS_CIO_name"),st_global("CTOOLS_CIO_excel_column"),st_local("case")))
            mata: st_local("__cio_fallback_`j'",st_global("CTOOLS_CIO_excel_column"))
        }
        mata: _cio_excel_names(`nvar')
    }
    forvalues j=1/`nvar' {
        plugin call ctools_plugin, "cio column `j'"
        if inlist("`kind'","xls","xlsx") {
            local vname "`__cio_resolved_`j''"
        }
        else mata: st_local("vname",ustrtoname(st_global("CTOOLS_CIO_name")))
        if "`case'" == "lower" & !inlist("`kind'","xls","xlsx") local vname = ustrlower("`vname'")
        if "`case'" == "upper" & !inlist("`kind'","xls","xlsx") local vname = ustrupper("`vname'")
        if inlist("`kind'","xls","xlsx") {
            if "`case'" == "lower" local vname = ustrlower("`vname'")
            if "`case'" == "upper" local vname = ustrupper("`vname'")
        }
        if !inlist("`kind'","xls","xlsx") {
            capture confirm new variable `vname'
            if _rc local vname "v`j'"
        }
        if $CTOOLS_CIO_external local external `external' `j'
        if $CTOOLS_CIO_numerictext local textnumbers `textnumbers' `j'
        local sourcecol`j' "${CTOOLS_CIO_excel_column}"
        local sourcevar`j' "`vname'"
        local type ${CTOOLS_CIO_type}
        if inlist("`kind'","xls","xlsx") & ("`vname'"=="_cons" | ustrregexm("`vname'","^str[0-9]+$")) {
            di as error "`vname' invalid varname"
            exit 198
        }
        mata: (void) st_addvar(st_local("type"),st_local("vname"))
        mata: st_varlabel(st_local("vname"),st_global("CTOOLS_CIO_labelnumber")=="" ? st_global("CTOOLS_CIO_label") : strtrim(strofreal(strtoreal(st_global("CTOOLS_CIO_labelnumber")),"%18.0g")))
        mata: st_varformat(st_local("vname"),st_global("CTOOLS_CIO_format"))
        if ("`novallabels'" == "" | "`kind'"=="sasxport5") & inlist("`kind'","spss","sasxport5") {
            * Label-set names come from the file and stay in Mata. Native
            * SPSS labels use their dictionary ordinal, beginning at 0.
            mata: st_local("labelname",_cio_labelset_name(st_local("kind"),st_global("CTOOLS_CIO_labelset"),1))
            * String variables take no value label: native SPSS notes it,
            * native XPORT fails unless novallabels is given.
            if "`labelname'" != "" & substr("`type'",1,3) == "str" {
                if "`kind'" == "sasxport5" & "`novallabels'" == "" {
                    mata: printf("{err}%s has string variable %s with format (value label) %s; value labels for string variables are not allowed in Stata\n", st_local("filename"), st_local("vname"), st_global("CTOOLS_CIO_labelset"))
                    exit 610
                }
                if "`kind'" == "spss" di as text "note: may not attach value label to string variable: no label applied"
                local labelname ""
            }
            if "`labelname'" != "" mata: st_varvaluelabel(st_local("vname"),st_local("labelname"))
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
    * Numbers in string columns: allstring()'s format, else native's %18.0g
    * (for which the plugin reserved 18 characters).
    local textfmt `"`allstringfmt'"'
    if `"`textfmt'"' == "" local textfmt "%18.0g"
    if `nobs' {
        foreach j of local textnumbers {
            tempvar rawnumber
            quietly gen double `rawnumber' = .
            plugin call ctools_plugin `rawnumber', "cio textnumbers `j'"
            local target `sourcevar`j''
            mata: _cio_format_text(st_local("target"),st_local("rawnumber"),st_local("textfmt"),st_local("allstringfmt")=="")
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
            mata: st_local("labelname",_cio_labelset_name(st_local("kind"),st_global("CTOOLS_CIO_labelset"),0))
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
    capture mata: mata drop __ctools_cio_sets
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
/* Lowercase hex pairs to bytes, decoded as whole vectors. */
string scalar _cio_hex_bytes(string scalar hex) {
    real matrix digits
    if (strlen(hex) < 2) return("")
    digits = ascii(hex)
    digits = colshape(digits :- 48 :- 39 :* (digits :> 57), 2)
    return(char((digits[., 1] * 16 + digits[., 2])'))
}
end

mata:
/* An empty .dta (format 118, no variables or observations) whose header
   holds label, the first 80 characters as label data keeps them: using it
   sets a data label that label data cannot quote (an unbalanced compound
   quote). */
void _cio_label_stub(string scalar path, string scalar label) {
    string scalar head, body, map, text
    real colvector at
    real scalar fh, i
    string rowvector tags
    text=usubstr(label,1,80)
    head="<stata_dta><header><release>118</release><byteorder>LSF</byteorder><K>"+_cio_le(0,2)+"</K><N>"+_cio_le(0,8)+"</N><label>"+_cio_le(strlen(text),2)+text+"</label><timestamp>"+char(0)+"</timestamp></header>"
    tags=("variable_types","varnames","sortlist","formats","value_label_names","variable_labels","characteristics","data","strls","value_labels")
    at=J(14,1,0)
    at[2]=strlen(head)
    body=""
    for (i=1;i<=10;i++) {
        at[i+2]=at[2]+strlen("<map>")+14*8+strlen("</map>")+strlen(body)
        body=body+"<"+tags[i]+">"+(tags[i]=="sortlist" ? _cio_le(0,2) : "")+"</"+tags[i]+">"
    }
    at[13]=at[2]+strlen("<map>")+14*8+strlen("</map>")+strlen(body)
    at[14]=at[13]+strlen("</stata_dta>")
    map="<map>"
    for (i=1;i<=14;i++) map=map+_cio_le(at[i],8)
    map=map+"</map>"
    if (fileexists(path)) unlink(path)
    fh=fopen(path,"w")
    fwrite(fh,head+map+body+"</stata_dta>")
    fclose(fh)
}
string scalar _cio_le(real scalar value, real scalar bytes) {
    string scalar out
    real scalar i, rest
    out=""
    rest=value
    for (i=1;i<=bytes;i++) {
        out=out+char(mod(rest,256))
        rest=floor(rest/256)
    }
    return(out)
}
void _cio_format_text(string scalar target, string scalar source, string scalar fmt, real scalar trim) {
    real colvector values, selected
    string colvector text
    real scalar width
    values=st_data(.,source)
    selected=selectindex(values:<.)
    if (!length(selected)) return
    text=trim ? strtrim(strofreal(values[selected],fmt)) : strofreal(values[selected],fmt)
    /* Widen the variable only as far as the text needs (strL past 2045) */
    width=max(strlen(text))
    if (st_vartype(target)!="strL" & width>strtoreal(substr(st_vartype(target),4,.))) stata("quietly recast "+(width>2045 ? "strL" : "str"+strofreal(width))+" "+target)
    st_sstore(selected,target,text)
}
end

mata:
/* Native cuts the header to 32 characters, keeps Unicode letters, decimal
   digits and underscores, drops leading digits, and uses the column letter
   when the result, in the requested case, is not a valid name. */
string scalar _cio_excel_name(string scalar source, string scalar fallback, string scalar casing) {
    string scalar name, cased
    name=ustrregexra(usubstr(source,1,32),"[^\p{L}\p{Nd}_]","")
    name=ustrregexrf(name,"^\p{Nd}+","")
    cased=casing=="lower" ? ustrlower(name) : (casing=="upper" ? ustrupper(name) : name)
    /* Native does not replace _cons or str#: creating them fails. */
    if (cased=="_cons" || ustrregexm(cased,"^str[0-9]+$")) return(name)
    if (name=="" || !st_isvarname(cased)) return(fallback)
    return(name)
}
/* Value-label name for a label set named in the file: for SPSS the reader's
   labels# (the set's dictionary ordinal, as native), the lowercased name for
   the other formats. "" when a variable has no label set (add=1). */
string scalar _cio_labelset_name(string scalar kind, string scalar set, real scalar add) {
    if (add & set=="") return("")
    if (kind!="spss") return(ustrlower(ustrtoname(set)))
    return(set)
}
void _cio_excel_names(real scalar n) {
    string rowvector names, fallback
    string scalar name
    real scalar j, at, previous
    transmorphic scalar used
    names=J(1,n,""); fallback=J(1,n,"")
    used=asarray_create("string",1)
    for(j=1;j<=n;j++) fallback[j]=st_local(sprintf("__cio_fallback_%g",j))
    for(j=1;j<=n;j++) {
        name=st_local(sprintf("__cio_candidate_%g",j))
        if(asarray_contains(used,name)) name=fallback[j]
        at=j
        while(asarray_contains(used,name)) {
            previous=asarray(used,name)
            names[at]=name
            asarray(used,name,at)
            at=previous
            name=fallback[at]
        }
        names[at]=name
        asarray(used,name,at)
    }
    for(j=1;j<=n;j++) st_local(sprintf("__cio_resolved_%g",j),names[j])
}
end
