/* Exact values (including binary strL), variable and dataset metadata, and
   label definitions. Native commands are reference fixtures, never fallback. */
clear all
set more off
set linesize 255
adopath ++ "build"
args pluginpath
if `"`pluginpath'"' == "" local pluginpath : environment CTOOLS_TEST_PLUGIN_DIR
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
        if(st_nobs()) {
            if(st_isnumvar(j)) s.values[,j]=strofreal(st_data(.,j),"%21x")
            else s.values[,j]=st_sdata(.,j)
        }
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
void cio_invalid_numeric_formats(struct cio_snapshot scalar s) {
    real scalar j
    for(j=1;j<=rows(s.metadata);j++) {
        if(substr(s.metadata[j,2],1,3)!="str" && !st_isnumfmt(s.metadata[j,3]))
            s.metadata[j,3]="%10.0g"
    }
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
    syntax using/, KIND(string) NAME(string) [OPTS(string asis) SELECT(string asis) LONGSTRING NATIVEBADFMT]
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
    if "`kind'"=="delimited" {
        local returned_encoding `"`r(encoding)'"'
        local returned_delimiters `"`r(delimiters)'"'
        local returned_N=r(N)
        local returned_k=r(k)
    }
    if inlist("`kind'","sas","spss","sasxport8","dbase") {
        local returned_N=r(N)
        local returned_k=r(k)
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
    if "`nativebadfmt'"!="" {
        mata: cio_invalid_numeric_formats(cio_ref)
        di "KNOWN GAP: `name': native assigns a string display format to numeric data; C uses %10.0g"
    }
    capture noisily cimport `kind' `select' using `"`filepath'"', clear `opts'
    local rc=_rc
    if !`rc' {
        if inlist("`kind'","sas","spss","sasxport8","dbase") {
            if r(N)!=`returned_N' | r(k)!=`returned_k' local rc=9
        }
        if "`kind'"=="delimited" {
            if `"`r(encoding)'"'!=`"`returned_encoding'"' | `"`r(delimiters)'"'!=`"`returned_delimiters'"' | r(N)!=`returned_N' | r(k)!=`returned_k' {
                di "DETAIL: native encoding=`returned_encoding', C encoding=`r(encoding)'"
                di `"DETAIL: native delimiters=[`macval(returned_delimiters)'], C delimiters=[`r(delimiters)']"'
                local rc=9
            }
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

program define cio_import_state_test
    syntax using/, KIND(string) STATE(string) [OPTS(string asis)]
    local filepath `"`using'"'
    local using ""
    local clause using
    if "`kind'"=="sasxport5" local clause ""
    tempfile initial
    clear
    if "`state'"!="empty" {
        quietly set obs 2
        quietly gen long sentinel=99
        label variable sentinel "Original data"
        label data "Original dataset"
        quietly save `initial', replace
        if "`state'"=="changed" quietly replace sentinel=100 in 1
        if "`state'"=="zeroobs" quietly drop in 1/2
        if "`state'"=="zerovars" drop sentinel
    }
    preserve
    capture noisily import `kind' `clause' "`filepath'", `opts'
    local nr=_rc
    local changed=c(changed)
    local datalabel : data label
    local sortedby : sortedby
    quietly label dir
    local labelsets `r(names)'
    mata: cio_ref=cio_snapshot_data()
    restore
    capture noisily cimport `kind' using "`filepath'", `opts'
    local cr=_rc
    local rc=cond(`nr'==`cr' & `changed'==c(changed),0,9)
    local datalabel : data label
    local sortedby : sortedby
    quietly label dir
    local labelsets `r(names)'
    mata: cio_compare(cio_ref)
    if !`same' local rc=9
    cio_record "`kind' `filepath' memory state `state' native rc=`nr' C rc=`cr'" `rc'
end

program define cio_companion_state_test
    syntax using/, COMPANION(string) [OPTS(string asis)]
    local filepath `"`using'"'
    local using ""
    clear
    quietly set obs 2
    quietly gen long sentinel=99
    label data "Caller data"
    capture noisily import sasxport8 using `"`filepath'"', vlabfile(`"`companion'"') clear `opts'
    local nr=_rc
    local rn=r(N)
    local rk=r(k)
    local rw=r(width)
    local rchanged=r(changed)
    local changed=c(changed)
    local datalabel : data label
    local sortedby : sortedby
    quietly label dir
    local labelsets `r(names)'
    mata: cio_ref=cio_snapshot_data()
    clear
    quietly set obs 2
    quietly gen long sentinel=99
    label data "Caller data"
    capture noisily cimport sasxport8 using `"`filepath'"', vlabfile(`"`companion'"') clear `opts'
    local cr=_rc
    local rc=cond(`nr'==`cr' & `changed'==c(changed),0,9)
    if r(N)!=`rn' | r(k)!=`rk' | r(width)!=`rw' | r(changed)!=`rchanged' local rc=9
    local datalabel : data label
    local sortedby : sortedby
    quietly label dir
    local labelsets `r(names)'
    mata: cio_compare(cio_ref)
    if !`same' local rc=9
    cio_record "XPORT8 companion partial state `companion' `opts' native rc=`nr' C rc=`cr'" `rc'
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

/* Cell formats of a workbook as Excel reads them: one sorted entry per cell,
   "sheet|row|column|format|kind", where format is a built-in number format id
   or a custom format code and kind is num, text, blank or other. XLSX parts
   are read after unzipfile; XLS through its compound file and BIFF8 records. */
mata:
string scalar cio_wb_text(string scalar path) {
    real scalar fh
    string scalar s, chunk
    fh = fopen(path, "r")
    s = ""
    while ((chunk = fread(fh, 1048576)) != J(0, 0, "")) s = s + chunk
    fclose(fh)
    return(s)
}
string scalar cio_wb_attr(string scalar tag, string scalar name) {
    if (regexm(tag, "[ :]" + name + `"="([^"]*)""')) return(regexs(1))
    return("")
}
string colvector cio_wb_tags(string scalar xml, string scalar name) {
    string colvector tags
    string scalar next
    real scalar p, q
    tags = J(0, 1, "")
    while ((p = strpos(xml, "<" + name))) {
        next = substr(xml, p + strlen(name) + 1, 1)
        xml = substr(xml, p + 1, .)
        if (next != " " & next != "/" & next != ">") continue
        q = strpos(xml, ">")
        tags = tags \ "<" + substr(xml, 1, q)
        xml = substr(xml, q + 1, .)
    }
    return(tags)
}
string scalar cio_wb_format(real scalar id, string colvector ids, string colvector codes) {
    real scalar i
    for (i = 1; i <= rows(ids); i++) if (strtoreal(ids[i]) == id) return("code:" + codes[i])
    return("id:" + strofreal(id))
}
string colvector cio_xlsx_styles(string scalar dir) {
    string scalar styles, book, rels, sheet, tag, ref, letters, kind, path
    string colvector fmts, fids, fcodes, xfs, sheets, relations, cells, out
    real scalar i, j, k, col, row, xf
    styles = cio_wb_text(dir + "/xl/styles.xml")
    fmts = cio_wb_tags(styles, "numFmt")
    fids = J(rows(fmts), 1, "")
    fcodes = fids
    for (i = 1; i <= rows(fmts); i++) {
        fids[i] = cio_wb_attr(fmts[i], "numFmtId")
        fcodes[i] = cio_wb_attr(fmts[i], "formatCode")
    }
    xfs = cio_wb_tags(substr(styles, strpos(styles, "<cellXfs"), strpos(styles, "</cellXfs>") - strpos(styles, "<cellXfs")), "xf")
    book = cio_wb_text(dir + "/xl/workbook.xml")
    rels = cio_wb_text(dir + "/xl/_rels/workbook.xml.rels")
    sheets = cio_wb_tags(book, "sheet")
    relations = cio_wb_tags(rels, "Relationship")
    out = J(0, 1, "")
    for (i = 1; i <= rows(sheets); i++) {
        path = ""
        for (j = 1; j <= rows(relations); j++) {
            if (cio_wb_attr(relations[j], "Id") == cio_wb_attr(sheets[i], "r:id")) path = cio_wb_attr(relations[j], "Target")
        }
        path = substr(path, 1, 1) == "/" ? dir + path : dir + "/xl/" + path
        sheet = cio_wb_text(path)
        cells = cio_wb_tags(sheet, "c")
        for (k = 1; k <= rows(cells); k++) {
            tag = cells[k]
            ref = cio_wb_attr(tag, "r")
            letters = regexr(ref, "[0-9]+$", "")
            row = strtoreal(substr(ref, strlen(letters) + 1, .))
            col = 0
            for (j = 1; j <= strlen(letters); j++) col = col * 26 + ascii(substr(letters, j, 1)) - 64
            xf = strtoreal(cio_wb_attr(tag, "s"))
            if (xf == .) xf = 0
            kind = cio_wb_attr(tag, "t")
            kind = substr(tag, -2, 2) == "/>" ? "blank" : (anyof(("s", "inlineStr", "str"), kind) ? "text" : (kind == "" | kind == "n" ? "num" : "other"))
            out = out \ cio_wb_attr(sheets[i], "name") + "|" + strofreal(row) + "|" + strofreal(col) + "|" + cio_wb_format(strtoreal(cio_wb_attr(xfs[xf + 1], "numFmtId")) == . ? 0 : strtoreal(cio_wb_attr(xfs[xf + 1], "numFmtId")), fids, fcodes) + "|" + kind
        }
    }
    return(sort(out, 1))
}
string colvector cio_wb_diff(string colvector a, string colvector b) {
    string colvector d
    real scalar i
    d = J(0, 1, "")
    for (i = 1; i <= rows(a); i++) if (!anyof(b, a[i])) d = d \ "native only: " + a[i]
    for (i = 1; i <= rows(b); i++) if (!anyof(a, b[i])) d = d \ "cexport only: " + b[i]
    return(rows(d) > 12 ? d[|1 \ 12|] : d)
}
real scalar cio_wb_u(real rowvector b, real scalar p, real scalar n) {
    real scalar v, i
    v = 0
    for (i = n; i >= 1; i--) v = v * 256 + b[p + i]
    return(v)
}
real rowvector cio_wb_chain(real rowvector b, real rowvector table, real scalar start, real scalar size, real scalar unit, real scalar base) {
    real rowvector s
    real scalar sector
    s = J(1, 0, .)
    sector = start
    while (sector < 4294967290 & (size < 0 | cols(s) < size)) {
        s = s, b[|base + sector * unit + 1 \ base + (sector + 1) * unit|]
        sector = table[sector + 1]
    }
    if (size >= 0 & cols(s) > size) s = s[|1 \ size|]
    return(s)
}
string colvector cio_xls_styles(string scalar path) {
    real rowvector b, fat, fatsectors, dir, mini, minifat, wb, xfs
    real scalar ss, ms, i, j, p, id, len, start, size, cutoff, sector, row, col, ixfe, n, ifmt
    string scalar name, kind
    string colvector fids, fcodes, names, out
    real colvector offsets
    b = ascii(cio_wb_text(path))
    ss = 2 ^ cio_wb_u(b, 30, 2)
    ms = 2 ^ cio_wb_u(b, 32, 2)
    cutoff = cio_wb_u(b, 56, 4)
    fatsectors = J(1, 0, .)
    for (i = 0; i < 109 & cols(fatsectors) < cio_wb_u(b, 44, 4); i++) fatsectors = fatsectors, cio_wb_u(b, 76 + 4 * i, 4)
    sector = cio_wb_u(b, 68, 4)
    while (cols(fatsectors) < cio_wb_u(b, 44, 4) & sector < 4294967290) {
        for (i = 0; i < ss / 4 - 1 & cols(fatsectors) < cio_wb_u(b, 44, 4); i++) fatsectors = fatsectors, cio_wb_u(b, (sector + 1) * ss + 4 * i, 4)
        sector = cio_wb_u(b, (sector + 1) * ss + ss - 4, 4)
    }
    fat = J(1, 0, .)
    for (j = 1; j <= cols(fatsectors); j++) {
        for (i = 0; i < ss / 4; i++) fat = fat, cio_wb_u(b, (fatsectors[j] + 1) * ss + 4 * i, 4)
    }
    dir = cio_wb_chain(b, fat, cio_wb_u(b, 48, 4), -1, ss, ss)
    minifat = J(1, 0, .)
    if (cio_wb_u(b, 64, 4)) {
        p = 0
        mini = cio_wb_chain(b, fat, cio_wb_u(b, 60, 4), -1, ss, ss)
        for (i = 0; i < cols(mini) / 4; i++) minifat = minifat, cio_wb_u(mini, 4 * i, 4)
    }
    mini = cio_wb_chain(b, fat, cio_wb_u(dir, 116, 4), cio_wb_u(dir, 120, 4), ss, ss)
    wb = J(1, 0, .)
    for (i = 0; i < cols(dir) / 128; i++) {
        name = ""
        for (j = 0; j < (cio_wb_u(dir, 128 * i + 64, 2) - 2) / 2; j++) name = name + char(cio_wb_u(dir, 128 * i + 2 * j, 2))
        if (name == "Workbook" | name == "Book") {
            start = cio_wb_u(dir, 128 * i + 116, 4)
            size = cio_wb_u(dir, 128 * i + 120, 4)
            wb = size < cutoff ? cio_wb_chain(mini, minifat, start, size, ms, 0) : cio_wb_chain(b, fat, start, size, ss, ss)
        }
    }
    fids = fcodes = names = J(0, 1, "")
    offsets = J(0, 1, .)
    xfs = J(1, 0, .)
    p = 0
    while (p + 4 <= cols(wb)) {
        id = cio_wb_u(wb, p, 2)
        len = cio_wb_u(wb, p + 2, 2)
        if (id == 1054) {
            n = cio_wb_u(wb, p + 6, 2)
            name = ""
            for (j = 0; j < n; j++) name = name + (mod(wb[p + 9], 2) ? uchar(cio_wb_u(wb, p + 9 + 2 * j, 2)) : char(wb[p + 10 + j]))
            fids = fids \ strofreal(cio_wb_u(wb, p + 4, 2))
            fcodes = fcodes \ name
        }
        if (id == 224) xfs = xfs, cio_wb_u(wb, p + 6, 2)
        if (id == 133) {
            offsets = offsets \ cio_wb_u(wb, p + 4, 4)
            name = ""
            for (j = 0; j < wb[p + 11]; j++) name = name + (mod(wb[p + 12], 2) ? uchar(cio_wb_u(wb, p + 12 + 2 * j, 2)) : char(wb[p + 13 + j]))
            names = names \ name
        }
        if (id == 10) break
        p = p + 4 + len
    }
    out = J(0, 1, "")
    for (i = 1; i <= rows(offsets); i++) {
        p = offsets[i]
        while (p + 4 <= cols(wb)) {
            id = cio_wb_u(wb, p, 2)
            len = cio_wb_u(wb, p + 2, 2)
            if (id == 10) break
            row = cio_wb_u(wb, p + 4, 2) + 1
            if (anyof((515, 638, 253, 516, 513, 517, 6), id)) {
                kind = anyof((515, 638), id) ? "num" : (anyof((253, 516), id) ? "text" : (id == 513 ? "blank" : "other"))
                ifmt = xfs[cio_wb_u(wb, p + 8, 2) + 1]
                out = out \ names[i] + "|" + strofreal(row) + "|" + strofreal(cio_wb_u(wb, p + 6, 2) + 1) + "|" + cio_wb_format(ifmt, ifmt >= 164 ? fids : J(0, 1, ""), fcodes) + "|" + kind
            }
            if (id == 189 | id == 190) {
                n = id == 189 ? (len - 6) / 6 : (len - 6) / 2
                for (j = 0; j < n; j++) {
                    ifmt = xfs[cio_wb_u(wb, p + 8 + j * (id == 189 ? 6 : 2), 2) + 1]
                    out = out \ names[i] + "|" + strofreal(row) + "|" + strofreal(cio_wb_u(wb, p + 6, 2) + 1 + j) + "|" + cio_wb_format(ifmt, ifmt >= 164 ? fids : J(0, 1, ""), fcodes) + "|" + (id == 189 ? "num" : "blank")
                }
            }
            p = p + 4 + len
        }
    }
    return(sort(out, 1))
}
end

/* cio_excel_styles using file, scratch(dir): cell formats into Mata
   cio_excel_styles_out (see cio_xlsx_styles). */
program define cio_excel_styles
    syntax using/, SCRATCH(string)
    mata: st_local("using", pathisabs(st_local("using")) ? st_local("using") : pathjoin(pwd(), st_local("using")))
    if lower(substr(`"`using'"', -5, .)) == ".xlsx" {
        capture mkdir `"`scratch'"'
        local here `"`c(pwd)'"'
        quietly cd `"`scratch'"'
        capture unzipfile `"`using'"', replace
        local rc = _rc
        quietly cd `"`here'"'
        if `rc' exit `rc'
        mata: cio_excel_styles_out = cio_xlsx_styles(st_local("scratch"))
    }
    else mata: cio_excel_styles_out = cio_xls_styles(st_local("using"))
end
