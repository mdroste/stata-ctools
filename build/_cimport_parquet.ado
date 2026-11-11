*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! cimport parquet: a native C Parquet reader with import parquet's syntax,
*! variables, messages and stored results.
program define _cimport_parquet, rclass
    version 14.1
    * As import parquet: a bare filename before the comma is the using file
    gettoken filename rest : 0, parse(" ,")
    gettoken comma : rest, parse(" ,")
    if `"`filename'"' != "" & (trim(`"`comma'"') == "," | trim(`"`comma'"') == "") {
        local 0 `"using `0'"'
    }
    local __desc 0
    capture syntax using/, DESCribe
    if !_rc local __desc 1
    else {
        syntax [anything(name=vlist id="vlist")] using/ [, clear ROWRange(string) ///
            favormemory SNAPshot FAST]
        if "`snapshot'" != "" & "`fast'" != "" {
            di as error "options snapshot and fast may not be combined"
            exit 198
        }
        if "`clear'" == "" & c(changed) error 4
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    * Inputs and results pass through CTOOLS_PQ_* globals, dropped at every exit
    macro drop CTOOLS_PQ_*
    local __pq_home : environment HOME
    if c(os) == "Windows" local __pq_home : environment USERPROFILE
    capture noisily mata: _cpq_path(st_local("using"))
    if _rc {
        local rc = _rc
        macro drop CTOOLS_PQ_*
        exit `rc'
    }
    if `__desc' {
        * Stata's describe decodes the data too; unreadable data show -1 rows
        capture noisily {
            _cimport_parquet_call describe
            mata: _cpq_describe()
        }
        local rc = _rc
        local k $CTOOLS_PQ_k
        local n $CTOOLS_PQ_N
        macro drop CTOOLS_PQ_*
        if `rc' exit `rc'
        return scalar n_columns = `k'
        return scalar N_rows = `n'
        exit
    }
    * Column selection, row range, supported types and variable names, in
    * the order import parquet checks them
    capture noisily {
        _cimport_parquet_call schema
        mata: _cpq_prepare(st_local("vlist"), st_local("rowrange"))
    }
    if _rc {
        local rc = _rc
        macro drop CTOOLS_PQ_*
        exit `rc'
    }
    local nout $CTOOLS_PQ_nout
    * The undo copy, as in _cio_import: a frame from Stata 16 when there are
    * data, otherwise preserve; snapshot always preserves; fast keeps none.
    local snap 0
    local snapframe ""
    if "`fast'" == "" {
        if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
            tempname snapframe
            capture noisily frame copy `c(frame)' `snapframe'
        }
        else {
            capture noisily preserve
            local snap 1
        }
        if _rc {
            local rc = _rc
            macro drop CTOOLS_PQ_*
            exit `rc'
        }
    }
    capture noisily {
        local nobs = $CTOOLS_PQ_last - $CTOOLS_PQ_first + 1
        if `nobs' > 2147483647 {
            di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
            di as error "This is a limitation of Stata's plugin API"
            exit 920
        }
        * As import parquet: an empty file, or an end row that is not an
        * integer, fails when the rows are loaded
        if $CTOOLS_PQ_badrows | `nobs' < 1 {
            di as error "out of memory"
            exit 3900
        }
        plugin call ctools_plugin, "cparquet scan"
        quietly drop _all
        quietly set obs `nobs'
        forvalues j = 1/`nout' {
            local type ${CTOOLS_PQ_vtype_`j'}
            local name ${CTOOLS_PQ_vname_`j'}
            * import parquet's own command, so that invalid names fail alike
            if substr("`type'", 1, 3) == "str" capture noisily quietly gen `type' `name' = ""
            else capture noisily quietly gen `type' `name' = .
            if _rc exit 3598
            if "${CTOOLS_PQ_fmt_`j'}" != "" format ${CTOOLS_PQ_fmt_`j'} `name'
        }
        mata: _cpq_labels(`nout')
        plugin call ctools_plugin *, "cparquet load"
        forvalues j = 1/`nout' {
            if ${CTOOLS_PQ_blob_`j'} {
                tempfile blobs
                global CTOOLS_PQ_blobpath `"`blobs'"'
                plugin call ctools_plugin, "cparquet blobfile `j'"
                mata: _cpq_read_blobs(st_global("CTOOLS_PQ_blobpath"), `j', `nobs')
                capture erase `"`blobs'"'
            }
        }
        * Dictionary (categorical) columns: value labels from their entries,
        * modifying a label of the same name as import parquet does
        forvalues j = 1/`nout' {
            if ${CTOOLS_PQ_dict_`j'} {
                tempfile entries
                global CTOOLS_PQ_blobpath `"`entries'"'
                plugin call ctools_plugin, "cparquet dictfile `j'"
                mata: _cpq_dict_labels(st_global("CTOOLS_PQ_blobpath"), st_global("CTOOLS_PQ_vname_`j'"), ${CTOOLS_PQ_dict_`j'})
                label values ${CTOOLS_PQ_vname_`j'} ${CTOOLS_PQ_vname_`j'}
                capture erase `"`entries'"'
            }
        }
        * Validation hook: fail as Break would once the data are loaded
        if "$CTOOLS_cimport_fault" == "load" error 1
    }
    local rc = _rc
    capture plugin call ctools_plugin, "cparquet clear"
    macro drop CTOOLS_PQ_*
    if `rc' {
        nobreak {
            if "`snapframe'" != "" frame copy `snapframe' `c(frame)', replace
            else if `snap' restore
        }
        exit `rc'
    }
    if "`snapframe'" != "" frame drop `snapframe'
    else if `snap' restore, not
    mata: printf("{txt}(%s, %s obs)\n", `nout' == 1 ? "1 var" : strofreal(`nout', "%16.0gc") + " vars", strofreal(st_nobs(), "%16.0gc"))
    return scalar k = `nout'
    return scalar N = _N
end

* Reads the footer; a file that cannot be opened fails as in import parquet.
program define _cimport_parquet_call
    version 14.1
    args op
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture noisily plugin call ctools_plugin, "cparquet `op'"
    local rc = _rc
    if `rc' == 603 di as error `"file {bf:${CTOOLS_PQ_hold}} could not be opened"'
    if `rc' exit `rc'
end

mata:
/* The file as import parquet finds it: no URLs, ~ for the home directory,
   relative to the working directory, .parquet added when there is no
   extension. CTOOLS_PQ_file is the absolute path. */
void _cpq_path(string scalar f)
{
    string scalar hold, home
    hold = f
    if (pathisurl(f)) {
        errprintf("invalid filename; filename may not be specified as a URL\n")
        exit(603)
    }
    if (substr(f, 1, 1) == "~" & (strlen(f) == 1 | anyof(("/", "\"), substr(f, 2, 1)))) {
        home = st_local("__pq_home")
        if (home != "") f = home + substr(f, 2, .)
    }
    if (!pathisabs(f)) f = _cpq_resolve(pwd(), f)
    if (fileexists(f) != 1) {
        if (pathsuffix(f) == "") f = f + ".parquet"
        if (fileexists(f) != 1) {
            errprintf("file {bf:%s} not found\n", hold)
            exit(601)
        }
    }
    st_global("CTOOLS_PQ_file", f)
    st_global("CTOOLS_PQ_hold", hold)
}

/* An absolute path from the working directory and a relative path, made as
   Stata 16's path resolution makes it (which Stata 14.1 lacks): forward
   slashes, no empty, . or .. components, no trailing slash. */
string scalar _cpq_resolve(string scalar base, string scalar f)
{
    string scalar s, prefix, part
    string rowvector stack
    real scalar at
    s = subinstr(base + "/" + f, "\", "/", .)
    prefix = ""
    if (substr(s, 1, 2) == "//") prefix = "//"
    else if (substr(s, 1, 1) == "/") prefix = "/"
    else if (regexm(s, "^[A-Za-z]:")) {
        prefix = substr(s, 1, 2) + "/"
        s = substr(s, 3, .)
    }
    stack = J(1, 0, "")
    while (s != "") {
        at = strpos(s, "/")
        if (at) {
            part = substr(s, 1, at - 1)
            s = substr(s, at + 1, .)
        }
        else {
            part = s
            s = ""
        }
        if (part == "" | part == ".") continue
        if (part == "..") {
            if (cols(stack)) stack = cols(stack) > 1 ? stack[1..cols(stack) - 1] : J(1, 0, "")
            continue
        }
        stack = (stack, part)
    }
    return(prefix + invtokens(stack, "/"))
}

void _cpq_describe()
{
    real scalar i, k
    k = strtoreal(st_global("CTOOLS_PQ_k"))
    printf("{txt}Contains data from {res}%s\n", st_global("CTOOLS_PQ_file"))
    printf("{txt} Observations:        {res}%s\n", strofreal(strtoreal(st_global("CTOOLS_PQ_N")), "%16.0gc"))
    printf("{txt}    Variables:        {res}%s\n", strofreal(k, "%16.0gc"))
    printf("{txt}{hline}\n")
    printf("Column                           Type\n")
    printf("{hline}\n")
    for (i = 1; i <= k; i++) {
        printf("{p 0 48}{res}{bind:%-32s }{txt}{res}{bind:%-15s}{p_end}", abbrev(st_global("CTOOLS_PQ_name_" + strofreal(i)), 32), st_global("CTOOLS_PQ_type_" + strofreal(i)))
    }
    printf("{txt}{hline}\n")
}

real scalar _cpq_range_error(string scalar why)
{
    errprintf("invalid range in option {bf:rowrange()}")
    errprintf("{p 4 4 2}%s{p_end}", why)
    return(198)
}

/* rowrange([start][:end]) in the forms import parquet accepts, with its
   messages: f and l name the first and last rows, and a negative number
   counts from the end. Returns (start, end), end missing for the last row. */
real rowvector _cpq_rowrange(string scalar r)
{
    real scalar colon, a, b
    string scalar left, right
    a = 1
    b = .
    colon = strpos(r, ":")
    if (colon > 1) {
        left = bsubstr(r, 1, colon - 1)
        a = anyof(("f", "F"), strtrim(left)) ? 1 : strtoreal(left)
        if (a == . | strpos(left, ".") > 0) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the starting row."))
        right = bsubstr(r, colon + 1, .)
        if (!anyof(("l", "L", ""), strtrim(right))) {
            b = strtoreal(right)
            if (b == .) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the ending row."))
        }
    }
    else if (colon == 1) {
        right = bsubstr(r, 2, .)
        b = strtoreal(right)
        if (b == . | strpos(right, ".") > 0) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the ending row."))
    }
    else {
        a = anyof(("f", "F"), strtrim(r)) ? 1 : strtoreal(r)
        if (a == . | strpos(r, ".") > 0) exit(_cpq_range_error("The option {bf:rowrange()} requires an integer for the starting row."))
    }
    return((a, b))
}

void _cpq_prepare(string scalar list, string scalar rowrange)
{
    real scalar k, nrows, i, j, first, last, found
    string colvector names, types, invalid
    string rowvector wanted, vnames, orig
    real rowvector sel, keep, range
    k = strtoreal(st_global("CTOOLS_PQ_k"))
    nrows = strtoreal(st_global("CTOOLS_PQ_N"))
    names = J(k, 1, "")
    types = J(k, 1, "")
    for (i = 1; i <= k; i++) {
        names[i] = st_global("CTOOLS_PQ_name_" + strofreal(i))
        types[i] = st_global("CTOOLS_PQ_type_" + strofreal(i))
    }
    /* Selection by exact name, the first column of that name; a name given
       twice is read once */
    sel = J(1, 0, .)
    if (strtrim(list) != "") {
        wanted = tokens(list)
        invalid = J(0, 1, "")
        for (i = 1; i <= cols(wanted); i++) {
            found = 0
            for (j = 1; j <= k; j++) {
                if (wanted[i] == names[j]) {
                    if (!anyof(sel, j)) sel = (sel, j)
                    found = 1
                    break
                }
            }
            if (!found) invalid = (invalid \ wanted[i])
        }
        if (cols(sel) == 0) {
            errprintf("invalid column names\n")
            exit(198)
        }
        for (i = 1; i <= rows(invalid); i++) printf("{txt}note: column {bf:%s} ignored because of invalid column name.\n", invalid[i])
    }
    else {
        if (k == 0) {
            errprintf("file contains 0 column\n")
            exit(198)
        }
        sel = 1..k
    }
    first = 1
    last = nrows
    if (strtrim(rowrange) != "") {
        range = _cpq_rowrange(rowrange)
        first = range[1]
        last = range[2]
        if (first < 0) first = nrows + first + 1
        if (last == .) last = nrows
        if (last < 0) last = nrows + last + 1
        if (first < 1 | first > nrows | last < 1 | last > nrows) exit(_cpq_range_error("The specified integers are outside of the valid range."))
        if (first > last) exit(_cpq_range_error("The ending row must be greater than the starting row."))
    }
    /* Types Stata cannot hold are skipped with a note */
    keep = J(1, 0, .)
    invalid = J(0, 1, "")
    for (i = 1; i <= cols(sel); i++) {
        if (anyof(("null", "struct", "list", "map", "duration"), types[sel[i]])) invalid = (invalid \ names[sel[i]])
        else keep = (keep, sel[i])
    }
    if (cols(keep) == 0) {
        errprintf("all columns are invalid due to unsupported data types\n")
        exit(198)
    }
    for (i = 1; i <= rows(invalid); i++) printf("{txt}note: column {bf:%s} ignored because of unsupported Parquet type.\n", invalid[i])
    orig = names[keep]'
    vnames = _cpq_varnames(orig)
    st_global("CTOOLS_PQ_cols", invtokens(strofreal(keep :- 1)))
    st_global("CTOOLS_PQ_first", strofreal(first, "%21.0g"))
    st_global("CTOOLS_PQ_last", strofreal(floor(last), "%21.0g"))
    st_global("CTOOLS_PQ_badrows", strofreal(last != floor(last)))
    st_global("CTOOLS_PQ_nout", strofreal(cols(keep)))
    for (i = 1; i <= cols(keep); i++) {
        st_global("CTOOLS_PQ_vname_" + strofreal(i), vnames[i])
        st_global("CTOOLS_PQ_orig_" + strofreal(i), orig[i])
    }
}

/* Variable names as import parquet makes them: letters, decimal digits and
   underscores are kept, leading digits are dropped, and the result is cut
   to 32 characters; an empty name becomes v#, # the column's 0-based
   position among the imported ones; and if any two names coincide, every
   column is named that way. */
string rowvector _cpq_varnames(string rowvector names)
{
    real scalar i, j
    string rowvector out
    out = J(1, cols(names), "")
    for (i = 1; i <= cols(names); i++) {
        out[i] = ustrregexra(names[i], "[^\p{L}\p{Nd}_]", "")
        out[i] = ustrregexrf(out[i], "^\p{Nd}+", "")
        out[i] = usubstr(out[i], 1, 32)
        if (out[i] == "") out[i] = "v" + strofreal(i - 1)
    }
    for (i = 1; i <= cols(out); i++) {
        for (j = i + 1; j <= cols(out); j++) {
            if (out[i] == out[j]) return("v" :+ strofreal(0..cols(out) - 1))
        }
    }
    return(out)
}

/* A variable whose name differs from its column is labeled with the column */
void _cpq_labels(real scalar n)
{
    real scalar i
    string scalar name, orig
    for (i = 1; i <= n; i++) {
        name = st_global("CTOOLS_PQ_vname_" + strofreal(i))
        orig = st_global("CTOOLS_PQ_orig_" + strofreal(i))
        if (name != orig) st_varlabel(name, usubstr(orig, 1, 80))
    }
}

void _cpq_dict_labels(string scalar path, string scalar name, real scalar n)
{
    real scalar fh, i, length
    real colvector io
    string colvector text
    fh = fopen(path, "r")
    io = bufio()
    text = J(n, 1, "")
    for (i = 1; i <= n; i++) {
        length = fbufget(io, fh, "%4bu")
        if (length) text[i] = fread(fh, length)
    }
    fclose(fh)
    st_vlmodify(name, (0::n - 1), text)
}

void _cpq_read_blobs(string scalar path, real scalar col, real scalar n)
{
    real scalar fh, i, length
    real colvector io
    fh = fopen(path, "r")
    io = bufio()
    for (i = 1; i <= n; i++) {
        length = fbufget(io, fh, "%4bu")
        st_sstore(i, col, length ? fread(fh, length) : "")
    }
    fclose(fh)
}
end
