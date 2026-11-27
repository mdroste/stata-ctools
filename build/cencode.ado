*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! original: Stata's encode (built in)
*! original: see help cencode for details

program define cencode, rclass
    version 14.1
    * Outputs are built under temporary names and committed only after every
    * variable is encoded, so a failure leaves the variables untouched. Value
    * labels are the only other state cencode changes: the labels it may write
    * are saved exactly (through Mata) and put back if anything fails.
    * encode leaves the caller's r() in place, so cencode keeps it too: it
    * is held while cencode works and returned along with r(N_unique) and
    * r(N_vars). It is in place while the if is evaluated, as it may name it.
    tempname __caller_r
    _return hold `__caller_r'
    _return restore `__caller_r', hold
    syntax varlist [if] [in], [Generate(string) replace Label(name) *]
    local __labels "`label'"
    if "`__labels'" == "" {
        local __labels = cond("`replace'" != "", "`varlist'", "`generate'")
    }
    capture mata: _cencode_labels_save(tokens(st_local("__labels")))
    local __backup_rc = _rc
    if `__backup_rc' {
        capture mata: mata drop __cencode_vl_saved __cencode_vl_names
        _return restore `__caller_r'
        return add
        exit `__backup_rc'
    }
    capture noisily _cencode_impl `0'
    local rc = _rc
    if `rc' {
        capture mata: _cencode_labels_restore()
        capture mata: mata drop __cencode_vl_saved __cencode_vl_names
        _return restore `__caller_r'
        return add
        exit `rc'
    }
    capture mata: mata drop __cencode_vl_saved __cencode_vl_names
    local __n_unique = r(N_unique)
    local __n_vars = r(N_vars)
    _return restore `__caller_r'
    return add
    return scalar N_unique = `__n_unique'
    return scalar N_vars = `__n_vars'
end

program define _cencode_impl, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    syntax varlist [if] [in], [Generate(string) replace Label(name) NOExtend Verbose THReads(integer 0)]
    if "`noextend'" != "" & "`label'" == "" {
        di as error "noextend option not allowed without label() option"
        exit 198
    }
    * The if is evaluated now, as encode does, before any helper replaces
    * the r() or e() results it may name.
    if `"`if'"' != "" {
        tempvar __iftouse
        mark `__iftouse' `if' `in'
        local if "if `__iftouse'"
    }

    * =========================================================================
    * UPFRONT VALIDATION
    * =========================================================================

    * Check that either generate or replace is specified, but not both
    if "`generate'" == "" & "`replace'" == "" {
        di as error "cencode: must specify either generate() or replace option"
        exit 198
    }
    if "`generate'" != "" & "`replace'" != "" {
        di as error "cencode: cannot specify both generate() and replace options"
        exit 198
    }

    * Count input variables
    local n_vars : word count `varlist'
    unab original_order : _all

    * Handle replace vs generate
    local __do_replace = 0
    if "`replace'" != "" {
        local __do_replace = 1
    }
    else {
        * Check generate list has same number of variables
        local n_gen : word count `generate'
        if `n_gen' != `n_vars' {
            di as error "cencode: generate() must specify `n_vars' variable(s) to match varlist"
            exit 198
        }
        _ctools_newvars `generate'
    }

    * Check that all source variables are string (encode's message, r(107))
    foreach v of local varlist {
        capture confirm string variable `v'
        if _rc != 0 {
            di as error "not possible with numeric variable"
            exit 107
        }
    }

    if _N == 0 {
        local staged ""
        forvalues i=1/`n_vars' {
            local srcvar : word `i' of `varlist'
            tempvar target
            local final_name = cond(`__do_replace', "`srcvar'", "`: word `i' of `generate''")
            local this_label "`label'"
            if "`this_label'"=="" local this_label "`final_name'"
            quietly encode `srcvar' `if' `in', generate(`target') label(`this_label') `noextend'
            local staged `staged' `target'
        }
        _cencode_commit "`varlist'" "`staged'" "`generate'" `__do_replace' "`original_order'"
        return scalar N_unique = 0
        return scalar N_vars = `n_vars'
        exit
    }

    * =========================================================================
    * END VALIDATION
    * =========================================================================

    * Start timing
    local __do_timing = ("`verbose'" != "")
    if `__do_timing' {
        local __ct_elapsed90 = 0
        local __ct_elapsed91 = 0
        local __ct_elapsed92 = 0
        local __ct_elapsed93 = 0
        if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
        quietly _ctools_clock __ct_start90
        if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
        quietly _ctools_clock __ct_start91
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0
    * Reset _rc (plugin load may leave _rc=110 for already-loaded plugin)
    capture confirm number 0

    * Build threads option
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Without if/in the plugin writes every observation, so the destination
    * can be created without first filling it with missing values.
    local __fill = cond(`"`if'`in'"' == "", "fillall", "")
    local __ifflag = cond(`"`if'"' != "", "ifcond", "")
    * An if was marked once after syntax; the marker is passed last to each
    * plugin call, which reads it in parallel (Stata would answer the if one
    * observation at a time).

    * =========================================================================
    * Process each variable
    * =========================================================================

    local staged ""
    forvalues i = 1/`n_vars' {
        local srcvar : word `i' of `varlist'

        * Every output is built under a temporary name; see _cencode_commit.
        tempvar __tempgen
        local destvar "`__tempgen'"
        local staged `staged' `destvar'
        if `__do_replace' {
            local __final_name "`srcvar'"
        }
        else {
            local __final_name : word `i' of `generate'
        }

        * Set label name for this variable
        if "`label'" != "" {
            local this_label "`label'"
        }
        else {
            local this_label "`__final_name'"
        }

        * An existing label is extended as encode does: the plugin receives
        * its mappings. Encode itself handles what the plugin hands back
        * (noextend with new strings, codes past the long range) and labels it
        * cannot pass (extended missing values, text over 2,045 bytes).
        local __existing = 0
        local __native = 0
        local __labelmap ""
        capture label list `this_label'
        if !_rc {
            local __existing = 1
            local __ce_ok = 0
            mata: _cencode_pack_label(st_local("this_label"))
            if `__ce_ok' local __labelmap "labelmap=`__ce_len' `noextend'"
            else local __native = 1
        }
        * The plugin reads strL text up to 2,045 bytes and no binary values;
        * encode itself takes the others
        if !`__native' & "`: type `srcvar''" == "strL" {
            mata: st_local("__native", strofreal(_cencode_strl_unsupported(st_local("srcvar"))))
        }
        if `__native' {
            _cencode_native `srcvar' `destvar' `this_label' `"`if' `in'"' `noextend'
            return scalar N_unique = r(N_unique)
            continue
        }

        * Create the destination variable (numeric, long type for sufficient range)
        if "`__fill'" != "" {
            mata: (void) st_addvar("long", st_local("destvar"), 1)
        }
        else {
            quietly generate long `destvar' = .
        }

        * Temp data file for C plugin to write length-prefixed label bytes
        tempfile __labelfile
        * The path can contain spaces, so the plugin reads it from local
        * __ce_labelfile and receives only its length in bytes.
        local __ce_labelfile `"`__labelfile'"'
        local __ce_lflen = strlen(`"`__labelfile'"')
        local labelfile_code "labelfile=`__ce_lflen'"

        * The plugin sees only the source (variable 1), the destination (2)
        * and any sample marker (3).
        if `__do_timing' {
            quietly _ctools_clock __ct_now
            local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
            if "`__ct_elapsed92'" == "" local __ct_elapsed92 = 0
            quietly _ctools_clock __ct_start92
        }
        capture scalar drop _cencode_max_len _cencode_protocol _cencode_fallback _cencode_n_new
        local __marker ""
        if `"`if'"' != "" {
            local __marker "`__iftouse'"
            local __ifflag "ifcond touse"
        }
        plugin call ctools_plugin `srcvar' `destvar' `__marker' `in', "cencode `threads_code' 1 2 `labelfile_code' stream `__fill' `__ifflag' `__labelmap'"
        local __ce_map ""
        local __ce_labelfile ""
        if `__do_timing' {
            quietly _ctools_clock __ct_now
            local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')
            if "`__ct_elapsed93'" == "" local __ct_elapsed93 = 0
            quietly _ctools_clock __ct_start93
        }

        * =====================================================================
        * Install the label vector without interpreting user text as source
        * =====================================================================

        local n_unique = scalar(_cencode_n_unique)

        if `__existing' {
            * A plugin that does not report protocol 2 numbered the strings
            * afresh; encode redoes those, as it does the handed-back cases.
            capture local __protocol = scalar(_cencode_protocol)
            if _rc local __protocol = 0
            capture local __fallback = scalar(_cencode_fallback)
            if _rc local __fallback = 1
            if `__protocol' != 2 | `__fallback' != 0 | `n_unique' == 0 {
                drop `destvar'
                _cencode_native `srcvar' `destvar' `this_label' `"`if' `in'"' `noextend'
                return scalar N_unique = r(N_unique)
                continue
            }
            local n_new = scalar(_cencode_n_new)
            if `n_new' > 0 {
                mata: _cencode_apply_labels(st_local("__labelfile"), st_local("this_label"), strtoreal(st_local("n_new")))
            }
        }
        else if `n_unique' > 0 {
            * Drop any pre-existing label with this name for a clean slate
            capture label drop `this_label'
            mata: _cencode_apply_labels(st_local("__labelfile"), st_local("this_label"), strtoreal(st_local("n_unique")))
        }

        * Apply label to variable if it exists
        capture label list `this_label'
        if _rc == 0 {
            label values `destvar' `this_label'
        }

        * encode sizes the display to the longest encoded value in bytes, with
        * a minimum of eight; the plugin reports that length. With no
        * categories it retains long's default format and empty label.
        if `n_unique' > 0 {
            capture local __max_len = scalar(_cencode_max_len)
            if _rc | "`__max_len'" == "" | "`__max_len'" == "." {
                * A plugin from before that report: measure in Mata.
                mata: st_local("__max_len",strofreal(max(strlen(select(st_sdata(.,st_local("srcvar")),st_data(.,st_local("destvar")):<.)))))
            }
            local display_width = max(8, `__max_len')
            format `destvar' %`display_width'.0g
        }
        else {
            mata: st_vlmodify(st_local("this_label"), J(0,1,.), J(0,1,""))
            label values `destvar' `this_label'
        }
        * Copied through Mata, so text such as $USD is not expanded
        mata: st_varlabel(st_local("destvar"), st_varlabel(st_local("srcvar")))

        * Store rclass results for this variable
        return scalar N_unique = `n_unique'

        if `__do_timing' {
            quietly _ctools_clock __ct_now
            local __ct_elapsed93 = `__ct_elapsed93' + (`__ct_now' - `__ct_start93')
            if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
            quietly _ctools_clock __ct_start91
        }
    }

    _cencode_commit "`varlist'" "`staged'" "`generate'" `__do_replace' "`original_order'"

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')
    }

    * Store rclass results
    return scalar N_vars = `n_vars'

    if `__do_timing' {
        quietly _ctools_clock __ct_now
        local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')

        * Extract timer values
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'

        * Extract sub-timer values
        * Elapsed seconds are held in caller-local timing state.
        local __time_preplugin = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_plugin = `__ct_elapsed92'
        * Elapsed seconds are held in caller-local timing state.
        local __time_postplugin = `__ct_elapsed93'

        * Calculate plugin call overhead
        capture local __plugin_time_total = scalar(_cencode_time_total)
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "cencode timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f scalar(_cencode_time_parse) " sec"
        di as text "    Read and code strings:  " as result %8.4f scalar(_cencode_time_load) " sec"
        di as text "    Merge unique values:    " as result %8.4f scalar(_cencode_time_collect) " sec"
        di as text "    Sort:                   " as result %8.4f scalar(_cencode_time_sort) " sec"
        di as text "    Encode:                 " as result %8.4f scalar(_cencode_time_encode) " sec"
        di as text "    Apply labels:           " as result %8.4f scalar(_cencode_time_labels) " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f scalar(_cencode_time_total) " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Pre-plugin setup:       " as result %8.4f `__time_preplugin' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f `__plugin_call_overhead' " sec"
        di as text "    Post-plugin cleanup:    " as result %8.4f `__time_postplugin' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_preplugin' + `__plugin_call_overhead' + `__time_postplugin'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Variables encoded:        " as result `n_vars'

        * Display thread diagnostics
        capture local __threads_max = scalar(_cencode_threads_max)
        if _rc == 0 {
            capture local __openmp_enabled = scalar(_cencode_openmp_enabled)
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clean up timing scalars
        capture scalar drop _cencode_time_parse _cencode_time_load _cencode_time_collect
        capture scalar drop _cencode_time_sort _cencode_time_encode _cencode_time_labels _cencode_time_total
        capture scalar drop _cencode_threads_max _cencode_openmp_enabled
    }

    * Clean up any remaining scalars
    capture scalar drop _cencode_n_unique _cencode_max_len _cencode_protocol
    capture scalar drop _cencode_fallback _cencode_n_new
end

* Encode one variable with encode itself into destvar; r(N_unique) counts the
* distinct nonmissing codes (in C). No encoding phases run in C for it.
program define _cencode_native, rclass
    args srcvar destvar label qualifiers noextend
    quietly encode `srcvar' `qualifiers', generate(`destvar') label(`label') `noextend'
    _ctools_group `destvar', count
    return scalar N_unique = r(N_groups)
    foreach phase in parse load collect sort encode labels total {
        scalar _cencode_time_`phase' = 0
    }
end

* Give the staged outputs their final names; a Break cannot split this step.
* staged lists the temporary outputs in varlist order. st_varrename() is built
* in; rename is an ado file whose first use in a session costs about 10 ms.
program define _cencode_commit
    args varlist staged generate do_replace original_order
    nobreak {
        local i = 0
        foreach srcvar of local varlist {
            local ++i
            local destvar : word `i' of `staged'
            if `do_replace' {
                * replace keeps the source's notes and characteristics
                char rename `srcvar' `destvar'
                drop `srcvar'
                mata: st_varrename("`destvar'", "`srcvar'")
            }
            else {
                mata: st_varrename("`destvar'", "`: word `i' of `generate''")
            }
        }
        if `do_replace' order `original_order'
    }
end

mata:
/* 1 if a strL variable holds a value the plugin cannot read: longer than
   2,045 bytes, or binary */
real scalar _cencode_strl_unsupported(string scalar var)
{
    string colvector text

    text = st_sdata(., var)
    if (!rows(text)) return(0)
    return(any(strlen(text) :> 2045) | any(strpos(text, char(0))))
}

struct _cencode_vl {
    real scalar exists
    real colvector values
    string colvector text
}

/* Record each named value label exactly, or that it does not exist. */
void _cencode_labels_save(string rowvector names)
{
    external struct _cencode_vl rowvector __cencode_vl_saved
    external string rowvector __cencode_vl_names
    real scalar i
    real colvector values
    string colvector text

    __cencode_vl_names = cols(names) ? uniqrows(names')' : J(1, 0, "")
    __cencode_vl_saved = _cencode_vl(cols(__cencode_vl_names))
    for (i = 1; i <= cols(__cencode_vl_names); i++) {
        __cencode_vl_saved[i].exists = st_vlexists(__cencode_vl_names[i])
        if (__cencode_vl_saved[i].exists) {
            st_vlload(__cencode_vl_names[i], values, text)
            __cencode_vl_saved[i].values = values
            __cencode_vl_saved[i].text = text
        }
    }
}

/* Put the recorded labels back: drop labels that did not exist, restore the rest. */
void _cencode_labels_restore()
{
    external struct _cencode_vl rowvector __cencode_vl_saved
    external string rowvector __cencode_vl_names
    real scalar i

    for (i = 1; i <= cols(__cencode_vl_names); i++) {
        if (st_vlexists(__cencode_vl_names[i])) st_vldrop(__cencode_vl_names[i])
        if (__cencode_vl_saved[i].exists) {
            st_vlmodify(__cencode_vl_names[i], __cencode_vl_saved[i].values,
                __cencode_vl_saved[i].text)
        }
    }
}

/* Hand an existing value label to the plugin as local __ce_map (format: see
   cencode_label_read in src/cencode/cencode_impl.c). __ce_ok is 0 when encode
   itself must handle the label: extended missing or non-integer values, more
   than 65,536 mappings, text over 2,045 bytes, or a mapping too large for a
   macro. */
void _cencode_pack_label(string scalar label)
{
    real colvector values, lengths
    string colvector text, entries
    string scalar packed

    st_local("__ce_ok", "0")
    st_vlload(label, values, text)
    if (rows(values) > 65536) return
    if (rows(values) == 0) {
        packed = "0 |"
    }
    else {
        if (missing(values) | any(values :!= floor(values))) return
        lengths = strlen(text)
        if (max(lengths) > 2045) return
        entries = strtrim(strofreal(values, "%12.0f")) :+ " " :+ strofreal(lengths) :+ " " :+ text
        packed = strofreal(rows(values)) + " " + invtokens(entries', " ") + " |"
        if (strlen(packed) != strlen(strofreal(rows(values))) + sum(strlen(entries)) + rows(values) + 2) return
    }
    if (strlen(packed) >= c("macrolen")) return
    st_local("__ce_map", packed)
    st_local("__ce_len", strofreal(strlen(packed), "%21.0f"))
    st_local("__ce_ok", "1")
}

void _cencode_apply_labels(string scalar path, string scalar name, real scalar n)
{
    real scalar fh, i, valid
    real rowvector header
    real colvector codes
    string colvector texts
    string matrix part

    codes = J(n, 1, .)
    texts = J(n, 1, "")
    fh = fopen(path, "r")
    part = _fget(fh)
    valid = 0
    if (rows(part) == 1 & cols(part) == 1) valid = (part == "CTOOLS_LABELS_1")
    for (i = 1; i <= n & valid; i++) {
        part = _fget(fh)
        if (rows(part) != 1 | cols(part) != 1) {
            valid = 0
            break
        }
        header = strtoreal(tokens(part))
        if (cols(header) != 2) {
            valid = 0
            break
        }
        if (any(missing(header)) | any(header :!= floor(header))) {
            valid = 0
            break
        }
        if (abs(header[1]) > 2147483647 | header[2] < 0 | header[2] > 2045) {
            valid = 0
            break
        }
        part = (header[2] == 0 ? "" : _fread(fh, header[2]))
        if (rows(part) != 1 | cols(part) != 1) {
            valid = 0
            break
        }
        if (strlen(part) != header[2]) {
            valid = 0
            break
        }
        codes[i] = header[1]
        texts[i] = part
    }
    if (valid) {
        part = _fread(fh, 1)
        valid = (rows(part) == 0 & fstatus(fh) == -1)
    }
    fclose(fh)
    if (!valid) _error(610)
    st_vlmodify(name, codes, texts)
}
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* Private monotonic timestamp in seconds; no public Stata timer IDs.
program define _ctools_clock
    version 14.1
    args target
    tempname held
    _return hold `held'
    capture noisily {
        * Once registered, this plugin stays loaded for the Stata session.
        * Resolve its path and check its protocol only on first use (or
        * after discard), rather than searching the ado-path at every lap.
        capture plugin call ctools_clock_plugin, "clock"
        if _rc {
            _ctools_load
            capture program ctools_clock_plugin, plugin using("`__ctools_plugin'")
            if _rc != 0 & _rc != 110 exit _rc
            plugin call ctools_clock_plugin, "clock"
        }
    }
    local rc = _rc
    _return restore `held'
    if `rc' exit `rc'
    c_local `target' `ctools_clock_value'
end

* _ctools_group — group identifiers from the ctools plugin, in place of egen group()
* 
* Usage:  _ctools_group varlist [if] [in], generate(newvar) [label]
* _ctools_group varlist [if] [in], count [strict]
* strict returns plugin errors instead of sorting in Stata.
* 
* Numbers the distinct values of varlist 1, 2, ... in sort order, as
* egen group() does: observations excluded by if/in, or with a missing value
* in any variable of varlist (string "" counts as missing), get missing. The
* new variable is a long labeled "group(varlist)". label also defines a value
* label named after the new variable, with egen group(, label)'s text for
* each group. count only counts the groups. Returns r(N_groups), and with
* label r(reps), one observation number of each group in group order.
* 
* The plugin reads, sorts and numbers the keys. Stata numbers them by sorting
* only when the plugin is unavailable or cannot read a key (strL text over
* 2,045 bytes).

program define _ctools_group, rclass
    version 14.1
    syntax varlist [if] [in], [GENerate(name) LABel COUNT STRICT]
    if ("`generate'" == "") == ("`count'" == "") | ("`label'" != "" & "`generate'" == "") {
        di as error "_ctools_group: specify generate() [label] or count"
        exit 198
    }
    if "`generate'" != "" {
        confirm new variable `generate'
        if "`label'" != "" {
            capture label list `generate'
            if !_rc {
                di as error "default label `generate' already defined"
                exit 110
            }
        }
    }
    local nkeys : word count `varlist'

    local rc = 199
    capture _ctools_load
    if !_rc {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if _rc == 0 | _rc == 110 local rc = 0
    }
    capture confirm number 0
    if !`rc' {
        _ctools_strw `varlist'
        local out ""
        if "`generate'" != "" {
            * The plugin writes every observation of the range, so without in
            * the variable is created without initial values.
            if `"`in'"' == "" mata: (void) st_addvar("long", st_local("generate"), 1)
            else quietly generate long `generate' = .
            local out "out"
        }
        local reps = cond("`label'" != "", "reps", "")
        * Stata answers an if one observation at a time. Without one the
        * plugin skips the check; an if that names a numeric variable (the
        * callers' touse) passes that variable last, as a sample marker the
        * plugin reads in parallel. Any other if is marked first, unless in
        * limits the range (marking covers every observation); that marker is
        * the last variable until it is dropped, which is then free.
        local select "skipif"
        local marker ""
        local marked ""
        local ifcond ""
        if `"`if'"' != "" {
            local select "touse"
            local marker = strtrim(substr(`"`if'"', 3, .))
            capture confirm numeric variable `marker', exact
            if _rc | wordcount(`"`marker'"') != 1 {
                local marker ""
                if `"`in'"' == "" {
                    tempvar marker
                    capture mark `marker' `if'
                    local rc = _rc
                    local marked `marker'
                }
                else {
                    local select ""
                    local ifcond `"`if'"'
                }
            }
        }
        local ctools_group_reps ""
        if !`rc' {
            capture plugin call ctools_plugin `varlist' `generate' `marker' `ifcond' `in', ///
                "group nkeys=`nkeys' `out' `reps' `select'"
            local rc = _rc
        }
        if "`marked'" != "" capture drop `marked'
        * Every failure, Break included, drops the unfilled output variable.
        if `rc' & "`generate'" != "" capture drop `generate'
        if `rc' == 1 exit 1
        if `rc' == 134 {
            di as error "too many values for {bf:label define}"
            exit 134
        }
        if !`rc' {
            local n_groups = scalar(_ctools_group_n)
            local group_reps `"`ctools_group_reps'"'
        }
        capture scalar drop _ctools_group_n
    }
    if `rc' & "`strict'" != "" exit `rc'
    if `rc' {
        _ctools_group_sort `varlist' `if' `in', generate(`generate') `count' `label'
        local n_groups = r(N_groups)
        local group_reps `"`r(reps)'"'
    }

    if "`generate'" != "" {
        if udstrlen("group(`varlist')") > 80 {
            note `generate' : group(`varlist')
            label variable `generate' "see notes"
        }
        else {
            label variable `generate' "group(`varlist')"
        }
        if "`label'" != "" {
            if `n_groups' > 0 {
                mata: _ctools_group_label(st_local("varlist"), st_local("group_reps"), ///
                    st_local("generate"))
            }
            capture label list `generate'
            if !_rc label values `generate' `generate'
        }
    }
    return scalar N_groups = `n_groups'
    * label: the observation number of one observation of each group
    if "`label'" != "" return local reps `"`group_reps'"'
end

* Number groups by sorting in Stata (the data and their sort order are restored).
program define _ctools_group_sort, rclass sortpreserve
    syntax varlist [if] [in], [GENerate(name) COUNT LABel]
    marksample touse, strok
    tempvar order g
    quietly generate long `order' = _n
    sort `touse' `varlist' `order'
    quietly by `touse' `varlist': generate long `g' = (_n == 1) if `touse'
    quietly replace `g' = sum(`g') if `touse'
    quietly summarize `g', meanonly
    return scalar N_groups = cond(r(N), r(max), 0)
    if "`label'" != "" {
        mata: _ctools_group_firsts(st_local("g"), st_local("order"))
        return local reps `"`__reps'"'
    }
    if "`generate'" != "" rename `g' `generate'
end

mata:
/* Original observation number of the first observation of each group. */
void _ctools_group_firsts(string scalar g, string scalar order)
{
    real colvector id, obs, previous, first

    id = st_data(., g)
    obs = st_data(., order)
    st_local("__reps", "")
    if (rows(id) == 0) return
    previous = (rows(id) > 1 ? (. \ id[|1 \ rows(id) - 1|]) : .)
    first = selectindex(id :< . :& (id :!= previous))
    if (rows(first)) st_local("__reps", invtokens(strofreal(obs[first], "%12.0g")'))
}

/*
    Value label of egen group(varlist), label for groups 1..G, from one
    observation of each group (reps, in group order). Each variable
    contributes, separated by spaces: a string's value; a numeric value's
    label text when its value label defines one, otherwise the value
    ("%18.0g" for integers, "%9.0g" or "%9,0g" otherwise). Text is cut to 32
    characters and trimmed, and a string's NUL bytes become \0, which is what
    egen group() does when called under version 16 or earlier, as the ctools
    ado files are.
*/
void _ctools_group_label(string scalar varlist, string scalar reps, string scalar lname)
{
    real rowvector v
    real colvector obs, s
    real scalar i, j
    string scalar fmt, vl
    string colvector x, y

    v = st_varindex(tokens(varlist))
    obs = strtoreal(tokens(reps))'
    if (rows(obs) > 65536) {
        errprintf("too many values for {bf:label define}\n")
        exit(134)
    }
    fmt = (st_global("c(dp)") == "comma" ? "%9,0g" : "%9.0g")
    for (j = 1; j <= cols(v); j++) {
        if (st_isnumvar(v[j])) {
            s = st_data(obs, v[j])
            x = J(rows(s), 1, "")
            vl = st_varvaluelabel(v[j])
            if (vl != "") {
                if (st_vlexists(vl)) x = ustrtrim(usubstr(st_vlmap(vl, s), 1, 32))
            }
            for (i = 1; i <= rows(s); i++) {
                if (x[i] == "") x[i] = strofreal(s[i], (s[i] != trunc(s[i]) ? fmt : "%18.0g"))
            }
        }
        else {
            x = st_sdata(obs, v[j])
            x = usubinstr(x, char(0), "\0", .)
            x = ustrtrim(usubstr(x, 1, 32))
        }
        y = (j == 1 ? x : y :+ " " :+ x)
    }
    st_vlmodify(lname, 1::rows(y), y)
}
end

* Shared plugin loader and identity contract
program define _ctools_load
    version 14.1
    * Not rclass, and findfile's r() is undone, so the caller's r() survives;
    * identity comes back through c_local.
    local __os = c(os)
    local __machine = c(machine_type)
    local __is_mac = 0
    if "`__os'" == "MacOSX" {
        local __is_mac = 1
    }
    else if strpos(lower("`__machine'"), "mac") > 0 {
        local __is_mac = 1
    }

    local __plugin = ""
    if "`__os'" == "Windows" {
        local __plugin "ctools_windows.plugin"
    }
    else if `__is_mac' {
        local __is_arm = 0
        if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
            local __is_arm = 1
        }
        if `__is_arm' == 0 & !strpos(lower("`__machine'"), "intel") & !strpos(lower("`__machine'"), "x86") {
            tempfile __archfile
            quietly shell uname -m > "`__archfile'" 2>&1
            tempname __fh
            file open `__fh' using "`__archfile'", read text
            file read `__fh' __archline
            file close `__fh'
            capture erase "`__archfile'"
            if strpos("`__archline'", "arm64") > 0 {
                local __is_arm = 1
            }
        }
        if `__is_arm' {
            local __plugin "ctools_mac_arm.plugin"
        }
        else {
            local __plugin "ctools_mac_x86.plugin"
        }
    }
    else if "`__os'" == "Unix" {
        local __plugin "ctools_linux.plugin"
    }
    else {
        local __plugin "ctools.plugin"
    }

    * The ado/plugin argument protocol, stamped from config/release.json.
    local __ado_protocol "4"
    local __reinstall `"net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace"'

    tempname __rhold
    _return hold `__rhold'
    capture findfile `__plugin'
    local __found = (_rc == 0)
    local __path `"`r(fn)'"'
    _return restore `__rhold'
    local __wanted "`__plugin'"
    if !`__found' & "`__plugin'" != "ctools.plugin" local __plugin "ctools.plugin"
    capture program ctools_plugin, plugin using("`__plugin'")
    local __load_rc = _rc
    if `__load_rc' != 0 & `__load_rc' != 110 {
        di as error "ctools: could not load the plugin `__plugin'"
        di as error `"  Stata reports c(os) "`__os'" and c(machine_type) "`__machine'"; ctools selected `__wanted'."'
        if `__found' di as error `"  The plugin file is `__path'."'
        else di as error "  `__wanted' was not found on the ado-path; reinstall ctools."
        di as error "  Stata's plugin loader reported:"
        capture noisily program ctools_plugin, plugin using("`__plugin'")
        if "`__os'" == "Unix" & !`__is_mac' {
            di as error "  The distribution Linux plugin needs x86-64 and glibc 2.17 or later; OpenMP is linked statically."
            di as error "  Run ldd -r on the plugin to inspect missing libraries or symbols."
        }
        else if `__is_mac' {
            di as error "  On macOS, check that the plugin is the notarized ctools release and is not quarantined."
        }
        di as error "  See https://github.com/mdroste/stata-ctools/blob/main/docs/INSTALLATION.md#troubleshooting"
        exit `__load_rc'
    }
    local ctools_plugin_hdfe_screen
    local ctools_plugin_hdfe_markout
    local ctools_plugin_hdfe_fvlevels
    capture plugin call ctools_plugin, "version"
    local __id_rc = _rc
    if `__id_rc' | "`ctools_plugin_protocol'" != "`__ado_protocol'" {
        local __has = cond("`ctools_plugin_protocol'" == "", "none", "`ctools_plugin_protocol'")
        di as error "ctools: the plugin `__plugin' does not match the installed ado files."
        local __ident = cond("`ctools_plugin_version'" == "", "", " (plugin version `ctools_plugin_version', revision `ctools_plugin_revision')")
        di as error "  The ado files need plugin protocol `__ado_protocol'; the loaded plugin reports `__has'`__ident'."
        di as error "  Reinstall ctools so its ado files and plugins come from the same release:"
        di as error `"    `__reinstall'"'
        di as error "  Then restart Stata, because a plugin stays loaded until Stata exits."
        exit 9
    }
    c_local __ctools_plugin "`__plugin'"
    c_local __ctools_plugin_version "`ctools_plugin_version'"
    c_local __ctools_plugin_revision "`ctools_plugin_revision'"
    c_local __ctools_plugin_protocol "`ctools_plugin_protocol'"
    c_local __ctools_plugin_hdfe_screen = ("`ctools_plugin_hdfe_screen'" == "1")
    c_local __ctools_plugin_hdfe_markout = ("`ctools_plugin_hdfe_markout'" == "1")
    c_local __ctools_plugin_hdfe_fvlevels = ("`ctools_plugin_hdfe_fvlevels'" == "1")
    capture confirm number 0
end

* ctools output-name validation
program define _ctools_newvars
    version 14.1
    syntax anything(name=outputs)
    local seen ""
    foreach name of local outputs {
        if `: list name in seen' {
            di as error "ctools: duplicate output variable `name'"
            exit 198
        }
        confirm new variable `name'
        local seen `seen' `name'
    }
end

* _ctools_strw — set string width metadata for ctools flat buffer optimization
* 
* Usage:  _ctools_strw varlist
* 
* Sets local __ctools_strw in CALLER's scope containing comma-separated
* string variable widths (0 = numeric, >0 = str width e.g. 17 for str17).
* The ctools C plugin reads this local automatically via SF_macro_use
* to enable flat buffer string I/O (eliminates per-string strlen, memcpy,
* and arena atomic CAS overhead).
* 
* Example:
* _ctools_strw `allvars'
* plugin call ctools_plugin `allvars', "mycommand args"
* * C plugin automatically detects __ctools_strw local

program define _ctools_strw
    syntax varlist

    local __strw ""
    local __sep ""
    foreach v of local varlist {
        local __type : type `v'
        if substr("`__type'", 1, 3) == "str" {
            local __w = real(substr("`__type'", 4, .))
            if missing(`__w') local __w = 0
        }
        else {
            local __w = 0
        }
        local __strw "`__strw'`__sep'`__w'"
        local __sep ","
    }
    c_local __ctools_strw `__strw'
end
* ctools-helpers:end
