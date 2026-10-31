*! version 1.2.3 20261030 github.com/mdroste/stata-ctools

program define cencode, rclass
    version 14.1
    * Outputs are built under temporary names and committed only after every
    * variable is encoded, so a failure leaves the variables untouched. Value
    * labels are the only other state cencode changes: the labels it may write
    * are saved exactly (through Mata) and put back if anything fails.
    syntax varlist [if] [in], [Generate(string) replace Label(name) *]
    local __labels "`label'"
    if "`__labels'" == "" {
        local __labels = cond("`replace'" != "", "`varlist'", "`generate'")
    }
    capture mata: _cencode_labels_save(tokens(st_local("__labels")))
    local __saved = (_rc == 0)
    capture noisily _cencode_impl `0'
    local rc = _rc
    if `rc' {
        if `__saved' capture mata: _cencode_labels_restore()
        capture mata: mata drop __cencode_vl_saved __cencode_vl_names
        exit `rc'
    }
    capture mata: mata drop __cencode_vl_saved __cencode_vl_names
    return add
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

    * Check that all source variables are string (match encode rc=107)
    foreach v of local varlist {
        capture confirm string variable `v'
        if _rc != 0 {
            di as error "`v' is not a string variable"
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
    * For each plugin call an if is marked and passed last, as a sample
    * marker the plugin reads in parallel; Stata would answer the if one
    * observation at a time. The marker is the last variable until it is
    * dropped after the call, which is then free. With in, which limits the
    * range, the plugin keeps the if (marking covers every observation).
    local __ifcond = cond(`"`in'"' == "", "", `"`if'"')

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
        if `"`if'"' != "" & `"`in'"' == "" {
            tempvar __marker
            mark `__marker' `if'
            local __ifflag "ifcond touse"
        }
        plugin call ctools_plugin `srcvar' `destvar' `__marker' `__ifcond' `in', "cencode `threads_code' 1 2 `labelfile_code' stream `__fill' `__ifflag' `__labelmap'"
        if "`__marker'" != "" drop `__marker'
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

        local n_unique = _cencode_n_unique

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
            local n_new = _cencode_n_new
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
        local source_label : variable label `srcvar'
        label variable `destvar' `"`source_label'"'

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
        capture local __plugin_time_total = _cencode_time_total
        if _rc != 0 local __plugin_time_total = 0
        local __plugin_call_overhead = `__time_plugin' - `__plugin_time_total'

        di as text ""
        di as text "{hline 55}"
        di as text "cencode timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Argument parsing:       " as result %8.4f _cencode_time_parse " sec"
        di as text "    Read and code strings:  " as result %8.4f _cencode_time_load " sec"
        di as text "    Merge unique values:    " as result %8.4f _cencode_time_collect " sec"
        di as text "    Sort:                   " as result %8.4f _cencode_time_sort " sec"
        di as text "    Encode:                 " as result %8.4f _cencode_time_encode " sec"
        di as text "    Apply labels:           " as result %8.4f _cencode_time_labels " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cencode_time_total " sec"
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
        capture local __threads_max = _cencode_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cencode_openmp_enabled
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
