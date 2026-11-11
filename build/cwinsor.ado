*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define cwinsor, rclass
    version 14.1
    * A failed or interrupted call leaves the data as they were. New
    * variables (suffix(), prefix()) are only
    * added at the end and dropped again. In place, the plugin keeps the
    * original values of the observations it rewrites until cwinsor finishes
    * and writes them back after a failure; _ctools_undo then restores the
    * variable labels, the sort marker and c(changed). Option snapshot
    * preserves the whole dataset instead (on disk above c(max_preservemem));
    * option fast keeps no copy and undoes nothing. Replacing existing
    * generated targets uses a snapshot to retain their original contents.
    _ctools_mode `0'
    local __command `"`0'"'
    if "`mode'" == "" {
        syntax varlist(numeric) [if] [in] [aw fw iw pw], [REPLACE ///
            P(string) PREfix(string) Suffix(string) GENerate(namelist) *]
        local __target_prefix `"`prefix'"'
        if `"`p'"' != "" {
            capture confirm number `p'
            if _rc local __target_prefix `"`p'"'
        }
        if "`replace'" != "" & "`__target_prefix'`suffix'`generate'" != "" {
            local mode snapshot
            local __command `"`__command' snapshot"'
        }
        local 0 `"`__command'"'
    }
    * Break is recognized only inside _cwinsor_impl.
    nobreak {
        macro drop CTOOLS_cwinsor_plugin CTOOLS_cwinsor_targets CTOOLS_cwinsor_recast
        if "`mode'" == "snapshot" preserve
        else if "`mode'" == "" _ctools_undo record CTOOLS_u_cwinsor_
        capture noisily break _cwinsor_impl `0'
        local rc = _rc
        * gstats winsor's results, kept before anything else runs
        if !`rc' {
            foreach r in J minJ maxJ cutlow cuthigh {
                local r_`r' = r(`r')
            }
        }
        if `rc' {
            if "`mode'" == "" {
                * Values first, then labels, sort marker and c(changed)
                if "$CTOOLS_cwinsor_targets" != "" {
                    capture program ctools_plugin, plugin using("$CTOOLS_cwinsor_plugin")
                    capture noisily plugin call ctools_plugin $CTOOLS_cwinsor_targets, "cwinsor rollback"
                    if _rc di as error "cwinsor: the original values could not be restored"
                }
                * Storage promoted for fractional cutoffs goes back too (recast
                * refuses if a value would change)
                local recast "$CTOOLS_cwinsor_recast"
                while "`recast'" != "" {
                    gettoken target recast : recast
                    gettoken vtype recast : recast
                    capture recast `vtype' `target'
                }
                _ctools_undo rollback CTOOLS_u_cwinsor_
            }
            else if "`mode'" == "snapshot" restore
            if "$CTOOLS_cwinsor_plugin" != "" {
                capture program ctools_plugin, plugin using("$CTOOLS_cwinsor_plugin")
                capture plugin call ctools_plugin, "cwinsor commit"
            }
            macro drop CTOOLS_cwinsor_plugin CTOOLS_cwinsor_targets CTOOLS_cwinsor_recast
            exit `rc'
        }
        if "$CTOOLS_cwinsor_plugin" != "" {
            capture program ctools_plugin, plugin using("$CTOOLS_cwinsor_plugin")
            capture plugin call ctools_plugin, "cwinsor commit"
        }
        macro drop CTOOLS_cwinsor_plugin CTOOLS_cwinsor_targets CTOOLS_cwinsor_recast
        if "`mode'" == "snapshot" restore, not
        else if "`mode'" == "" _ctools_undo clear CTOOLS_u_cwinsor_
        foreach r in J minJ maxJ cutlow cuthigh {
            return scalar `r' = `r_`r''
        }
    }
end

program define _cwinsor_impl, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start total timing FIRST
    local __ct_elapsed90 = 0
    local __ct_elapsed91 = 0
    local __ct_elapsed92 = 0
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * The abbreviations of winsor2 and gstats winsor: s() suffix, p() prefix,
    * t trim, c() cuts, l label, gen() generate. p() holding a number is
    * cwinsor's lower percentile (a prefix cannot be a number).
    syntax varlist(numeric) [if] [in] [aw fw iw pw], [Cuts(numlist min=2 max=2 >=0 <=100) P(string) PREfix(string) Q(real 99) Trim BY(varlist) Suffix(string) GENerate(namelist) REPLACE Label noMISSing Verbose THReads(integer 0) SNAPshot FAST ///
        compress forcestrl BENCHmark BENCHmarklevel(integer 0) HASHmethod(string) oncollision(string)]

    * Freeze incoming-r() qualifiers before parsing helpers replace r().
    if `"`if'"' != "" {
        tempvar __incoming_sample
        mark `__incoming_sample' `if' `in'
        local if "if `__incoming_sample'"
    }

    local __p_arg `"`p'"'
    local p = 1
    if `"`prefix'"' != "" {
        capture confirm number `prefix'
        if !_rc {
            local p = `prefix'
            local prefix
        }
    }
    if `"`__p_arg'"' != "" {
        capture confirm number `__p_arg'
        if !_rc local p = `__p_arg'
        else {
            if `"`prefix'"' != "" error 198
            local prefix `"`__p_arg'"'
        }
    }

    * gstats winsor's gtools options. compress and forcestrl concern reading
    * strL variables, which cwinsor reads without them; hashmethod() and
    * oncollision() choose gtools' hashing of the by() groups, which cwinsor
    * sorts instead. They are checked as gtools checks them and have no
    * effect; benchmark and benchmarklevel() show cwinsor's timing.
    if `"`hashmethod'"' != "" & !inlist(`"`hashmethod'"', "0", "1", "2", "default", "biject", "spooky") {
        di as error `"`hashmethod' not found"'
        exit 111
    }
    if `"`oncollision'"' != "" & !inlist(`"`oncollision'"', "fallback", "error") {
        di as error "option {opt oncollision()} must be 'fallback' or 'error'"
        exit 198
    }
    if "`benchmark'" != "" | `benchmarklevel' > 0 local verbose verbose

    local __do_timing = ("`verbose'" != "")

    * Handle percentile options
    if "`cuts'" != "" {
        local p : word 1 of `cuts'
        local q : word 2 of `cuts'
    }

    * Validate percentiles
    if `p' < 0 | `p' > 100 {
        di as error "cwinsor: lower percentile must be between 0 and 100"
        exit 198
    }
    if `q' < 0 | `q' > 100 {
        di as error "cwinsor: upper percentile must be between 0 and 100"
        exit 198
    }
    * Equal percentiles are allowed, as in winsor2
    if `p' > `q' {
        di as error "cwinsor: lower percentile (`p') must not exceed upper (`q')"
        exit 198
    }

    * Validate suffix/prefix/generate. As in gstats winsor, replace with them
    * overwrites targets that exist.
    * prefix() and suffix() together name prefix+source+suffix, as in gstats
    if "`generate'" != "" & "`suffix'`prefix'" != "" {
        di as error "cwinsor: specify only one of prefix()/suffix() or generate()"
        exit 198
    }
    if "`generate'" != "" {
        if `: word count `generate'' != `: word count `varlist'' {
            di as error "cwinsor: specify the same number of targets as sources with generate()"
            exit 198
        }
    }

    * Check numeric variables
    foreach v of local varlist {
        capture confirm numeric variable `v'
        if _rc != 0 {
            di as error "cwinsor: `v' is not numeric"
            exit 109
        }
    }

    local nvars : word count `varlist'

    * The plugin's stores into a frame alias variable have no effect; refuse
    * them as replace and winsor2 do
    local inplace = "`suffix'`prefix'`generate'" == ""
    if `inplace' & c(stata_version) >= 18 {
        foreach v of local varlist {
            mata: st_local("__cw_alias", strofreal(st_isalias(st_local("v"))))
            if `__cw_alias' {
                di as error "cannot change values of alias variables"
                di as error "    You tried to change values in alias variable `v'."
                exit 109
            }
        }
    }

    * Mark sample. With weights, as gstats winsor marks it: missing and zero
    * weights are excluded, negative ones refused (also iweights, which
    * gstats accepts; the weighted percentile needs positive weights), and
    * fweights must be integers.
    marksample touse, novarlist
    local weighted = ("`weight'" != "")
    if `weighted' {
        tempvar wvar
        quietly generate double `wvar' `exp' if `touse'
        quietly count if `touse' & `wvar' < 0
        if r(N) {
            di as error "negative weights encountered"
            exit 402
        }
        markout `touse' `wvar'
        quietly replace `touse' = 0 if `wvar' == 0
    }
    * nomissing: with by(), groups with a missing by-value are left out
    if "`missing'" == "nomissing" & "`by'" != "" markout `touse' `by', strok

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Count observations
    quietly count if `touse'
    local nobs_sample = r(N)
    if `nobs_sample' == 0 {
        di as error "cwinsor: no observations"
        exit 2000
    }

    * Count by-variables (sorting now done in C for speed)
    if "`__ct_elapsed91'" == "" local __ct_elapsed91 = 0
    quietly _ctools_clock __ct_start91
    local nby = 0
    if "`by'" != "" {
        local nby : word count `by'
    }
    quietly _ctools_clock __ct_now
    local __ct_elapsed91 = `__ct_elapsed91' + (`__ct_now' - `__ct_start91')

    * Floating-point outputs can be written directly from the plugin's double
    * working columns. Native recast otherwise copies and scans each column
    * again, then moves it back through the entire dataset.
    local direct_outputs = !`inplace' & "`trim'" == ""
    local overlap : list varlist & by
    if "`overlap'" != "" local direct_outputs = 0
    if "`replace'" != "" local direct_outputs = 0
    foreach v of local varlist {
        local output_type : type `v'
        if "`by'" != "" local output_type "`c(type)'"
        if !inlist("`output_type'", "float", "double") local direct_outputs = 0
    }

    * Handle generate new variables
    local target_varlist "`varlist'"
    if !`inplace' {
        local target_varlist ""
        local gen_index = 0
        * The names first, so nothing is created before they are checked
        foreach v of local varlist {
            local ++gen_index
            if "`suffix'`prefix'" != "" local newvar "`prefix'`v'`suffix'"
            else local newvar : word `gen_index' of `generate'
            confirm name `newvar'
            capture confirm variable `newvar', exact
            if _rc == 0 {
                if "`replace'" == "" {
                    di as error "cwinsor: variable `newvar' already exists"
                    exit 110
                }
                capture confirm numeric variable `newvar', exact
                if _rc {
                    di as error "cwinsor: variable `newvar' is not numeric"
                    exit 109
                }
                * A source as its own target is winsorized in place, which
                * an if/in or weighted sample would turn into missing values
                * outside the sample (gstats winsor refuses it with if/in)
                local source_pos : list posof "`newvar'" in varlist
                if `source_pos' & `source_pos' != `gen_index' {
                    di as error "cwinsor: target `newvar' is another source"
                    exit 198
                }
                if `source_pos' & (`"`if'`in'"' != "" | `weighted') {
                    di as error "cwinsor: source `newvar' as its own target is not allowed with if/in or weights"
                    exit 198
                }
            }
            local target_varlist "`target_varlist' `newvar'"
        }
        local dups : list dups target_varlist
        if "`dups'" != "" {
            di as error "cwinsor: target `: word 1 of `dups'' named twice"
            exit 198
        }
        local gen_index = 0
        foreach newvar of local target_varlist {
            local ++gen_index
            local v : word `gen_index' of `varlist'
            capture confirm variable `newvar', exact
            if _rc == 0 {
                * replace: the existing target is overwritten, keeping its
                * place in the data
                if "`newvar'" != "`v'" {
                    quietly recast double `newvar'
                    quietly replace `newvar' = cond(`touse' & !missing(`v'), `v', .)
                }
                continue
            }
            if `direct_outputs' {
                local output_type : type `v'
                if "`by'" != "" local output_type "`c(type)'"
                * Unfilled storage only when the plugin writes every
                * observation (no if/in, and no weights or nomissing that
                * leave observations out)
                local nofill = (`nobs_sample' == _N)
                mata: (void) st_addvar(st_local("output_type"), st_local("newvar"), `nofill')
            }
            else quietly gen double `newvar' = `v' if `touse' & !missing(`v')
        }
    }

    * In place, the default keeps an undo copy of the targets in the plugin,
    * and their variable labels in the record
    local undo_opt ""
    if `inplace' & "`snapshot'`fast'" == "" {
        local undo_opt "undo"
        _ctools_undo varlabels CTOOLS_u_cwinsor_ `varlist'
    }
    * Validation hooks: fail as Break would once new variables exist, or have
    * the plugin fail after writing half the observations
    if "$CTOOLS_cwinsor_fault" == "outputs" error 1
    if "$CTOOLS_cwinsor_fault" == "store" local undo_opt "`undo_opt' fault_store"

    * Checked SPI stores use positions in the plugin varlist, just like reads.
    * Only target variables need store indices; by-variables are read-only.
    * The targets are consecutive, so the plugin receives only the position
    * before the first one (a per-variable list overflows the plugin argument
    * on wide data).
    local store_indices = "storebase=" + string(cond(`direct_outputs', `nvars' + `nby' + `weighted', 0))

    * Build options
    local opts "p=`p' q=`q'"
    if `direct_outputs' local opts "`opts' newoutputs"
    * Every observation is in the sample without if/in (and without weights
    * or nomissing leaving some out): the plugin then reads no marker
    if `nobs_sample' == _N local opts "`opts' skipif"
    if "`trim'" != "" local opts "`opts' trim"
    if `weighted' local opts "`opts' weighted"
    * winsor2 expands non-grouped bounds through numeric local macros. This
    * rounding affects both clipping and trimming, including float midpoints.
    * Keep C's loaded columns and quickselect across this scalar-only step.
    if "`by'"=="" {
        tempname cwinsor_bounds
        matrix `cwinsor_bounds' = J(`nvars', 2, .)
        local opts "`opts' deferbounds"
    }
    else {
        * winsor2 stores by() cutoffs with egen in the source's type: the
        * plugin truncates them for byte, int and long (t) and rounds them to
        * float for float (f); double cutoffs are used as computed (d).
        local cwinsor_rounding
        foreach v of local varlist {
            local vtype : type `v'
            * With weights the cutoffs are gstats winsor's, used as computed
            if `weighted' local cwinsor_rounding = "`cwinsor_rounding'd"
            else local cwinsor_rounding = "`cwinsor_rounding'" + cond(inlist("`vtype'", "byte", "int", "long"), "t", cond("`vtype'" == "float", "f", "d"))
        }
        local opts "`opts' rounding"
    }

    if "`verbose'" != "" local opts "`opts' verbose"

    local threads_opt ""
    if `threads' > 0 {
        local threads_opt "threads(`threads')"
    }

    * Record pre-plugin time
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    * Elapsed seconds are held in caller-local timing state.
    local __time_preplugin = `__ct_elapsed90'
    if "`__ct_elapsed90'" == "" local __ct_elapsed90 = 0
    quietly _ctools_clock __ct_start90

    * Load source + by columns; separate final-type outputs follow them when
    * requested. C uses plugin-local indices for both loading and writing.
    local nvars_target : word count `target_varlist'

    local plugin_varlist `target_varlist' `by' `wvar'
    if `direct_outputs' local plugin_varlist `varlist' `by' `wvar' `target_varlist'
    * Pass the sample marker directly, avoiding serial SF_ifobs evaluation.
    * Keep the same varlist for the deferred-bound application below.
    local plugin_varlist `plugin_varlist' `touse'
    local opts `opts' touse
    plugin call ctools_plugin `plugin_varlist', ///
        "cwinsor `threads_opt' `nvars_target' `nby' `store_indices' `opts' `undo_opt'"
    global CTOOLS_cwinsor_plugin "`__ctools_plugin'"
    if strpos("`undo_opt'", "undo") {
        * The plugin holds the original values until cwinsor commits
        global CTOOLS_cwinsor_targets "`target_varlist'"
    }

    * Include metadata, conversion and cleanup in the post-plugin measurement.
    quietly _ctools_clock __ct_start92

    if "`by'"=="" {
        if "$CTOOLS_cwinsor_fault" == "bounds" error 1
        local bound_index=0
        foreach target of local target_varlist {
            local ++bound_index
            local lower = `cwinsor_bounds'[`bound_index',1]
            local upper = `cwinsor_bounds'[`bound_index',2]
            * winsor2's macro expansion; gstats winsor (weights) uses the
            * cutoffs as computed
            if !`weighted' {
                matrix `cwinsor_bounds'[`bound_index',1] = `lower'
                matrix `cwinsor_bounds'[`bound_index',2] = `upper'
            }
            * winsor2's replace promotes an integer variable that receives a
            * fractional cutoff: byte and int to float, long to double.
            local vtype : type `target'
            if `inplace' & "`trim'" == "" & inlist("`vtype'", "byte", "int", "long") & ///
                (`lower' != floor(`lower') | `upper' != floor(`upper')) {
                local promoted = cond("`vtype'" == "long", "double", "float")
                global CTOOLS_cwinsor_recast "$CTOOLS_cwinsor_recast `target' `vtype'"
                recast `promoted' `target'
            }
        }
        capture noisily plugin call ctools_plugin `plugin_varlist', "cwinsor `threads_opt' apply"
        local apply_rc = _rc
        if `apply_rc' {
            * The plugin restores a failed store before releasing its cache.
            macro drop CTOOLS_cwinsor_targets
            exit `apply_rc'
        }
    }

    * winsor2's label text: -Winsor(p#,p#) or -Trim(p#,p#), a 0 before a
    * lower percentile below 1
    if "`label'" != "" {
        local label_low = `p'
        local label_high = `q'
        if `label_low' < 1 local label_low "0`label_low'"
        local label_note = cond("`trim'" != "", "-Trim", "-Winsor") + "(p`label_low',p`label_high')"
    }

    * gstats winsor's results: groups, the smallest and largest, and the cuts
    return scalar J = _cwinsor_ngroups
    return scalar minJ = _cwinsor_minJ
    return scalar maxJ = _cwinsor_maxJ
    return scalar cutlow = `p'
    return scalar cuthigh = `q'

    * Match winsor2's output metadata and final storage. Keep the working
    * copy double until after computing the bounds, so creation does not round
    * the input before percentile calculation.
    local target_index=0
    foreach v of local varlist {
        local ++target_index
        local target : word `target_index' of `target_varlist'
        local source_label : variable label `v'
        if `"`source_label'"'=="" local source_label "`v'"
        * label: winsor2's note of the cuts
        if "`label'" != "" local source_label `"`source_label'`label_note'"'
        label variable `target' `"`source_label'"'
        if "`target'"!="`v'" & !`direct_outputs' {
            local output_type : type `v'
            if "`by'"!="" & "`trim'"=="" local output_type "`c(type)'"
            recast `output_type' `target', force
            * recast retains the double working copy's format. gen/egen use
            * the default format for the final storage type instead.
            tempvar format_probe
            quietly gen `output_type' `format_probe' = .
            local output_format : format `format_probe'
            format `target' `output_format'
            drop `format_probe'
        }
        * prefix(), generate() and weights are gstats winsor's, whose
        * targets keep the source's display format (winsor2's gen gives the
        * default one)
        if "`target'" != "`v'" & ("`prefix'`generate'" != "" | `weighted') {
            format `target' `: format `v''
        }
    }

    * The plugin changed the variables in place: as with replace, drop Stata's
    * sort marker when it involves one of them (trimmed values sort last, and
    * clamping only the if/in or by() subsets can break the order)
    if "`suffix'" == "" & "`prefix'" == "" {
        local marked : sortedby
        local hit : list marked & target_varlist
        if "`hit'" != "" _cwinsor_clear_sortedby
    }
    * Validation hook: fail after the targets were written and relabeled
    if "$CTOOLS_cwinsor_fault" == "post" error 1
    * Plugin stores do not mark the data changed; replace and winsor2 do
    if `inplace' mata: st_updata(1)

    * Finish timing
    quietly _ctools_clock __ct_now
    local __ct_elapsed90 = `__ct_elapsed90' + (`__ct_now' - `__ct_start90')
    quietly _ctools_clock __ct_now
    local __ct_elapsed92 = `__ct_elapsed92' + (`__ct_now' - `__ct_start92')

    if `__do_timing' {
        * Elapsed seconds are held in caller-local timing state.
        local __time_total = `__ct_elapsed90'
        * Elapsed seconds are held in caller-local timing state.
        local __time_sort = `__ct_elapsed91'
        * Elapsed seconds are held in caller-local timing state.
        local __time_postplugin = `__ct_elapsed92'

        * Calculate pre-plugin time (excluding sort which is timed separately)
        local __time_preplugin_other = `__time_preplugin' - `__time_sort'

        di as text ""
        di as text "{hline 60}"
        di as text "cwinsor timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata pre-plugin:"
        di as text "    Syntax/validation:      " as result %8.4f `__time_preplugin_other' " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f `__time_sort' " sec"
        }
        di as text "  {hline 58}"
        di as text "    Pre-plugin total:       " as result %8.4f `__time_preplugin' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f _cwinsor_time_load " sec"
        if `nby' > 0 {
            di as text "    Sort by-variables:      " as result %8.4f _cwinsor_time_sort " sec"
        }
        di as text "    Group detection:        " as result %8.4f _cwinsor_time_groups " sec"
        di as text "    Winsorization:          " as result %8.4f _cwinsor_time_winsor " sec"
        di as text "    Data store:             " as result %8.4f _cwinsor_time_store " sec"
        di as text "  {hline 58}"
        di as text "    C plugin total:         " as result %8.4f _cwinsor_time_total " sec"
        di as text "  {hline 58}"
        di as text "  Stata post-plugin:        " as result %8.4f `__time_postplugin' " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Dataset: " as result `nobs_sample' as text " obs x " as result `nvars_target' as text " vars"
        di as text "  Groups:  " as result %8.0f _cwinsor_ngroups
        di as text "  Mode:    " as result cond("`trim'"!="", "trim", "winsorize") as text " at [" as result `p' as text ", " as result `q' as text "] percentiles"
        di as text "  Threads: " as result %8.0f _cwinsor_threads
        di as text "{hline 60}"
    }

    * Summary
    if "`verbose'" == "" {
        local mode_str = cond("`trim'" != "", "trimmed", "winsorized")
        di as text "(`nvars_target' variables `mode_str' at [`p', `q'] percentiles)"
    }
end

program define _cwinsor_clear_sortedby
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    * Stata drops the sort marker when a sort variable changes: change one
    * value with replace, then put the original back through Mata (no copy
    * of the variable is needed)
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_sb_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_sb_hold)
    }
    else {
        mata: __ctools_sb_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_sb_hold)
    }
    mata: mata drop __ctools_sb_hold
end
