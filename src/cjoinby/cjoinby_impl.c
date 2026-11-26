/*
    cjoinby_impl.c
    C-accelerated joinby: plugin calls, argument parsing and state

    Plugin calls, in order:

    using <nkeys> <npay> <master_nobs> [verbose]
        Run in the using data. The varlist is the key variables followed by
        the payload variables: the variables new to the master and, under
        update, the ones it shares. Stages every observation and indexes the
        keys (cjoinby_index.c).

    plan <nkeys> <unmatched> <nsort> [verbose]
        Run in the untouched master. The varlist is the master's key
        variables, in the using call's order, followed by any other variables
        of the master's sort order; local __cj_sortcols holds the nsort
        positions (one-based, in this varlist) of the sort variables.
        unmatched is a mask: 1 keeps unmatched master observations, 2
        unmatched using observations. Matches the keys and lays out the output
        (cjoinby_plan.c); nothing has changed if this call fails.

    write <nmaster> <nnew> [merge] [mobs] [uobs] [update] [replace] [undo]
          [fault_store] [verbose]
        Run after the ado added the new variables and observations. The
        varlist is the master variables the plugin moves (nmaster of them, the
        master's own variables in order), the nnew new variables in payload
        order, then the variables the flags name: merge (_merge), mobs and uobs
        (the master and using observation numbers of each output row, for
        strL values the ado copies afterwards). Local __cj_keypos holds the
        position (one-based, among the master variables) of each key, and
        __cj_paymap one entry per payload variable: the position of the
        master variable it shares, or 0 for a new variable. Stages the master
        observations that move and writes the output (cjoinby_write.c). With
        undo the staged observations are kept until clear or rollback;
        fault_store is a validation hook that fails the write half way.

    clear
        Release everything; the join stands.

    rollback
        The varlist is the write call's master variables. Write the staged
        master observations back, then release everything.

    Every call saves its counts and timings in scalars named _cjoinby_*, which
    the ado reads; verbose also prints one line per call (for calls made by
    hand: the ado prints its own timing table).
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_spi.h"
#include "cjoinby_impl.h"
#include "cjoinby_internal.h"

#define CJ_MODULE "cjoinby"

static struct
{
    int phase; /* 1 using indexed, 2 planned, 3 written */
    int nkeys, npay;
    cj_table utab; /* using observations: keys, then payload */
    cj_index ix;
    cj_plan plan;
    cj_table mtab; /* master observations mfirst .. nm - 1, kept for rollback */
    int mheld;
    size_t mfirst;
} cj;

/* Idempotent. The staged data and the index go to the release helper when a
 * release session is open (the plugin call's final cleanup). */
void cjoinby_cleanup_cache(void)
{
    cj_table_release(&cj.utab);
    cj_table_release(&cj.mtab);
    cj_index_release(&cj.ix);
    cj_plan_free(&cj.plan);
    memset(&cj, 0, sizeof(cj));
    cj_table_init(&cj.utab);
    cj_table_init(&cj.mtab);
    cj_index_init(&cj.ix);
    cj_plan_init(&cj.plan);
}

/* Read n integers from the ado's local `name`. Per-variable lists travel
 * through locals: in the plugin argument they would exceed its length on wide
 * data. Fails on truncation, a missing value or trailing tokens. */
static int cj_macro_ints(const char *name, int *out, size_t n)
{
    char local[48];
    if (snprintf(local, sizeof(local), "_%s", name) >= (int)sizeof(local))
        return 0;
    if (n > (size_t)INT_MAX / 16)
        return 0;
    const size_t size = n * 12 + 16;
    char *buf = (char *)malloc(size);
    if (!buf)
        return 0;
    buf[0] = '\0';
    int ok = SF_macro_use(local, buf, (int)size) == 0 && strlen(buf) < size - 1;
    const char *p = buf;
    ok = ok && ctools_parse_int_array(out, n, &p) == 0;
    while (ok && (*p == ' ' || *p == '\t'))
        p++;
    ok = ok && *p == '\0';
    free(buf);
    return ok;
}

static int *cj_sequence(size_t n)
{
    int *idx = (int *)ctools_safe_malloc2(n ? n : 1, sizeof(int));
    if (idx)
        for (size_t j = 0; j < n; j++)
            idx[j] = (int)j + 1;
    return idx;
}

static ST_retcode cj_fail(ST_retcode rc)
{
    cjoinby_cleanup_cache();
    return rc;
}

/* ----------------------------------------------------------------------------
 * using
 * -------------------------------------------------------------------------- */

static ST_retcode cj_using(const char *args)
{
    int nkeys, npay, master_nobs;
    cjoinby_cleanup_cache();
    if (ctools_parse_next_int(&args, &nkeys) || ctools_parse_next_int(&args, &npay) ||
        ctools_parse_next_int(&args, &master_nobs) || nkeys < 1 || npay < 0 || master_nobs < 0 ||
        SF_nvars() != nkeys + npay)
        return cj_fail(198);
    const int verbose = ctools_parse_bool_option(args, "verbose");
    const size_t nvars = (size_t)(nkeys + npay);
    const size_t nobs = (size_t)SF_nobs();

    int *idx = cj_sequence(nvars);
    if (!idx)
        return cj_fail(920);
    const double t0 = ctools_timer_seconds();
    stata_retcode rc = cj_table_load(&cj.utab, idx, nvars, 1, nobs);
    free(idx);
    if (rc)
    {
        ctools_error(CJ_MODULE, "could not load the using data");
        return cj_fail(ctools_stata_rc(rc));
    }
    const double t1 = ctools_timer_seconds();
    rc = cj_index_build(&cj.ix, &cj.utab, nkeys, (size_t)master_nobs);
    if (rc)
    {
        if (rc == STATA_ERR_MEMORY)
            ctools_error_alloc(CJ_MODULE);
        return cj_fail(ctools_stata_rc(rc));
    }
    const double t2 = ctools_timer_seconds();

    cj.nkeys = nkeys;
    cj.npay = npay;
    cj.phase = 1;
    ST_retcode src = 0;
    src = src ? src : ctools_scal_save("_cjoinby_using_N", (double)nobs);
    src = src ? src : ctools_scal_save("_cjoinby_groups", (double)cj.ix.ngroups);
    src = src ? src : ctools_scal_save("_cjoinby_direct", (double)cj.ix.use_direct);
    src = src ? src : ctools_scal_save("_cjoinby_rowmajor_using", (double)cj.utab.rowmajor);
    src = src ? src : ctools_scal_save("_cjoinby_t_using_load", t1 - t0);
    src = src ? src : ctools_scal_save("_cjoinby_t_using_index", t2 - t1);
    if (src)
        return cj_fail(src);
    ctools_verbose(CJ_MODULE, verbose,
                   "using: %zu obs, %zu key groups, %s lookup; load %.4fs, index %.4fs", nobs,
                   cj.ix.ngroups, cj.ix.use_direct ? "direct" : "hashed", t1 - t0, t2 - t1);
    return 0;
}

/* ----------------------------------------------------------------------------
 * plan
 * -------------------------------------------------------------------------- */

static ST_retcode cj_plan_call(const char *args)
{
    int nkeys, unmatched, nsort;
    if (cj.phase != 1 || ctools_parse_next_int(&args, &nkeys) ||
        ctools_parse_next_int(&args, &unmatched) || ctools_parse_next_int(&args, &nsort) ||
        nkeys != cj.nkeys || unmatched < 0 || unmatched > 3 || nsort < 0 || SF_nvars() < nkeys ||
        nsort > SF_nvars())
        return cj_fail(198);
    const int verbose = ctools_parse_bool_option(args, "verbose");
    const size_t nm = (size_t)SF_nobs();
    const int nvars = SF_nvars();

    int *sortcols = NULL;
    cj_sortspec spec = {0, NULL};
    if (nsort > 0)
    {
        sortcols = (int *)malloc((size_t)nsort * sizeof(int));
        if (!sortcols)
            return cj_fail(920);
        if (!cj_macro_ints("__cj_sortcols", sortcols, (size_t)nsort))
        {
            free(sortcols);
            return cj_fail(198);
        }
        for (int s = 0; s < nsort; s++)
        {
            if (sortcols[s] < 1 || sortcols[s] > nvars)
            {
                free(sortcols);
                return cj_fail(198);
            }
            sortcols[s]--;
        }
        spec.nsort = nsort;
        spec.col = sortcols;
    }

    const double t0 = ctools_timer_seconds();
    ctools_filtered_data mkeys;
    ctools_filtered_data_init(&mkeys);
    if (nm > 0)
    {
        stata_retcode lrc = ctools_data_load_ex(
            &mkeys, NULL, 0, 0, 0, CTOOLS_LOAD_SKIP_IF | CTOOLS_LOAD_NO_SORT_ORDER, NULL);
        if (lrc == STATA_OK && (mkeys.data.nobs != nm || mkeys.data.nvars != (size_t)nvars))
            lrc = STATA_ERR_STATA_READ;
        if (lrc)
        {
            ctools_filtered_data_free(&mkeys);
            free(sortcols);
            ctools_error(CJ_MODULE, "could not load the master's key variables");
            return cj_fail(ctools_stata_rc(lrc));
        }
    }
    const double t1 = ctools_timer_seconds();
    int rc = cj_plan_build(&cj.plan, &cj.ix, nm ? &mkeys.data : NULL, nm, unmatched,
                           nsort > 0 ? &spec : NULL);
    ctools_filtered_data_free(&mkeys);
    free(sortcols);
    if (rc)
    {
        if (rc == 920)
            ctools_error(CJ_MODULE, "the joined data would exceed Stata's limit of "
                                    "2,147,483,647 observations, or memory ran out");
        else if (rc == 106)
            ctools_error(CJ_MODULE, "key variables must be of the same generic type "
                                    "(both numeric or both string) in the master and using data");
        return cj_fail((ST_retcode)rc);
    }
    const double t2 = ctools_timer_seconds();

    cj.phase = 2;
    ST_retcode src = 0;
    src = src ? src : ctools_scal_save("_cjoinby_N", (double)cj.plan.nout);
    src = src ? src : ctools_scal_save("_cjoinby_N_1", (double)cj.plan.n_master);
    src = src ? src : ctools_scal_save("_cjoinby_N_2", (double)cj.plan.n_using);
    src = src ? src : ctools_scal_save("_cjoinby_N_3", (double)cj.plan.n_both);
    src = src ? src : ctools_scal_save("_cjoinby_stage_from", (double)cj.plan.stage_from);
    src = src ? src : ctools_scal_save("_cjoinby_placed", (double)cj.plan.placed);
    src = src ? src : ctools_scal_save("_cjoinby_t_master_load", t1 - t0);
    src = src ? src : ctools_scal_save("_cjoinby_t_match", t2 - t1);
    if (src)
        return cj_fail(src);
    ctools_verbose(CJ_MODULE, verbose,
                   "plan: %zu master obs, %zu output obs (%zu matched, %zu master only, "
                   "%zu using only); load %.4fs, match %.4fs",
                   nm, cj.plan.nout, cj.plan.n_both, cj.plan.n_master, cj.plan.n_using, t1 - t0,
                   t2 - t1);
    return 0;
}

/* ----------------------------------------------------------------------------
 * write
 * -------------------------------------------------------------------------- */

static ST_retcode cj_write_call(const char *args)
{
    int nmaster, nnew;
    if (cj.phase != 2 || ctools_parse_next_int(&args, &nmaster) ||
        ctools_parse_next_int(&args, &nnew) || nmaster < cj.nkeys || nnew < 0 || nnew > cj.npay)
        return cj_fail(198);
    const int has_merge = ctools_parse_bool_option(args, "merge");
    const int has_mobs = ctools_parse_bool_option(args, "mobs");
    const int has_uobs = ctools_parse_bool_option(args, "uobs");
    const int update = ctools_parse_bool_option(args, "update");
    const int replace = ctools_parse_bool_option(args, "replace");
    const int undo = ctools_parse_bool_option(args, "undo");
    const int fault_store = ctools_parse_bool_option(args, "fault_store");
    const int verbose = ctools_parse_bool_option(args, "verbose");
    const int nkeys = cj.nkeys, npay = cj.npay;
    const size_t nm = cj.plan.nm, nout = cj.plan.nout;

    if (SF_nvars() != nmaster + nnew + has_merge + has_mobs + has_uobs ||
        (size_t)SF_nobs() < nout || (size_t)SF_nobs() < nm)
        return cj_fail(198);

    int *key_mcol = (int *)malloc((size_t)nkeys * sizeof(int));
    int *pay_mcol = (int *)malloc((size_t)(npay ? npay : 1) * sizeof(int));
    int *pay_var = (int *)malloc((size_t)(npay ? npay : 1) * sizeof(int));
    int *idx = cj_sequence((size_t)nmaster);
    ST_retcode rc = 0;
    if (!key_mcol || !pay_mcol || !pay_var || !idx)
        rc = 920;
    if (!rc && !cj_macro_ints("__cj_keypos", key_mcol, (size_t)nkeys))
        rc = 198;
    if (!rc && npay > 0 && !cj_macro_ints("__cj_paymap", pay_mcol, (size_t)npay))
        rc = 198;
    int n_shared = 0;
    for (int k = 0; !rc && k < nkeys; k++)
    {
        if (key_mcol[k] < 1 || key_mcol[k] > nmaster)
            rc = 198;
        key_mcol[k]--;
    }
    for (int c = 0, next = nmaster; !rc && c < npay; c++)
    {
        if (pay_mcol[c] < 0 || pay_mcol[c] > nmaster)
            rc = 198;
        pay_var[c] = pay_mcol[c] == 0 ? ++next : 0;
        n_shared += pay_mcol[c] > 0;
        pay_mcol[c]--;
        if (c == npay - 1 && next != nmaster + nnew)
            rc = 198;
    }
    if (!rc && npay == 0 && nnew != 0)
        rc = 198;

    /* Master observations before stage_from stay where they are, unless
     * update may change their shared variables. */
    const size_t mfirst = (update && n_shared > 0) ? 0 : cj.plan.stage_from;
    const double t0 = ctools_timer_seconds();
    if (!rc && nm > mfirst)
    {
        stata_retcode lrc = cj_table_load(&cj.mtab, idx, (size_t)nmaster, mfirst + 1, nm - mfirst);
        if (lrc)
        {
            ctools_error(CJ_MODULE, "could not load the master data");
            rc = ctools_stata_rc(lrc);
        }
        else
        {
            cj.mheld = 1;
            cj.mfirst = mfirst;
        }
    }
    const double t1 = ctools_timer_seconds();

    if (!rc)
    {
        cj_writer w;
        memset(&w, 0, sizeof(w));
        w.ix = &cj.ix;
        w.plan = &cj.plan;
        w.utab = &cj.utab;
        w.mtab = cj.mheld ? &cj.mtab : NULL;
        w.stage_from = mfirst;
        w.nkeys = nkeys;
        w.key_mcol = key_mcol;
        w.npay = npay;
        w.pay_mcol = pay_mcol;
        w.pay_var = pay_var;
        w.update = update;
        w.replace = replace;
        w.merge_var = has_merge ? nmaster + nnew + 1 : 0;
        w.mobs_var = has_mobs ? nmaster + nnew + has_merge + 1 : 0;
        w.uobs_var = has_uobs ? nmaster + nnew + has_merge + has_mobs + 1 : 0;
        w.limit = fault_store ? nout / 2 : nout;
        stata_retcode wrc = cj_write_output(&w);
        if (!wrc && fault_store)
            wrc = STATA_ERR_MEMORY;
        if (wrc)
        {
            /* The runtime releases this command's state after a failed call,
             * so the master observations go back now. */
            if (cj.mheld && undo && cj_table_restore(&cj.mtab, cj.mfirst + 1) != STATA_OK)
                ctools_error(CJ_MODULE, "the master observations could not be restored");
            rc = ctools_stata_rc(wrc);
        }
    }
    const double t2 = ctools_timer_seconds();
    free(key_mcol);
    free(pay_mcol);
    free(pay_var);
    free(idx);
    if (rc)
        return cj_fail(rc);

    /* The using data and the plan are done with; the master copy stays for a
     * rollback only. */
    cj_table_release(&cj.utab);
    cj_index_release(&cj.ix);
    cj_plan_free(&cj.plan);
    if (!undo)
    {
        cj_table_release(&cj.mtab);
        cj.mheld = 0;
    }
    cj.phase = 3;
    ST_retcode src = 0;
    src = src ? src : ctools_scal_save("_cjoinby_t_stage", t1 - t0);
    src = src ? src : ctools_scal_save("_cjoinby_t_write", t2 - t1);
    src = src ? src : ctools_scal_save("_cjoinby_staged", (double)(nm > mfirst ? nm - mfirst : 0));
    src = src ? src : ctools_scal_save("_cjoinby_threads_max", (double)ctools_get_max_threads());
    src = src ? src : ctools_scal_save("_cjoinby_openmp_enabled", (double)CTOOLS_OPENMP_ENABLED);
    if (src)
        return cj_fail(src);
    ctools_verbose(CJ_MODULE, verbose,
                   "write: staged %zu master obs in %.4fs, wrote %zu obs in %.4fs",
                   nm > mfirst ? nm - mfirst : 0, t1 - t0, nout, t2 - t1);
    return 0;
}

/* ----------------------------------------------------------------------------
 * rollback
 * -------------------------------------------------------------------------- */

static ST_retcode cj_rollback(void)
{
    ST_retcode rc = 0;
    /* Only a completed write left changed master observations; a failed one
     * restored them itself. */
    if (cj.phase == 3 && cj.mheld)
    {
        if ((size_t)SF_nobs() < cj.mfirst + cj.mtab.nobs || (size_t)SF_nvars() < cj.mtab.nvars)
            rc = 198;
        else
            rc = ctools_stata_rc(cj_table_restore(&cj.mtab, cj.mfirst + 1));
        if (rc)
            ctools_error(CJ_MODULE, "the master observations could not be restored");
    }
    cjoinby_cleanup_cache();
    return rc;
}

ST_retcode cjoinby_main(const char *args)
{
    if (!args)
        return cj_fail(198);
    while (*args == ' ' || *args == '\t')
        args++;
    if (!strcmp(args, "clear"))
    {
        cjoinby_cleanup_cache();
        return 0;
    }
    if (!strcmp(args, "rollback"))
        return cj_rollback();
    if (!strncmp(args, "using ", 6))
        return cj_using(args + 6);
    if (!strncmp(args, "plan ", 5))
        return cj_plan_call(args + 5);
    if (!strncmp(args, "write ", 6))
        return cj_write_call(args + 6);
    return cj_fail(198);
}
