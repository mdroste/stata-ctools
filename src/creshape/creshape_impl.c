/*
 * creshape_impl.c - plugin calls of creshape.
 *
 * The dataset is rebuilt in four steps, three of them here:
 *
 *   plan      Every variable of every observation is staged in C memory
 *             (ctools_rowstage). The request is checked on the staged copy,
 *             as reshape checks it, before anything in Stata changes.
 *   (ado)     The ado file drops the observations and saves the empty
 *             dataset: a template holding types, formats, labels and
 *             characteristics.
 *   schema    The template is rewritten into the description of the result,
 *             which the ado file loads with -use- and extends with -set obs-.
 *   store     The values are written in parallel.
 *
 * Until the store succeeds the staged copy is the only copy of the user's
 * data. A failure in schema or store is therefore reported through the local
 * __creshape_rc with a zero return code, so that the dispatcher keeps the
 * state (ctools_runtime.c frees it on any other return code) and the ado file
 * can load the template again and call restore.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "creshape_impl.h"
#include "creshape_internal.h"
#include "ctools_config.h"
#include "ctools_file.h"
#include "ctools_parse.h"
#include "ctools_release.h"
#include "ctools_runtime.h"

static creshape_state state;

void creshape_cleanup_state(void)
{
    creshape_state *st = &state;
    /* the large blocks go to the release helper at the end of a call, and
     * are freed at once outside one */
    ctools_rowstage_release(&st->staged);
    ctools_release_later_aligned(st->order);
    ctools_release_later(st->cell);
    free(st->group_start);
    free(st->first);
    free(st->names);
    free(st->names_buf);
    free(st->name_slots);
    free(st->is_str);
    free(st->icols);
    free(st->atwl);
    free(st->xicols);
    for (size_t k = 0; k < st->nstubs; k++)
    {
        creshape_stub *sb = &st->stubs[k];
        for (size_t l = 0; sb->names && l < st->nlevels; l++)
            free(sb->names[l]);
        free(sb->names);
        free(sb->cols);
        free(sb->spec);
        free(sb->prefix);
        free(sb->suffix);
        free(sb->long_name);
    }
    free(st->stubs);
    for (size_t l = 0; l < st->nlevels; l++)
        free(st->levels[l].text);
    free(st->levels);
    free(st->report_levels);
    free(st->out);
    memset(st, 0, sizeof(*st));
}

static void creshape_set_local(const char *name, const char *text)
{
    char spi[96];
    snprintf(spi, sizeof(spi), "_%s", name);
    SF_macro_save(spi, (char *)(uintptr_t)text);
}

static void creshape_set_local_number(const char *name, double value)
{
    char text[40];
    snprintf(text, sizeof(text), "%.0f", value);
    creshape_set_local(name, text);
}

/* The names of the i variables, for messages. */
static void creshape_i_names(const creshape_state *st, char *buf, size_t cap)
{
    size_t at = 0;
    buf[0] = '\0';
    for (size_t k = 0; k < st->ni && at + 1 < cap; k++)
        at += (size_t)snprintf(buf + at, cap - at, "%s%s", k ? " " : "", st->names[st->icols[k]]);
}

static void creshape_error_i_unique(const creshape_state *st)
{
    char i[1024], text[4096];
    creshape_i_names(st, i, sizeof(i));
    creshape_error("variable {bf:%s} does not uniquely identify the observations", i);
    snprintf(text, sizeof(text),
             "{p 4 4 2}Your data are currently wide. You are performing a {bf:creshape long}. "
             "You specified {bf:i(%s)} and {bf:j(%s)}. In the current wide form, variable "
             "{bf:%s} should uniquely identify the observations. Remember this picture:{p_end}\n",
             i, st->jname, i);
    creshape_error_text(text);
    creshape_error_text(creshape_picture);
    creshape_error_text(
        "{p 4 4 2}Type {bf:creshape error} for a list of the problem observations.{p_end}");
}

static void creshape_error_j_unique(const creshape_state *st)
{
    char i[1024], text[4096];
    creshape_i_names(st, i, sizeof(i));
    creshape_error("values of variable {bf:%s} not unique within {bf:%s}", st->jname, i);
    snprintf(text, sizeof(text),
             "{p 4 4 2}Your data are currently long. You are performing a {bf:creshape wide}. "
             "You specified {bf:i(%s)} and {bf:j(%s)}. There are observations within {bf:i(%s)} "
             "with the same value of {bf:j(%s)}. In the long data, variables {bf:i()} and "
             "{bf:j()} together must uniquely identify the observations.{p_end}\n",
             i, st->jname, i, st->jname);
    creshape_error_text(text);
    creshape_error_text(creshape_picture);
    creshape_error_text(
        "{p 4 4 2}Type {bf:creshape error} for a list of the problem variables.{p_end}");
}

static void creshape_error_xi_varies(const creshape_state *st, const unsigned char *varying)
{
    char i[1024], text[4096];
    creshape_i_names(st, i, sizeof(i));
    for (size_t x = 0; x < st->nxi; x++)
        if (varying[x])
            creshape_error("variable {bf:%s} not constant within {bf:%s}", st->names[st->xicols[x]],
                           i);
    snprintf(text, sizeof(text),
             "{p 4 4 2}Your data are currently long. You are performing a {bf:creshape wide}. "
             "You typed something like{p_end}\n\n"
             "{p 8 8 2}{bf:. creshape wide a b, i(%s) j(%s)}{p_end}\n\n"
             "{p 4 4 2}There are variables other than {bf:a}, {bf:b}, {bf:%s}, {bf:%s} in your "
             "data. They must be constant within {bf:%s} because that is the only way they can "
             "fit into wide data without loss of information.{p_end}\n\n",
             i, st->jname, i, st->jname, i);
    creshape_error_text(text);
    snprintf(text, sizeof(text),
             "{p 4 4 2}The variable or variables listed above are not constant within {bf:%s}. "
             "Perhaps the values are in error. Type {bf:creshape error} for a list of the "
             "problem observations.{p_end}\n\n"
             "{p 4 4 2}Either that, or the values vary because they should vary, in which case "
             "you must either add the variables to the list of xij variables to be reshaped, or "
             "{bf:drop} them.{p_end}\n",
             i);
    creshape_error_text(text);
}

/* ------------------------------------------------------------------------
 * plan
 * ------------------------------------------------------------------------ */

/* 1 when the ado file must use reshape itself: the reason is in the local
 * __creshape_fallback. Nothing stays allocated. */
static int creshape_fall_back(const char *reason)
{
    creshape_cleanup_state();
    creshape_set_local("__creshape_fallback", reason);
    return 0;
}

static ST_retcode creshape_plan(const char *args)
{
    creshape_state *st = &state;
    creshape_cleanup_state();
    creshape_set_local("__creshape_fallback", "");
    const double t0 = ctools_timer_seconds();

    if (!strncmp(args, "long", 4))
        st->dir = CRESHAPE_LONG;
    else if (!strncmp(args, "wide", 4))
        st->dir = CRESHAPE_WIDE;
    else
        return 198;

    int rc = creshape_plan_read(st, args);
    if (rc)
        return rc;

    /* The SPI cannot write strL values; reshape itself handles such data. */
    for (size_t j = 0; j < st->nvars; j++)
        if (st->is_str[j] && SF_var_is_strl((ST_int)j + 1))
            return creshape_fall_back("strL");

    /* what the names alone decide */
    if (st->dir == CRESHAPE_LONG)
    {
        if ((rc = creshape_plan_long(st)) || (rc = creshape_plan_outputs(st)))
            return rc;
    }
    else if ((rc = creshape_plan_wide(st)))
        return rc;

    /* stage every variable: the result needs them, and so does a restore */
    if (st->nobs)
    {
        int *vars = (int *)malloc(st->nvars * sizeof(int));
        if (!vars)
        {
            ctools_error_alloc("creshape");
            return 920;
        }
        for (size_t j = 0; j < st->nvars; j++)
            vars[j] = (int)j + 1;
        if (!ctools_rowstage_supported(vars, st->nvars, NULL, st->nobs))
        {
            free(vars);
            return creshape_fall_back("strings too wide to stage");
        }
        stata_retcode loaded =
            ctools_rowstage_load(&st->staged, vars, st->nvars, NULL, 1, st->nobs);
        free(vars);
        if (loaded != STATA_OK)
        {
            if (loaded == STATA_ERR_MEMORY)
                creshape_error("creshape: not enough memory to copy the dataset");
            else
                creshape_error("creshape: could not read the dataset");
            return ctools_stata_rc(loaded);
        }
    }
    st->t_load = ctools_timer_seconds() - t0;
    const double t1 = ctools_timer_seconds();

    if (st->dir == CRESHAPE_LONG)
    {
        int unique = 1;
        if ((rc = creshape_order_rows(st, &unique)))
        {
            ctools_error_alloc("creshape");
            return rc;
        }
        if (!unique)
        {
            creshape_error_i_unique(st);
            return 9;
        }
    }
    else
    {
        if ((rc = creshape_wide_levels(st)))
        {
            if (rc == 920)
                ctools_error_alloc("creshape");
            return rc;
        }
        st->t_levels = ctools_timer_seconds() - t1;
        const double t2 = ctools_timer_seconds();
        int unique = 1;
        if ((rc = creshape_order_rows(st, &unique)))
        {
            ctools_error_alloc("creshape");
            return rc;
        }
        if (!unique)
        {
            creshape_error_j_unique(st);
            return 9;
        }
        st->t_sort = ctools_timer_seconds() - t2;
        const double t3 = ctools_timer_seconds();
        unsigned char *varying = NULL;
        if ((rc = creshape_wide_groups(st, &varying)))
        {
            creshape_error("creshape: not enough memory to lay out the wide data");
            return rc;
        }
        st->t_groups = ctools_timer_seconds() - t3;
        int any = 0;
        for (size_t x = 0; x < st->nxi; x++)
            any |= varying[x];
        if (any)
            creshape_error_xi_varies(st, varying);
        free(varying);
        if (any)
            return 9;
        if ((rc = creshape_plan_outputs(st)))
            return rc;
    }
    st->t_order = ctools_timer_seconds() - t1;

    creshape_set_local_number("__creshape_N", (double)st->out_nobs);
    creshape_set_local_number("__creshape_k", (double)st->out_nvars);
    SF_scal_save("_creshape_time_load", st->t_load);
    SF_scal_save("_creshape_time_plan", st->t_order);
    SF_scal_save("_creshape_time_levels", st->t_levels);
    SF_scal_save("_creshape_time_sort", st->dir == CRESHAPE_LONG ? st->t_order : st->t_sort);
    SF_scal_save("_creshape_time_groups", st->t_groups);
    CTOOLS_SAVE_THREAD_INFO("_creshape");
    st->stage = 1;
    return 0;
}

/* ------------------------------------------------------------------------
 * schema, store, restore
 * ------------------------------------------------------------------------ */

/* Report a step's outcome in the local __creshape_rc and keep the state. */
static ST_retcode creshape_status(int rc)
{
    creshape_set_local_number("__creshape_rc", (double)rc);
    return 0;
}

static ST_retcode creshape_schema(void)
{
    creshape_state *st = &state;
    if (st->stage != 1)
        return creshape_status(198);
    const double t0 = ctools_timer_seconds();
    char *template_path = creshape_local("__creshape_template");
    char *out_path = creshape_local("__creshape_schema");
    char msg[256] = "";
    int rc = 920;
    if (template_path && out_path)
        rc = creshape_schema_write(st, template_path, out_path, msg, sizeof(msg));
    free(template_path);
    free(out_path);
    if (rc)
    {
        creshape_error("creshape: %s", msg[0] ? msg : "out of memory");
        return creshape_status(rc);
    }
    st->t_schema = ctools_timer_seconds() - t0;
    SF_scal_save("_creshape_time_schema", st->t_schema);
    st->stage = 2;
    return creshape_status(0);
}

/* 1,234,567 when it fits in 12 columns, as %12.0gc shows it */
static void creshape_commas(size_t value, char *buf, size_t cap)
{
    char digits[32], out[48];
    int n = snprintf(digits, sizeof(digits), "%zu", value), at = 0;
    for (int i = 0; i < n; i++)
    {
        if (i && (n - i) % 3 == 0)
            out[at++] = ',';
        out[at++] = digits[i];
    }
    out[at] = '\0';
    snprintf(buf, cap, "%s", at <= 12 ? out : digits);
}

/* reshape's summary of what changed */
static void creshape_report(const creshape_state *st)
{
    const int to_long = st->dir == CRESHAPE_LONG;
    char line[512], a[48], b[48];
    SF_display("{txt}\n");
    snprintf(line, sizeof(line), "{txt}Data{col 36}%s{col 43}->{col 48}%s\n",
             to_long ? "Wide" : "Long", to_long ? "Long" : "Wide");
    SF_display(line);
    SF_display("{txt}{hline 77}\n");
    creshape_commas(st->old_nobs, a, sizeof(a));
    creshape_commas(st->out_nobs, b, sizeof(b));
    snprintf(line, sizeof(line),
             "{txt}Number of observations{col 28}{res}%12s{txt}{col 43}->   {res}%-12s\n", a, b);
    SF_display(line);
    creshape_commas(st->old_nvars, a, sizeof(a));
    creshape_commas(st->out_nvars, b, sizeof(b));
    snprintf(line, sizeof(line),
             "{txt}Number of variables{col 28}{res}%12s{txt}{col 43}->   {res}%-12s\n", a, b);
    SF_display(line);
    if (to_long)
        snprintf(line, sizeof(line), "{txt}j variable (%zu values){col 43}->{col 48}{res}%s\n",
                 st->nreport, st->jname);
    else
    {
        int col = 31 + (9 - (int)strlen(st->jname));
        snprintf(line, sizeof(line),
                 "{txt}j variable (%zu values){col %d}{res}%s{txt}{col 43}->{col 48}(dropped)\n",
                 st->nreport, col < 1 ? 1 : col, st->jname);
    }
    SF_display(line);
    SF_display("{txt}xij variables:\n");
    for (size_t k = 0; k < st->nstubs; k++)
    {
        /* the first two wide variables and the last, as reshape lists them */
        const creshape_stub *sb = &st->stubs[k];
        const size_t n = st->nreport;
        char list[4 * CRESHAPE_NAME_BYTES + 16] = "";
        size_t at = 0;
        if (n >= 1)
            at += (size_t)snprintf(list + at, sizeof(list) - at, "%s",
                                   sb->names[st->report_levels[0]]);
        if (n >= 2)
            at += (size_t)snprintf(list + at, sizeof(list) - at, " %s",
                                   sb->names[st->report_levels[1]]);
        if (n >= 3)
            snprintf(list + at, sizeof(list) - at, "%s %s", n > 3 ? " ..." : "",
                     sb->names[st->report_levels[n - 1]]);
        const char *left = to_long ? list : sb->long_name;
        const char *right = to_long ? sb->long_name : list;
        int skip = 39 - (int)strlen(left);
        char *row = (char *)malloc(strlen(left) + strlen(right) + 96);
        if (!row)
            continue;
        sprintf(row, "{res}{space %d}%s{col 43}{txt}->{res}{col 48}%s\n", skip < 0 ? 0 : skip, left,
                right);
        SF_display(row);
        free(row);
    }
    SF_display("{txt}{hline 77}\n");
}

static ST_retcode creshape_store_call(void)
{
    creshape_state *st = &state;
    if (st->stage != 2)
        return creshape_status(198);
    const double t0 = ctools_timer_seconds();
    stata_retcode rc = creshape_store(st);
    if (rc != STATA_OK)
    {
        creshape_error("creshape: could not write the reshaped data");
        return creshape_status(ctools_stata_rc(rc));
    }
    st->t_store = ctools_timer_seconds() - t0;
    SF_scal_save("_creshape_time_store", st->t_store);
    if (st->report)
        creshape_report(st);
    /* the result stands: the staged copy is no longer needed */
    creshape_cleanup_state();
    return creshape_status(0);
}

static ST_retcode creshape_restore(void)
{
    creshape_state *st = &state;
    stata_retcode rc = STATA_OK;
    if (!st->stage)
        return 198;
    if (st->nobs)
    {
        /* the ado file has loaded the template and set the observations */
        if ((size_t)SF_nvars() != st->nvars || (size_t)SF_nobs() != st->nobs)
            rc = STATA_ERR_INVALID_INPUT;
        else
            rc = ctools_rowstage_store(&st->staged, NULL, st->nobs, 1);
    }
    creshape_cleanup_state();
    if (rc != STATA_OK)
        creshape_error("creshape: could not restore the original data");
    return ctools_stata_rc(rc);
}

ST_retcode creshape_main(const char *args)
{
    while (*args == ' ')
        args++;
    if (!strncmp(args, "plan ", 5))
    {
        ST_retcode rc = creshape_plan(args + 5);
        if (rc)
            creshape_cleanup_state();
        return rc;
    }
    if (!strcmp(args, "schema"))
        return creshape_schema();
    if (!strcmp(args, "store"))
        return creshape_store_call();
    if (!strcmp(args, "restore"))
        return creshape_restore();
    if (!strcmp(args, "clear"))
    {
        creshape_cleanup_state();
        return 0;
    }
    creshape_cleanup_state();
    return 198;
}
