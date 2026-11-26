/*
    cegen_impl.c
    Extensions to generate: the functions of egen in C

    Two plugin calls per command (see cegen_impl.h):

    compute  parse the request, load the sample (sources, by-variables and
             the sample marker in one parallel pass), run the function's
             kernel and keep the result with the sample's observation map;
    store    write the kept result into the new variable, which the ado
             created in between with its final storage type.

    The loaded sources are released as soon as the result exists. Memory is
    returned on the helper thread after the call, so Stata does not wait for
    it.

    License: MIT
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stplugin.h"
#include "ctools_alloc.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_release.h"
#include "ctools_runtime.h"
#include "ctools_types.h"
#include "cegen_impl.h"
#include "cegen_internal.h"

/* ---------------------------------------------------------------------------
   The result kept between compute and store
   --------------------------------------------------------------------------- */

static struct
{
    int active;
    cegen_result result;
    ctools_filtered_data fd; /* the observation map; sources a string result borrows */
    size_t n;                /* rows of the sample */
    ST_int host_nobs;        /* _N when the result was computed */
} pending;

static void result_free(cegen_result *result, int later)
{
    if (later)
    {
        ctools_release_later_aligned(result->num);
        ctools_release_later(result->str);
        ctools_release_later(result->arena);
    }
    else
    {
        ctools_aligned_free(result->num);
        free(result->str);
        free(result->arena);
    }
    memset(result, 0, sizeof(*result));
}

void cegen_cleanup_state(void)
{
    result_free(&pending.result, 0);
    ctools_filtered_data_free(&pending.fd);
    memset(&pending, 0, sizeof(pending));
}

/* The same at the end of a store: the helper thread returns the memory */
static void pending_release(void)
{
    result_free(&pending.result, 1);
    ctools_filtered_data_release(&pending.fd);
    memset(&pending, 0, sizeof(pending));
}

/* ---------------------------------------------------------------------------
   Request parsing
   --------------------------------------------------------------------------- */

static const struct
{
    const char *name;
    cegen_fn fn;
} cegen_functions[] = {
    {"count", CEGEN_FN_COUNT},         {"mean", CEGEN_FN_MEAN},
    {"total", CEGEN_FN_TOTAL},         {"sd", CEGEN_FN_SD},
    {"min", CEGEN_FN_MIN},             {"max", CEGEN_FN_MAX},
    {"mdev", CEGEN_FN_MDEV},           {"skew", CEGEN_FN_SKEW},
    {"kurt", CEGEN_FN_KURT},           {"std", CEGEN_FN_STD},
    {"pc", CEGEN_FN_PC},               {"median", CEGEN_FN_MEDIAN},
    {"pctile", CEGEN_FN_PCTILE},       {"iqr", CEGEN_FN_IQR},
    {"mad", CEGEN_FN_MAD},             {"seq", CEGEN_FN_SEQ},
    {"rank", CEGEN_FN_RANK},           {"mode", CEGEN_FN_MODE},
    {"rowtotal", CEGEN_FN_ROWTOTAL},   {"rowmean", CEGEN_FN_ROWMEAN},
    {"rowsd", CEGEN_FN_ROWSD},         {"rowmin", CEGEN_FN_ROWMIN},
    {"rowmax", CEGEN_FN_ROWMAX},       {"rowfirst", CEGEN_FN_ROWFIRST},
    {"rowlast", CEGEN_FN_ROWLAST},     {"rowmiss", CEGEN_FN_ROWMISS},
    {"rownonmiss", CEGEN_FN_ROWNONMISS}, {"rowmedian", CEGEN_FN_ROWMEDIAN},
    {"rowpctile", CEGEN_FN_ROWPCTILE}, {"anycount", CEGEN_FN_ANYCOUNT},
    {"anymatch", CEGEN_FN_ANYMATCH},   {"anyvalue", CEGEN_FN_ANYVALUE},
    {"diff", CEGEN_FN_DIFF},           {"group", CEGEN_FN_GROUP},
    {"tag", CEGEN_FN_TAG},             {"cut", CEGEN_FN_CUT},
    {"fill", CEGEN_FN_FILL},           {"ma", CEGEN_FN_MA},
    {"concat", CEGEN_FN_CONCAT},       {"ends", CEGEN_FN_ENDS},
};

static void job_free(cegen_job *job)
{
    free(job->list);
    free(job->punct);
    job->list = NULL;
    job->punct = NULL;
}

static void report_free(cegen_report *report)
{
    free(report->needfmt);
    free(report->cuts);
    free(report->reps);
    memset(report, 0, sizeof(*report));
}

/* A real option that may be absent; 0 when it is present but not a number */
static int parse_real(const char *args, const char *key, double *value)
{
    return ctools_parse_double_checked(args, key, value) != CTOOLS_PARSE_INVALID;
}

/* Local macro `name' of the calling ado, of the stated byte length (owned
 * copy, terminated); NULL when it cannot be read */
static char *read_local(const char *name, size_t length)
{
    char spi_name[64];
    char *text = (char *)malloc(length + 2);
    if (!text)
        return NULL;
    snprintf(spi_name, sizeof(spi_name), "_%s", name);
    text[0] = '\0';
    if (length && SF_macro_use(spi_name, text, (ST_int)(length + 1)) != 0)
    {
        free(text);
        return NULL;
    }
    text[length] = '\0';
    return text;
}

/* The numbers of local cegen_list, as Stata wrote them */
static int parse_list(cegen_job *job, size_t length)
{
    char *text = read_local("cegen_list", length);
    if (!text)
        return 0;
    size_t cap = 16, n = 0;
    double *list = (double *)malloc(cap * sizeof(double));
    char *p = text;
    int ok = list != NULL;
    while (ok)
    {
        while (*p == ' ')
            p++;
        if (!*p)
            break;
        char *end;
        const double v = strtod(p, &end);
        if (end == p || (*end && *end != ' ') || !isfinite(v))
        {
            ok = 0;
            break;
        }
        if (n == cap)
        {
            cap *= 2;
            double *grown = (double *)realloc(list, cap * sizeof(double));
            if (!grown)
            {
                ok = 0;
                break;
            }
            list = grown;
        }
        list[n++] = v;
        p = end;
    }
    free(text);
    if (!ok)
    {
        free(list);
        return 0;
    }
    job->list = list;
    job->nlist = n;
    return 1;
}

static ST_retcode job_parse(const char *args, cegen_job *job)
{
    char name[32] = "", type[16] = "";
    int value = 0;
    memset(job, 0, sizeof(*job));

    if (ctools_parse_string_option(args, "fn", name, sizeof(name)) != CTOOLS_PARSE_OK)
        return 198;
    for (size_t i = 0; i < sizeof(cegen_functions) / sizeof(cegen_functions[0]); i++)
        if (!strcmp(name, cegen_functions[i].name))
            job->fn = cegen_functions[i].fn;
    if (job->fn == CEGEN_FN_NONE)
        return 198;

    if (ctools_parse_int_option(args, "nsrc", &value) != CTOOLS_PARSE_OK || value < 0)
        return 198;
    job->nsrc = (size_t)value;
    if (ctools_parse_int_option(args, "nby", &value) != CTOOLS_PARSE_OK || value < 0)
        return 198;
    job->nby = (size_t)value;
    job->type = CEGEN_T_DOUBLE;
    if (ctools_parse_string_option(args, "type", type, sizeof(type)) == CTOOLS_PARSE_OK &&
        !cegen_type_parse(type, &job->type))
        return 198;

    job->xfloat = ctools_parse_bool_option(args, "xfloat");
    job->stats = ctools_parse_bool_option(args, "stats");
    job->missing = ctools_parse_bool_option(args, "missing");
    job->prop = ctools_parse_bool_option(args, "prop");
    job->field = ctools_parse_bool_option(args, "field");
    job->track = ctools_parse_bool_option(args, "track");
    job->unique = ctools_parse_bool_option(args, "unique");
    job->minmode = ctools_parse_bool_option(args, "minmode");
    job->maxmode = ctools_parse_bool_option(args, "maxmode");
    job->strok = ctools_parse_bool_option(args, "strok");
    job->icodes = ctools_parse_bool_option(args, "icodes");
    job->label = ctools_parse_bool_option(args, "label");
    job->nomiss = ctools_parse_bool_option(args, "nomiss");
    job->head = ctools_parse_bool_option(args, "head");
    job->last = ctools_parse_bool_option(args, "last");
    job->tail = ctools_parse_bool_option(args, "tail");
    job->trim = ctools_parse_bool_option(args, "trim");

    job->p = 50.0;
    job->std_mean = 0.0;
    job->std_sd = 1.0;
    job->seq_from = 1.0;
    job->seq_block = 1.0;
    job->nobs_total = (double)SF_nobs();
    job->seq_to = job->nobs_total;
    double real = 0.0;
    if (!parse_real(args, "p", &job->p) || !parse_real(args, "smean", &job->std_mean) ||
        !parse_real(args, "ssd", &job->std_sd) || !parse_real(args, "from", &job->seq_from) ||
        !parse_real(args, "to", &job->seq_to) || !parse_real(args, "block", &job->seq_block))
        return 198;
    if (!parse_real(args, "nummode", &real))
        return 198;
    job->nummode = (long long)real;
    real = 0.0;
    if (!parse_real(args, "nq", &real))
        return 198;
    job->nq = (long long)real;
    real = 0.0;
    if (!parse_real(args, "window", &real))
        return 198;
    job->window = (long long)real;
    if (job->seq_block < 1.0 || job->nummode < 0 || job->nq < 0 || job->window < 0 ||
        job->nq > 2147483647LL)
        return 198;

    if (ctools_parse_int_option(args, "listlen", &value) == CTOOLS_PARSE_OK)
    {
        if (value < 0 || !parse_list(job, (size_t)value))
            return 198;
    }
    value = 0;
    if (ctools_parse_int_option(args, "punctlen", &value) == CTOOLS_PARSE_INVALID || value < 0)
        return 198;
    job->punct = read_local("cegen_punct", (size_t)value);
    if (!job->punct)
        return 920;
    job->punct_len = strlen(job->punct);
    return 0;
}

/* ---------------------------------------------------------------------------
   Reports: locals of the calling ado
   --------------------------------------------------------------------------- */

static ST_retcode save_text(const char *name, const char *text)
{
    char spi_name[64];
    snprintf(spi_name, sizeof(spi_name), "_%s", name);
    return SF_macro_save(spi_name, (char *)text) ? 920 : 0;
}

static ST_retcode save_count(const char *name, size_t value)
{
    char text[32];
    snprintf(text, sizeof(text), "%llu", (unsigned long long)value);
    return save_text(name, text);
}

static ST_retcode save_real(const char *name, double value)
{
    char text[40];
    snprintf(text, sizeof(text), "%.17g", value);
    return save_text(name, text);
}

static ST_retcode report_save(const cegen_job *job, const cegen_report *report,
                              const cegen_result *result, size_t n)
{
    ST_retcode rc = save_count("cegen_n", n);
    rc = rc ? rc : save_count("cegen_ngroups", report->ngroups);
    rc = rc ? rc : save_text("cegen_kind", result->str ? "str" : "num");
    rc = rc ? rc : save_count("cegen_strlen", report->strlen_max);
    rc = rc ? rc : save_count("cegen_allmiss", (size_t)report->all_missing);
    rc = rc ? rc : save_count("cegen_groupmiss", (size_t)report->group_missing);
    rc = rc ? rc : save_count("cegen_multiple", (size_t)report->multiple);
    rc = rc ? rc : save_count("cegen_maxmodes", report->max_modes);
    rc = rc ? rc : save_count("cegen_badgroups", report->bad_groups);
    rc = rc ? rc : save_text("cegen_needfmt", report->needfmt ? report->needfmt : "");
    rc = rc ? rc : save_count("cegen_overflow", (size_t)report->reps_overflow);
    rc = rc ? rc
            : save_text("cegen_ftype",
                        report->has_final_type ? cegen_type_name(report->final_type) : "");
    if (rc)
        return rc;

    /* The range and kind of a numeric result, for promotion and autotype */
    cegen_range range = {0, 1, 1, 0.0, 0.0};
    if (job->stats && result->num)
        cegen_result_range(result->num, n, &range);
    rc = save_count("cegen_nonmissing", range.nonmissing);
    rc = rc ? rc : save_count("cegen_isint", (size_t)range.all_integer);
    rc = rc ? rc : save_count("cegen_isfloat", (size_t)range.all_float);
    rc = rc ? rc : save_real("cegen_min", range.min);
    rc = rc ? rc : save_real("cegen_max", range.max);
    if (rc)
        return rc;

    /* cut(): the cut points; group(, label): one observation per group */
    size_t items = report->ncuts > report->nreps ? report->ncuts : report->nreps;
    char *text = (char *)malloc(items * 26 + 1);
    if (!text)
        return 920;
    size_t used = 0;
    text[0] = '\0';
    for (size_t i = 0; i < report->ncuts; i++)
        used += (size_t)snprintf(text + used, 27, i ? " %.17g" : "%.17g", report->cuts[i]);
    rc = save_text("cegen_cuts", text);
    used = 0;
    text[0] = '\0';
    for (size_t i = 0; i < report->nreps; i++)
        used += (size_t)snprintf(text + used, 13, i ? " %u" : "%u", (unsigned)report->reps[i]);
    rc = rc ? rc : save_text("cegen_reps", text);
    free(text);
    return rc;
}

/* ---------------------------------------------------------------------------
   compute
   --------------------------------------------------------------------------- */

static int fn_sorts(cegen_fn fn)
{
    return fn == CEGEN_FN_RANK || fn == CEGEN_FN_MODE || fn == CEGEN_FN_GROUP ||
           fn == CEGEN_FN_TAG;
}

static stata_retcode run_kernel(const cegen_job *job, ctools_filtered_data *fd,
                                cegen_result *result, cegen_report *report)
{
    switch (job->fn)
    {
    case CEGEN_FN_RANK:
    case CEGEN_FN_MODE:
        return cegen_order_run(job, &fd->data, result, report);
    case CEGEN_FN_GROUP:
    case CEGEN_FN_TAG:
    case CEGEN_FN_CUT:
    case CEGEN_FN_FILL:
    case CEGEN_FN_MA:
        return cegen_keys_run(job, &fd->data, fd->obs_map, result, report);
    case CEGEN_FN_CONCAT:
    case CEGEN_FN_ENDS:
        return cegen_strings_run(job, &fd->data, result, report);
    default:
        if (job->fn >= CEGEN_FN_ROWTOTAL && job->fn <= CEGEN_FN_DIFF)
            return cegen_rows_run(job, &fd->data, result, report);
        return cegen_stats_run(job, &fd->data, result, report);
    }
}

static ST_retcode cegen_compute(const char *args)
{
    const double t_start = ctools_timer_seconds();
    cegen_job job;
    cegen_report report;
    cegen_result result;
    ctools_filtered_data fd;
    memset(&report, 0, sizeof(report));
    memset(&result, 0, sizeof(result));
    ctools_filtered_data_init(&fd);
    cegen_cleanup_state();

    ST_retcode rc = job_parse(args, &job);
    if (rc)
    {
        job_free(&job);
        ctools_error(CEGEN_MODULE, "invalid request");
        return rc;
    }

    /* Sources, by-variables, then the sample marker */
    const int touse = ctools_parse_bool_option(args, "touse");
    const int skipif = ctools_parse_bool_option(args, "skipif");
    const size_t nload = job.nsrc + job.nby;
    if ((size_t)SF_nvars() != nload + (touse ? 1 : 0))
    {
        job_free(&job);
        ctools_error(CEGEN_MODULE, "unexpected number of variables");
        return 198;
    }
    int flags = skipif ? CTOOLS_LOAD_SKIP_IF : touse ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_CHECK_IF;
    if (!job.nby && !fn_sorts(job.fn))
        flags |= CTOOLS_LOAD_NO_SORT_ORDER;

    stata_retcode status;
    if (nload)
    {
        int *indices = (int *)malloc(nload * sizeof(int));
        if (!indices)
        {
            job_free(&job);
            return 920;
        }
        for (size_t j = 0; j < nload; j++)
            indices[j] = (int)j + 1;
        status = ctools_data_load(&fd, indices, nload, 0, 0, flags);
        free(indices);
    }
    else
        status = ctools_data_select(&fd, 0, 0, flags);
    const double t_load = ctools_timer_seconds();

    const size_t n = fd.data.nobs;
    result.fill = SV_missval;
    if (job.fn == CEGEN_FN_ANYCOUNT || job.fn == CEGEN_FN_ANYMATCH || job.fn == CEGEN_FN_TAG)
        result.fill = 0.0;
    if (status == STATA_OK)
        status = run_kernel(&job, &fd, &result, &report);
    const double t_compute = ctools_timer_seconds();

    if (status == STATA_OK)
        rc = report_save(&job, &report, &result, n);
    else
        rc = ctools_stata_rc(status);
    if (!rc)
    {
        char text[32];
        snprintf(text, sizeof(text), "%.6f", t_load - t_start);
        rc = save_text("cegen_t_load", text);
        snprintf(text, sizeof(text), "%.6f", t_compute - t_load);
        rc = rc ? rc : save_text("cegen_t_compute", text);
        rc = rc ? rc : save_count("cegen_threads", (size_t)ctools_get_max_threads());
    }

    if (!rc && (result.num || result.str))
    {
        /* Keep the result and the observation map. A string result that
         * points into the sources (mode() of a string) keeps them too. */
        const int borrows = result.str && !result.arena;
        pending.active = 1;
        pending.result = result;
        pending.n = n;
        pending.host_nobs = SF_nobs();
        if (borrows)
        {
            pending.fd = fd;
            ctools_filtered_data_init(&fd);
        }
        else
        {
            pending.fd.obs_map = fd.obs_map;
            fd.obs_map = NULL;
        }
        memset(&result, 0, sizeof(result));
    }
    result_free(&result, 1);
    ctools_filtered_data_release(&fd);
    report_free(&report);
    job_free(&job);
    return rc;
}

/* ---------------------------------------------------------------------------
   store
   --------------------------------------------------------------------------- */

static ST_retcode cegen_store(const char *args)
{
    const double t_start = ctools_timer_seconds();
    char type_name[16] = "";
    cegen_type type = CEGEN_T_DOUBLE;
    if (!pending.active || pending.host_nobs != SF_nobs() || SF_nvars() != 1)
    {
        cegen_cleanup_state();
        ctools_error(CEGEN_MODULE, "no result to store");
        return 198;
    }
    if (ctools_parse_string_option(args, "type", type_name, sizeof(type_name)) !=
            CTOOLS_PARSE_OK ||
        !cegen_type_parse(type_name, &type))
    {
        cegen_cleanup_state();
        return 198;
    }

    const size_t nobs_total = (size_t)SF_nobs();
    size_t nmissing = 0;
    stata_retcode status;
    if (pending.result.str)
        status = cegen_store_strings(pending.result.str, pending.n, pending.fd.obs_map, 1,
                                     nobs_total, &nmissing);
    else
        status = cegen_store_numeric(pending.result.num, pending.n, pending.fd.obs_map, 1,
                                     nobs_total, type, pending.result.fill, &nmissing);
    pending_release();
    if (status != STATA_OK)
    {
        if (status == STATA_ERR_STATA_WRITE)
            ctools_error(CEGEN_MODULE, "failed to write Stata data");
        return ctools_stata_rc(status);
    }
    char text[32];
    snprintf(text, sizeof(text), "%.6f", ctools_timer_seconds() - t_start);
    ST_retcode rc = save_text("cegen_t_store", text);
    return rc ? rc : save_count("cegen_nmiss", nmissing);
}

/* ---------------------------------------------------------------------------
   Entry point
   --------------------------------------------------------------------------- */

ST_retcode cegen_main(const char *args)
{
    if (!args)
        return 198;
    while (*args == ' ')
        args++;
    if (!strncmp(args, "compute ", 8))
        return cegen_compute(args + 8);
    if (!strncmp(args, "store ", 6))
        return cegen_store(args + 6);
    if (!strcmp(args, "discard"))
    {
        cegen_cleanup_state();
        return 0;
    }
    ctools_error(CEGEN_MODULE, "unknown request");
    return 198;
}
