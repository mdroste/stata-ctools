/*
 * ccollapse_impl.c
 *
 * The plugin side of ccollapse: the compute and store calls.
 *
 * compute reads the selected observations of the by(), source and weight
 * variables in one tiled pass, applies collapse's sample rules, finds the
 * groups (ccollapse_group.c) and fills one column per target
 * (ccollapse_stats.c, ccollapse_quantile.c). It writes nothing to the
 * dataset: the group keys and the columns stay here, as the collapsed
 * dataset in memory, until the ado has cut the dataset down to one
 * observation per group and store writes them.
 *
 * The sample, as collapse marks it: observations with a zero or missing
 * weight are left out; a negative weight is an error unless the weights are
 * iweights, and so is a noninteger fweight; cw leaves out observations with
 * a missing value in any source variable.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ccollapse_impl.h"
#include "ccollapse_internal.h"
#include "ctools_alloc.h"
#include "ctools_arena.h"
#include "ctools_cancel.h"
#include "ctools_config.h"
#include "ctools_data_io.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"
#include "ctools_spi.h"

/* ===========================================================================
   What compute keeps for store
   =========================================================================== */

static struct
{
    int ready;
    stata_data out; /* one row per group: the by() columns, then the targets */
} cc_state;

void ccollapse_cleanup_state(void)
{
    stata_data_free(&cc_state.out);
    cc_state.ready = 0;
}

/* ===========================================================================
   Arguments
   =========================================================================== */

typedef struct
{
    size_t nby, nsources, ntargets;
    cc_wtype wtype;
    cc_target *targets;
    int marker, cw, verbose;
} cc_config;

static int cc_parse_stat(const char *name, cc_target *target)
{
    static const struct
    {
        const char *name;
        cc_stat stat;
    } names[] = {
        {"mean", CC_STAT_MEAN},             {"sum", CC_STAT_SUM},
        {"rawsum", CC_STAT_RAWSUM},         {"count", CC_STAT_COUNT},
        {"percent", CC_STAT_PERCENT},       {"sd", CC_STAT_SD},
        {"semean", CC_STAT_SEMEAN},         {"sebinomial", CC_STAT_SEBINOMIAL},
        {"sepoisson", CC_STAT_SEPOISSON},   {"min", CC_STAT_MIN},
        {"max", CC_STAT_MAX},               {"first", CC_STAT_FIRST},
        {"last", CC_STAT_LAST},             {"firstnm", CC_STAT_FIRSTNM},
        {"lastnm", CC_STAT_LASTNM},         {"iqr", CC_STAT_IQR},
    };
    target->pct = 0;
    for (size_t k = 0; k < sizeof(names) / sizeof(names[0]); k++)
    {
        if (!strcmp(name, names[k].name))
        {
            target->stat = names[k].stat;
            return 0;
        }
    }
    if (!strcmp(name, "median"))
    {
        target->stat = CC_STAT_PCTILE;
        target->pct = 50;
        return 0;
    }
    if (name[0] == 'p' && name[1] >= '1' && name[1] <= '9' &&
        (name[2] == '\0' || (name[2] >= '0' && name[2] <= '9' && name[3] == '\0')))
    {
        target->stat = CC_STAT_PCTILE;
        target->pct = atoi(name + 1);
        return 0;
    }
    return -1;
}

/* targets=<source>:<stat>[,<source>:<stat>...] */
static int cc_parse_targets(char *list, cc_config *config)
{
    size_t count = 1;
    for (const char *p = list; *p; p++)
        count += *p == ',';
    config->targets = malloc(count * sizeof(*config->targets));
    if (!config->targets)
        return 920;
    size_t n = 0;
    for (char *item = list; item && *item;)
    {
        char *next = strchr(item, ',');
        if (next)
            *next++ = '\0';
        char *colon = strchr(item, ':');
        char *end = NULL;
        if (!colon)
            return 198;
        *colon = '\0';
        long source = strtol(item, &end, 10);
        if (end != colon || source < 1 || (size_t)source > config->nsources ||
            cc_parse_stat(colon + 1, &config->targets[n]))
            return 198;
        config->targets[n++].source = (int)(source - 1);
        item = next;
    }
    if (n == 0)
        return 198;
    config->ntargets = n;
    return 0;
}

static int cc_parse_config(const char *args, cc_config *config)
{
    int nby = -1, nsources = -1;
    char wtype[16] = "";
    memset(config, 0, sizeof(*config));
    if (ctools_parse_int_option(args, "nby", &nby) != CTOOLS_PARSE_OK || nby < 0 ||
        ctools_parse_int_option(args, "nsrc", &nsources) != CTOOLS_PARSE_OK || nsources < 1 ||
        ctools_parse_string_option(args, "wtype", wtype, sizeof(wtype)) != CTOOLS_PARSE_OK)
        return 198;
    config->nby = (size_t)nby;
    config->nsources = (size_t)nsources;
    if (!strcmp(wtype, "none"))
        config->wtype = CC_W_NONE;
    else if (!strcmp(wtype, "fweight"))
        config->wtype = CC_W_FWEIGHT;
    else if (!strcmp(wtype, "aweight"))
        config->wtype = CC_W_AWEIGHT;
    else if (!strcmp(wtype, "iweight"))
        config->wtype = CC_W_IWEIGHT;
    else if (!strcmp(wtype, "pweight"))
        config->wtype = CC_W_PWEIGHT;
    else
        return 198;
    config->marker = ctools_parse_bool_option(args, "marker");
    config->cw = ctools_parse_bool_option(args, "cw");
    config->verbose = ctools_parse_bool_option(args, "verbose");

    size_t cap = strlen(args) + 1;
    char *list = malloc(cap);
    if (!list)
        return 920;
    int rc = ctools_parse_string_option(args, "targets", list, cap) == CTOOLS_PARSE_OK
                 ? cc_parse_targets(list, config)
                 : 198;
    free(list);
    return rc;
}

/* ===========================================================================
   The sample
   =========================================================================== */

/* Keep the rows with keep[i] != 0, in order. Returns the number kept. */
static size_t cc_compact_rows(stata_data *data, const unsigned char *keep, int threads)
{
    size_t n = data->nobs, kept = 0;
    for (size_t i = 0; i < n; i++)
        kept += keep[i] != 0;
    if (kept == n)
        return n;
    int team = (size_t)threads < data->nvars ? threads : (int)data->nvars;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(team) if (team > 1 && n >= 100000)
    for (size_t v = 0; v < data->nvars; v++)
    {
        stata_variable *var = &data->vars[v];
        size_t at = 0;
        if (var->type == STATA_TYPE_DOUBLE)
        {
            double *x = var->data.dbl;
            for (size_t i = 0; i < n; i++)
                if (keep[i])
                    x[at++] = x[i];
        }
        else
        {
            /* strings outside the arena are owned one by one */
            ctools_string_arena *arena = var->_arena;
            int owned = var->string_storage == CTOOLS_STRINGS_OWNED &&
                        (!arena || arena->has_fallback);
            char **s = var->data.str;
            for (size_t i = 0; i < n; i++)
            {
                if (keep[i])
                    s[at++] = s[i];
                else if (owned && s[i] && (!arena || !ctools_string_arena_owns(arena, s[i])))
                    free(s[i]);
            }
        }
        var->nobs = kept;
    }
    data->nobs = kept;
    return kept;
}

/*
    Apply the weight rules and cw to the loaded rows. Returns 0, or 402
    (negative weight), 401 (noninteger frequency weight), 2000 (nothing
    left), 920.
*/
static int cc_apply_sample(stata_data *data, const cc_config *config, int threads)
{
    size_t n = data->nobs;
    const double missval = SV_missval;
    const int weighted = config->wtype != CC_W_NONE;
    int scan_cw = config->cw, string_source = 0;
    for (size_t s = 0; s < config->nsources; s++)
        string_source |= data->vars[config->nby + s].type != STATA_TYPE_DOUBLE;
    if (!weighted && !scan_cw)
        return 0;
    if (!weighted && string_source)
        return 2000; /* markout without strok marks out every observation */

    unsigned char *keep = malloc(n ? n : 1);
    if (!keep)
        return 920;
    memset(keep, 1, n);
    int team = n >= 200000 ? threads : 1;
    int dropped = 0;

    if (weighted)
    {
        const double *w = data->vars[config->nby + config->nsources].data.dbl;
        const int allow_negative = config->wtype == CC_W_IWEIGHT;
        const int whole = config->wtype == CC_W_FWEIGHT;
        int negative = 0, fractional = 0;
        #pragma omp parallel for schedule(static) num_threads(team) if (team > 1) \
            reduction(|: negative, fractional, dropped)
        for (size_t i = 0; i < n; i++)
        {
            double v = w[i];
            if (!(v < missval) || v == 0.0)
            {
                keep[i] = 0;
                dropped |= 1;
            }
            else
            {
                negative |= v < 0.0;
                fractional |= whole && v != floor(v);
            }
        }
        if ((negative && !allow_negative) || fractional)
        {
            /* a noninteger frequency weight is reported first */
            free(keep);
            return fractional ? 401 : 402;
        }
    }
    if (scan_cw && string_source)
    {
        /* after the weights, whose errors come first: markout without strok
         * marks out every observation */
        free(keep);
        return 2000;
    }
    if (scan_cw)
    {
        for (size_t s = 0; s < config->nsources; s++)
        {
            const double *x = data->vars[config->nby + s].data.dbl;
            #pragma omp parallel for schedule(static) num_threads(team) if (team > 1) \
                reduction(|: dropped)
            for (size_t i = 0; i < n; i++)
            {
                if (!(x[i] < missval))
                {
                    keep[i] = 0;
                    dropped |= 1;
                }
            }
        }
    }
    size_t kept = dropped ? cc_compact_rows(data, keep, threads) : n;
    free(keep);
    return kept ? 0 : 2000;
}

/* ===========================================================================
   The collapsed dataset
   =========================================================================== */

/* dst: the strings of src at rows[0 .. count), in an arena of their own */
static stata_retcode cc_gather_strings(const stata_variable *src, const perm_idx_t *rows,
                                       size_t count, stata_variable *dst)
{
    size_t bytes = 0, longest = 0;
    for (size_t g = 0; g < count; g++)
    {
        const char *s = src->data.str[rows[g]];
        size_t length = s ? strlen(s) : 0;
        if (length > longest)
            longest = length;
        bytes += length + 1;
    }
    dst->type = STATA_TYPE_STRING;
    dst->nobs = count;
    dst->str_maxlen = longest;
    dst->data.str = ctools_safe_cacheline_alloc2(count ? count : 1, sizeof(char *));
    ctools_string_arena *arena =
        ctools_string_arena_create(bytes ? bytes : 1, CTOOLS_STRING_ARENA_NO_FALLBACK);
    dst->_arena = arena;
    /* a column freed before it is filled must hold no stray pointers */
    if (dst->data.str)
        memset(dst->data.str, 0, (count ? count : 1) * sizeof(char *));
    if (!dst->data.str || !arena)
        return STATA_ERR_MEMORY;
    for (size_t g = 0; g < count; g++)
    {
        const char *s = src->data.str[rows[g]];
        dst->data.str[g] = ctools_string_arena_strdup(arena, s ? s : "");
        if (!dst->data.str[g])
            return STATA_ERR_MEMORY;
    }
    return STATA_OK;
}

static stata_retcode cc_gather_numbers(const stata_variable *src, const perm_idx_t *rows,
                                       size_t count, stata_variable *dst)
{
    dst->type = STATA_TYPE_DOUBLE;
    dst->nobs = count;
    dst->data.dbl = ctools_safe_cacheline_alloc2(count ? count : 1, sizeof(double));
    if (!dst->data.dbl)
        return STATA_ERR_MEMORY;
    for (size_t g = 0; g < count; g++)
        dst->data.dbl[g] = src->data.dbl[rows[g]];
    return STATA_OK;
}

/* ===========================================================================
   compute
   =========================================================================== */

static ST_retcode cc_save_results(const cc_groups *groups, size_t nobs, int threads,
                                  const double *seconds)
{
    static const char *const names[] = {"_ccollapse_time_load", "_ccollapse_time_sample",
                                        "_ccollapse_time_group", "_ccollapse_time_stats",
                                        "_ccollapse_time_total"};
    ST_retcode rc = ctools_scal_save("_ccollapse_ngroups", (double)groups->ngroups);
    rc = rc ? rc : ctools_scal_save("_ccollapse_nobs", (double)nobs);
    rc = rc ? rc : ctools_scal_save("_ccollapse_threads", (double)threads);
    rc = rc ? rc : ctools_scal_save("_ccollapse_method", (double)groups->method_code);
    for (size_t k = 0; k < sizeof(names) / sizeof(names[0]) && !rc; k++)
        rc = ctools_scal_save(names[k], seconds[k]);
    return rc ? rc : SF_macro_save("_ccollapse_method", (char *)groups->method);
}

static ST_retcode cc_compute(const char *args)
{
    cc_config config;
    ctools_filtered_data loaded;
    cc_groups groups;
    stata_data out;
    perm_idx_t **string_rows = NULL;
    double **values = NULL;
    int *columns = NULL;
    ST_retcode rc = 0;
    stata_retcode status = STATA_OK;
    double seconds[5] = {0};

    ccollapse_cleanup_state();
    ctools_filtered_data_init(&loaded);
    stata_data_init(&out);
    memset(&groups, 0, sizeof(groups));

    double started = ctools_timer_seconds();
    rc = cc_parse_config(args, &config);
    if (rc)
    {
        if (rc == 198)
            ctools_error(CCOLLAPSE_MODULE, "invalid arguments");
        else
            ctools_error_alloc(CCOLLAPSE_MODULE);
        goto done;
    }
    const int weighted = config.wtype != CC_W_NONE;
    const size_t nload = config.nby + config.nsources + (size_t)weighted;
    if ((size_t)SF_nvars() != nload + (size_t)config.marker)
    {
        ctools_error(CCOLLAPSE_MODULE, "the variables do not match the arguments");
        rc = 198;
        goto done;
    }
    for (size_t t = 0; t < config.ntargets; t++)
    {
        cc_stat stat = config.targets[t].stat;
        int column = (int)(config.nby + (size_t)config.targets[t].source) + 1;
        if (SF_var_is_string(column) && stat != CC_STAT_FIRST && stat != CC_STAT_LAST &&
            stat != CC_STAT_FIRSTNM && stat != CC_STAT_LASTNM)
        {
            rc = 109;
            goto done;
        }
    }
    if (weighted && SF_var_is_string((int)nload))
    {
        rc = 109;
        goto done;
    }
    /* A strL value can be read but not written back: refuse before any work,
     * so that the ado never cuts the dataset down for a store that fails */
    for (size_t v = 1; v <= nload; v++)
    {
        if (SF_var_is_strl((int)v))
        {
            ctools_error(CCOLLAPSE_MODULE, "strL variables are not supported by the plugin");
            rc = 109;
            goto done;
        }
    }

    int threads = ctools_get_max_threads();
    if (threads < 1)
        threads = 1;

    /* Load: every variable of a tile of observations in one visit */
    columns = malloc(nload * sizeof(*columns));
    if (!columns)
    {
        ctools_error_alloc(CCOLLAPSE_MODULE);
        rc = 920;
        goto done;
    }
    for (size_t v = 0; v < nload; v++)
        columns[v] = (int)v + 1;
    status = ctools_data_load(&loaded, columns, nload, 0, 0,
                              (config.marker ? CTOOLS_LOAD_MARKER : CTOOLS_LOAD_SKIP_IF) |
                                  CTOOLS_LOAD_NO_SORT_ORDER | CTOOLS_LOAD_ROW_TILES);
    if (status != STATA_OK)
    {
        if (status != STATA_ERR_CANCELLED)
            ctools_error(CCOLLAPSE_MODULE, "failed to load data");
        rc = ctools_stata_rc(status);
        goto done;
    }
    seconds[0] = ctools_timer_seconds() - started;
    if (loaded.data.nobs == 0)
    {
        rc = 2000;
        goto done;
    }

    /* Sample: weights and cw */
    double mark = ctools_timer_seconds();
    rc = cc_apply_sample(&loaded.data, &config, threads);
    if (rc)
    {
        if (rc == 920)
            ctools_error_alloc(CCOLLAPSE_MODULE);
        goto done;
    }
    const size_t nobs = loaded.data.nobs;
    seconds[1] = ctools_timer_seconds() - mark;
    if (ctools_cancel_requested())
    {
        rc = 1;
        goto done;
    }

    /* Groups */
    mark = ctools_timer_seconds();
    status = cc_group_rows(&loaded.data, config.nby, threads, &groups);
    int quantiles = 0;
    for (size_t t = 0; t < config.ntargets; t++)
        quantiles |= cc_stat_is_quantile(config.targets[t].stat);
    if (status == STATA_OK && quantiles)
        status = cc_groups_order(&groups, threads);
    if (status != STATA_OK)
        goto failed;
    seconds[2] = ctools_timer_seconds() - mark;
    if (ctools_cancel_requested())
    {
        rc = 1;
        goto done;
    }

    /* The collapsed dataset: by() columns, then one column per target */
    mark = ctools_timer_seconds();
    const size_t ngroups = groups.ngroups;
    out.nvars = config.nby + config.ntargets;
    out.nobs = ngroups;
    out.vars = ctools_safe_cacheline_alloc2(out.nvars, sizeof(*out.vars));
    values = calloc(config.ntargets, sizeof(*values));
    string_rows = calloc(config.ntargets, sizeof(*string_rows));
    if (out.vars)
        memset(out.vars, 0, out.nvars * sizeof(*out.vars));
    else
        out.nvars = 0;
    if (!out.vars || !values || !string_rows)
    {
        status = STATA_ERR_MEMORY;
        goto failed;
    }
    for (size_t k = 0; k < config.nby && status == STATA_OK; k++)
    {
        const stata_variable *src = &loaded.data.vars[k];
        status = src->type == STATA_TYPE_DOUBLE
                     ? cc_gather_numbers(src, groups.first, ngroups, &out.vars[k])
                     : cc_gather_strings(src, groups.first, ngroups, &out.vars[k]);
    }
    for (size_t t = 0; t < config.ntargets && status == STATA_OK; t++)
    {
        const stata_variable *src = &loaded.data.vars[config.nby + (size_t)config.targets[t].source];
        stata_variable *dst = &out.vars[config.nby + t];
        if (src->type == STATA_TYPE_DOUBLE)
        {
            dst->type = STATA_TYPE_DOUBLE;
            dst->nobs = ngroups;
            dst->data.dbl = ctools_safe_cacheline_alloc2(ngroups, sizeof(double));
            values[t] = dst->data.dbl;
            if (!values[t])
                status = STATA_ERR_MEMORY;
        }
        else
        {
            string_rows[t] = malloc(ngroups * sizeof(**string_rows));
            if (!string_rows[t])
                status = STATA_ERR_MEMORY;
        }
    }
    if (status != STATA_OK)
        goto failed;

    cc_job job;
    memset(&job, 0, sizeof(job));
    job.data = &loaded.data;
    job.nby = config.nby;
    job.nsources = config.nsources;
    job.weight = weighted ? loaded.data.vars[config.nby + config.nsources].data.dbl : NULL;
    job.wtype = config.wtype;
    job.groups = &groups;
    job.targets = config.targets;
    job.ntargets = config.ntargets;
    job.values = values;
    job.string_rows = string_rows;
    job.threads = threads;

    status = cc_compute_stats(&job);
    if (status == STATA_OK && quantiles)
        status = cc_compute_quantiles(&job);
    for (size_t t = 0; t < config.ntargets && status == STATA_OK; t++)
    {
        if (!string_rows[t])
            continue;
        const stata_variable *src = &loaded.data.vars[config.nby + (size_t)config.targets[t].source];
        status = cc_gather_strings(src, string_rows[t], ngroups, &out.vars[config.nby + t]);
    }
    if (status != STATA_OK)
        goto failed;
    seconds[3] = ctools_timer_seconds() - mark;
    seconds[4] = ctools_timer_seconds() - started;

    rc = cc_save_results(&groups, nobs, threads, seconds);
    if (!rc)
    {
        cc_state.out = out;
        cc_state.ready = 1;
        stata_data_init(&out);
    }
    goto done;

failed:
    if (status == STATA_ERR_MEMORY)
        ctools_error_alloc(CCOLLAPSE_MODULE);
    else if (status == STATA_ERR_UNSUPPORTED_TYPE)
        ctools_error(CCOLLAPSE_MODULE, "too many observations");
    else if (status != STATA_ERR_CANCELLED)
        ctools_error(CCOLLAPSE_MODULE, "internal error");
    rc = ctools_stata_rc(status);

done:
    if (string_rows)
        for (size_t t = 0; t < config.ntargets; t++)
            free(string_rows[t]);
    free(string_rows);
    free(values);
    free(columns);
    free(config.targets);
    cc_groups_free(&groups);
    stata_data_free(&out);
    /* the loaded columns go back on the release helper, after the call */
    ctools_filtered_data_release(&loaded);
    return rc;
}

/* ===========================================================================
   store
   =========================================================================== */

static ST_retcode cc_store(void)
{
    if (!cc_state.ready)
    {
        ctools_error(CCOLLAPSE_MODULE, "nothing to store");
        return 198;
    }
    stata_data *out = &cc_state.out;
    if ((size_t)SF_nvars() != out->nvars || (size_t)SF_nobs() < out->nobs)
    {
        ctools_error(CCOLLAPSE_MODULE, "the dataset does not match the results");
        return 198;
    }
    for (size_t v = 0; v < out->nvars; v++)
    {
        int is_string = SF_var_is_string((int)v + 1) != 0;
        if (is_string != (out->vars[v].type == STATA_TYPE_STRING))
        {
            ctools_error(CCOLLAPSE_MODULE, "the dataset does not match the results");
            return 198;
        }
    }
    double started = ctools_timer_seconds();
    stata_retcode status = ctools_data_store(out, 1);
    if (status != STATA_OK)
    {
        ctools_error(CCOLLAPSE_MODULE, "failed to store the results");
        return ctools_stata_rc(status);
    }
    ST_retcode rc = ctools_scal_save("_ccollapse_time_store", ctools_timer_seconds() - started);
    ccollapse_cleanup_state();
    return rc;
}

/* ===========================================================================
   Entry point
   =========================================================================== */

ST_retcode ccollapse_main(const char *args)
{
    char op[16] = "";
    if (args)
        sscanf(args, "%15s", op);
    if (!strcmp(op, "compute"))
        return cc_compute(args);
    if (!strcmp(op, "store"))
        return cc_store();
    ccollapse_cleanup_state();
    ctools_error(CCOLLAPSE_MODULE, "unknown operation");
    return 198;
}
