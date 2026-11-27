#include <stddef.h>
#include "cpsmatch_internal.h"

int cpsmatch_parse_options(const char *args, cpsmatch_options *opts)
{
    const char *p = args;
    char *endptr;
    int val;

    /* Initialize defaults */
    memset(opts, 0, sizeof(cpsmatch_options));
    opts->n_neighbors = 1;
    opts->caliper = 0.0; /* 0 means no caliper */
    opts->with_replace = 1;
    opts->ties = 0;
    opts->descending = 0;
    opts->verbose = 0;
    opts->kernel = KERNEL_EPAN;
    opts->bandwidth = 0.06;

    /* Parse space-separated arguments */
    /* Expected format: treat_idx pscore_idx outcome_idx out_weight_idx out_match_idx
     *                  out_support_idx nobs method neighbor caliper common
     *                  kernel bwidth with_replace ties descending verbose */

    /* Parse treat_idx */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->treat_idx = val;
    p = endptr;

    /* Parse pscore_idx */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->pscore_idx = val;
    p = endptr;

    /* Parse outcome_idx */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->outcome_idx = val;
    p = endptr;

    /* Parse out_weight_idx */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->out_weight_idx = val;
    p = endptr;

    /* Parse out_match_idx */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->out_match_idx = val;
    p = endptr;

    /* Parse out_support_idx */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->out_support_idx = val;
    p = endptr;

    /* Parse nobs */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    size_t nobs_val = strtoul(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->nobs = nobs_val;
    p = endptr;

    /* Parse method */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->method = (match_method_t)val;
    p = endptr;

    /* Parse n_neighbors */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->n_neighbors = val;
    p = endptr;

    /* Parse caliper */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    double dval = strtod(p, &endptr);
    if (endptr == p)
        return -1;
    opts->caliper = dval;
    p = endptr;

    /* Parse common_support */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->common_support = val;
    p = endptr;

    /* Parse kernel */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->kernel = (kernel_type_t)val;
    p = endptr;

    /* Parse bandwidth */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    dval = strtod(p, &endptr);
    if (endptr == p)
        return -1;
    opts->bandwidth = dval;
    p = endptr;

    /* Parse with_replace */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->with_replace = val;
    p = endptr;

    /* Parse ties */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->ties = val;
    p = endptr;

    /* Parse descending */
    while (*p == ' ')
        p++;
    if (*p == '\0')
        return -1;
    val = strtol(p, &endptr, 10);
    if (endptr == p)
        return -1;
    opts->descending = val;
    p = endptr;

    /* Parse verbose (optional) */
    while (*p == ' ')
        p++;
    if (*p != '\0')
    {
        val = strtol(p, &endptr, 10);
        if (endptr != p)
        {
            opts->verbose = val;
            p = endptr;
        }
    }

    /* psmatch2 accepts descending only when matching without replacement */
    if (opts->descending && opts->with_replace)
        return -1;

    /* Optional wire flags, the sample marker always the final plugin
     * variable: a bare name or name=integer. */
    static const struct
    {
        const char *name;
        size_t offset;
    } keys[] = {
        {"presup", offsetof(cpsmatch_options, presupport_idx)},
        {"arm", offsetof(cpsmatch_options, arm_idx)},
        {"nnvar", offsetof(cpsmatch_options, out_nn_idx)},
        {"ids", offsetof(cpsmatch_options, neighbor_ids_idx)},
        {"pdif", offsetof(cpsmatch_options, pdif_idx)},
        {"matched", offsetof(cpsmatch_options, matched_idx)},
        {"self", offsetof(cpsmatch_options, self_idx)},
        {"ai", offsetof(cpsmatch_options, ai)},
        {"ny", offsetof(cpsmatch_options, ny)},
        {"yeff", offsetof(cpsmatch_options, yeff_idx)},
        {"nx", offsetof(cpsmatch_options, nx)},
        {"xcol", offsetof(cpsmatch_options, xcol_idx)},
        {"xself", offsetof(cpsmatch_options, xself_idx)},
        {"supscore", offsetof(cpsmatch_options, supscore_idx)},
    };
    static const struct
    {
        const char *name;
        size_t offset;
    } flags[] = {
        {"touse", offsetof(cpsmatch_options, touse)},
        {"stats", offsetof(cpsmatch_options, statistics)},
        {"ate", offsetof(cpsmatch_options, ate)},
        {"altvar", offsetof(cpsmatch_options, altvariance)},
        {"w2float", offsetof(cpsmatch_options, float_squares)},
        {"mkarm", offsetof(cpsmatch_options, make_arm)},
        {"fill", offsetof(cpsmatch_options, fill)},
        {"outids", offsetof(cpsmatch_options, outside_ids)},
        {"single", offsetof(cpsmatch_options, single)},
        {"effects", offsetof(cpsmatch_options, effects)},
        {"samplevar", offsetof(cpsmatch_options, samplevar)},
        {"warn", offsetof(cpsmatch_options, warn)},
        {"kdirect", offsetof(cpsmatch_options, kernel_direct)},
        {"radall", offsetof(cpsmatch_options, radius_all)},
    };
    while (*p)
    {
        while (*p == ' ') p++;
        size_t length = strcspn(p, " ");
        int flagged = 0;
        for (size_t k = 0; k < sizeof(flags) / sizeof(flags[0]); k++)
        {
            if (length == strlen(flags[k].name) && !strncmp(p, flags[k].name, length))
            {
                *(int *)((char *)opts + flags[k].offset) = 1;
                flagged = 1;
            }
        }
        if (flagged)
            ;
        else
        {
            for (size_t k = 0; k < sizeof(keys) / sizeof(keys[0]); k++)
            {
                size_t n = strlen(keys[k].name);
                if (length > n + 1 && !strncmp(p, keys[k].name, n) && p[n] == '=')
                {
                    long value = strtol(p + n + 1, &endptr, 10);
                    if (endptr != p + length || value < 0 || value > 1000000) return -1;
                    *(int *)((char *)opts + keys[k].offset) = (int)value;
                }
            }
        }
        p += length;
    }
    if (opts->n_neighbors < 1) return -1;

    return 0;
}
