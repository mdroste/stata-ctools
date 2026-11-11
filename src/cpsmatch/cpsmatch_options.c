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
    };
    while (*p)
    {
        while (*p == ' ') p++;
        size_t length = strcspn(p, " ");
        if (length == 5 && !strncmp(p, "touse", 5)) opts->touse = 1;
        else if (length == 5 && !strncmp(p, "stats", 5)) opts->statistics = 1;
        else if (length == 3 && !strncmp(p, "ate", 3)) opts->ate = 1;
        else if (length == 6 && !strncmp(p, "altvar", 6)) opts->altvariance = 1;
        else if (length == 7 && !strncmp(p, "w2float", 7)) opts->float_squares = 1;
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
