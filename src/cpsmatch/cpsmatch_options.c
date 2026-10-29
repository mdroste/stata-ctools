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

    /* Optional wire flags. The stored _treated column (stats) precedes the
     * sample marker, which always remains the final plugin variable. */
    while (*p)
    {
        while (*p == ' ') p++;
        size_t length = strcspn(p, " ");
        if (length == 5 && !strncmp(p, "touse", 5)) opts->touse = 1;
        if (length == 5 && !strncmp(p, "stats", 5)) opts->statistics = 1;
        p += length;
    }

    return 0;
}
