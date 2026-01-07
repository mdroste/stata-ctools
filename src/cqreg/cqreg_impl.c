/*
 * cqreg_impl.c
 *
 * Main dispatcher for C-accelerated quantile regression.
 * Part of the ctools suite.
 */

#include "cqreg_impl.h"
#include "cqreg_main.h"
#include "cqreg_types.h"
#include "../ctools_error.h"

#include <string.h>
#include <ctype.h>

/* Skip leading whitespace */
static const char *skip_whitespace(const char *s)
{
    while (*s && isspace((unsigned char)*s)) s++;
    return s;
}

/* Extract next word from string */
static const char *get_word(const char *s, char *buf, size_t bufsize)
{
    size_t i = 0;
    s = skip_whitespace(s);
    while (*s && !isspace((unsigned char)*s) && i < bufsize - 1) {
        buf[i++] = *s++;
    }
    buf[i] = '\0';
    return s;
}

/*
 * Main dispatcher for cqreg subcommands.
 */
ST_retcode cqreg_main(const char *args)
{
    char subcmd[64];
    const char *rest;

    if (args == NULL || *args == '\0') {
        ctools_error("cqreg", "No subcommand specified");
        return 198;
    }

    /* Parse subcommand */
    rest = get_word(args, subcmd, sizeof(subcmd));

    if (strcmp(subcmd, "full_regression") == 0) {
        return cqreg_full_regression(rest);
    }
    else if (strcmp(subcmd, "ipm_only") == 0) {
        /* For testing: run IPM solver without full VCE computation */
        return cqreg_full_regression(rest);  /* Same for now */
    }
    else {
        ctools_error("cqreg", "Unknown subcommand: %s", subcmd);
        return 198;
    }
}
