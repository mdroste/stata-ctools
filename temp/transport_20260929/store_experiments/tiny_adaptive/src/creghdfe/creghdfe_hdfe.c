/*
 * creghdfe_hdfe.c
 *
 * Global HDFE state (g_state) and its cleanup
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_hdfe.h"
#include "creghdfe_utils.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_config.h"
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */

/* Define the global state pointer */
HDFE_State *g_state = NULL;

/*
 * Clean up global state
 */
void cleanup_state(void)
{
    if (g_state == NULL) return;

    /* Free weights (ctools_hdfe_state_cleanup does not free weights) */
    if (g_state->weights != NULL) free(g_state->weights);
    g_state->weights = NULL;

    /* Free factors, thread buffers, etc. */
    ctools_hdfe_state_cleanup(g_state);

    free(g_state);
    g_state = NULL;
}

/* ============================================================================
 * Cleanup function for ctools_cleanup system
 * ============================================================================ */

void creghdfe_cleanup_state(void)
{
    /* Reuse existing cleanup function */
    cleanup_state();
}
