/*
 * creghdfe_hdfe.h
 *
 * Global HDFE state (g_state) and its cleanup
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_HDFE_H
#define CREGHDFE_HDFE_H

#include "creghdfe_types.h"

/* Private to creghdfe; IV and PPML own independent state. */
extern HDFE_State *g_state;

/* Clean up global state */
void cleanup_state(void);

#endif /* CREGHDFE_HDFE_H */
