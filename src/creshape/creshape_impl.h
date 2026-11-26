#ifndef CRESHAPE_IMPL_H
#define CRESHAPE_IMPL_H
#include "stplugin.h"

/*
 * creshape: reshape between wide and long form.
 *
 * Plugin calls, in order (see creshape_impl.c):
 *   creshape plan long|wide ...   stage the dataset, check it, plan the result
 *   creshape schema               description of the result from the ado's template
 *   creshape store                fill the result's variables
 *   creshape restore              put the staged observations back (after a failure)
 *   creshape clear                free everything
 */
ST_retcode creshape_main(const char *args);

/* Idempotent: frees the staged dataset and the plan. */
void creshape_cleanup_state(void);

#endif /* CRESHAPE_IMPL_H */
