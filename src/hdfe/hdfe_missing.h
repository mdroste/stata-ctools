#ifndef CTOOLS_HDFE_MISSING_H
#define CTOOLS_HDFE_MISSING_H

#include "stplugin.h"
#include "../ctools_data_io.h"

/* Drop rows missing in any loaded numeric input, before FE remapping.
 * The caller loads inputs only, with NO_SORT_ORDER. marker is the private
 * Stata sample variable: clear dropped rows, preserving already-excluded rows.
 * Keep column storage/ownership and original row order; compact obs_map too.
 * No N-sized allocation or compaction when every input is nonmissing.
 * Returns a native Stata rc (0, 198, 920, or an SPI write error).
 */
ST_retcode hdfe_screen_missing(ctools_filtered_data *fd, ST_int marker);

/* Pre-expansion screening: numeric inputs followed by the private sample
 * marker. Stream bounded row blocks through SPI, updating only that marker;
 * do not materialize/copy the data or allocate an observation-sized mask. */
ST_retcode hdfe_markout_main(const char *args);

/* The levels of factor variables f1 ... fk on the rows of the marker (the
 * last plugin variable) as locals ctools_fvl_<j>, with ctools_fvl_ok_<j>
 * (hdfe_fvlevels.c) */
ST_retcode hdfe_fvlevels_main(const char *args);

#endif
