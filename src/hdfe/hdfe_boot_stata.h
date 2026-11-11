/*
 * hdfe_boot_stata.h
 *
 * The Stata side of the pairs bootstrap shared by the HDFE estimators: a
 * chunk's draws read from the plugin call's variables, and the replicates'
 * coefficients stored for saving(). Design:
 * docs/development/BOOTSTRAP_PLAN.md (sections 5.8, 6.9 and 11.1).
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_HDFE_BOOT_STATA_H
#define CTOOLS_HDFE_BOOT_STATA_H

#include "hdfe_boot.h"

/*
 * Unit multiplicities of count replicates (count x d->U, row-major,
 * allocated into *out) from the variables of the plugin call.
 *
 * supplied == 0 (draws(stata)): the uniforms that `_ctools_hdfe_boot
 * uniforms` packed, replicate r's hdfe_boot_draw_units(d) draws and then
 * as many tie-breaks, mapped by hdfe_boot_map_stata.
 *
 * supplied != 0 (the hidden _bsdraws()): variable r + 1 of the call holds
 * replicate r's multiplicity on each row, read at observation
 * obs_of_row[i] of compact row i, the same on every row of a unit, a
 * nonnegative integer. With frequency-weighted rows resampled, a row's
 * value is the number of its observations drawn.
 *
 * cmd names the command in error messages. Returns 0, 198 or 920.
 */
ST_retcode hdfe_boot_read_draws(const hdfe_boot_design *d, int cluster_mode, int supplied,
                                ST_int count, ST_int N, const ST_int *unit_of_row,
                                const ST_int *obs_of_row, const char *cmd, ST_int **out);

/* Replicates first .. first + count - 1 of acc into rows 1 .. count of the
 * Stata matrix name (at least count x K_b), missing where a replicate has
 * no complete coefficients. Returns 0 or the SPI's error. */
ST_retcode hdfe_boot_store_rows(const hdfe_boot_acc *acc, ST_int first, ST_int count,
                                const char *name);

/*
 * For draws(stata) with observations resampled within strata: bsample, as
 * the bootstrap prefix calls it, takes hdfe_boot_draw_units(d) uniforms more
 * after a replicate's draws only when expand added rows, that is when some
 * stratum drew one of its units twice (BOOTSTRAP_REVIEW.md, M4). The ado
 * decides this from the draws; it needs the strata that can draw a unit
 * twice, those of two or more units in the drawn frame (fweight-expanded
 * when fweighted rows are resampled). Saves scalar __<cmd>_bs_dupK, their
 * number K, and for j = 1 .. K __<cmd>_bs_dupo<j> (0-based offset of the
 * stratum's first unit in the draw order) and __<cmd>_bs_dupn<j> (its
 * units). With more than HDFE_BOOT_DUP_MAX such strata, dupK is missing:
 * a replicate in which none of them draws a unit twice then has probability
 * at most 2^-HDFE_BOOT_DUP_MAX, and the extra uniforms are always taken.
 * Returns 0 or the SPI's error.
 */
#define HDFE_BOOT_DUP_MAX 64
ST_retcode hdfe_boot_save_dup_strata(const hdfe_boot_design *d, const char *cmd);

#endif /* CTOOLS_HDFE_BOOT_STATA_H */
