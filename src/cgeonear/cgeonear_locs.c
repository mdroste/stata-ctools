/* Base and neighbor locations as geonear_prep_locations makes them. */
#include <stdlib.h>
#include <string.h>
#include "cgeonear_internal.h"
#include "ctools_types.h"
#include "ctools_order.h"
#include "ctools_runtime.h"

void cgeo_locs_init(cgeo_locs *l)
{
    memset(l, 0, sizeof(*l));
}

void cgeo_locs_free(cgeo_locs *l)
{
    free(l->id_num);
    free(l->id_str);
    free(l->id_buf);
    free(l->lat);
    free(l->lon);
    cgeo_locs_init(l);
}

int cgeo_same_id(const cgeo_locs *a, size_t i, const cgeo_locs *b, size_t j)
{
    if (a->id_is_str)
        return !strcmp(a->id_str[i], b->id_str[j]);
    return a->id_num[i] == b->id_num[j];
}

int cgeo_cmp_id(const cgeo_locs *a, size_t i, const cgeo_locs *b, size_t j)
{
    if (a->id_is_str)
        return strcmp(a->id_str[i], b->id_str[j]);
    double x = a->id_num[i], y = b->id_num[j];
    return x < y ? -1 : (x > y ? 1 : 0);
}

/* Stata's missing(): a missing code for numbers, "" for strings */
static int value_missing(const stata_variable *v, size_t row)
{
    if (v->type == STATA_TYPE_DOUBLE)
        return v->data.dbl[row] >= SV_missval;
    const char *s = v->data.str[row];
    return !s || !*s;
}

static int value_equal(const stata_variable *v, size_t r1, size_t r2)
{
    if (v->type == STATA_TYPE_DOUBLE)
        return v->data.dbl[r1] == v->data.dbl[r2];
    return !strcmp(v->data.str[r1], v->data.str[r2]);
}

int cgeo_locs_read(cgeo_locs *l, uint32_t **obs_loc_out, size_t *nobs_out, int *in_order,
                   cgeo_fail *fail)
{
    *fail = CGEO_FAIL_NONE;
    *in_order = 0;
    cgeo_locs_free(l);
    if (obs_loc_out)
        *obs_loc_out = NULL;
    size_t nobs = (size_t)SF_nobs();
    *nobs_out = nobs;
    if (SF_nvars() < 3)
        return 198;
    if (nobs >= CGEO_NONE)
        return 920;
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);
    uint32_t *obs_loc = NULL;
    size_t *starts = NULL;
    int rc = 0;
    if (!nobs)
    {
        *fail = CGEO_FAIL_NOOBS;
        return 459;
    }
    int vars[3] = {1, 2, 3};
    rc = ctools_stata_rc(ctools_data_load(&fd, vars, 3, 1, nobs, CTOOLS_LOAD_SKIP_IF));
    if (rc)
        goto done;
    if (fd.data.nobs != nobs)
    {
        rc = 459;
        goto done;
    }
    const stata_variable *id = &fd.data.vars[0], *lat = &fd.data.vars[1], *lon = &fd.data.vars[2];
    /* geonear sorts by the identifier and keeps one observation of each */
    int key = 0;
    rc = ctools_stata_rc(ctools_order_stable(&fd.data, &key, 1));
    if (rc)
        goto done;
    /* a stable sort leaves rows that are already in order where they are */
    {
        size_t i = 0;
        while (i < nobs && fd.data.sort_order[i] == i)
            i++;
        *in_order = i == nobs;
    }
    starts = malloc((nobs + 1) * sizeof(*starts));
    if (obs_loc_out)
        obs_loc = malloc(nobs * sizeof(*obs_loc));
    if (!starts || (obs_loc_out && !obs_loc))
    {
        rc = 920;
        goto done;
    }
    size_t nruns = ctools_group_starts_parallel(&fd.data, fd.data.sort_order, &key, 1, starts);
    const perm_idx_t *order = fd.data.sort_order;
    /* First pass: which runs have coordinates, and are they constant */
    size_t nloc = 0, idbytes = 0;
    int bad = 0;
    for (size_t g = 0; g < nruns && !bad; g++)
    {
        size_t first = order[starts[g]];
        if (value_missing(id, first))
            continue;
        size_t ref = SIZE_MAX;
        for (size_t p = starts[g]; p < starts[g + 1]; p++)
        {
            size_t r = order[p];
            if (value_missing(lat, r) || value_missing(lon, r))
                continue;
            if (ref == SIZE_MAX)
                ref = r;
            else if (!value_equal(lat, r, ref) || !value_equal(lon, r, ref))
            {
                bad = 1;
                break;
            }
        }
        if (ref != SIZE_MAX)
        {
            nloc++;
            if (id->type != STATA_TYPE_DOUBLE)
                idbytes += strlen(id->data.str[first]) + 1;
        }
    }
    if (bad)
    {
        *fail = CGEO_FAIL_NOTCONST;
        rc = 459;
        goto done;
    }
    if (!nloc)
    {
        *fail = CGEO_FAIL_NOOBS;
        rc = 459;
        goto done;
    }
    if (lat->type != STATA_TYPE_DOUBLE || lon->type != STATA_TYPE_DOUBLE)
    {
        *fail = CGEO_FAIL_STRING;
        rc = 7;
        goto done;
    }
    l->n = nloc;
    l->id_is_str = id->type != STATA_TYPE_DOUBLE;
    l->lat = malloc(nloc * sizeof(double));
    l->lon = malloc(nloc * sizeof(double));
    if (l->id_is_str)
    {
        l->id_str = malloc(nloc * sizeof(char *));
        l->id_buf = malloc(idbytes ? idbytes : 1);
    }
    else
        l->id_num = malloc(nloc * sizeof(double));
    if (!l->lat || !l->lon || (l->id_is_str ? !l->id_str || !l->id_buf : !l->id_num))
    {
        rc = 920;
        goto done;
    }
    /* Second pass: the locations, and each observation's location */
    size_t k = 0, used = 0;
    for (size_t g = 0; g < nruns; g++)
    {
        size_t first = order[starts[g]];
        uint32_t loc = CGEO_NONE;
        if (!value_missing(id, first))
        {
            for (size_t p = starts[g]; p < starts[g + 1]; p++)
            {
                size_t r = order[p];
                if (value_missing(lat, r) || value_missing(lon, r))
                    continue;
                loc = (uint32_t)k;
                l->lat[k] = lat->data.dbl[r];
                l->lon[k] = lon->data.dbl[r];
                if (l->id_is_str)
                {
                    size_t len = strlen(id->data.str[first]) + 1;
                    memcpy(l->id_buf + used, id->data.str[first], len);
                    l->id_str[k] = l->id_buf + used;
                    used += len;
                }
                else
                    l->id_num[k] = id->data.dbl[first];
                k++;
                break;
            }
        }
        if (obs_loc)
            for (size_t p = starts[g]; p < starts[g + 1]; p++)
                obs_loc[fd.obs_map[order[p]] - 1] = loc;
    }
    /* geonear checks the range of the latitudes, then of the longitudes */
    for (size_t i = 0; i < nloc; i++)
        if (l->lat[i] > 90 || l->lat[i] < -90)
        {
            *fail = CGEO_FAIL_LAT;
            rc = 198;
            goto done;
        }
    for (size_t i = 0; i < nloc; i++)
        if (l->lon[i] > 180 || l->lon[i] < -180)
        {
            *fail = CGEO_FAIL_LON;
            rc = 198;
            goto done;
        }
done:
    if (rc)
    {
        cgeo_locs_free(l);
        free(obs_loc);
        obs_loc = NULL;
    }
    if (obs_loc_out)
        *obs_loc_out = obs_loc;
    free(starts);
    ctools_filtered_data_free(&fd);
    return rc;
}
