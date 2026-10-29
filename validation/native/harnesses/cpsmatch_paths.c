/* Independent scan/enumeration oracles for the optimized matching paths.
 * The tracked allocator is prepended by the Python runner. */
#include "cpsmatch/cpsmatch_internal.h"
#include "cpsmatch/cpsmatch_index.h"
#include "cpsmatch/cpsmatch_kernels.c"
#include "cpsmatch/cpsmatch_radius.c"

static ST_plugin host;
ST_plugin *_stata_ = &host;
int ctools_get_max_threads(void)
{
    return 4;
}

static void availability(void)
{
    enum
    {
        N = 1024
    };
    perm_idx_t previous[N + 1], next[N + 1];
    int live[N];
    for (size_t i = 0; i <= N; i++)
        previous[i] = next[i] = (perm_idx_t)i;
    for (size_t i = 0; i < N; i++)
        live[i] = 1;
    for (size_t step = 0; step <= N; step++)
    {
        for (size_t q = 0; q <= N; q++)
        {
            size_t r = q;
            while (r < N && !live[r])
                r++;
            size_t l = q;
            while (l > 0 && !live[l - 1])
                l--;
            assert(cpsmatch_available_find(next, (perm_idx_t)q) == r);
            assert(cpsmatch_available_find(previous, (perm_idx_t)q) == l);
        }
        if (step == N)
            break;
        size_t deleted = step * 733 % N;
        live[deleted] = 0;
        previous[deleted + 1] = cpsmatch_available_find(previous, (perm_idx_t)deleted);
        next[deleted] = cpsmatch_available_find(next, (perm_idx_t)(deleted + 1));
    }
}

static void radius(int shape, int inject)
{
    enum
    {
        NC = 513,
        NT = 777,
        N = NC + NT
    };
    sorted_control_t controls[NC];
    double ps[NT], weights[N] = {0}, ids[N], nn[N] = {0}, support[N];
    double expected[N] = {0};
    size_t ci[NC], ti[NT];
    perm_idx_t obs[N];
    for (size_t i = 0; i < N; i++)
    {
        obs[i] = (perm_idx_t)(3 * i + 1);
        ids[i] = support[i] = host.missval;
    }
    for (size_t j = 0; j < NC; j++)
    {
        ci[j] = j;
        controls[j] = (sorted_control_t){0};
        controls[j].orig_idx = NC - 1 - j;
        controls[j].pscore = shape == 0 ? .5 : (double)(j / 3) / 171;
    }
    for (size_t t = 0; t < NT; t++)
    {
        ti[t] = NC + t;
        ps[t] = shape == 0 ? .5 : shape == 4 ? 2.0 : (double)((t * 173) % 1000) / 1000;
    }
    double width = shape == 2 ? 1e-12 : shape == 3 ? .03 : .25;
    int common = shape == 3;
    fail_at = inject;
    calls = 0;
    int matched = 0, off = 0;
    int rc = cpsmatch_match_radius(controls, NC, ps, ti, NT, ci, obs, width, common, .2, .8,
                                   weights, ids, nn, support, &matched, &off);
    if (inject)
    {
        assert(rc == 920 && matched == 0 && off == 0);
        for (size_t i = 0; i < N; i++)
            assert(weights[i] == 0 && ids[i] == host.missval);
        assert(alive == 0);
        return;
    }
    assert(rc == 0);
    int expected_matched = 0;
    for (size_t t = 0; t < NT; t++)
    {
        size_t start = NC, count = 0;
        if (!common || (ps[t] >= .2 && ps[t] <= .8))
        {
            for (size_t j = 0; j < NC; j++)
                if (controls[j].pscore >= ps[t] - width && controls[j].pscore <= ps[t] + width)
                {
                    if (!count)
                        start = j;
                    count++;
                }
        }
        assert(nn[ti[t]] == (double)count);
        assert(support[ti[t]] == (count ? 1 : 0));
        if (!count)
        {
            assert(ids[ti[t]] == host.missval);
            continue;
        }
        expected_matched++;
        assert(ids[ti[t]] == obs[ci[controls[start].orig_idx]]);
        expected[ti[t]] = 1;
        for (size_t j = start; j < start + count; j++)
            expected[ci[controls[j].orig_idx]] += 1.0 / count;
    }
    assert(matched == expected_matched && off == NT - expected_matched);
    for (size_t i = 0; i < N; i++)
    {
        if (expected[i] == 0)
            assert(weights[i] == 0);
        else
            assert(weights[i] > 0 && fabs(weights[i] - expected[i]) / expected[i] < 1e-7);
    }
    assert(alive == 0);
}

int main(void)
{
    host.missval = 8.98846567431158e307;
    int counts[] = {1, 2, 3, 17, 1023, 2147483647};
    for (size_t i = 0; i < sizeof(counts) / sizeof(counts[0]); i++)
    {
        double weight = 1.0 / counts[i];
        cpsmatch_radius_sum small = cpsmatch_radius_encode(weight);
        assert(cpsmatch_radius_value(small) == weight);
        cpsmatch_radius_sum big = cpsmatch_radius_encode(1.0);
        for (int bit = 0; bit < 31; bit++)
            big = cpsmatch_radius_add(big, big);
        cpsmatch_radius_sum recovered = cpsmatch_radius_sub(cpsmatch_radius_add(big, small), big);
        assert(recovered.hi == small.hi && recovered.lo == small.lo);
        recovered = cpsmatch_radius_add(cpsmatch_radius_sub(small, big), big);
        assert(recovered.hi == small.hi && recovered.lo == small.lo);
    }
    availability();
    for (int shape = 0; shape < 5; shape++)
        for (int injection = 0; injection <= 2; injection++)
            radius(shape, injection);
    puts("CPSMATCH_PATHS_NATIVE_COMPLETE");
    return 0;
}
