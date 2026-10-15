#include "crangestat/crangestat_internal.h"
#include <assert.h>

static ST_plugin host;
ST_plugin *_stata_ = &host;

int main(void)
{
    host.missval = 8.98846567431158e307;
    double x[] = {1.0, 2.0, 3.0, host.missval};
    double constant[] = {2.0, 2.0};
    for (int stat = STAT_SKEWNESS; stat <= STAT_KURTOSIS; stat++)
    {
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 0, 0, 0) == host.missval);
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 1, 0, 0) == host.missval);
        assert(crangestat_compute_window_stat_simple(stat, constant, 0, 2, 0, 0) == host.missval);
        double pair = stat == STAT_SKEWNESS ? 0.0 : 1.0;
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 2, 0, 0) == pair);
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 4, 1, 1) == pair);
        double triple = stat == STAT_SKEWNESS ? 0.0 : 1.5;
        assert(crangestat_compute_window_stat_simple(stat, x, 0, 4, 0, 0) == triple);
    }
    return 0;
}
