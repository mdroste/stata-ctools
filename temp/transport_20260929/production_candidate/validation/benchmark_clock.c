/* Benchmark-only clock for monotonic elapsed-time measurement.
 * This plugin is not part of the distributed ctools plugin. */
#include <string.h>
#include <time.h>
#include "stplugin.h"

STDLL stata_call(int argc, char *argv[])
{
    static struct timespec start;
    struct timespec now;
    if (argc != 1 || clock_gettime(CLOCK_MONOTONIC, &now)) return 198;
    if (!strcmp(argv[0], "start")) { start = now; return 0; }
    if (strcmp(argv[0], "stop")) return 198;
    double elapsed = (double)(now.tv_sec-start.tv_sec) + (now.tv_nsec-start.tv_nsec)*1e-9;
    return SF_scal_save("__ctools_perf_seconds", elapsed);
}
