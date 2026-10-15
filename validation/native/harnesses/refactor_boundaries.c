#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <string.h>
#include "stplugin.h"
#include "ctools_rng.h"
#include "ctools_group.h"
#include "ctools_platform.h"
#include "ctools_runtime.h"
#include "formats/excel_address.h"

static ST_plugin mock;
ST_plugin *_stata_ = &mock;
static char message[1024];
static ST_retcode capture(char *text)
{
    assert(strlen(text) < sizeof(message));
    strcpy(message, text);
    return 0;
}

static atomic_int ready, go;
static ctools_prefetch_distances snapshots[16];
static void *worker(void *argument)
{
    size_t i = (size_t)argument;
    atomic_fetch_add(&ready, 1);
    while (!atomic_load(&go))
    {
    }
    snapshots[i] = ctools_get_prefetch_distances();
    double first = ctools_timer_seconds();
    for (int j = 0; j < 1000; j++)
        assert(ctools_timer_seconds() >= first);
    return NULL;
}

int main(void)
{
    /* Simultaneous first use must publish fully initialized process-wide state. */
    pthread_t threads[16];
    for (size_t i = 0; i < 16; i++)
        assert(!pthread_create(&threads[i], NULL, worker, (void *)i));
    while (atomic_load(&ready) != 16)
    {
    }
    atomic_store(&go, 1);
    for (int i = 0; i < 16; i++)
        assert(!pthread_join(threads[i], NULL));
    for (int i = 0; i < 16; i++)
    {
        assert(snapshots[i].near_dist > 0);
        assert(snapshots[i].near_dist == snapshots[0].near_dist);
        assert(snapshots[i].far_dist == snapshots[0].far_dist);
        assert(snapshots[i].stream == snapshots[0].stream);
        assert(snapshots[i].general == snapshots[0].general);
    }

    /* Published xoshiro256** transition for state {1,2,3,4}; no seeding ambiguity. */
    ctools_rng rng = {{1, 2, 3, 4}};
    assert(ctools_rng_uniform(&rng, 0) == 0);
    assert(ctools_rng_next(&rng) == UINT64_C(11520));
    assert(ctools_rng_next(&rng) == UINT64_C(0));
    assert(ctools_rng_next(&rng) == UINT64_C(1509978240));

    int col = -1, row = -1;
    assert(ctools_excel_address("XFD1048576", 16384, 1048576, false, &col, &row, NULL));
    assert(col == 16384 && row == 1048576);
    const char *bad[] = {"XFE1",
                         "A1048577",
                         "A0",
                         "A",
                         "1",
                         "$A$1",
                         "A1suffix",
                         "ZZZZZZZZZZZZZZ1",
                         "A9999999999999999999"};
    for (size_t i = 0; i < sizeof(bad) / sizeof(*bad); i++)
    {
        col = row = -1;
        assert(!ctools_excel_address(bad[i], 16384, 1048576, false, &col, &row, NULL));
        assert(col == -1 && row == -1);
    }
    const char *end;
    assert(ctools_excel_address("$az$12:B20", 16384, 1048576, true, &col, &row, &end));
    assert(col == 52 && row == 12 && !strcmp(end, ":B20"));

    double values[] = {0.0, -0.0, 8.98846567431158e307, 8.98846567431158e307,
                       nextafter(8.98846567431158e307, INFINITY)};
    char *strings[] = {NULL, "", "a", "a", "a"};
    stata_variable vars[2] = {{.type = STATA_TYPE_DOUBLE, .data.dbl = values},
                              {.type = STATA_TYPE_STRING, .data.str = strings}};
    stata_data data = {.nobs = 5, .nvars = 2, .vars = vars};
    ctools_group_range groups[5];
    size_t count;
    ctools_group_ranges(&data, NULL, NULL, 2, groups, &count);
    assert(count == 3 && groups[0].count == 2 && groups[1].count == 2);
    perm_idx_t order[] = {4, 3, 2, 1, 0};
    ctools_group_ranges(&data, order, NULL, 2, groups, &count);
    assert(count == 3 && groups[0].count == 1 && groups[1].start == 1 && groups[2].count == 2);

    mock.spouterr = capture;
    mock.spoutsml = capture;
    char oversized[4096];
    memset(oversized, 'x', sizeof(oversized) - 1);
    oversized[sizeof(oversized) - 1] = 0;
    ctools_error(oversized, "%s", oversized);
    assert(strlen(message) == 1023 && message[1022] == '\n');
    ctools_msg("unit", "%s", oversized);
    assert(strlen(message) == 1023 && !strncmp(message, "unit: ", 6) && message[1022] == '\n');
    ctools_verbose("unit", 0, "hidden");
    assert(strlen(message) == 1023);
    ctools_msg(NULL, "complete\n");
    assert(!strcmp(message, "complete\n"));
    return 0;
}
