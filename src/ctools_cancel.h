/* Cooperative cancellation. Only the command's calling thread polls the host;
 * workers observe a latched atomic flag. Begin/end bracket a joined command. */
#ifndef CTOOLS_CANCEL_H
#define CTOOLS_CANCEL_H
#include <stdatomic.h>

typedef int (*ctools_cancel_poll)(void);
void ctools_cancel_begin(ctools_cancel_poll poll);
void ctools_cancel_end(void);
int ctools_cancel_requested(void);

/* Use after an OpenMP for nowait, before the parallel region's final barrier.
 * The caller keeps polling even if its own columns/replicates finished first.
 * Each team (including a nested team) owns a separate zero-initialized counter. */
typedef struct { atomic_int arrived; } ctools_cancel_team;
void ctools_cancel_wait(ctools_cancel_team *team);
#endif
