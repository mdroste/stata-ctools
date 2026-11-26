/* Cooperative cancellation. Only the command's calling thread polls the host;
 * workers observe a latched atomic flag. Begin/end bracket a joined command. */
#ifndef CTOOLS_CANCEL_H
#define CTOOLS_CANCEL_H
#include <stdatomic.h>

typedef int (*ctools_cancel_poll)(void);
void ctools_cancel_begin(ctools_cancel_poll poll);
void ctools_cancel_end(void);
int ctools_cancel_requested(void);

/* Use after an OpenMP for nowait, as the last statement of the parallel
 * region. The command's calling thread stays here, polling the host about
 * once a millisecond, until every member of the team has arrived, also when
 * its own columns/replicates finished first; the member that arrives last
 * wakes it at once. Every other member returns immediately and waits at the
 * region's barrier, so nothing after this call may rely on the team being
 * done. Each team (including a nested team) owns a separate zero-initialized
 * counter. */
typedef struct { atomic_int arrived; } ctools_cancel_team;
void ctools_cancel_wait(ctools_cancel_team *team);
#endif
