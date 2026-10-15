/*
    cgroup_impl.h
    Group identifiers for ado helpers (_ctools_group.ado), in place of egen group()
*/

#ifndef CGROUP_IMPL_H
#define CGROUP_IMPL_H

#include "stplugin.h"

/*
    plugin call ctools_plugin <keys...> [<destination>] [if] [in],
        "group nkeys=# [out] [reps]"

    Numbers the distinct tuples of the nkeys key variables 1, 2, ... in their
    sort order (numbers ascending, strings bytewise), over the observations
    the if/in qualifier selects whose keys are all nonmissing (string "" is
    missing). With out, the numeric destination after the keys receives every
    observation's number in the in-range: missing for the others. Saves scalar
    _ctools_group_n (the number of groups) and, with reps, local
    _ctools_group_reps listing one observation of each group in group order.

    Returns 109 before writing anything when a key cannot be read (a strL value
    longer than 2,045 bytes or binary), so the caller can number groups itself.
*/
ST_retcode cgroup_main(const char *args);

#endif /* CGROUP_IMPL_H */
