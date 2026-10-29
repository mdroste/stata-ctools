"""Statistics contracts, scalar publication failures, and parallel grouping fallback."""

import tempfile
from pathlib import Path

from native.support.runner import harness, run_test


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-match-winsor-") as tmp:
        for omp in (False, True):
            run_test(
                Path(tmp),
                "match_winsor_" + str(omp),
                harness("match_winsor_optimizations.c"),
                openmp=omp,
            )


if __name__ == "__main__":
    main()
