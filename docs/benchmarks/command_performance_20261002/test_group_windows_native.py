"""Sanitized checks of the population moments used by rangestat's tiny windows."""

import tempfile
from pathlib import Path

from native.support.runner import harness, run_test


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-group-moments-") as tmp:
        run_test(
            Path(tmp),
            "group_moments",
            harness("group_moments.c"),
            ["crangestat/crangestat_kernels.c"],
        )


if __name__ == "__main__":
    main()
