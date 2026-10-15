"""Allocation/store failures and sparse row restoration in production cwinsor."""

import tempfile
from pathlib import Path

from native.support.runner import harness, run_test


def main():
    allocator = Path(__file__).with_name("native") / "memory_safety/alloc.h"
    source = allocator.read_text() + "\n" + harness("winsor_transport.c")
    with tempfile.TemporaryDirectory(prefix="ctools-winsor-transport-") as tmp:
        run_test(Path(tmp), "winsor_transport", source, openmp=True)


if __name__ == "__main__":
    main()
