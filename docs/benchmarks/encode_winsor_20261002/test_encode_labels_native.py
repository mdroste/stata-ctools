"""Sanitized byte-framing and I/O failure checks for the cencode label writer."""

import tempfile
from pathlib import Path

from native.support.runner import harness, run_test


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-encode-labels-") as tmp:
        run_test(
            Path(tmp),
            "encode_label_file",
            harness("encode_label_file.c"),
            ["ctools_arena.c"],
        )


if __name__ == "__main__":
    main()
