"""Repository validation package; run modules from the repository root."""

from pathlib import Path

if not __debug__:
    raise RuntimeError("Validation requires Python assertions; remove -O/PYTHONOPTIMIZE")

ROOT = Path(__file__).resolve().parents[1]
