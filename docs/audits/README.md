# Retained audit evidence

These are historical snapshots, not descriptions of the current implementation.
Their assertions, paths and counts refer to the revision inspected at the time.

| Evidence | Purpose |
| --- | --- |
| [C pre-release audit, October 14](C_RELEASE_AUDIT_2026-10-14.md) | Memory-safety, crash and leak findings in the C sources, with fix status |
| [Proposed patches for that audit](C_RELEASE_AUDIT_2026-10-14_patches.md) | Per-finding diffs as proposed; the applied fixes differ where the report says so |
| [Compatibility fixes, October 2](compatibility_fixes_20261002.json) | Before/after test counts for the cencode, cwinsor and crangestat repairs |
| [Option coverage, October 2](option_coverage_20261002.json) | Scope and results of the parsed-option coverage run |

The earlier September audit reports and evidence archives, and the September 28
to October 2 refactoring notes, were removed from the repository on 15 October
2026. The September reports remain in git history at the repository root.

Current validation is defined by `validation/suites.json`. Active performance
reports remain under `docs/`, with data under `docs/benchmarks/`. Historical
transport experiments already tracked under `temp/transport_20260929/` remain in
place so their provenance and internal links are preserved.
