# Agent Instructions

These rules apply to every coding agent working in this repository. Read `CLAUDE.md` for the full project conventions (building, validation, project structure); they apply to all agents, not only Claude.

## Reports

**Write every agent-generated report to `reports/`, never to the repository root or `docs/`.** Audits, bug hunts, performance write-ups, investigation notes and similar reports go in the top-level `reports/` folder, which is gitignored and never committed. Use `reports/audits/`, `reports/performance/` or a new descriptive subfolder, and put the date in the file name (for example `reports/performance/PERFORMANCE_CSORT_20261028.md`). Logs, generated data and other scratch output go in the untracked `temp/`, never the repository root.

**Tracked files describe the current code, not its history.** `docs/` holds only user- and developer-facing documentation of how things work now.

- Never add dated files or directories, run logs, driver do-files, patches or build identities under `docs/`. `python3 validation/check_docs_tree.py` enforces this in CI.
- `docs/benchmarks/` holds only evidence that a script regenerates in place at a fixed, undated path: `readme_profile/` (from `validation/profile_readme.do`) and `hdfe_datasets_*.csv` (from `validation/run_hdfe_datasets.py --write-doc`). A new entry needs such a script and a tracked doc that cites it.
- Experiments, A/B runs, profiles, tile sweeps and investigations go in `reports/performance/<topic>_<YYYYMMDD>/`. Do not cite their files from tracked docs or code comments, because they cannot be found from a clean clone.
- When a measurement justifies a threshold or design choice, record the conclusion where the choice lives (a sentence in `DEVELOPERS.md` or a code comment), state the numbers that matter, and name the `validation/` script that reproduces them. Put the full story in the commit message.
- Do not write dated narrative such as "the October 28 investigation" or "the earlier benchmark is retained as historical evidence" in tracked docs. When results are superseded, overwrite them; Git keeps the old ones.

## Versioning

**`release.json` is the only place to set the release version, date and minimum Stata version.** Never hand-edit a version string anywhere else; every copy must state the same version.

1. Edit `release.json`.
2. Run `python3 scripts/sync_release.py`. It stamps the README badge, line 1 (`*! version ...`) of every `build/*.ado`, the `{* *! version ...}` line of every `build/*.sthlp`, `build/ctools.pkg`, `build/stata.toc` and `src/ctools_version.h`, and copies `LICENSE` and `THIRD_PARTY_NOTICES` into `build/`.
3. Run it again with `--check` and commit every regenerated file together. CI fails on stale copies.

New `.ado` and `.sthlp` files get the stamp on the next sync. The plugin's embedded version changes only when it is rebuilt (CI rebuilds on pushes to `main`).
