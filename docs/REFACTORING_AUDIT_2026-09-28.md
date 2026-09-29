ctools organization and refactoring audit — September 28, 2026

The six P1 items are addressed in [the implementation and validation report](P1_REFACTORING_FIXES_2026-09-29.md). The user subsequently specified `stata` as the current launcher; item 4 reflects that correction.

The highest-return work is to make lifecycle, parsing, ownership, error handling, and validation contracts explicit before moving large bodies of code. The repository already has useful shared infrastructure; the problem is incomplete adoption and shared functionality that still lives inside individual commands.

This review covers the current, extensively modified working tree on `dev`, including untracked source files. It inventories all 20 command directories, shared infrastructure, Stata wrappers, validation scripts, build and packaging code, and documentation. There are 185 C/header/include files outside the vendored directories and Stata SDK, totaling 74,194 lines, plus 31 ado files totaling 12,341 lines. These counts include generated tables and the locally adapted charset recognizer. Vendored implementations were reviewed for integration and organization, not audited line by line. Source inspection and targeted searches establFix alish the organizational findings; this is not a complete behavioral or numerical correctness certification.

Two temporary native probes compiled the actual runtime and parser sources and reproduced the behaviors in findings 1 and 2. No Stata process, full validation suite, performance benchmark, or platform plugin build was run. No implementation or validation files were changed for this review.

Priority means refactoring order: **P1** addresses correctness-sensitive contracts or the ability to validate subsequent changes; **P2** addresses substantial maintenance cost; **P3** addresses hygiene and consistency. Effort is qualitative: **S** is a narrow change, **M** spans several related files, and **L** crosses major module boundaries. Each item gives a proposed verification requirement, not a claim that that verification has already passed.

1. **P1 — Replace implicit command-history cleanup with explicit lifecycle ownership. Effort: M.**

   Evidence: [ctools_runtime.c:226](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_runtime.c:226), [ctools_plugin.c:296](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_plugin.c:296). Cleanup compares the incoming command with `g_previous_command`, but the dispatcher updates history only afterward. A native probe of `cmerge → cimport → cmerge` recorded zero import-cleanup callbacks on the final switch, where one is expected. Command classification, handler dispatch, and cleanup hooks also live in separate lists.

   Solution: first repair the transition logic in a small change. Then introduce a command descriptor containing handler, state owner, and cleanup hook. Give multi-phase operations explicit phase states and preserve only the active owner's valid session. Keep wire aliases such as `cio` explicit. Do not infer continuation solely from an earlier command name.

   Verify: same-command continuation; A→B→A and A→B→C; failed scan/prepare; repeated cleanup; and recovery after interrupted operations. Mock cleanup callbacks can cover the dispatcher without a licensed Stata process.

2. **P1 — Complete and harden the shared argument parser before migrating callers. Effort: M.**

   Evidence: [ctools_parse.c:20](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_parse.c:20), [csample_impl.c:214](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c:214), [crangestat_impl.c:1156](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:1156). Sampling and range statistics retain local parsers; merge/import use additional tokenization schemes. The shared helpers themselves need repair: the native probe returned false for `notverbose verbose`, selected `wrong` from `nolabel=wrong label=right`, and accepted `2147483648` as an `int`, producing `-2147483648` on this machine.

   Solution: provide a bounded token cursor with whole-token matching, checked signed/unsigned conversions, and distinct absent/invalid/present statuses. Add explicit quoting support only where the wire protocol needs it. Move seed-half parsing into the shared layer. Migrate one command family at a time; merely replacing local functions with today's shared functions would retain defects.

   Verify: token prefixes, repeated options, whitespace, quoted values, trailing junk, missing values, integer limits, and oversized input. These are helper-level defects; the probes do not establish that every example is reachable through public ado syntax.

3. **P1 — Establish one authoritative validation-suite registry. Effort: M.**

   Evidence: [validate_all.do:29](/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_all.do:29), [run_stata_audit.py:11](/Users/Mike/Documents/GitHub/stata-ctools/validation/run_stata_audit.py:11), [build.yml:76](/Users/Mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:76). `memory_safety` is in the master Stata list but absent from the Python completion list. The separate I/O parity/options/formats suites, transport-strL checks, command-optimization checks, and big-data checks are outside that master runner. Native CI is more comprehensive: it explicitly runs several special cases and discovers the remaining `test_*_native.py` files, but its policy is encoded in shell lists.

   Solution: define a registry with suite name, runner, fixture preparation, completion marker, dependency requirements, and category. Generate or validate the Stata and Python lists against it. Make correctness, optional large-data validation, and performance runs explicit categories; do not silently include costly benchmarks in every CI run. Require every suite to be registered or deliberately excluded with a reason.

   Verify: adding a suite cannot leave it undiscovered; missing markers fail the gate; fixture requirements are prepared; CI and local runners select the same correctness suites. The existing master still counts `memory_safety` failures, so the list mismatch alone is not proof that its failures pass publication.

4. **P1 — Centralize Stata launching and eliminate machine-specific validation paths. Effort: S–M.**

   Evidence: [run_stata_audit.py:56](/Users/Mike/Documents/GitHub/stata-ctools/validation/run_stata_audit.py:56), [run_io_parity.py:68](/Users/Mike/Documents/GitHub/stata-ctools/validation/run_io_parity.py:68), [run_io_excel_options.py:33](/Users/Mike/Documents/GitHub/stata-ctools/validation/run_io_excel_options.py:33), [validate_io_excel_options.do:3](/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_io_excel_options.do:3). Five Python runners independently construct `stata -q -b do ...`; the user has confirmed that `stata` is the current launcher. Several do-files embed `/Users/Mike/Documents/GitHub/stata-ctools/temp/...`.

   Solution: one runner utility should build `/bin/zsh -lic 'stata ...'` on this machine, quote the absolute driver path, capture nonce-tagged completion, and guarantee a driver ending with `exit, clear`. Pass fixture/output directories as arguments or configured globals. A different machine must use an explicitly confirmed launcher configuration rather than assuming this wrapper exists there.

   Verify: launcher construction with spaces in paths, nonzero Stata return codes, incomplete logs, and temporary-directory cleanup. Use the confirmed `stata` alias through interactive login zsh.

5. **P1 — Separate internal status codes from Stata return codes. Effort: M.**

   Evidence: [ctools_types.h:26](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_types.h:26), [cipolate_impl.c:37](/Users/Mike/Documents/GitHub/stata-ctools/src/cipolate/cipolate_impl.c:37), [cipolate_impl.c:106](/Users/Mike/Documents/GitHub/stata-ctools/src/cipolate/cipolate_impl.c:106), [crangejoin_impl.c:78](/Users/Mike/Documents/GitHub/stata-ctools/src/crangejoin/crangejoin_impl.c:78). Internal statuses are 0–5, while command failures use Stata codes such as 109, 198, 459, and 920. Several callers convert any load/store failure to 920; others distinguish allocation and write failures. Low-level numerical routines also use `-1` for different causes.

   Solution: define one explicit boundary conversion helper, and retain the original SPI failure where useful. Numerical routines should distinguish allocation failure, invalid input, singularity, and nonconvergence. Convert statuses once at command boundaries, with command-specific context in the message.

   Verify: injected allocation, read, write, cancellation, and unsupported-type failures produce the intended status and leave no successful output publication. This improves diagnosis as well as organization.

6. **P1 — Encode owned versus borrowed data instead of relying on cleanup conventions. Effort: L, with small initial extractions.**

   Evidence: [ctools_types.h:59](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_types.h:59), [ctools_types.c:68](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_types.c:68), [cmerge_impl.c:1392](/Users/Mike/Documents/GitHub/stata-ctools/src/cmerge/cmerge_impl.c:1392), [cmerge_impl.c:1493](/Users/Mike/Documents/GitHub/stata-ctools/src/cmerge/cmerge_impl.c:1493). The generic destructor treats a string column without an arena as individually owned strings. Merge deliberately uses `_arena == NULL` for borrowed strings and therefore repeats custom cleanup instead. Allocation provenance also matters on Windows because aligned buffers need a different deallocator.

   Solution: introduce explicit column ownership/storage modes or separate owning buffers from read-only views. First extract merge's current output cleanup into one helper, preserving its borrowing optimization. Then use typed workspaces with idempotent destructors across the largest commands. Preserve owner lifetimes until all parallel stores finish.

   Verify: allocation failure at each stage, partial output, repeated destruction, and borrowed strings with both normal and error exits. Cover Windows allocator pairing rather than assuming POSIX `free` behavior.

7. **P2 — Split estimator orchestration into explicit phases. Effort: L.**

   Evidence: [cpplmhdfe_irls.c:453](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:453), [creghdfe_regress.c:115](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_regress.c:115), [civreghdfe_impl.c:297](/Users/Mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.c:297). The main PPML, OLS-HDFE, and IV functions span approximately 1,909, 1,529, and 1,482 physical lines respectively. Each combines input preparation, sample changes, numerical estimation, covariance, output, and cleanup. IV's `ivest_compute_2sls` is another 1,075-line function.

   Solution: introduce command-owned options, workspace, sample map, and result structures. Extract load/validate, sample construction, FE preparation, estimation, inference, and posting phases. Keep estimator-specific statistical choices explicit. Extract cohesive phases with unchanged arithmetic before attempting algorithmic unification.

   Verify: coefficients, covariance, retained sample, residuals, FE output, convergence status, and posted metadata against existing references. Preserve reduction order and tolerances during the structural change.

8. **P2 — Finish extracting HDFE infrastructure from `creghdfe`. Effort: L.**

   Evidence: [civreghdfe_impl.h:14](/Users/Mike/Documents/GitHub/stata-ctools/src/civreghdfe/civreghdfe_impl.h:14), [cpplmhdfe_irls.c:29](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:29), [creghdfe_types.h:37](/Users/Mike/Documents/GitHub/stata-ctools/src/creghdfe/creghdfe_types.h:37). Shared FE types already live in `ctools_hdfe_utils`, but IV imports five `creghdfe` headers and PPML imports its solver/utilities. `g_state` remains exported through a command type header. Binscatter imports remapping utilities from the same command directory.

   Solution: put generic FE remapping, factor lifecycle, projection interfaces, and solver contracts in a neutral HDFE module. Pass state explicitly. Keep creghdfe posting and OLS-specific adjustments in its adapter. Preserve PPML's reference-sensitive projection controls and binscatter's simpler implementation until equivalence is established.

   Verify: all consumers build against the neutral interface; no generic HDFE header depends on a command entry-point header; projection failure and convergence propagate unchanged.

9. **P2 — Make the ado-to-C protocol an explicit boundary. Effort: M–L.**

   Evidence: [cpplmhdfe_irls.c:45](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:45), [cpplmhdfe_irls.c:2266](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_irls.c:2266), [cio.c:72](/Users/Mike/Documents/GitHub/stata-ctools/src/io/cio.c:72), [_ctools_load.ado:57](/Users/Mike/Documents/GitHub/stata-ctools/build/_ctools_load.ado:57). Options/results move through positional strings, local macros, global macros, and individually named scalars. PPML reads global options into a static struct; I/O has special macro-name translation to survive ado helper scopes. A plugin version query exists, but it is not a command schema.

   Solution: give each command one protocol reader and one result writer, backed by typed C options/results and centralized field definitions. Add a protocol version/capability mechanism only as needed for incompatible changes. Keep Stata scope and tempvar ownership explicit. Initially preserve the existing wire representation to avoid a simultaneous protocol migration.

   Verify: missing or malformed fields, stale globals, wrapper/plugin mismatch, failed result writes, and helper-program scope changes. Numerical kernels should not reach back into Stata for configuration.

10. **P2 — Turn the format I/O include chain into real modules. Effort: L.**

   Evidence: [cio.c:20](/Users/Mike/Documents/GitHub/stata-ctools/src/io/cio.c:20), [cio.c:686](/Users/Mike/Documents/GitHub/stata-ctools/src/io/cio.c:686), [cimport_impl.c:2100](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_impl.c:2100). `cio.c` defines shared dataset types and a static cache, then includes six implementation fragments for formats. Together these files comprise 3,064 lines in one translation unit. Generic I/O cleanup is invoked through the CSV import cleanup function, and I/O imports an export command's parsing header.

   Solution: separate dataset/metadata types, the Stata adapter, format readers/writers, and session ownership. Give each format a narrow scan/load/export interface receiving an explicit context. Convert handwritten implementation `.inc` files to `.c` modules once their dependencies are explicit; generated data tables can remain textual includes.

   Verify: format-specific metadata, dates, labels, tagged missing values, long strings, malformed inputs, and scan/load/clear transitions. Use the existing format fixtures and differential suites.

11. **P2 — Move shared format utilities to their actual owners. Effort: M.**

   Evidence: [cexport_parse.h:16](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_parse.h:16), [cexport_xlsx.c:11](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_xlsx.c:11), [cimport_xlsx_xml.c:528](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_xlsx_xml.c:528), [cexport_xlsx.c:352](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_xlsx.c:352). Output-file prepare/commit/cleanup lives in an export argument-parsing header and is used by other formats. XLSX export obtains ZIP/XML utilities from the import directory. Cell references have multiple parsers, including bounded and NUL-terminated forms.

   Solution: extract neutral file-output transactions, XML/ZIP support, and Excel address/format primitives. Implement a bounded cell-reference parser with small adapters where accepted syntax differs. Separate XLSX writing, workbook editing, and formula adjustment behind explicit interfaces. Keep CSV and Excel number-formatting policies distinct.

   Verify: atomic replacement behavior and cleanup; Unicode and XML escaping; address boundaries; existing workbook/formula preservation; byte/structural comparisons of exported files.

12. **P2 — Give long-text transport a neutral, policy-based home. Effort: M.**

   Evidence: [cexport_xlsx.c:1547](/Users/Mike/Documents/GitHub/stata-ctools/src/cexport/cexport_xlsx.c:1547), [ctools_data_io.c:107](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_data_io.c:107), [test_xlsx_strl_native.py:33](/Users/Mike/Documents/GitHub/stata-ctools/validation/test_xlsx_strl_native.py:33). The generally named `cexport_load_text_data` sits inside the XLSX writer and duplicates selection, allocation, and numeric/string loading. Its long-string policy differs intentionally from the shared loader's bounded str2045 handling.

   Solution: extract serial, length-aware text loading into transport with explicit byte/character limits and a binary-string policy. Let the Excel adapter enforce its UTF-16 cell limit. Reuse observation selection and ownership helpers. Keep strL SPI calls on the calling thread as required by the existing implementation evidence.

   Verify: long strings, supplementary Unicode, binary strL rejection, filtered observations, cancellation, and read/allocation failures. Do not replace this with the bounded shared loader unchanged.

13. **P2 — Separate transport planning from transport kernels. Effort: M–L.**

   Evidence: [ctools_data_io.c:77](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_data_io.c:77), [ctools_data_io.c:1267](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_data_io.c:1267), [ctools_data_io.c:1987](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_data_io.c:1987). The 2,335-line file contains validation, observation selection, multiple string storage paths, row/column/tile scheduling, worker kernels, permutation stores, and public dispatch.

   Solution: introduce a validated transfer plan containing columns, row mapping, widths, operation, and selected schedule. Separate selection, allocation, load/store kernels, and scheduling while retaining the public facade. Let tuned paths share validation and ownership without forcing their inner loops through costly generic callbacks.

   Verify: transport and scheduling native suites, checked-store invariants, thread-count limits, strL behavior, and representative narrow/wide numeric/string benchmarks. Preserve the documented checked write callbacks.

14. **P2 — Separate statistics and matching kernels from command orchestration. Effort: M–L.**

   Evidence: [crangestat_impl.c:1187](/Users/Mike/Documents/GitHub/stata-ctools/src/crangestat/crangestat_impl.c:1187), [cpsmatch_impl.c:280](/Users/Mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_impl.c:280), [cpsmatch_impl.c:701](/Users/Mike/Documents/GitHub/stata-ctools/src/cpsmatch/cpsmatch_impl.c:701). Range statistics occupy a 2,399-line implementation with a 1,213-line main function. Matching combines option parsing, candidate search, kernels, workspaces, treatment-effect accumulation, and SPI output in one file.

   Solution: split range parsing/planning, sorted window traversal, statistic kernels, and output. Split matching configuration, neighbor/radius/kernel search, effect accumulation, and the Stata adapter. Use command-specific workspace structs; reuse shared ordering where semantics match.

   Verify: empty windows, excludeself, ties, missing values, by-groups, collinear window regressions, calipers, replacement rules, and documented standard-error method differences. Preserve the statistical definitions.

15. **P2 — Deduplicate sampling RNG and seed handling. Effort: S.**

   Evidence: [csample_impl.c:40](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c:40), [cbsample_impl.c:42](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c:42). Both commands carry xoshiro256**, SplitMix64, bounded rejection sampling, and the same two-half seed protocol.

   Solution: extract a small private RNG utility and shared seed decoder. Keep sampling with replacement and Fisher–Yates sampling as distinct algorithms. Preserve seed expansion, draw order, and thread/group assignment.

   Verify: exact seeded outputs before/after, including grouped and clustered cases, zero/one-sized populations, and rejection-sampling boundaries. Distributional checks alone are insufficient for a behavior-preserving extraction.

16. **P2 — Consolidate grouping and clarify index conventions. Effort: M.**

   Evidence: [csample_impl.c:98](/Users/Mike/Documents/GitHub/stata-ctools/src/csample/csample_impl.c:98), [cbsample_impl.c:98](/Users/Mike/Documents/GitHub/stata-ctools/src/cbsample/cbsample_impl.c:98), [ctools_order.h:4](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_order.h:4), [ctools_types.h:209](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_types.h:209). Sampling duplicates group-boundary detection despite `ctools_group_starts`. The newer ordering interface uses zero-based columns; sort dispatch takes one-based variables. Observation maps and permutations have further distinct meanings.

   Solution: provide one grouping contract with adapters for physically sorted data and index-sorted data. Use explicit names or lightweight wrapper types for plugin variable indices, C column indices, Stata observations, and permutations. Convert at documented boundaries rather than throughout kernels.

   Verify: strings, empty strings, signed zero, distinct extended missings, stable ties, empty datasets, and filtered/permuted observations. Preserve a parallel group-detection path if measurements justify it.

17. **P2 — Share sort support code without collapsing specialized algorithms. Effort: M.**

   Evidence: [ctools_sort_merge.c:85](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_merge.c:85), [ctools_sort_sample.c:226](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_sample.c:226), [ctools_sort_radix_lsd.c:86](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_sort_radix_lsd.c:86). Search found matching radix histogram/scatter sequences in merge/sample support and similar scratch setup in LSD/MSD. The seven algorithms also carry local key, allocation, and prefetch scaffolding.

   Solution: extract proven-common key transforms, stable block-radix primitives, checked scratch allocation, and partition helpers into private sorting support. Keep algorithm selection and specialized hot loops separate. Organize sorting under a dedicated directory only after interfaces stabilize.

   Verify: algorithm-specific parity, stability, missing-key semantics, constrained OpenMP teams, allocation-failure fallbacks, and existing performance workloads. Similar code alone does not prove every implementation should be merged.

18. **P2 — Give linear algebra explicit layout, mutation, and precision contracts. Effort: M–L.**

   Evidence: [ctools_ols.c:25](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_ols.c:25), [cqreg_linalg.c:140](/Users/Mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_linalg.c:140), [cqreg_blas.c:27](/Users/Mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_blas.c:27), [ctools_matrix.c:29](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_matrix.c:29). Cholesky and solve logic overlap, while BLAS/fallback and matrix routines have separate interfaces. There are real behavioral distinctions: shared Cholesky zeros the upper triangle; cqreg's version leaves it alone, and solve APIs differ in whether they refactor or consume a factor.

   Solution: document storage layout, leading dimensions, overwrite behavior, factorization state, scratch requirements, precision, and failure codes. Extract exact common triangular/factorization primitives first, using thin adapters. Keep compensated accumulations and reference-sensitive solver choices explicit.

   Verify: singular/near-singular inputs, ill-conditioning, weighted cases, backend equivalence, and existing regression tolerances. Do not substitute a faster BLAS path as part of a structural refactor without separate numerical evidence.

19. **P2 — Decompose the largest ado wrappers around Stata scope boundaries. Effort: M–L.**

   Evidence: [civreghdfe.ado:3](/Users/Mike/Documents/GitHub/stata-ctools/build/civreghdfe.ado:3), [creghdfe.ado:141](/Users/Mike/Documents/GitHub/stata-ctools/build/creghdfe.ado:141), [cpplmhdfe.ado:175](/Users/Mike/Documents/GitHub/stata-ctools/build/cpplmhdfe.ado:175), [cmerge.ado:59](/Users/Mike/Documents/GitHub/stata-ctools/build/cmerge.ado:59). The IV wrapper is 2,111 lines, merge 1,414, and PPML 1,094. Syntax handling, variable expansion, protocol fields, result reconstruction, and display are interleaved. Shared loader, weight, and new-variable helpers already provide a useful starting point.

   Solution: extract pure validation, factor/name bookkeeping, and result-display helpers first. Keep tempvar allocation and plugin registration in their required caller scopes. Separate estimator-specific sample/DOF rules from reusable syntax mechanics; avoid a universal regression wrapper.

   Verify: factor and time-series expansion, coefficient names, `e(sample)`, prediction, caller data preservation, weight expressions, and plugin registration after `clear all`. Repeated registration lines are not automatically removable: Stata scopes them to the ado caller.

20. **P2 — Separate editable package source from compiled/staged output. Effort: M.**

   Evidence: [Makefile:35](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:35), [scripts/stage_package.py:65](/Users/Mike/Documents/GitHub/stata-ctools/scripts/stage_package.py:65), [build/ctools.pkg](/Users/Mike/Documents/GitHub/stata-ctools/build/ctools.pkg). Handwritten ado/help files, generated identity/notices, and platform binaries share `build/`. This makes the directory simultaneously authoritative source and output. Binary publication there is intentional for installation, so simply deleting or ignoring it is inappropriate.

   Solution: introduce canonical ado/help/package-source directories and a deterministic staging step. Keep the published `build/` layout if existing Stata installation URLs require it, but document generated files and prevent hand editing of staged copies. Use an unambiguous development object/output directory.

   Verify: an installed package contains every helper and help file; package manifests remain complete; staging from a clean checkout reproduces the intended install tree; installation URLs continue to work.

21. **P2 — Refactor platform builds into reusable, incremental object rules. Effort: M.**

   Evidence: [Makefile:68](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:68), [Makefile:406](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:406), [Makefile:443](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:443), [Makefile:509](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:509). Phony platform targets recompile the full source set; four near-identical shell recipes manage private objects, compile libdeflate, link all other sources, and publish the result. All source directories become global include paths.

   Solution: use per-platform/toolchain object directories, compiler-generated dependency files, and a shared compile/link recipe parameterized by platform flags. Track flag and source-inventory changes so old objects cannot survive a configuration change or source deletion. Narrow include paths and require explicit module-qualified includes. Preserve atomic final publication and fail-fast behavior.

   Verify: no-op builds do no compilation; changing one header rebuilds the correct consumers, including `.inc` dependencies; source deletion and flag changes invalidate the right objects; an intentionally failed compile/link leaves the previous plugin intact. Retain the existing build-contract tests.

22. **P2 — Consolidate vendored dependencies and isolate their build flags. Effort: M.**

   Evidence: [cimport/miniz/miniz.h](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/miniz/miniz.h), [io/vendor/README.md](/Users/Mike/Documents/GitHub/stata-ctools/src/io/vendor/README.md), [Makefile:69](/Users/Mike/Documents/GitHub/stata-ctools/Makefile:69). Compression libraries are under the import command; ReadStat/libxls/shims are under generic I/O. Only libdeflate gets a separate warning policy. The vendor README explicitly records header renames needed to avoid shadowing system headers under recursive include discovery.

   Solution: put third-party packages under a single vendor root with scoped include directories and per-package source/flag declarations. Maintain upstream version, license, local patches, and regeneration details together. Keep the zlib-to-miniz shim clearly named and scoped; do not expose its substitute header globally by accident.

   Verify: all platform dependency checks and notices pass, no extra dynamic dependency appears, and first-party warnings remain visible. Directory movement should not become an upstream-library upgrade.

23. **P2 — Use one package staging implementation locally and in CI. Effort: S–M.**

   Evidence: [stage_package.py:30](/Users/Mike/Documents/GitHub/stata-ctools/scripts/stage_package.py:30), [build.yml:263](/Users/Mike/Documents/GitHub/stata-ctools/.github/workflows/build.yml:263), [check_package.py:9](/Users/Mike/Documents/GitHub/stata-ctools/validation/check_package.py:9). Local staging validates entries, declares platform scope, stages atomically, and writes hashes. CI independently copies globs and parses the package manifest with shell commands, then writes separate provenance. Packaging code imports its checker by adding `validation/` to `sys.path`.

   Solution: use the same staging library/CLI for host and all-platform packages, parameterized by binary input locations and revision. Move package invariants into a small reusable tooling module consumed by staging and tests. Extend the existing `release.json` source of truth rather than replacing it with another version list.

   Verify: equivalent manifests/content for the same scope, deliberate failure on a missing helper or binary, correct provenance/hashes, and unchanged publication gates. Preserve the licensed-Stata release requirement.

24. **P2 — Extract native test infrastructure from individual test scripts. Effort: M.**

   Evidence: [test_transport_native.py:11](/Users/Mike/Documents/GitHub/stata-ctools/validation/test_transport_native.py:11), [test_io_parity_native.py:4](/Users/Mike/Documents/GitHub/stata-ctools/validation/test_io_parity_native.py:4), [test_split_tokens_native.py:46](/Users/Mike/Documents/GitHub/stata-ctools/validation/test_split_tokens_native.py:46), [test_newcommands_native.py:15](/Users/Mike/Documents/GitHub/stata-ctools/validation/test_newcommands_native.py:15). Tests import mocks and compilers from other tests, embed large C programs in Python strings, and sometimes rename another harness's `main` using string replacement. Many include production `.c` files to access static internals. The newer `validation/native/memory_safety/` directory demonstrates a cleaner separation.

   Solution: create native support for SPI mocks, fault injection, compiler/backend flags, and sanitizer profiles. Move harnesses to ordinary C files, leaving Python to orchestrate. Prefer linking production objects through stable internal test interfaces; retain direct inclusion only for deliberate white-box tests.

   Verify: preserve every existing assertion and allocation-failure iteration; compare discovered test counts; cover both serial and OpenMP configurations. Refactoring tests must not relax tolerances or reinterpret failures.

25. **P2 — Split oversized core headers and move mutable initialization into implementation files. Effort: M.**

   Evidence: [ctools_types.h:107](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_types.h:107), [ctools_types.h:218](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_types.h:218), [ctools_config.h:214](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_config.h:214), [ctools_config.h:350](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_config.h:350). The types header exposes transport, sorting, conversion, and data structures. The config header contains allocation helpers, OS cache detection, threading APIs, and mutable static prefetch state, giving including translation units their own copies. Its initialization comments permit simultaneous writes to a non-atomic structure.

   Solution: separate data types, transport API, sort API, allocation helpers, and platform/thread configuration. Put runtime cache discovery behind one implementation-owned initialization path, ideally initialized before workers or protected by an appropriate once primitive. Keep only pure, justified hot-path inline helpers in headers. Provide an umbrella compatibility header during migration.

   Verify: public headers compile independently, clients include their actual dependencies, initialization runs safely under parallel access, and optimized kernels retain performance. The initialization concern is based on source inspection, not a demonstrated production race.

26. **P3 — Give generated tables an explicit generation workflow. Effort: S–M.**

   Evidence: [generate_cimport_locales.py:1](/Users/Mike/Documents/GitHub/stata-ctools/validation/generate_cimport_locales.py:1), [cimport_charset.inc:11](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_charset.inc:11), [cimport_locale_profiles.inc](/Users/Mike/Documents/GitHub/stata-ctools/src/cimport/cimport_locale_profiles.inc). Generated locale/charset tables and handwritten implementation fragments all use `.inc`, while their generators live under validation.

   Solution: move generators to a generation-tool directory; place generated outputs in an identifiable location or naming scheme; record exact source/runtime versions and checksums. Add a check mode that compares regenerated output without rewriting it. Keep generated data distinct from handwritten `.inc` implementations proposed for modularization.

   Verify: deterministic generation with pinned inputs and an explicit check for drift. Ordinary builds should continue to use checked-in tables without Java or external downloads.

27. **P3 — Standardize diagnostics and timing without adding work to hot loops. Effort: S–M.**

   Evidence: [ctools_runtime.c:24](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_runtime.c:24), [ctools_runtime.h:119](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_runtime.h:119), [cqreg_fn.c:1091](/Users/Mike/Documents/GitHub/stata-ctools/src/cqreg/cqreg_fn.c:1091), [civreghdfe.ado:14](/Users/Mike/Documents/GitHub/stata-ctools/build/civreghdfe.ado:14). Message functions repeat formatting logic. Multiple timing macro pairs, command-local timing fields, debug instrumentation, and hard-coded ado timer IDs coexist; the IDs can overwrite users' timer state.

   Solution: share a bounded `va_list` message formatter and one phase-timing convention with explicit units. Use a per-command timing result and serialize it centrally. Avoid fixed public timer IDs where possible, or define/document a deliberate timer policy. Keep expensive diagnostics opt-in.

   Verify: existing verbose output and benchmark consumers, phase units, truncated messages, and caller timer state. Preserve useful phase granularity rather than replacing it with one total time.

28. **P3 — Separate scratch outputs, fixtures, and historical audit evidence. Effort: S.**

   Evidence: [.gitignore](/Users/Mike/Documents/GitHub/stata-ctools/.gitignore), [validation/temp_cppl_pathology_debug.do](/Users/Mike/Documents/GitHub/stata-ctools/validation/temp_cppl_pathology_debug.do), [validation/fixtures/README.md](/Users/Mike/Documents/GitHub/stata-ctools/validation/fixtures/README.md). This working tree contains many ignored logs, an untracked `temp/` tree, temporary-named do-files in validation, root-level dated audits, and intentionally retained benchmark evidence. These categories are difficult to distinguish during discovery.

   Solution: route ephemeral runs to a documented scratch directory; keep fixtures with provenance under fixtures; promote useful probes into named regressions; archive historical reports under an audit directory with an index. Review untracked work before moving or deleting anything. Narrow broad ignore patterns where they obscure legitimate fixture data.

   Verify: a clean test run leaves only declared outputs, fixture discovery still works, and retained benchmark/audit provenance remains accessible. Working-tree clutter is not evidence that all those files are tracked.

29. **P3 — Update architecture documentation to reflect actual ownership and contracts. Effort: S.**

   Evidence: [DEVELOPERS.md](/Users/Mike/Documents/GitHub/stata-ctools/DEVELOPERS.md), [CLAUDE.md](/Users/Mike/Documents/GitHub/stata-ctools/CLAUDE.md), [ctools_plugin.c:8](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_plugin.c:8), [ctools_runtime.h:167](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_runtime.h:167). Guidance still refers to the `stata` alias, parts of the dispatcher header describe retired cdecode behavior, and lifecycle comments refer to cleanup at plugin exit while dispatch cleans stale command state on entry and destroys the worker pool on exit. Thread-default descriptions also coexist with code that captures OS CPU capacity.

   Solution: document the real module boundaries, phase lifetimes, index bases, memory ownership, thread policy, launcher, and validation entry points. Derive command inventory where practical and make the add-command checklist cover descriptor registration, package entries, and suite registration. Keep SDK files untouched.

   Verify: walkthrough for adding one small command and reproducing a validation run matches current code. Retain the deliberate legacy cdecode rejection shim while correcting its description.

30. **P3 — Adopt incremental formatting and naming conventions for maintained code. Effort: S to establish, ongoing thereafter.**

   Evidence: [cpplmhdfe_solvers.h:23](/Users/Mike/Documents/GitHub/stata-ctools/src/cpplmhdfe/cpplmhdfe_solvers.h:23), [cio.c:44](/Users/Mike/Documents/GitHub/stata-ctools/src/io/cio.c:44), [ctools_parse.c:20](/Users/Mike/Documents/GitHub/stata-ctools/src/ctools_parse.c:20). Dense multi-statement lines, varying indentation, generic private names, mixed include styles, and implementation-heavy `.h` files coexist. No repository formatting configuration was found in the inspected tree.

   Solution: add a small formatting/style configuration for first-party C and Python; format extracted or modified modules in separate mechanical changes. Use consistent names for context init/free, checked statuses, and private helpers. Exclude vendored sources, generated tables, and the Stata SDK. Move PPML's implementation-heavy solver headers to private source modules as part of its phase extraction rather than renaming them in isolation.

   Verify: mechanical diffs preserve code behavior, generation checks remain clean, and substantive numerical changes are reviewed separately. Avoid a repository-wide formatting sweep across the current ongoing work.

Recommended implementation order: first repair and characterize lifecycle/parser behavior and make the validation launcher/registry reliable (1–4). Then standardize statuses and ownership (5–6) and extract native test support (24). Next make small shared-library extractions—RNG, grouping adapters, format utilities, and header boundaries—before splitting estimator/I/O orchestration. Build, package, and vendor reorganization should be separate changes with install/build-contract checks. Apply documentation, generation, and formatting improvements alongside the modules they describe.

Existing work worth preserving: shared matrix/VCE and HDFE types, the shared loader/weight/new-variable ado helpers, atomic output staging, atomic plugin publication, package completeness checks, pinned fixture provenance, native fault-injection tests, constrained-team tests, explicit legacy cdecode rejection, and performance reports. The goal is to complete these abstractions and clarify their boundaries. It is not to replace specialized algorithms, alter numerical tolerances, change statistical definitions, or remove compatibility behavior merely to reduce line counts.
