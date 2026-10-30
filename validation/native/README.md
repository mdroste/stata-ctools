# C tests

These tests call C functions directly. Python scripts compile and run the test
programs. From the repository root, run an individual test with:

```sh
python3 -m validation.native.tests.test_sort_keys_native
```

Run all C tests and repository checks with
`python3 -m validation.suite_registry --native`.
See [build requirements](../../docs/development/PLATFORMS.md) for compiler and
OpenMP setup.

| Directory | Contents |
| --- | --- |
| `tests/` | Python scripts that compile tests and check their results |
| `harnesses/` | C test programs for individual components |
| `memory_safety/` | C tests for memory handling and allocation failures |
| `support/` | Compiler helpers and a test substitute for the Stata Plugin Interface (SPI) |

Some tests include the production `.c` files directly so they can call internal
functions and simulate allocation failures or Stata responses.
`support/compiler.py` includes any other source files those functions need in
the same compilation unit. This lets the tests simulate allocation failures in
those files too. The plugin build compiles the files separately.

For new independent tests, prefer calling functions through small headers and
linking the compiled files. Put shared helpers in `support/` so each test script
can run independently. Add each test to `validation/suites.json`; see
[adding a test](../README.md#add-or-move-a-test).

When reorganizing tests, preserve their checks, numerical tolerances, the range
of allocation failures they simulate, the implementations they exercise, and
runs with and without OpenMP.
