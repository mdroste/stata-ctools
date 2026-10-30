# Native tests

These tests exercise C code directly, using Python drivers to compile and run
small test programs. Run a driver from the repository root:

```sh
python3 -m validation.native.tests.test_sort_keys_native
```

To run all registered native tests and repository checks, use
`python3 -m validation.suite_registry --native`.
See [platform dependencies](../../docs/development/PLATFORMS.md) for compiler
and OpenMP setup.

| Directory | Contents |
| --- | --- |
| `tests/` | Python drivers that compile tests and check their results |
| `harnesses/` | C test programs for individual components |
| `memory_safety/` | C tests for memory handling and allocation failures |
| `support/` | Compiler helpers and a mock Stata Plugin Interface (SPI) |

Some tests include production `.c` files directly to inspect internal functions
and substitute allocation or SPI calls. `support/compiler.py` adds the companion
source files needed by those implementations to the same compilation unit. This
lets the tests inject allocation failures throughout the code they exercise.
Plugin builds compile those sources separately.

For new independent tests, prefer small headers and separately compiled objects.
Put shared helpers in `support/` so each driver can run independently. Register
every correctness driver in `validation/suites.json`; see
[adding a test](../README.md#add-or-move-a-test).

When reorganizing tests, preserve assertions, numerical tolerances, allocation
failure loop bounds, backend selections and serial/OpenMP variants.
