# Native test organization

Python drivers orchestrate compilation and execution. C fixtures live in
`harnesses/` or `memory_safety/`; shared SPI and compiler support live in
`support/`. Drivers never import executable tests from other tests.

The existing suites intentionally include production `.c` files for white-box
assertions and allocator/SPI interposition. `support/compiler.py` declares the
companion modules extracted from those implementations and includes them in the
same test translation unit, preserving every injected allocation. Ordinary
plugin builds compile separate objects. This explicit exception avoids exposing
private numerical functions as a public testing API.

Assertions, numerical tolerances, fault-loop bounds, backend selections, serial
and OpenMP variants are part of each suite's contract. Refactoring infrastructure
must preserve them. New independent tests should prefer linking narrow headers
and objects. Register every correctness driver in `validation/suites.json`.
