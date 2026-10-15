# Maintained-code conventions

Use `.clang-format` for new or extracted first-party C modules and Ruff's
`pyproject.toml` settings for Python. Apply formatting only to code being worked
on; do not sweep unrelated numerical routines. Vendor code, the Stata SDK,
generated tables and pinned generator inputs are excluded.

Private cross-file helpers have an owning module prefix and a private header.
Use `*_init` and `*_free` for zero initialization and idempotent cleanup. State
who owns each allocation and whether indices refer to zero-based C columns,
zero-based permutations or one-based Stata columns/observations. Preserve the
existing failure/status domain at adapters, with checked SPI writes.

Numerical refactoring preserves arithmetic order, precision, stopping rules,
missing-value semantics and method-specific inference. Factorization and solve
APIs must state layout, mutation and scratch requirements. Similar-looking
algorithms are not necessarily interchangeable.
