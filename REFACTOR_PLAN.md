# ctools Infrastructure Refactoring Plan

## Overview
Standardize all ctools commands to use common infrastructure for better performance and maintainability.

---

## Phase 1: Replace `qsort()` with `ctools_sort_*` (6 commands)

### 1.1 cwinsor
- [ ] Find qsort usage for by-variable sorting
- [ ] Replace with `ctools_sort_dispatch()` or appropriate sort function
- [ ] Test with validation script

### 1.2 csample
- [ ] Find qsort usage for by-group sorting
- [ ] Replace with `ctools_sort_dispatch()`
- [ ] Test with validation script

### 1.3 cbsample
- [ ] Find qsort usage for cluster grouping
- [ ] Replace with `ctools_sort_dispatch()`
- [ ] Test with validation script

### 1.4 cencode
- [ ] Find qsort usage for sorting unique values
- [ ] Replace with appropriate ctools sort
- [ ] Test with validation script

### 1.5 cbinscatter
- [ ] Find qsort usage in binning logic
- [ ] Replace with `ctools_sort_dispatch()`
- [ ] Test with validation script

### 1.6 cpsmatch
- [ ] Find qsort usage for match sorting
- [ ] Replace with `ctools_sort_dispatch()`
- [ ] Test with validation script

---

## Phase 2: Standardize Data Loading (4 commands)

### 2.1 cbinscatter
- [ ] Replace direct SF_vdata calls with `ctools_data_load_selective()`
- [ ] Update variable access patterns
- [ ] Test with validation script

### 2.2 cpsmatch
- [ ] Replace custom variable loading with `ctools_data_load_selective()`
- [ ] Update variable access patterns
- [ ] Test with validation script

### 2.3 csample
- [ ] Evaluate if `ctools_data_load_selective()` is appropriate
- [ ] Update if beneficial
- [ ] Test with validation script

### 2.4 cbsample
- [ ] Evaluate if `ctools_data_load_selective()` is appropriate
- [ ] Update if beneficial
- [ ] Test with validation script

---

## Phase 3: Use `ctools_persistent_pool` (2 commands)

### 3.1 cimport
- [ ] Replace custom pthread pool with `ctools_get_global_pool()`
- [ ] Update task submission to use `ctools_persistent_pool_submit_batch()`
- [ ] Test with validation script

### 3.2 cpsmatch
- [ ] Utilize `ctools_persistent_pool` for parallel NN search
- [ ] Test with validation script

---

## Phase 4: Centralize Numeric Parsing

### 4.1 Audit current usage
- [ ] Find all numeric parsing functions across commands
- [ ] Identify common patterns

### 4.2 Consolidate in ctools_types
- [ ] Ensure `ctools_parse_double_with_separators()` is in ctools_types.h
- [ ] Update all commands to use centralized version

---

## Testing Strategy
- Run validation scripts after each change
- Run `make` to ensure compilation
- Run full validation suite at end of each phase

---

## Status: Not Started
