# Proposed patches — C pre-release audit, 2026-10-14

Companion to [C_RELEASE_AUDIT_2026-10-14.md](C_RELEASE_AUDIT_2026-10-14.md). These are the verifiers' and auditors' proposed changes, **not applied and not compiled** unless the main report says otherwise. They were written against the working tree of 14 October 2026 while other sessions were editing it, so hunks may need re-basing. Review each one before applying.

## F004 [P0, confirmed] Pool worker threads get the 512 KiB macOS default stack, but the cimport chunk parser needs 524,416 bytes of stack frame: stack overflow on Intel macOS

Location: The large frame comes from src/cimport/cimport_impl.c, cimport_parse_chunk_parallel (about line 732): `CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];`. Two thread-creation sites give that function too small a stack. One is the persistent pool in src/ctools_threads.c, ctools_persistent_pool_init (about line 233): `pthread_create(&pool->workers[i], NULL, ...)`. The other is cimport's raw-pthread fallback in cimport_parse_csv (cimport_impl.c about line 1255). The dispatch is in cimport_parse_csv: num_chunks = cpu_count at about lines 1093-1110, and submit_batch at about line 1243.

Give pool workers an explicit stack size in ctools_persistent_pool_init:
```c
#define CTOOLS_POOL_STACK_SIZE ((size_t)8 * 1024 * 1024)
pthread_attr_t attr, *attrp = NULL;
if (pthread_attr_init(&attr) == 0) {
    size_t cur = 0;
    if (pthread_attr_getstacksize(&attr, &cur) != 0 || cur < CTOOLS_POOL_STACK_SIZE)
        (void)pthread_attr_setstacksize(&attr, CTOOLS_POOL_STACK_SIZE);
    attrp = &attr;
}
... pthread_create(&pool->workers[i], attrp, persistent_worker_thread, pool) ...
if (attrp) pthread_attr_destroy(attrp);   /* also on the failure path */
```
Also move field_buf off the stack in cimport_parse_chunk_parallel: `CImportFieldRef *field_buf = malloc(CTOOLS_MAX_COLUMNS * sizeof *field_buf);`, set ctx->error_code on failure and free it at exit. Use the same attr for the raw-pthread fallbacks at cimport_impl.c:1255/1563. Add a large-CSV (more than 2 MiB) cimport regression to validation; the current suite's largest file has 10k rows.

## F009 [P0, confirmed] Streaming csort writes numeric columns from OpenMP threads through the unchecked SF_vstore store (SD_FASTMODE), which the project bans for concurrent writes

**Verified minimal fix:** add `#define SD_SAFEMODE` on the line before `#include "stplugin.h"` in `src/csort/csort_stream.c`. A scratch plugin built with only that change passed the post-estimation reproduction (`cf` rc 0; the current build gives 71 mismatches in `yhat`). The verifier's fuller diff below also propagates store errors; it has not been compiled.

Location: src/csort/csort_stream.c:556-580 (csort_stream_apply_permutation, OpenMP Phase 2 writeback); also the serial fallback at :690-700

```diff
--- a/src/csort/csort_stream.c
+++ b/src/csort/csort_stream.c
@@ -414,6 +414,7 @@
        vars_per_batch controls how many variables are loaded at a time */
     #ifdef _OPENMP
     int success = 1;
+    int write_err = 0;
     double **batch_bufs = NULL;
     int nthreads = 0;
     int actual_batch_size = vars_per_batch;
@@ -550,33 +551,27 @@
             batch_scatter_time = ctools_timer_seconds() - t_batch_scatter;
             t_scatter += batch_scatter_time;
 
-            /* Phase 2: Write back all variables in batch (parallel) */
+            /* Phase 2: Write back all variables in batch (parallel). Checked
+             * store (IO_VSTORE_FN): under SD_FASTMODE SF_vstore is the unchecked
+             * callback, which must not write different variables at once. */
             double t_batch_writeback = ctools_timer_seconds();
+            const ST_IIID store_fn = (_stata_)->safestore;
 
             #pragma omp parallel for schedule(static) num_threads(batch_threads)
             for (int b = 0; b < this_batch; b++) {
                 size_t var_idx = batch_start + b;
                 double *num_buf = batch_bufs[b];
                 ST_int stata_var = (ST_int)numeric_vars[var_idx];
-                size_t i;
+                int err = 0;
 
-                size_t nobs_batch8 = (nobs / 8) * 8;
-                for (i = 0; i < nobs_batch8; i += 8) {
-                    /* Prefetch next chunk with non-temporal hint (won't be reused) */
-                    CTOOLS_PREFETCH_NTA(&num_buf[i + 64]);
-
-                    /* Store 8 values to Stata */
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 0], num_buf[i + 0]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 1], num_buf[i + 1]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 2], num_buf[i + 2]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 3], num_buf[i + 3]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 4], num_buf[i + 4]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 5], num_buf[i + 5]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 6], num_buf[i + 6]);
-                    SF_vstore(stata_var, (ST_int)obs_map[i + 7], num_buf[i + 7]);
+                for (size_t i = 0; i < nobs && !err; i++) {
+                    /* Prefetch ahead with non-temporal hint (won't be reused) */
+                    if ((i & 7) == 0) CTOOLS_PREFETCH_NTA(&num_buf[i + 64]);
+                    err = store_fn(stata_var, (ST_int)obs_map[i], num_buf[i]);
                 }
-                for (i = nobs_batch8; i < nobs; i++) {
-                    SF_vstore(stata_var, (ST_int)obs_map[i], num_buf[i]);
+                if (err) {
+                    #pragma omp atomic write
+                    write_err = err;
                 }
             }
 
@@ -584,6 +579,7 @@
             t_writeback += batch_writeback_time;
 
             batch_start += this_batch;
+            if (write_err) break;
         }
 
         /* Free batch buffers */
@@ -687,17 +683,17 @@
                 CTOOLS_PREFETCH_NTA(&num_buf[i + 64]);
 
                 /* Store 8 values to Stata */
-                SF_vstore(stata_var, (ST_int)obs_map[i + 0], num_buf[i + 0]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 1], num_buf[i + 1]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 2], num_buf[i + 2]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 3], num_buf[i + 3]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 4], num_buf[i + 4]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 5], num_buf[i + 5]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 6], num_buf[i + 6]);
-                SF_vstore(stata_var, (ST_int)obs_map[i + 7], num_buf[i + 7]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 0], num_buf[i + 0]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 1], num_buf[i + 1]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 2], num_buf[i + 2]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 3], num_buf[i + 3]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 4], num_buf[i + 4]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 5], num_buf[i + 5]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 6], num_buf[i + 6]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i + 7], num_buf[i + 7]);
             }
             for (i = nobs_batch8; i < nobs; i++) {
-                SF_vstore(stata_var, (ST_int)obs_map[i], num_buf[i]);
+                (_stata_)->safestore(stata_var, (ST_int)obs_map[i], num_buf[i]);
             }
 
             t_var_writeback = ctools_timer_seconds() - t_var_writeback;
@@ -713,6 +709,7 @@
     #ifdef _OPENMP
     {
         int string_threads = ctools_get_max_threads();
+        if (write_err) n_string = 0;
         if (string_threads > (int)n_string) string_threads = (int)n_string;
         if (string_threads < 1) string_threads = 1;
 
@@ -750,6 +747,9 @@
     ctools_aligned_free(inv_perm);
 
     #ifdef _OPENMP
+    if (write_err) {
+        return ctools_spi_status(write_err, STATA_ERR_STATA_WRITE);
+    }
     if (!success) {
         return STATA_ERR_MEMORY;
     }
```

Auditor's original suggestion: In both writeback loops (569-579, and 690-700 for consistency), replace `SF_vstore(...)` with `(_stata_)->safestore(...)`, as cmerge_io.c:72 does. Check each return code: on nonzero, set an `atomic_int` or `#pragma omp atomic write` error flag, break, and return STATA_ERR_STATA_WRITE after the region. Alternatively, call ctools_store_filtered(num_buf, nobs, stata_var, obs_map), which uses IO_VSTORE_FN and validates its arguments. Add a validation case with K >= 8 numeric payload variables, N >= 2M and stream(8).

## F033 [P0, confirmed] cexport sas/spss/sasxport5/sasxport8: the strL width scan writes into a NULL or undersized buffer because stringcap is raised without reallocating s

Location: src/io/cio_readstat.c:751-763 (cio_do_export; the bad call is in the strL scan at 632-649)

```diff
--- a/src/io/cio_readstat.c
+++ b/src/io/cio_readstat.c
@@ -569,7 +569,7 @@ int cio_do_export(const char *kind)
     int *dates = calloc((size_t)nvar, sizeof(int));
     char *s = NULL;
-    size_t stringcap = 0;
+    size_t stringcap = 0, stringneed = 0; /* stringcap is always the real size of s */
     FILE *fp = NULL;
     if (!vars || !dates || !names || !sets)
     {
@@ -748,17 +748,18 @@ int cio_do_export(const char *kind)
             }
         }
-        if (isstr && (size_t)width + 1 > stringcap)
-            stringcap = (size_t)width + 1;
+        if (isstr && (size_t)width + 1 > stringneed)
+            stringneed = (size_t)width + 1;
     }
-    if (stringcap)
+    if (stringneed > stringcap)
     {
-        char *tmp = realloc(s, stringcap);
+        char *tmp = realloc(s, stringneed);
         if (!tmp)
         {
             rc = 909;
             goto done;
         }
         s = tmp;
+        stringcap = stringneed;
     }
```

Auditor's original suggestion: Keep `stringcap` as the true capacity and track the requirement separately: size_t stringneed = 0; /* next to stringcap */ ... if (isstr && (size_t)width + 1 > stringneed) stringneed = (size_t)width + 1; /* replaces lines 751-752 */ } if (stringneed > stringcap) { char *tmp = realloc(s, stringneed); if (!tmp) { rc = 909; goto done; } s = tmp; stringcap = stringneed; } Add a validation case with a str# variable before a strL variable, and one with two strL variables for sasxport8.

## F005 [P1, confirmed] cexport with any strL variable allocates 2,046 bytes for every str# cell, regardless of its width

Location: src/ctools_text.c, ctools_load_text_data(): `int length = islong ? SF_sdatalen(j + 1, obs) : 2045;` followed by `v->data.str[i] = malloc((size_t)length + 1);` (lines 58-62 in the current tree). Callers: cexport_main in src/cexport/cexport_impl.c (the longtext detection loop and the `longtext ? ctools_load_text_data(&g_ctx.filtered, SF_nvars(), (ctools_text_policy){0, 0})` dispatch, around lines 771-786), and cexport_xlsx_main -> xlsx_load_long_text in src/cexport/cexport_xlsx.c (around lines 1545 and 1606-1610).

Size the allocation from the value actually read. For non-strL variables, read into a caller-thread stack buffer, then copy:
```c
char tmp[2046];
...
if (islong) { /* unchanged: malloc(length+1) + SF_strldata */ }
else {
    int rc = SF_sdata(j + 1, obs, tmp);
    if (rc) return rc;
    size_t len = strlen(tmp);
    if (!(v->data.str[i] = malloc(len + 1))) return 920;
    memcpy(v->data.str[i], tmp, len + 1);
}
```
Move the strL `malloc(length + 1)` into the islong branch. For better results, copy into a per-column ctools_string_arena (CTOOLS_STRING_ARENA_STRDUP_FALLBACK) and set v->_arena, so cleanup stays O(1).

## F010 [P1, confirmed] Repeated sort key in streaming mode overflows two heap arrays and leaves keys sorted but all other columns unsorted (SORT-8 still present)

Location: src/csort/csort_stream.c, csort_stream_sort: line 927 (`nvars_nonkey = nvars - nkeys`), the allocations at 930-931, and the unchecked writes at 947 and 951 inside the fill loop (938-955). The sanity check at 958 runs after the overflow and after Phase 4 (875-883) has already stored the sorted keys. Upstream: package/ado/csort.ado:13 (`syntax varlist` keeps duplicates) and :116-126 (each occurrence gets the same index), and csort_impl.c parse_sort_vars (137-210) and csort_main (408), which pass the list to the stream path without de-duplicating.

Remove duplicate keys before either path. In csort_main after parse_sort_vars, compact sort_vars to its first occurrences; sorting by (x, x) equals sorting by x. In csort.ado, `local sortkeys : list uniq varlist` for the plugin indices, keeping `sort `varlist', stable` unchanged. In csort_stream_sort, build `unsigned char *is_key = calloc(nvars + 1, 1)`, set nvars_nonkey to the count of !is_key, and guard the write with `if (nk >= nvars_nonkey) { rc = STATA_ERR_INVALID_INPUT; goto cleanup; }`. Do this before Phase 4 writes anything.

## F011 [P1, confirmed] Streaming csort writes sorted keys before allocating the payload buffers; a realistic allocation failure leaves the dataset's rows scrambled (SORT-9 still present)

Location: src/csort/csort_stream.c. In csort_stream_sort, the Phase 4 key stores (about lines 875-883: ctools_store_filtered / ctools_store_filtered_str for each key) run before every Phase 5 allocation. The allocations that come after the first write are: (a) csort_stream_sort's own nonkey_var_indices/nonkey_str_widths malloc/calloc (about 930-935), which the finding does not list; (b) in csort_stream_apply_permutation, inv_perm (360), numeric_vars/string_vars/string_widths (383-385), batch_bufs and the per-batch nobs-double buffers (432-451, OpenMP) or num_buf (601, non-OpenMP); (c) stream_permute_string_var_flat's calloc(nobs, width+1) (222), called from the parallel string loop (715-725). The non-OpenMP string loop (728-732) ignores the failure, and the function returns STATA_OK. The caller, the streaming branch of csort_main in src/csort/csort_impl.c (about 370-420), just maps rc to a Stata error with no recovery. package/ado/csort.ado runs `plugin call ctools_plugin ...` (about line 166) with no preserve/snapshot.

Allocate everything Phase 5 needs before Phase 4 writes anything: inv_perm, the numeric batch buffers, and a fixed pool of string_threads flat buffers of nobs x (max_string_width+1) bytes, reused for successive string variables instead of one calloc per variable. If the pool does not fit, lower string_threads (down to 1) and retry. If even that fails, return STATA_ERR_MEMORY before the first write. Pass the preallocated resources into csort_stream_apply_permutation, so the only failures possible after the first write are SPI errors; report those as data-modified errors. In the non-OpenMP branch, propagate a stream_permute_string_var failure.

## F021 [P1, confirmed] crangestat: a multi-variable source after '=' desynchronises the positional plugin arguments and overwrites an existing variable

Location: The root cause is in package/ado/crangestat.ado, in _crangestat_impl's stat-spec loop. At lines 205-209, `source_var` takes everything after '='. At lines 217-226, `confirm variable` and `confirm numeric variable` accept a varlist. At lines 240-244, `source_var_list` gains every token while `nsource` goes up by one. `_crangestat_varpos source_indices` at line 286 then emits one position per token.

The C side is src/crangestat/crangestat_impl.c, crangestat_main. Its fixed-count positional parser has no consistency checks: the result-index loop near line 111 only checks the range [1, SF_nvars()], the stat type at line 156 is never validated, and the source map at line 172 is never validated. The compute loops skip a stat whose source map is out of range (lines ~557, 678, 796, 936, 1072), and the store loop (~lines 1165-1181) writes into specs[s].result_var_idx through SF_vstore.

Ado: after `confirm variable`, run `unab source_var : `source_var'` and reject more than one name, e.g. `if `: word count `source_var'' != 1 { di as error "crangestat: one variable per statistic"; exit 103 }`.
C (defence in depth): at lines 147-172, return 198 when `stat_int < 0 || stat_int >= STAT_COUNT_ALL` or `src_idx < 0 || src_idx >= (int)nsource`. Reject (198) any result index equal to key_idx, a source index or a by index. Optionally require that the next remaining token is not a bare integer, so a count mismatch fails instead of shifting fields.

## F026 [P1, confirmed] cimport delimited stringcols()/numericcols() lists are cut at 4 KiB: columns past roughly 1,040 silently lose their forced type

Location: src/cimport/cimport_impl.c:1587-1592 cimport_parse_col_list (char buf[4096]); callers cimport_do_scan:1650-1652 and cimport_do_load:1856-1859; package/ado/cimport.ado:68-74 (numlist expansion), 320-326 (globals), 373 (drop)

```diff
--- a/package/ado/cimport.ado
+++ b/package/ado/cimport.ado
@@ -319,9 +319,13 @@
-    * Pass numericcols/stringcols via global macros (space-separated column numbers)
+    * Pass numericcols/stringcols via global macros (space-separated column
+    * numbers) with their byte lengths, so a truncated read fails in C
+    capture macro drop CIMPORT_NUMCOLS CIMPORT_STRCOLS CIMPORT_NUMCOLS_LEN CIMPORT_STRCOLS_LEN
     if "`numcols_str'" != "" {
         global CIMPORT_NUMCOLS "`numcols_str'"
+        global CIMPORT_NUMCOLS_LEN = strlen("`numcols_str'")
     }
     if "`strcols_str'" != "" {
         global CIMPORT_STRCOLS "`strcols_str'"
+        global CIMPORT_STRCOLS_LEN = strlen("`strcols_str'")
     }
@@ -373,1 +377,1 @@
-    capture macro drop CIMPORT_NUMCOLS CIMPORT_STRCOLS
+    capture macro drop CIMPORT_NUMCOLS CIMPORT_STRCOLS CIMPORT_NUMCOLS_LEN CIMPORT_STRCOLS_LEN
--- a/src/cimport/cimport_impl.c
+++ b/src/cimport/cimport_impl.c
@@ -1586,15 +1586,30 @@
-/* Parse column list from Stata global macro */
-static int *cimport_parse_col_list(const char *macro_name, int *out_count) {
-    char buf[4096];
-    *out_count = 0;
-
-    if (SF_macro_use((char *)macro_name, buf, sizeof(buf)) != 0 || strlen(buf) == 0) {
-        return NULL;
-    }
+/* Parse a column list from a Stata global macro. The ado passes its byte
+ * length in <name>_LEN so truncation fails instead of dropping columns.
+ * Returns 0 with *out_cols NULL when no list was given. */
+static ST_retcode cimport_parse_col_list(const char *macro_name, int **out_cols, int *out_count) {
+    char name[64], length_text[32] = "", *stop;
+    *out_cols = NULL;
+    *out_count = 0;
+    snprintf(name, sizeof(name), "%s_LEN", macro_name);
+    if (SF_macro_use(name, length_text, sizeof(length_text) - 1) != 0 || !*length_text) return 0;
+    long length = strtol(length_text, &stop, 10);
+    if (stop == length_text || *stop || length <= 0 || length >= INT_MAX) return 198;
+    char *buf = calloc((size_t)length + 1, 1);
+    if (!buf) return 920;
+    if (SF_macro_use((char *)macro_name, buf, (int)length + 1) != 0 || strlen(buf) != (size_t)length) {
+        free(buf);
+        return 198;
+    }
 
     if (strcmp(buf, "_all") == 0) {
-        int *all = calloc(1, sizeof(int));
-        if (all) *out_count = 1;
-        return all;
+        free(buf);
+        if (!(*out_cols = calloc(1, sizeof(int)))) return 920;
+        *out_count = 1;
+        return 0;
     }
@@ -1610,5 +1625,5 @@
-    if (count == 0) return NULL;
+    if (count == 0) { free(buf); return 0; }
 
     int *cols = ctools_safe_malloc2(count, sizeof(int));
-    if (!cols) return NULL;
+    if (!cols) { free(buf); return 920; }
@@ -1624,4 +1639,5 @@
         if (errno == ERANGE || endptr == p || val <= 0 || val > INT_MAX) {
             free(cols);
-            return NULL;  /* Invalid column number */
+            free(buf);
+            return 198;  /* Invalid column number */
         }
@@ -1632,4 +1648,6 @@
-    *out_count = idx;
-    return cols;
+    free(buf);
+    *out_cols = cols;
+    *out_count = idx;
+    return 0;
 }
@@ -1649,4 +1667,6 @@ static ST_retcode cimport_do_scan(
-    int num_force_num = 0, num_force_str = 0;
-    int *force_num = cimport_parse_col_list("CIMPORT_NUMCOLS", &num_force_num);
-    int *force_str = cimport_parse_col_list("CIMPORT_STRCOLS", &num_force_str);
+    int num_force_num = 0, num_force_str = 0, *force_num = NULL, *force_str = NULL;
+    ST_retcode list_rc = cimport_parse_col_list("CIMPORT_NUMCOLS", &force_num, &num_force_num);
+    if (!list_rc) list_rc = cimport_parse_col_list("CIMPORT_STRCOLS", &force_str, &num_force_str);
+    if (list_rc) { free(force_num); return list_rc; }
   (apply the same three-line change to the two calls in cimport_do_load ~1857-1859)
```

Auditor's original suggestion: Transfer the list with its length, as cexport_read_local does. The ado sets __cimport_strcols_len/__cimport_strcols locals; C reads the length, mallocs length+1, calls SF_macro_use(name, buf, length+1), and returns 198 when strlen(buf) != length. Alternatively pass ranges (a/b pairs) instead of an expanded numlist. At minimum use sizeof(buf)-1 and fail on truncation instead of silently ignoring the list.

## F028 [P1, confirmed] XLSX writer ignores chunk-formatting and task-submit failures: rows silently dropped and DEFLATE stream left unflushed, export reported as saved

Location: src/cexport/cexport_xlsx.c: xlsx_format_chunk (~1042-1153; escape malloc at 1050-1051, strbuf grow-failure semantics at 529-588), xlsx_gen_worksheet_segments (~1321-1365: submit/wait/serial results ignored, chunks[c].sb.data never checked, chunk_bufs NULL not handled), xlsx_comp_worker (686-705: empty last segment skips SYNC_FLUSH), xlsx_compress_parallel (781-785: submit result ignored). The edit path reaches the same code via cexport_xlsx_edit.inc:586-603.

```diff
--- a/src/cexport/cexport_xlsx.c
+++ b/src/cexport/cexport_xlsx.c
@@ -683,15 +683,29 @@ static void *xlsx_comp_worker(void *arg) {
         return NULL;
     }
 
-    for (size_t i = 0; i < task->num_segs; i++) {
+    /* The task's closing flush must go on its last non-empty segment, or the
+     * buffered tail is lost and the stream is not byte-aligned. */
+    size_t last = task->num_segs;
+    for (size_t i = task->num_segs; i-- > 0;) {
+        if (task->segments[task->first_seg + i].len > 0) { last = i; break; }
+    }
+    if (last == task->num_segs) {
+        if (task->is_last &&
+            tdefl_compress_buffer(task->comp, NULL, 0, TDEFL_FINISH) != TDEFL_STATUS_DONE)
+            task->failed = true;
+        return NULL;
+    }
+
+    for (size_t i = 0; i <= last; i++) {
         const xlsx_segment *seg = &task->segments[task->first_seg + i];
         if (seg->len == 0) continue;
 
         tdefl_flush flush;
-        if (i < task->num_segs - 1) {
+        if (i < last) {
             flush = TDEFL_NO_FLUSH;
         } else {
             flush = task->is_last ? TDEFL_FINISH : TDEFL_SYNC_FLUSH;
@@ -780,7 +794,8 @@ static bool xlsx_compress_parallel(const xlsx_worksheet_segments *segs,
     if (use_parallel) {
         for (int t = 0; t < num_threads; t++) {
-            ctools_persistent_pool_submit(pool, xlsx_comp_worker, &tasks[t]);
+            if (ctools_persistent_pool_submit(pool, xlsx_comp_worker, &tasks[t]) != 0)
+                xlsx_comp_worker(&tasks[t]);   /* each task owns its compressor */
         }
         ctools_persistent_pool_wait(pool);
     } else {
@@ -1036,6 +1051,7 @@ typedef struct {
     const char *escaped_missing;
     size_t escaped_missing_len;
     StringBuffer sb;            /* output buffer */
+    bool failed;                /* escape buffer or sb growth failed */
 } xlsx_chunk_args_t;
@@ -1048,7 +1064,10 @@ static void *xlsx_format_chunk(void *arg) {
     /* Per-thread escape buffer for strings (heap-allocated for thread safety) */
     char *escaped = malloc(MAX_CELL_CONTENT * 6 + 1);
-    if (!escaped) return (void *)(intptr_t)-1;
+    if (!escaped) {
+        a->failed = true;
+        return (void *)(intptr_t)-1;
+    }
@@ -1149,5 +1168,9 @@ static void *xlsx_format_chunk(void *arg) {
     }
 
     free(escaped);
+    if (!sb->data) {            /* a strbuf grow failed and dropped the chunk */
+        a->failed = true;
+        return (void *)(intptr_t)-1;
+    }
     return NULL;
 }
@@ -1321,19 +1344,27 @@ static bool xlsx_gen_worksheet_segments(XLSXExportContext *ctx,
     if (use_parallel) {
         for (size_t c = 0; c < num_chunks; c++) {
-            ctools_persistent_pool_submit(pool, xlsx_format_chunk, &chunks[c]);
+            if (ctools_persistent_pool_submit(pool, xlsx_format_chunk, &chunks[c]) != 0)
+                xlsx_format_chunk(&chunks[c]);
         }
         ctools_persistent_pool_wait(pool);
     } else {
         for (size_t c = 0; c < num_chunks; c++) {
             xlsx_format_chunk(&chunks[c]);
         }
     }
 
+    bool format_ok = preamble.data != NULL;
+    for (size_t c = 0; c < num_chunks; c++)
+        if (chunks[c].failed || !chunks[c].sb.data) format_ok = false;
+
     /* Build segment array: preamble + num_chunks + footer */
     size_t total_segments = 1 + num_chunks + 1;
-    segs->segments = (xlsx_segment *)malloc(total_segments * sizeof(xlsx_segment));
-    if (!segs->segments) {
+    segs->segments = format_ok ? (xlsx_segment *)malloc(total_segments * sizeof(xlsx_segment)) : NULL;
+    segs->chunk_bufs = format_ok ? (StringBuffer *)malloc(num_chunks * sizeof(StringBuffer)) : NULL;
+    if (!segs->segments || !segs->chunk_bufs) {
+        free(segs->segments);
+        free(segs->chunk_bufs);
+        memset(segs, 0, sizeof(*segs));
         for (size_t c = 0; c < num_chunks; c++) {
             free(chunks[c].sb.data);
         }
@@ -1353,14 +1384,11 @@ static bool xlsx_gen_worksheet_segments(XLSXExportContext *ctx,
     /* Transfer chunk StringBuffers to segments */
-    segs->chunk_bufs = (StringBuffer *)malloc(num_chunks * sizeof(StringBuffer));
     segs->num_chunks = num_chunks;
     size_t total_size = preamble.len + footer_len;
 
     for (size_t c = 0; c < num_chunks; c++) {
         segs->segments[1 + c].data = chunks[c].sb.data;
         segs->segments[1 + c].len = chunks[c].sb.len;
-        if (segs->chunk_bufs) {
-            segs->chunk_bufs[c] = chunks[c].sb;
-        }
+        segs->chunk_bufs[c] = chunks[c].sb;
         total_size += chunks[c].sb.len;
     }
```

Auditor's original suggestion: Add `bool failed;` to xlsx_chunk_args_t. In xlsx_format_chunk set a->failed = true on escape-buffer malloc failure and also when `sb->data == NULL` before returning. In xlsx_gen_worksheet_segments, run a chunk inline when ctools_persistent_pool_submit returns nonzero. Treat a nonzero ctools_persistent_pool_wait, any chunks[c].failed, a NULL chunks[c].sb.data, a NULL preamble.data, or a NULL segs->chunk_bufs as failure: free every chunk buffer, the preamble and the column arrays, and return false so the callers return 920.

## F045 [P1, confirmed] fweights summing above 2^31-1: (ST_int)N_eff conversion is UB, and the dof arithmetic that follows overflows (prior IV-23, still present in civreghdfe)

Location: src/civreghdfe/civreghdfe_impl.c:1441 (do_iv_regression, call to ivest_compute_2sls with (ST_int)N_eff), plus impl.c:1527/1535 (df_r_val). Callees take ST_int N_eff: civreghdfe_estimate.c:663 (ivest_compute_2sls; int arithmetic at 788, 1326, 1536), civreghdfe_vce.c:93/124 (ivvce_compute_twoway), vce.c:189/220 (ivvce_compute_full), civreghdfe_tests.c:569/594 (civreghdfe_compute_underid_test).

```diff
--- a/src/civreghdfe/civreghdfe_impl.c
+++ b/src/civreghdfe/civreghdfe_impl.c
@@ -1438,7 +1438,7 @@
     if (!rc) rc = ivest_compute_2sls(
         y_dem, X_exog_dem, X_endog_dem, Z_dem,
         weights_c, weight_type,
-        N, (ST_int)N_eff, K_exog, K_endog, K_iv,
+        N, N_eff, K_exog, K_endog, K_iv,
         beta, V, first_stage_F,
@@ -1523,9 +1523,12 @@
     ST_int df_r_val;
+    /* N_eff (sum of fweights) may exceed INT_MAX: subtract in double, cap */
+    ST_double df_r_obs = N_eff - K_for_dof - df_a - dofminus;
+    if (df_r_obs > 2147483647.0) df_r_obs = 2147483647.0;
     if (kiefer) {
         /* Kiefer VCE: df_r = N_eff - K - df_a - dofminus (not cluster-based) */
-        df_r_val = (ST_int)(N_eff - K_for_dof - df_a - dofminus);
+        df_r_val = (ST_int)df_r_obs;
     } else if (has_cluster2) {
@@ -1534,5 +1537,5 @@
     } else {
-        df_r_val = (ST_int)(N_eff - K_for_dof - df_a - dofminus);
+        df_r_val = (ST_int)df_r_obs;
     }
--- a/src/civreghdfe/civreghdfe_estimate.h
+++ b/src/civreghdfe/civreghdfe_estimate.h
@@ -188,1 +188,1 @@
-    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
+    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
--- a/src/civreghdfe/civreghdfe_estimate.c
+++ b/src/civreghdfe/civreghdfe_estimate.c
@@ -663,1 +663,1 @@
-    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
+    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
@@ -1326,3 +1326,3 @@
-        ST_int df_r_s = N_eff - K_dof - df_a - dofminus;
+        ST_double df_r_s = N_eff - K_dof - df_a - dofminus;
         if (df_r_s <= 0) df_r_s = 1;
         dof_adj_s = (ST_double)(N_eff - dofminus) / (ST_double)df_r_s;
@@ -1536,1 +1536,1 @@
-            ST_int denom_df = N_eff - Kiv_dof - df_a - dofminus;
+            ST_double denom_df = N_eff - Kiv_dof - df_a - dofminus;
--- a/src/civreghdfe/civreghdfe_vce.h
+++ b/src/civreghdfe/civreghdfe_vce.h
@@ -66,1 +66,1 @@
-    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
+    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
@@ -130,1 +130,1 @@
-    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
+    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
--- a/src/civreghdfe/civreghdfe_vce.c
+++ b/src/civreghdfe/civreghdfe_vce.c
@@ -93,1 +93,1 @@
-    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
+    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
@@ -124,2 +124,2 @@
-    ST_int df_r = N_eff - K_total - df_a;
+    ST_double df_r = N_eff - K_total - df_a;
     if (df_r <= 0) df_r = 1;
@@ -189,1 +189,1 @@
-    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
+    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
@@ -220,2 +220,2 @@
-    ST_int df_r = N_eff - K_total - df_a - dofminus;
+    ST_double df_r = N_eff - K_total - df_a - dofminus;
     if (df_r <= 0) df_r = 1;
--- a/src/civreghdfe/civreghdfe_tests.h
+++ b/src/civreghdfe/civreghdfe_tests.h
@@ -46,1 +46,1 @@
-    ST_int N_eff,
+    ST_double N_eff,
--- a/src/civreghdfe/civreghdfe_tests.c
+++ b/src/civreghdfe/civreghdfe_tests.c
@@ -569,1 +569,1 @@
-    ST_int N_eff,
+    ST_double N_eff,
@@ -594,2 +594,2 @@
-    ST_int df_resid = N_eff - K_iv - df_a;
+    ST_double df_resid = N_eff - K_iv - df_a;
     if (df_resid <= 0) df_resid = 1;
```

Auditor's original suggestion: Minimal release fix, right after N_eff = sum_fw in impl.c: if (!(sum_fw <= (ST_double)INT_MAX - 1.0)) { SF_error("civreghdfe: the sum of fweights exceeds 2,147,483,646\n"); iv_rc = 198; goto iv_cleanup; } Proper fix, mirroring the creghdfe repair: - Change the `ST_int N_eff` parameters of ivest_compute_2sls, ivvce_compute_full/_twoway and civreghdfe_compute_underid_test to ST_double. - Do the dof arithmetic in double. - Compute df_r_val in double and clamp it to [1, INT_MAX] before converting.

## F001 [P2, uncertain] Row-staged direct string reads can silently overwrite already-staged cells when the width hint is narrower than the variable's real width

Location: src/ctools_rowstage.c:224-228 (ctools_rowstage_load, string branch of the tile loop); hint source src/ctools_transport.c ctools_io_ctools_read_strw_from_stata (16 KiB buffer); caller src/cmerge/cmerge_impl.c:1184 and :1214 (str_widths = NULL)

```diff
--- a/src/ctools_rowstage.c
+++ b/src/ctools_rowstage.c
@@ -195,8 +195,6 @@
         tid = (size_t)omp_get_thread_num();
 #endif
         rowstage_share(nobs, team, tid, &start, &end);
-        /* Each worker owns its rows; only the last one may read into the guard. */
-        const char *safe_end = rs->rows + end * stride + (end == nobs ? IO_SLOT_GUARD : 0);
         for (size_t tile = start; tile < end && !atomic_load(&error); tile += tile_rows)
         {
             size_t count = end - tile < tile_rows ? end - tile : tile_rows;
@@ -221,11 +219,15 @@
                 }
                 else
                 {
+                    /* Bounced reads only: tiles are staged column by column, so
+                     * an overrun past this slot can reach cells of the next row
+                     * that are already staged. A direct read under a narrow
+                     * width hint could overwrite them undetected. */
                     ST_IIIS read_fn = (_stata_)->sdata;
                     size_t width = (size_t)rs->width[j];
                     for (size_t i = 0; i < count && !io_rc; i++, slot += stride)
-                        io_rc = ctools_io_read_string_slot_direct(read_fn, idx, obs + (ST_int)i,
-                                                                  slot, width, safe_end);
+                        io_rc = ctools_io_read_string_slot(read_fn, idx, obs + (ST_int)i,
+                                                           slot, width);
                 }
                 if (io_rc)
                 {
```

Auditor's original suggestion: Minimal fix: in the row-major stager, always use the bounced read, so the SPI never writes outside the slot's own width+1 bytes. Replace lines 226-228 with: for (size_t i = 0; i < count && !io_rc; i++, slot += stride) io_rc = ctools_io_read_string_slot(read_fn, idx, obs + (ST_int)i, slot, width); This drops the speedup from direct reads, for row staging only. To keep direct reads, fill each tile in address order instead: rows as the outer loop, and inside each row the columns in offset order (numerics first, then strings, as rowstage_plan lays them out).

## F002 [P2, confirmed] __ctools_strw width hints are still read into a fixed 16 KiB stack buffer; a cut entry becomes a wrong (too narrow) width

Location: src/ctools_transport.c, ctools_io_ctools_read_strw_from_stata (around lines 109-136): the fixed `char buf[16384]`, the call `SF_macro_use("___ctools_strw", buf, sizeof(buf) - 1)`, the check `if (rc != 0 || buf[0] == '\0')` that cannot detect a cut value, and the strtol parse loop. The cut widths reach callers through ctools_io_plan_load (src/ctools_data_plan.c, around line 267, the auto-detect branch when str_widths == NULL). From there they go to the column loaders in src/ctools_data_io.c and to rowstage_plan (src/ctools_rowstage.c, around lines 116-119, when cmerge passes NULL widths). The producer is package/ado/_ctools_strw.ado, which writes one comma-separated entry per variable, "0" for numeric.

Treat possibly cut metadata as no hint. After the SF_macro_use call, add: if (rc != 0 || buf[0] == '\0' || strlen(buf) + 2 >= sizeof(buf)) return NULL; so every variable takes the arena path. The better fix follows cexport_read_local: have _ctools_strw also c_local __ctools_strw_len (the byte length), malloc length+1, read with SF_macro_use(name, buf, length + 1), and return NULL unless strlen(buf) == length.

## F006 [P2, confirmed] submit_batch partially runs a batch before failing; the pairs-sort serial retry then re-runs scatter tasks with advanced offsets (heap overflow)

Location: The root cause is in src/ctools_threads.c, function ctools_persistent_pool_submit_batch (about lines 300-322; the failing loop is about 311-319). It enqueues items one at a time. When item k fails it calls ctools_persistent_pool_wait, so items 0..k-1 run to completion, and then it returns -1. The defect shows up in src/ctools_sort_pairs.c, function ctools_sort_order_pairs, at the scatter fallback (about lines 206-215). That fallback serially reruns every scatter task without rebuilding all_offsets. The out-of-bounds write is in ctools_scatter_thread at line 79 (`dst[offsets[digit]++] = src[i];`). The caller that can reach it is the cmerge preserve_order path at src/cmerge/cmerge_impl.c:849. A secondary leak-only site is src/cimport/cimport_impl.c, in the fallbacks at about lines 1243-1264 (cimport_parse_chunk_parallel) and 1551-1572 (cimport_build_cache_worker). None of src/ctools_threads.c, src/ctools_threads.h or src/ctools_sort_pairs.c differs from HEAD.

Make batch submission all-or-nothing. Allocate and link all `count` work items first; on any allocation failure free the partial list and return -1 without enqueuing anything. Then append the whole list under one queue_mutex acquisition (queue_size += count; pending_items += count; pthread_cond_broadcast(&work_available)). With that contract every existing 'serial fallback after failed submit' path becomes correct. Document in ctools_threads.h that -1 means no item ran.

## F012 [P2, confirmed] cmerge's hard-coded `_merge_temp_filter` plus `merge_var_idx = c(k)` makes the plugin write _merge codes into a user variable and skip reordering it (MRG-10 still present)

Location: package/ado/cmerge.ado, _cmerge_impl: line 1101 (fixed name `_merge_temp_filter` with no existence check), line 1152 (_cmerge_template builds it), the append of the template, line 1258 (`local merge_var_idx = c(k)`), and lines 1380-1389 (`keep if inlist(`filter_var', ...)` at 1384, then `qui drop `temp_merge_var'` at 1388). src/cmerge/cmerge_impl.c, cmerge_execute: argument parsing (merge_var_idx is accepted with no range check), the stage_vars loop (line 1174), and the _merge write loop (lines 1393-1405, safestore at 1404). build/cmerge.ado is identical to the package copy.

In cmerge.ado, use `tempvar __mfilter` for the scratch variable, and after the append pass `merge_var_idx` from `st_varindex("`generate'`__mfilter'")`, not `c(k)`. In cmerge_execute, also reject `merge_var_idx` (and every non-shared keepusing placeholder) that is <= the original master variable count or > SF_nvars(), returning 198 before any write. The original count could be a new `master_k <K>` argument.

## F016 [P2, confirmed] cencode label-file path is cut at the first space: labels go to a truncated path and cencode fails (DSTR-8 still open)

Location: src/cencode/cencode_impl.c, cencode_main: lines 430-431 read the path (`ctools_parse_string_option(args, "labelfile", ...)`). The truncated path is used by cencode_stream (the `ctools_label_write_file` call near line 350) and by the strL/bulk path in cencode_main (near line 727). Root cause on the ado side: package/ado/cencode.ado lines 192-193 build `labelfile=`__labelfile'` with no quotes, and line 203 passes it in the plugin call. Tokenizer: src/ctools_parse.c, `ctools_args_next` (about lines 15-33) and `option()`/`ctools_parse_string_option` (about lines 35-86). The file is opened by `ctools_label_write_file` in src/ctools_hash.c (about lines 231-253, `fopen(filepath, "wb")`). The failure appears in the Mata function `_cencode_apply_labels` (package/ado/cencode.ado about line 468, `fopen(path, "r")`).

Pass the path outside the argument string, the same way __ce_map is passed. In cencode.ado, replace `labelfile=...` with `local __ce_labelfile `"`__labelfile'"'` before the plugin call, and clear it afterwards. In cencode_main:
    char labelfile[1024] = "";
    if (SF_macro_use("___ce_labelfile", labelfile, sizeof(labelfile)) != 0 || !labelfile[0] ||
        strlen(labelfile) >= sizeof(labelfile) - 1) { SF_error("cencode: invalid label file path\n"); return 198; }
A smaller alternative is to quote the value: local labelfile_code `"labelfile="`__labelfile'""' and put the plugin-call string in compound quotes. ctools_args_next/option() already honour and strip double quotes.

## F017 [P2, confirmed] cwinsor allocates one nobs-sized work buffer per CPU even when only nvars threads use them (SMP-2 still open)

Location: src/cwinsor/cwinsor_impl.c, cwinsor_main: num_threads is set at line 398 (`num_threads = ctools_get_max_threads();`). The per-thread buffers are allocated at lines 560-586 (`malloc(max_group_size * sizeof(double))` for each t < num_threads). The no-by() loop is at lines 601-639 (`#pragma omp parallel for schedule(static)` over v < nvars, using work_arrays[omp_get_thread_num()]). The buffers are freed at lines 733-735. These line numbers match the report.

Compute the worker count from the available work:
    size_t tasks = ngroups == 1 ? nvars : nvars * ngroups;
    int nwork = (int)(tasks < (size_t)num_threads ? tasks : (size_t)num_threads);
Allocate nwork buffers, and free them with the same count in both cleanup loops. Add num_threads(nwork) to both `#pragma omp parallel for` lines (603 and 644), so omp_get_thread_num() < nwork still holds. Optionally allocate each buffer lazily on the thread's first use.

## F019 [P2, not reproduced] crangestat stores result columns from parallel OpenMP threads through the unchecked SPI store (SF_vstore is (_stata_)->store under SD_FASTMODE)

Location: src/crangestat/crangestat_impl.c:1161-1184 (crangestat_main store phase; SF_vstore at :1174 and :1179)

```diff
--- a/src/crangestat/crangestat_impl.c
+++ b/src/crangestat/crangestat_impl.c
@@ -1159,6 +1159,10 @@
      * fwd_perm[sorted_idx] = original_idx, so data[sorted_idx] is sequential.
      */
     int store_rc = 0;
+    /* Checked store, as IO_VSTORE_FN: under SD_FASTMODE SF_vstore is the
+     * unchecked callback, which must not write different variables at once
+     * and does not resolve plugin-varlist indices past hidden variables. */
+    const ST_IIID store_fn = (_stata_)->safestore;
     #ifdef _OPENMP
     #pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads) reduction(max:store_rc) if(nstats > 1)
     #endif
@@ -1171,12 +1175,12 @@
             /* Sequential reads from data[sorted_idx], scattered writes via obs_map */
             for (size_t sorted_idx = 0; sorted_idx < nobs; sorted_idx++) {
                 perm_idx_t orig_idx = fwd_perm[sorted_idx];
-                local_rc = SF_vstore(result_idx, (ST_int)obs_map[orig_idx], data[sorted_idx]);
+                local_rc = store_fn(result_idx, (ST_int)obs_map[orig_idx], data[sorted_idx]);
                 if (local_rc) break;
             }
         } else {
             for (size_t i = 0; i < nobs; i++) {
-                local_rc = SF_vstore(result_idx, (ST_int)obs_map[i], data[inv_perm[i]]);
+                local_rc = store_fn(result_idx, (ST_int)obs_map[i], data[inv_perm[i]]);
                 if (local_rc) break;
             }
         }
```

Auditor's original suggestion: Use the checked callback in the parallel loop, as the transport does: before the region add `const ST_IIID store_fn = (_stata_)->safestore;` and replace both SF_vstore calls with `store_fn(result_idx, ...)`. Alternatively call ctools_store_filtered()/IO_VSTORE_FN. Audit the remaining SF_vstore call sites in parallel code (see the crangejoin and csort stream findings).

## F020 [P2, not reproduced] crangejoin write and rollback store numeric columns concurrently through the unchecked SF_vstore

Location: src/crangejoin/crangejoin_impl.c:74-86 (restore_master) and :228-251 (crangejoin_main "write")

```diff
--- a/src/crangejoin/crangejoin_impl.c
+++ b/src/crangejoin/crangejoin_impl.c
@@ -71,6 +71,8 @@
     if ((size_t)SF_nobs() < n || SF_nvars() < state.nm) return 198;
     int rc = 0;
     const size_t tile = 4096;
+    /* Checked store: workers write different variables at once. */
+    const ST_IIID store_fn = (_stata_)->safestore;
     #pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if(n > 10000) reduction(max:rc)
     for (size_t begin = 0; begin < n; begin += tile) {
         size_t end = begin + tile < n ? begin + tile : n;
@@ -79,7 +81,7 @@
             for (size_t i = begin; i < end; i++) {
                 int err = v->type == STATA_TYPE_STRING
                     ? SF_sstore(col+1, (ST_int)(i+1), v->data.str[i])
-                    : SF_vstore(col+1, (ST_int)(i+1), v->data.dbl[i]);
+                    : store_fn(col+1, (ST_int)(i+1), v->data.dbl[i]);
                 if (err > rc) rc = err;
             }
         }
@@ -225,6 +227,8 @@
          * Workers own disjoint row ranges, avoiding cross-column false sharing. */
         const size_t tile = 4096, nm = state.master.data.nobs;
         const size_t nwrite = fault_store ? state.nout / 2 : state.nout;
+        /* Checked store: workers on different tiles write different variables. */
+        const ST_IIID store_fn = (_stata_)->safestore;
         #pragma omp parallel for schedule(static) num_threads(ctools_get_max_threads()) if(state.nout > 10000) reduction(max:rc)
         for (size_t begin = 0; begin < nwrite; begin += tile) {
             size_t end = begin + tile < nwrite ? begin + tile : nwrite;
@@ -248,7 +252,7 @@
                         if (v->type == STATA_TYPE_STRING)
                             err = SF_sstore(dest, (ST_int)(j+1), (ismaster || matched) ? v->data.str[source] : "");
                         else
-                            err = SF_vstore(dest, (ST_int)(j+1), (ismaster || matched) ? v->data.dbl[source] : SV_missval);
+                            err = store_fn(dest, (ST_int)(j+1), (ismaster || matched) ? v->data.dbl[source] : SV_missval);
                         if (err > rc) rc = err;
                     }
                 }
```

Auditor's original suggestion: Hoist `const ST_IIID store_fn = (_stata_)->safestore;` and use it at lines 82 and 251 (SF_sstore may stay; the transport already writes strings concurrently through sstore).

## F022 [P2, confirmed] cpsmatch nearest-neighbour search can index sorted_controls[-1] (or loop forever) when a distance equals the DBL_MAX sentinel

Location: src/cpsmatch/cpsmatch_kernels.c, cpsmatch_find_nearest_neighbors_fast, first two-pointer loop: the side choice `if (dist_left <= dist_right)` (~line 224), the caliper reset `left = -1` (~line 228) and the read `sorted_controls[left].orig_idx` (~line 232). The garbage index is then used at `ws->match_ids[i] = control_indices[ws->candidates[i].idx]` (~line 314). In the caller cpsmatch_main (src/cpsmatch/cpsmatch_impl.c, the NN/kernel `omp parallel` block), the bad value reaches `obs_map[result.match_ids[0]]` (~line 722) and the write `local_weights[match_filtered_idx] += result.match_weights[m]` (~line 728), or the atomic `out_weight[...]` write (~731). It is not at impl:597 as the finding says. The ties-loop part (~line 286) cannot be reached from the shipped ado.

Choose the side explicitly: `int take_left = left >= 0 && (right >= n_controls || dist_left <= dist_right);`. Use it in the first loop. In the ties loop, use `if (take_left && dist_left <= last_dist + 1e-12) ... else if (right < n_controls && dist_right <= last_dist + 1e-12) ... else break;`. When the caliper rejects one side, set that side exhausted; when both are exhausted, break.

## F024 [P2, confirmed] strL cache builder writes one pointer past string_data (and rowrange() reads an uninitialized header pointer) when chunk 0 holds no rows

Location: src/cimport/cimport_impl.c: cimport_parse_csv multi-chunk branch (~1209-1221, chunking starts at skip_rows result), cimport_parse_chunk_parallel:755 (is_header_row only for chunk_id 0), total_rows-- at ~1281, cimport_build_cache_worker:1455 (unbounded string_data[row_idx++]); cimport_delimiters.inc:137/182 (cimport_select_rows write=header); cimport_infer_column_types:522-524 (reads chunks[0].rows[0])

```diff
--- a/src/cimport/cimport_context.h
+++ b/src/cimport/cimport_context.h
@@ -124,4 +124,5 @@ typedef struct CImportContext {
     int header_row;         /* 1-based row number for header (1 = first row, 0 = no header) */
     int skip_rows;          /* Number of rows to skip at start of file */
+    size_t first_row_offset; /* Offset of the first row read (header or first data row) */
     CImportBindQuotesMode bindquotes;
--- a/src/cimport/cimport_impl.c
+++ b/src/cimport/cimport_impl.c
@@ -325,4 +325,5 @@ static int cimport_parse_header(CImportContext *ctx) {
     ctx->num_columns = cimport_parse_header_row(ctx, fields, &header_start, &next);
+    ctx->first_row_offset = (size_t)(header_start - ctx->file_data);
     /* Only empty rows: native imports nothing, like an empty file. */
     if (ctx->num_columns == 0) return 0;
@@ -1209,12 +1210,9 @@ static CImportContext *cimport_parse_csv(
     } else {
-        /* Skip rows if requested (for varnames(N) where N > 1) */
-        const char *data_start = ctx->file_data;
-        size_t data_size = ctx->file_size;
-        if (ctx->skip_rows > 0) {
-            data_start = cimport_skip_rows(ctx->file_data, ctx->file_data + ctx->file_size,
-                                           ctx->skip_rows, ctx->quote_char, ctx->bindquotes);
-            data_size = ctx->file_size - (data_start - ctx->file_data);
-        }
+        /* Start at the first row read (after varnames(#) skips and the empty
+         * rows emptylines(skip) drops), so chunk 0 always holds the header
+         * that every consumer skips as chunks[0].rows[0]. */
+        const char *data_start = ctx->file_data + ctx->first_row_offset;
+        size_t data_size = ctx->file_size - ctx->first_row_offset;
 
         size_t chunk_size = data_size / num_chunks;
@@ -1453,4 +1451,7 @@ static void *cimport_build_cache_worker(void *arg) {
                     }
 
+                    if (row_idx >= ctx->total_rows) {
+                        atomic_store(&ctx->error_code,459);free(field_buf);return (void*)1;
+                    }
                     cache->string_data[row_idx++] = str;
--- a/src/cimport/cimport_delimiters.inc
+++ b/src/cimport/cimport_delimiters.inc
@@ -136,3 +136,3 @@ static void cimport_select_rows(CImportContext *ctx) {
     CImportParsedChunk *chunk = &ctx->chunks[c];
-    size_t header = c == 0 && ctx->has_header ? 1 : 0, write = header;
+    size_t header = c == 0 && ctx->has_header && chunk->num_rows ? 1 : 0, write = header;
     cimport_init_col_stats(chunk->col_stats, chunk->num_col_stats);
```

Auditor's original suggestion: Start chunking at the header row so that chunk 0 always begins with it. Have cimport_parse_header store `header_start - ctx->file_data` in the context. In the multi-chunk branch, use that offset as start_offset instead of the skip_rows result; only empty rows, which are dropped anyway, lie before it. Then add defensive guards: - In cimport_build_cache_worker: `if (row_idx >= ctx->total_rows) { atomic_store(&ctx->error_code, 459); free(field_buf); return (void *)1; }` before line 1455. - In cimport_select_rows: `write = (header && chunk->num_rows) ? 1 : 0`.

## F025 [P2, uncertain] Scan call chain keeps about 2.2 MB of CTOOLS_MAX_COLUMNS stack arrays on the Stata thread

Location: src/cimport/cimport_impl.c: cimport_parse_header:321/347-348 (fields, next_fields, candidates), cimport_parse_csv:1133 (single-chunk field_buf), cimport_do_scan:1645 (macro_val[65536]) and :1743 (headings[CTOOLS_MAX_COLUMNS]); sibling cimport_parse_chunk_parallel:732 (field_buf on pool worker threads)

```diff
--- a/src/cimport/cimport_impl.c
+++ b/src/cimport/cimport_impl.c
@@ -320,5 +320,7 @@
 static int cimport_parse_header(CImportContext *ctx) {
-    CImportFieldRef fields[CTOOLS_MAX_COLUMNS];
+    /* fields and next_fields share one heap block (off the Stata thread's stack). */
+    CImportFieldRef *fields = ctools_safe_malloc2(2 * (size_t)CTOOLS_MAX_COLUMNS, sizeof(CImportFieldRef));
     const char *end = ctx->file_data + ctx->file_size;
     const char *header_start, *next;
+    if (!fields) return -1;
 
     ctx->num_columns = cimport_parse_header_row(ctx, fields, &header_start, &next);
     /* Only empty rows: native imports nothing, like an empty file. */
-    if (ctx->num_columns == 0) return 0;
+    if (ctx->num_columns == 0) { free(fields); return 0; }
@@ -346,6 +348,7 @@
         ctx->has_header = false;
-        CImportFieldRef next_fields[CTOOLS_MAX_COLUMNS];
-        bool candidates[CTOOLS_MAX_COLUMNS];
+        CImportFieldRef *next_fields = fields + CTOOLS_MAX_COLUMNS;
+        bool *candidates = calloc((size_t)ctx->num_columns, sizeof(bool));
+        if (!candidates) { free(fields); return -1; }
         for (int i = 0; i < ctx->num_columns; i++) {
@@ -378,7 +381,8 @@
         for (int i = 0; i < ctx->num_columns; i++) {
             if (candidates[i]) { ctx->has_header = true; break; }
         }
+        free(candidates);
     }
 
     ctx->columns = ctools_safe_calloc2(ctx->num_columns, sizeof(CImportColumnInfo));
-    if (!ctx->columns) return -1;
+    if (!ctx->columns) { free(fields); return -1; }
@@ -396,4 +400,5 @@
         ctx->columns[i].max_strlen = 0;
     }
 
+    free(fields);
     return 0;
@@ -731,5 +736,9 @@ static void *cimport_parse_chunk_parallel(void *arg) {
     cimport_init_col_stats(chunk->col_stats, ctx->num_columns);
 
-    CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
+    CImportFieldRef *field_buf = ctools_safe_malloc2(CTOOLS_MAX_COLUMNS, sizeof(CImportFieldRef));
+    if (!field_buf) {
+        atomic_store(&ctx->error_code, 1);
+        return NULL;
+    }
     const char *ptr = task->start;
@@ -805,5 +814,6 @@
         ptr = next_row;
     }
 
+    free(field_buf);
     return NULL;
 }
@@ -1131,5 +1141,9 @@ (cimport_parse_csv, single chunk)
         cimport_init_col_stats(chunk->col_stats, ctx->num_columns);
 
-        CImportFieldRef field_buf[CTOOLS_MAX_COLUMNS];
+        CImportFieldRef *field_buf = ctools_safe_malloc2(CTOOLS_MAX_COLUMNS, sizeof(CImportFieldRef));
+        if (!field_buf) {
+            cimport_free_context(ctx);
+            return NULL;
+        }
         const char *ptr = ctx->file_data;
   (in each of the four `cimport_free_context(ctx); return NULL;` exits inside the row loop, add `free(field_buf);` first)
@@ -1205,3 +1219,4 @@
             ptr = next_row;
         }
+        free(field_buf);
     } else {
@@ -1645,3 +1660,3 @@ static ST_retcode cimport_do_scan(
-    char macro_val[65536];
+    static char macro_val[65536];  /* off the Stata thread's stack; plugin calls are serial */
@@ -1742,6 +1757,9 @@
     /* Transport original headings independently, including quotes and spaces. */
-    CImportFieldRef headings[CTOOLS_MAX_COLUMNS];
+    CImportFieldRef *headings = NULL;
     int count = 0;
     if (ctx->has_header) {
+        headings = ctools_safe_malloc2(CTOOLS_MAX_COLUMNS, sizeof(CImportFieldRef));
+        if (!headings) return 920;
         count = cimport_parse_header_row(ctx, headings, NULL, NULL);
     }
@@ -1753,4 +1771,5 @@
         SF_macro_save(key, value);
     }
+    free(headings);
```

Auditor's original suggestion: Heap-allocate these buffers, as in the P0 fix. In cimport_parse_header, allocate one `CImportFieldRef *` buffer of 2*CTOOLS_MAX_COLUMNS entries for fields and next_fields, and a `bool *candidates = calloc(ctx->num_columns, 1)`. Free them on every return path; the function's two returns are easy to cover. In cimport_parse_csv, malloc field_buf before the single-chunk loop and free it after (or on each early `cimport_free_context` return). In cimport_do_scan, malloc headings sized ctx->num_columns and macro_val. All frames then drop to a few KB.

## F027 [P2, confirmed] Windows iconv shim reports EINVAL for odd-length UTF-16 without converting the complete prefix, so cimport shifts the whole file by one byte and imports garbage

Location: src/io/cio_iconv_win.h:96-103 and 153-156 cio_iconv (all-or-nothing consumption); caller src/cimport/cimport_encoding.c:619-626 cimport_convert_named_to_utf8 (skip-one-byte retry); reached from cimport_parse_csv ~956-1003

```diff
--- a/src/io/cio_iconv_win.h
+++ b/src/io/cio_iconv_win.h
@@ -95,15 +95,18 @@ static inline size_t cio_iconv(iconv_t c, char **in, size_t *inleft, char **out, size_t *outleft)
     int n;
+    /* Like iconv, convert every complete unit and leave an incomplete
+       trailing UTF-16 byte in *inleft with EINVAL. */
+    size_t usable = *inleft;
+    int partial = 0;
     if (c->source == 1200 || c->source == 1201)
     {
-        if (*inleft % 2)
-        {
-            errno = EINVAL;
-            return (size_t)-1;
-        }
-        n = (int)(*inleft / 2);
+        partial = (int)(usable % 2);
+        usable -= (size_t)partial;
+        n = (int)(usable / 2);
     }
     else
-        n = MultiByteToWideChar(c->source, 0, *in, (int)*inleft, NULL, 0);
-    if (!n && *inleft)
+        n = MultiByteToWideChar(c->source, 0, *in, (int)usable, NULL, 0);
+    if (!n && usable)
     {
@@ -126,8 +129,9 @@
     else
-        MultiByteToWideChar(c->source, 0, *in, (int)*inleft, wide, n);
+        MultiByteToWideChar(c->source, 0, *in, (int)usable, wide, n);
     int bytes = c->destination == 1200 || c->destination == 1201
                     ? n * 2
                     : WideCharToMultiByte(c->destination, 0, wide, n, NULL, 0, NULL, NULL);
-    if (bytes < 0 || (size_t)bytes > *outleft)
+    /* A failed size query (e.g. output over INT_MAX) is E2BIG, not EILSEQ,
+       so callers that skip bad bytes do not retry the whole buffer. */
+    if (bytes < 0 || (!bytes && n) || (size_t)bytes > *outleft)
     {
@@ -152,7 +156,11 @@
     free(wide);
-    *in += *inleft;
-    *inleft = 0;
+    *in += usable;
+    *inleft -= usable;
     *out += bytes;
     *outleft -= (size_t)bytes;
+    if (partial)
+    {
+        errno = EINVAL;
+        return (size_t)-1;
+    }
     return 0;
 }
```

Auditor's original suggestion: Match iconv's partial-progress semantics. Convert the even prefix and report EINVAL for the leftover byte: size_t usable = *inleft; int partial = 0; if (c->source == 1200 \|\| c->source == 1201) { partial = (int)(usable % 2); usable -= (size_t)partial; n = (int)(usable / 2); } else n = MultiByteToWideChar(c->source, 0, *in, (int)usable, NULL, 0); ... (use usable in place of *inleft) ... *in += usable; *inleft -= usable; *out += bytes; *outleft -= (size_t)bytes; if (partial) { errno = EINVAL; return (size_t)-1; } return 0;

## F029 [P2, confirmed] Windows batched async writes return on a synchronous WriteFile error while earlier writes are still in flight; the caller then frees their OVERLAPPED blocks, events and buffers (kernel use-after-free)

Location: src/cexport/cexport_io.c: cexport_io_write_batch, Windows branch (~409-442: synchronous WriteFile failure return at ~421-426, WAIT_FAILED return at ~437-442). The memory is freed by the caller export_pwrite (src/cexport/cexport_impl.c ~662-676) and by win_cleanup_batch (~121-136) through cexport_io_close (~496).

```diff
--- a/src/cexport/cexport_io.c
+++ b/src/cexport/cexport_io.c
@@ -136,6 +136,21 @@ static void win_cleanup_batch(cexport_io_file *file)
     file->batch_capacity = 0;
 }
 
+/*
+ * Retire batch writes already issued before reporting an error: the caller
+ * frees their source buffers and cexport_io_close frees their OVERLAPPED
+ * blocks and events. Batch writes are issued only from the calling thread,
+ * so CancelIo reaches all of them.
+ */
+static void win_drain_batch(cexport_io_file *file, size_t issued)
+{
+    CancelIo(file->handle);
+    for (size_t j = 0; j < issued; j++) {
+        DWORD ignored = 0;
+        GetOverlappedResult(file->handle, &file->batch_overlapped[j], &ignored, TRUE);
+    }
+}
+
 int cexport_io_open(cexport_io_file *file, const char *filename,
@@ -418,11 +433,15 @@ ssize_t cexport_io_write_batch(cexport_io_file *file,
         BOOL success = WriteFile(file->handle, entries[i].buf,
                                   (DWORD)entries[i].len, &written, ov);
 
-        if (!success && GetLastError() != ERROR_IO_PENDING) {
-            snprintf(file->error_message, sizeof(file->error_message),
-                     "WriteFile batch[%zu] failed: Windows error %lu",
-                     i, GetLastError());
-            return -1;
+        if (!success) {
+            DWORD err = GetLastError();
+            if (err != ERROR_IO_PENDING) {
+                win_drain_batch(file, i);
+                snprintf(file->error_message, sizeof(file->error_message),
+                         "WriteFile batch[%zu] failed: Windows error %lu", i, err);
+                return -1;
+            }
         }
     }
@@ -436,9 +455,11 @@ ssize_t cexport_io_write_batch(cexport_io_file *file,
     if (wait_result == WAIT_FAILED) {
+        DWORD err = GetLastError();
+        win_drain_batch(file, count);
         snprintf(file->error_message, sizeof(file->error_message),
-                 "WaitForMultipleObjects failed: Windows error %lu",
-                 GetLastError());
+                 "WaitForMultipleObjects failed: Windows error %lu", err);
         return -1;
     }
```

Auditor's original suggestion: Before any early return in the async section, retire the writes already issued. Capture `DWORD err = GetLastError();` first, then call `CancelIoEx(file->handle, NULL);` and, for each issued entry j, `DWORD w; GetOverlappedResult(file->handle, &file->batch_overlapped[j], &w, TRUE);` (j < i for the synchronous-failure branch, j < count for WAIT_FAILED). Only then return -1. As defense in depth, have cexport_io_close call CancelIoEx and wait on any batch entries still marked in flight before win_cleanup_batch, and free the wave buffers only after cexport_io_close in export_pwrite's error paths.

## F030 [P2, confirmed] sheet(..., modify) shared-formula anchor relocation is quadratic: every anchor rescans the merged sheet from the start and each relocation copies the whole document

Location: src/cexport/cexport_xlsx_formula.inc: xedit_shared_anchors (87-171; the per-anchor rescan `q = result, limit = q + strlen(q)` at ~110-124, and the full-document xedit_replace copy at ~162). Called from cexport_xlsx_edit.inc:618-626 when sheet_action == 1 (sheetmodify).

```diff
--- a/src/cexport/cexport_xlsx_formula.inc
+++ b/src/cexport/cexport_xlsx_formula.inc
@@ -78,94 +78,140 @@
   return out.data;
 }
-static int xedit_shared_id(const xedit_node *f, const char *id) {
-  const char *a, *b, *v;
-  size_t length;
-  return xedit_attr_span(f, "si", &a, &b, &v, &length) &&
-         strlen(id) == length && !memcmp(id, v, length);
-}
+/* Shared <f> nodes of the merged sheet, indexed once by si then position. */
+typedef struct {
+  char *id;
+  const char *start, *end;
+  int col, row, has_cell, anchor;
+} xedit_shared_ref;
+typedef struct {
+  const char *start, *end;
+  char *text;
+} xedit_shared_edit;
+static int xedit_shared_ref_cmp(const void *x, const void *y) {
+  const xedit_shared_ref *a = x, *b = y;
+  int c = strcmp(a->id, b->id);
+  return c ? c : (a->start > b->start) - (a->start < b->start);
+}
+static int xedit_shared_edit_cmp(const void *x, const void *y) {
+  const xedit_shared_edit *a = x, *b = y;
+  return (a->start > b->start) - (a->start < b->start);
+}
 static char *xedit_shared_anchors(const char *oldxml, const char *merged) {
-  char *result = strdup(merged);
-  if (!result)
-    return NULL;
+  const char *mend = merged + strlen(merged);
+  xedit_shared_ref *refs = NULL;
+  xedit_shared_edit *edits = NULL;
+  size_t nrefs = 0, cap = 0, nedits = 0, ecap = 0, extra = 0;
+  char *result = NULL;
+  xedit_node cell, anchor, f;
+  for (const char *q = merged; xedit_find(q, mend, "c", &cell);) {
+    q = cell.end;
+    if (!xedit_find(cell.open_end, cell.close, "f", &f))
+      continue;
+    char *si = xedit_attr(&f, "si");
+    if (!si)
+      continue;
+    if (nrefs == cap) {
+      size_t grown = cap ? cap * 2 : 64;
+      xedit_shared_ref *r = realloc(refs, grown * sizeof(*refs));
+      if (!r) {
+        free(si);
+        goto done;
+      }
+      refs = r;
+      cap = grown;
+    }
+    xedit_shared_ref *r = &refs[nrefs++];
+    r->id = si;
+    r->start = f.start;
+    r->end = f.end;
+    r->anchor = !f.empty && f.close > f.open_end;
+    r->has_cell = xedit_cell(&cell, &r->col, &r->row);
+  }
+  if (nrefs)
+    qsort(refs, nrefs, sizeof(*refs), xedit_shared_ref_cmp);
   const char *p = oldxml, *end = p + strlen(p);
-  xedit_node cell, anchor;
   while (xedit_find(p, end, "c", &cell)) {
     p = cell.end;
@@
     int col, row;
     if (!xedit_cell(&cell, &col, &row)) {
       free(id);
       continue;
     }
-    const char *q = result, *limit = q + strlen(q);
-    xedit_node next, f, first = {0};
-    int nc = 0, nr = 0, surviving_anchor = 0;
-    while (xedit_find(q, limit, "c", &next)) {
-      q = next.end;
-      if (!xedit_find(next.open_end, next.close, "f", &f) ||
-          !xedit_shared_id(&f, id))
-        continue;
-      if (!f.empty && f.close > f.open_end) {
-        surviving_anchor = 1;
-        break;
-      }
-      if (!first.start && xedit_cell(&next, &nc, &nr))
-        first = f;
-    }
+    size_t lo = 0, hi = nrefs, first = nrefs;
+    while (lo < hi) {
+      size_t mid = lo + (hi - lo) / 2;
+      if (strcmp(refs[mid].id, id) < 0)
+        lo = mid + 1;
+      else
+        hi = mid;
+    }
+    int surviving_anchor = 0;
+    for (size_t g = lo; g < nrefs && !strcmp(refs[g].id, id); g++) {
+      if (refs[g].anchor)
+        surviving_anchor = 1;
+      else if (first == nrefs && refs[g].has_cell)
+        first = g;
+    }
     free(id);
-    if (surviving_anchor || !first.start)
+    if (surviving_anchor || first == nrefs)
       continue;
+    int nc = refs[first].col, nr = refs[first].row;
     size_t length = (size_t)(anchor.close - anchor.open_end);
     char *formula = malloc(length + 1);
     if (!formula)
-      goto failed;
+      goto done;
@@  (remaining `goto failed` in the node-building code become `goto done`)
     strbuf_append_n(&node, anchor.close, (size_t)(anchor.end - anchor.close));
-    char *updated = node.data ? xedit_replace(result, &first, node.data) : NULL;
-    free(node.data);
-    if (!updated)
-      goto failed;
-    free(result);
-    result = updated;
+    if (!node.data)
+      goto done;
+    if (nedits == ecap) {
+      size_t grown = ecap ? ecap * 2 : 16;
+      xedit_shared_edit *e = realloc(edits, grown * sizeof(*edits));
+      if (!e) {
+        free(node.data);
+        goto done;
+      }
+      edits = e;
+      ecap = grown;
+    }
+    edits[nedits++] = (xedit_shared_edit){refs[first].start, refs[first].end, node.data};
+    extra += node.len;
+    refs[first].anchor = 1; /* later anchors with this si now see a survivor */
   }
-  return result;
-failed:
-  free(result);
-  return NULL;
+  if (nedits)
+    qsort(edits, nedits, sizeof(*edits), xedit_shared_edit_cmp);
+  StringBuffer out;
+  if (!strbuf_init(&out, (size_t)(mend - merged) + extra + 1))
+    goto done;
+  const char *keep = merged;
+  for (size_t e = 0; e < nedits; e++) {
+    strbuf_append_n(&out, keep, (size_t)(edits[e].start - keep));
+    strbuf_append(&out, edits[e].text);
+    keep = edits[e].end;
+  }
+  strbuf_append(&out, keep);
+  result = out.data; /* NULL if a grow failed; caller maps NULL to r(920) */
+done:
+  for (size_t i = 0; i < nrefs; i++)
+    free(refs[i].id);
+  free(refs);
+  for (size_t e = 0; e < nedits; e++)
+    free(edits[e].text);
+  free(edits);
+  return result;
 }
```

Auditor's original suggestion: Make two linear passes. First, scan `result` once and build a small table keyed by si: whether a non-empty anchor survives, and the first dependent <f> node with its cell coordinates. Next, scan the old sheet once, look up each shared anchor's si in that table, and collect the needed replacements (position, new node text). Finally, rebuild the document once by copying the spans between the sorted replacement positions. This makes the work O(document + F log F) and removes the per-anchor full copies.

## F031 [P2, confirmed] CSV row formatters return the row length as int; rows over INT_MAX bytes are misreported, and rows of 4 GiB or more wrap to a small positive length (silent truncation and overlap)

Location: src/cexport/cexport_format.c: cexport_format_row (371-430: `return (int)pos` at ~429) and cexport_write_quoted_string (311-330: `return (int)pos`). Consumers in src/cexport/cexport_impl.c: write_rows_buffered (~268-287, `buf_pos += row_len`), format_chunk_thread (~184-203), format_mmap_thread (~221-228). cexport_format_row_numeric (<= nvars*~34 bytes) and cexport_format_row_stream (strings capped at 2045 bytes, <= ~0.5 GB per row) cannot exceed INT_MAX, so the defect is latent in those two.

```diff
--- a/src/cexport/cexport_format.c
+++ b/src/cexport/cexport_format.c
@@ -9,6 +9,7 @@
 #include <stdlib.h>
 #include <string.h>
 #include <stdio.h>
+#include <limits.h>
 #include <stdbool.h>
@@ -313,6 +314,7 @@ static int cexport_write_quoted_string(const char *str, char *buf, size_t buf_size)
     size_t need = 3; /* outer quotes and terminator */
     for (const char *p = str; *p; p++) {
         size_t add = *p == '"' ? 2 : 1;
+        if (need + add > (size_t)INT_MAX) return -3; /* length must fit int */
         if (need > buf_size || add > buf_size - need) return -1;
         need += add;
     }
@@ -407,6 +409,7 @@ int cexport_format_row(const cexport_context *ctx, size_t row_idx,
             if (need_quote) {
                 field_len = cexport_write_quoted_string(str, buf + pos, buf_size - pos);
+                if (field_len == -3) return -3;
                 if (field_len < 0 || pos + (size_t)field_len + 2 + (size_t)line_end_len > buf_size) return -1;
             } else {
@@ -424,6 +427,8 @@ int cexport_format_row(const cexport_context *ctx, size_t row_idx,
             memcpy(buf + pos, line_end, line_end_len);
             pos += line_end_len;
         }
+        /* -3: the row is too long to report through the int return */
+        if (pos > (size_t)INT_MAX) return -3;
     }
 
     return (int)pos;
--- a/src/cexport/cexport_impl.c
+++ b/src/cexport/cexport_impl.c
@@ -196,7 +196,7 @@ static void *format_chunk_thread(void *arg)
         if (row_len < 0) {
-            args->read_error = 1;
+            args->read_error = row_len == -2;   /* -3: row over 2 GiB */
             args->success = 0;
             args->bytes_written = pos;
             return NULL;
@@ -281,7 +281,7 @@ static int write_rows_buffered(FILE *fp, size_t nobs, size_t avg_row_size, size_t *rows_written_out)
         if (row_len < 0) {
             free(buffer);
-            return -1;
+            return row_len == -3 ? -3 : -1;
         }
@@ -333,8 +333,10 @@
     size_t rows_written = 0;
-    if (write_rows_buffered(g_ctx.fp, nobs, avg_row_size, &rows_written) != 0) {
-        SF_error("cexport: write error during buffered write\n");
+    int wr = write_rows_buffered(g_ctx.fp, nobs, avg_row_size, &rows_written);
+    if (wr != 0) {
+        SF_error(wr == -3 ? "cexport: an observation exceeds 2 GiB of delimited text\n"
+                          : "cexport: write error during buffered write\n");
         fclose(g_ctx.fp);
         g_ctx.fp = NULL;
         return 693;
```

Auditor's original suggestion: Before returning, reject rows that do not fit the return type: `if (pos > INT_MAX) return -3;` in all three formatters, documented as 'row too large'. format_chunk_thread's `== -1` growth loop then stops immediately. In format_chunk_thread and write_rows_buffered, report -3 as "row exceeds 2 GiB" instead of a read/write error. Alternatively, change the formatter return type to ptrdiff_t throughout.

## F032 [P2, confirmed] xedit_styles ignores a failed re-find of cellXfs and passes a node from a freed-later, different buffer to xedit_replace (cross-object pointer subtraction)

Location: src/cexport/cexport_xlsx_edit.inc: xedit_styles (362-447). The second `xedit_find(result, ..., "cellXfs", &xfs)` at ~424 has an unchecked return, and the stale xfs then reaches xedit_set_attr and xedit_replace (~425-430, with the pointer subtraction in xedit_replace at ~173). It is reached from xlsx_edit_file:580 on every export into an existing workbook (no `replace`).

```diff
--- a/src/cexport/cexport_xlsx_edit.inc
+++ b/src/cexport/cexport_xlsx_edit.inc
@@ -421,7 +421,10 @@ static char *xedit_styles(char *xml, unsigned *base) {
   if (!result)
     return NULL;
-  xedit_find(result, result + strlen(result), "cellXfs", &xfs);
+  if (!xedit_find(result, result + strlen(result), "cellXfs", &xfs)) {
+    free(result); /* never apply a node from xml to result */
+    return NULL;
+  }
   char *node = xedit_set_attr(&xfs, "count", *base + 3);
```

Auditor's original suggestion: `if (!xedit_find(result, result + strlen(result), "cellXfs", &xfs)) { free(result); return NULL; }` The caller already maps NULL to r(920).

## F034 [P2, confirmed] NULL libxls sheet names are dereferenced: crash in `cimport excel .xls, describe`, with `sheet()`, and when exporting into an existing .xls

Location: src/io/cio_excel.c:62 (cio_excel_workbook) -> cio_excel_sheet_info:16 -> src/io/cio_context.c:52 (cio_save); src/io/cio_xls.c:157 (cio_scan_xls); src/io/cexport_xls_edit.inc:390 (xe_edit)

```diff
--- a/src/io/cio_context.c
+++ b/src/io/cio_context.c
@@ -49,5 +49,5 @@ int cio_save(const char *key, const char *value)
     else
         cio_copy_text(name, sizeof(name), key);
-    return SF_macro_save(name, (char *)value);
+    return SF_macro_save(name, (char *)(value ? value : ""));
 }
--- a/src/io/cio_xls.c
+++ b/src/io/cio_xls.c
@@ -155,5 +155,6 @@ int cio_scan_xls(cio_dataset *dataset, const char *path)
         which = -1;
         for (unsigned j = 0; j < w->sheets.count; j++)
-            if (!strcmp(sheet, w->sheets.sheet[j].name))
+            if (w->sheets.sheet[j].name && !strcmp(sheet, w->sheets.sheet[j].name))
             {
--- a/src/io/cexport_xls_edit.inc
+++ b/src/io/cexport_xls_edit.inc
@@ -389,3 +389,4 @@ static int xe_edit(const char *path, const char *sheet, int action, int keep,
   for (unsigned i = 0; i < wb->sheets.count; i++)
-    if (!*sheet ? i == 0 : !strcasecmp(sheet, wb->sheets.sheet[i].name))
+    if (!*sheet ? i == 0
+                : wb->sheets.sheet[i].name && !strcasecmp(sheet, wb->sheets.sheet[i].name))
       which = (int)i;
```

Auditor's original suggestion: Treat a NULL name as "" everywhere it is read. Add a helper to cio_internal.h, `static inline const char *cio_xls_sheet_name(xlsWorkBook *w, unsigned j) { const char *s = w->sheets.sheet[j].name; return s ? s : ""; }`, and use it at cio_excel.c:62, cio_xls.c:157 and cexport_xls_edit.inc:390. Also harden `cio_save`: `return SF_macro_save(name, (char *)(value ? value : ""));`. Optionally add the Mac code pages (10004-10082) and BIG5HKSCS (951) to `cio_codepage()` so that BIFF5 strings decode on Windows.

## F036 [P2, confirmed] cexport dbase converts out-of-range %td values with (int)x (UB), then overflows cio_days_civil

Location: src/io/cio_dbf.c:413-418 (cio_export_dbf, the date-field branch)

```diff
--- a/src/io/cio_dbf.c
+++ b/src/io/cio_dbf.c
@@ -413,7 +413,14 @@ int cio_export_dbf(void)
                 else if (kinds[j] == 1)
                 {
                     int year;
                     unsigned month, day;
-                    cio_days_civil((int)x, &year, &month, &day);
+                    double d = floor(x); /* %td displays the floor of fractional days */
+                    /* The 8-byte D field holds only 0000-01-01 .. 9999-12-31 */
+                    if (d < -715875.0 || d > 2936549.0)
+                    {
+                        rc = 108;
+                        goto done;
+                    }
+                    cio_days_civil((int)d, &year, &month, &day);
                     snprintf(value, sizeof(value), "%04d%02u%02u", year, month, day);
                 }
```

Auditor's original suggestion: `double d = floor(x); if (d < -2000000000.0 \|\| d > 2000000000.0) { rc = 108; goto done; } cio_days_civil((int)d, &year, &month, &day);`. floor() also fixes the off-by-one for negative fractional dates.

## F039 [P2, confirmed] xe_globals/xe_sheet read an uninitialized xe_record when the first record cannot be parsed

Location: src/io/cexport_xls_edit.inc:86 (xe_globals, `xe_record r, nr` uninitialized, read at ~152 `if (r.id != 0x0a)`) and :177 (xe_sheet, `xe_record r;`, read at ~229)

```diff
--- a/src/io/cexport_xls_edit.inc
+++ b/src/io/cexport_xls_edit.inc
@@ -85,3 +85,3 @@ static int xe_globals(FILE *f, const unsigned char *old, size_t oldlen,
   size_t p = 0, np = 0;
-  xe_record r, nr, sst = {0}, bound = {0};
+  xe_record r = {0}, nr = {0}, sst = {0}, bound = {0};
   *strings = 0;
@@ -176,3 +176,3 @@ static int xe_sheet(FILE *f, const unsigned char *old, size_t oldlen,
   size_t p = start;
-  xe_record r;
+  xe_record r = {0};
   int rc, skip_string = 0;
```

Auditor's original suggestion: Zero-initialize: `xe_record r = {0}, nr = {0}, ...` in xe_globals and `xe_record r = {0};` in xe_sheet, so a failed first read returns 610.

## F041 [P2, confirmed] Row-block count overflows for N within 2047 of INT_MAX, so ctools_matvec_add/ctools_matmul_xb skip every row and the results come out silently wrong

Location: src/ctools_matrix.c:81 matmul_xb_blocks (callers ctools_matvec_add / ctools_matmul_xb); siblings with the same pattern: src/civreghdfe/civreghdfe_impl.c:520, src/civreghdfe/civreghdfe_tests.c:77, src/cbinscatter/cbinscatter_impl.c:351, src/cbinscatter/cbinscatter_resid.c:302

```diff
--- a/src/ctools_matrix.c
+++ b/src/ctools_matrix.c
@@ -78,7 +78,8 @@
 static void matmul_xb_blocks(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                              ST_double *C, int zero_first)
 {
-    const ST_int num_blocks = (N + CTOOLS_ROW_BLOCK - 1) / CTOOLS_ROW_BLOCK;
+    /* N + CTOOLS_ROW_BLOCK - 1 overflows int when N is within a block of INT_MAX */
+    const ST_int num_blocks = N / CTOOLS_ROW_BLOCK + (N % CTOOLS_ROW_BLOCK != 0);
     ST_int blk;
 
     #pragma omp parallel for schedule(static) if((size_t)N * (size_t)K * (size_t)M > 400000)
--- a/src/civreghdfe/civreghdfe_impl.c
+++ b/src/civreghdfe/civreghdfe_impl.c
@@ -519,3 +519,3 @@
         const ST_int block = 16384;
-        const ST_int num_blocks = (N_total + block - 1) / block;
+        const ST_int num_blocks = N_total / block + (N_total % block != 0);
         ST_int blk;
--- a/src/civreghdfe/civreghdfe_tests.c
+++ b/src/civreghdfe/civreghdfe_tests.c
@@ -76,3 +76,3 @@
         const ST_int block = 2048;
-        const ST_int num_blocks = (N + block - 1) / block;
+        const ST_int num_blocks = N / block + (N % block != 0);
         ST_int blk;
--- a/src/cbinscatter/cbinscatter_impl.c
+++ b/src/cbinscatter/cbinscatter_impl.c
@@ -350,3 +350,3 @@
         const ST_int block = 16384;
-        const ST_int num_blocks = (N + block - 1) / block;
+        const ST_int num_blocks = N / block + (N % block != 0);
         ST_int blk;
--- a/src/cbinscatter/cbinscatter_resid.c
+++ b/src/cbinscatter/cbinscatter_resid.c
@@ -301,3 +301,3 @@
         enum { BLOCK = 2048 };
-        const ST_int num_blocks = (N + BLOCK - 1) / BLOCK;
+        const ST_int num_blocks = N / BLOCK + (N % BLOCK != 0);
         ST_int blk;
```

Auditor's original suggestion: Compute the block count without overflowing: `const ST_int num_blocks = N / CTOOLS_ROW_BLOCK + (N % CTOOLS_ROW_BLOCK != 0);` (or do the addition in size_t/long long). The same pattern appears in civreghdfe_impl.c:520 (`(N_total + block - 1) / block` with block = 16384). There it leaves valid_mask uninitialized and N_valid = 0, so that site should be fixed the same way.

## F042 [P2, confirmed] K x K sizes and indices are computed in int in the shared OLS/matrix kernels; with 46,341 or more regressors this overflows and corrupts the heap

Location: src/ctools_matrix.c:35,44,61,66 (ctools_matmul_atb/ata), :163,173,184 (ctools_matmul_atdb), :264-265,323-355 (ctools_vce_robust/cluster meat/temp); src/ctools_ols.c:302-330,364-413 (compute_xtx_xty), :437-479 (detect_collinearity), :509-522 (detect_collinearity_scaled). Reached from src/civreghdfe/civreghdfe_impl.c:1290-1292.

```diff
--- a/src/civreghdfe/civreghdfe_impl.c
+++ b/src/civreghdfe/civreghdfe_impl.c
@@ -370,6 +370,12 @@
     if (center) {
         SF_error("civreghdfe: center is not supported\n");
         return 198;
     }
+    /* K x K workspaces (shared kernels and the IV code) are indexed in int:
+       j*K + i overflows once K exceeds 46340 */
+    if (K_iv > 46340 || K_endog + K_exog > 46340) {
+        SF_error("civreghdfe: too many regressors or instruments (limit 46340)\n");
+        return 198;
+    }
 
     /* Calculate variable positions */
--- a/src/ctools_matrix.c
+++ b/src/ctools_matrix.c
@@ -35 +35 @@
-    memset(C, 0, K1 * K2 * sizeof(ST_double));
+    memset(C, 0, (size_t)K1 * K2 * sizeof(ST_double));
@@ -44 +44 @@
-            C[j * K1 + i] = ctools_dot_unrolled(a_col, b_col, N);
+            C[(size_t)j * K1 + i] = ctools_dot_unrolled(a_col, b_col, N);
@@ -61 +61 @@
-            C[j * K + i] = ctools_dot_unrolled(A + (size_t)i * N, a_j, N);
+            C[(size_t)j * K + i] = ctools_dot_unrolled(A + (size_t)i * N, a_j, N);
@@ -66 +66 @@
-            C[i * K + j] = C[j * K + i];
+            C[(size_t)i * K + j] = C[(size_t)j * K + i];
@@ -163 +163 @@
-    memset(C, 0, K1 * K2 * sizeof(ST_double));
+    memset(C, 0, (size_t)K1 * K2 * sizeof(ST_double));
@@ -173 +173 @@
-                C[j * K1 + i] = ctools_simd_dot_weighted(w, a_col, b_col, (size_t)N);
+                C[(size_t)j * K1 + i] = ctools_simd_dot_weighted(w, a_col, b_col, (size_t)N);
@@ -184 +184 @@
-                C[j * K1 + i] = ctools_simd_dot(a_col, b_col, (size_t)N);
+                C[(size_t)j * K1 + i] = ctools_simd_dot(a_col, b_col, (size_t)N);
@@ -264,2 +264,2 @@
-    long double *meat = (long double *)calloc(K * K, sizeof(long double));
-    long double *temp = (long double *)malloc(K * K * sizeof(long double));
+    long double *meat = (long double *)calloc((size_t)K * K, sizeof(long double));
+    long double *temp = (long double *)malloc((size_t)K * K * sizeof(long double));
--- a/src/ctools_ols.c
+++ b/src/ctools_ols.c
@@ -437,4 +437,4 @@
-    L = (ST_double *)malloc(K * K * sizeof(ST_double));
+    L = (ST_double *)ctools_safe_malloc3((size_t)K, (size_t)K, sizeof(ST_double));
     if (!L) return -1;
 
-    memcpy(L, xx, K * K * sizeof(ST_double));
+    memcpy(L, xx, (size_t)K * K * sizeof(ST_double));
```

Auditor's original suggestion: Index and size the K x K arrays in size_t, for example `C[(size_t)j * K + i]`, `memset(C, 0, (size_t)K1 * K2 * sizeof(ST_double))`, `calloc((size_t)K * K, ...)`, and `ctools_safe_malloc3((size_t)K, (size_t)K, sizeof(ST_double))` in detect_collinearity. Alternatively, reject K > 46340 with r(198)/r(920) at the civreghdfe/cpplmhdfe entry points, as creghdfe already does with its K > 1000 cap.

## F043 [P2, confirmed] cbsample with strata() and cluster(): finding each stratum's clusters takes time proportional to strata x clusters (near-hang on large data)

Location: src/cbsample/cbsample_impl.c:432-455 (cbsample_main, cluster-bootstrap stratum loop)

```diff
--- a/src/cbsample/cbsample_impl.c
+++ b/src/cbsample/cbsample_impl.c
@@ -429,30 +429,22 @@
         /* Cluster bootstrap: resample clusters with replacement within each stratum */
         /* For each stratum, identify which clusters are in it and sample from those */
 
+        /* Clusters nest in strata and both are in sort order, so one cursor
+         * walks the clusters once across all strata */
+        size_t cl_cursor = 0;
         for (size_t st = 0; st < nstrata_groups; st++) {
             size_t st_start = strata[st].start;
             size_t st_count = strata[st].count;
+            size_t st_end = st_start + st_count;
 
             /* Find clusters within this stratum */
-            /* Clusters are sorted by (strata, cluster), so find range */
-            size_t first_cluster = 0;
-            size_t last_cluster = nclusters;
-
-            /* Binary search for first cluster in stratum */
-            for (size_t cl = 0; cl < nclusters; cl++) {
-                if (clusters[cl].start >= st_start &&
-                    clusters[cl].start < st_start + st_count) {
-                    first_cluster = cl;
-                    break;
-                }
-            }
-            /* Find last cluster in stratum */
-            for (size_t cl = first_cluster; cl < nclusters; cl++) {
-                if (clusters[cl].start >= st_start + st_count) {
-                    last_cluster = cl;
-                    break;
-                }
-            }
+            while (cl_cursor < nclusters && clusters[cl_cursor].start < st_start) cl_cursor++;
+            size_t first_cluster = cl_cursor;
+            while (cl_cursor < nclusters && clusters[cl_cursor].start < st_end) cl_cursor++;
+            size_t last_cluster = cl_cursor;
 
             size_t n_clusters_in_stratum = last_cluster - first_cluster;
             if (n_clusters_in_stratum == 0) continue;
```

Auditor's original suggestion: Keep a running cursor across strata, since both strata and clusters are in sorted order and clusters nest in strata. Declare `size_t cl_cursor = 0;` before the `for (st ...)` loop. Then for each stratum: `size_t st_end = st_start + st_count; while (cl_cursor < nclusters && clusters[cl_cursor].start < st_start) cl_cursor++; size_t first_cluster = cl_cursor; while (cl_cursor < nclusters && clusters[cl_cursor].start < st_end) cl_cursor++; size_t last_cluster = cl_cursor;`. Total work becomes O(nclusters + nstrata) and the draws stay identical (same first/last indices, same RNG sequence).

## F046 [P2, confirmed] tsset time index truncated to int: spans of 2^31 or more delta units wrap tidx, which mis-pairs kernel lags, drops DK periods, and can make pair_next read past om->order

Location: src/civreghdfe/civreghdfe_impl.c:1393-1397 (do_iv_regression, tidx = (ST_int)llround((tv-tmin)/tdelta) with no span check). Consumers: civreghdfe_omega.c:201-206 (pair_next), 233-236 and 359/368 (omega_groups with G = tmax_idx+1), 42-49 (omega_max_lag), estimate.c:734 (kiefer s_bw = tmax_idx).

```diff
--- a/src/civreghdfe/civreghdfe_impl.c
+++ b/src/civreghdfe/civreghdfe_impl.c
@@ -1391,9 +1391,16 @@
                 if (tv[i] > tmax) tmax = tv[i];
             }
+            /* The time index is an int: refuse spans of 2^31-2 or more delta
+               units (e.g. a %tc clock with delta 1 ms over more than 24.8
+               days), which would wrap and break the kernel lag pairing */
+            if (!((tmax - tmin) / tdelta < 2147483646.0)) {
+                SF_error("civreghdfe: tsset time span too long for kernel estimators; tsset with a coarser delta()\n");
+                rc = 198;
+            }
             ST_int tmax_idx = 0;
-            for (ST_int i = 0; i < N; i++) {
+            for (ST_int i = 0; i < N && !rc; i++) {
                 tidx[i] = (ST_int)llround((tv[i] - tmin) / tdelta);
                 if (tidx[i] > tmax_idx) tmax_idx = tidx[i];
             }
--- a/src/civreghdfe/civreghdfe_omega.c
+++ b/src/civreghdfe/civreghdfe_omega.c
@@ -46,5 +46,5 @@ static ST_int omega_max_lag(ST_int kernel, ST_int bw, ST_int tmax_idx, ST_double tdelta)
         if (tdelta <= 0.0) tdelta = 1.0;
         ST_double tau = floor((ST_double)tmax_idx - 1.0 + 1.0 / tdelta + 1e-9);
-        return tau > 0.0 ? (ST_int)tau : 0;
+        return tau > 0.0 ? (tau < 2147483646.0 ? (ST_int)tau : 2147483646) : 0;
     }
```

Auditor's original suggestion: In impl.c, after the tmin/tmax loop and before computing tidx, reject spans that do not fit: if (!((tmax - tmin) / tdelta < (ST_double)INT_MAX - 1.0)) { SF_error("civreghdfe: the tsset time span is too long for kernel estimators (more than 2^31-2 periods of delta)\n"); rc = 198; } else { /* existing tidx loop */ } This keeps every tidx in [0, INT_MAX-2], so tmax_idx + 1 and t - tau cannot overflow. In omega_max_lag, clamp tau in double before the (ST_int) cast. Optionally, bound `a` by `b` in pair_next as a defensive check.

## F047 [P2, confirmed] Unchecked small allocations in GMM2S and the k-class path (NULL dereference); a failed CUE or GMM Hessian allocation silently falls back

Location: src/civreghdfe/civreghdfe_estimate.c:82-122 (ivest_compute_gmm2s: WZtX, XZW, XZWZX, XZWZy, XZWZX_L, beta_temp unchecked; 161-175 hessian inverse) and estimate.c:1014-1016 (ivest_compute_2sls k-class XtX/Xty), 1165 and 1219 (gmm_hessian_inv), 1187-1190 and 1251-1253 (silent 2SLS fallback).

```diff
--- a/src/civreghdfe/civreghdfe_estimate.c
+++ b/src/civreghdfe/civreghdfe_estimate.c
@@ -79,10 +79,18 @@ ST_retcode ivest_compute_gmm2s(
     /* Compute GMM estimator: β = (X'ZWZ'X)^-1 X'ZWZ'y */
-    /* 1. WZtX = W * Z'X */
-    ST_double *WZtX = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
-    ctools_matmul_ab(ZOmegaZ_inv, ctx->ZtX, K_iv, K_iv, K_total, WZtX);
-
-    /* 2. XZW = (Z'X)' * W (K_total x K_iv) */
-    ST_double *XZW = (ST_double *)calloc(K_total * K_iv, sizeof(ST_double));
+    ST_double *XZW = (ST_double *)calloc((size_t)K_total * K_iv, sizeof(ST_double));
+    ST_double *XZWZX = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
+    ST_double *XZWZy = (ST_double *)calloc((size_t)K_total, sizeof(ST_double));
+    ST_double *XZWZX_L = (ST_double *)malloc((size_t)K_total * K_total * sizeof(ST_double));
+    ST_double *beta_temp = (ST_double *)calloc((size_t)K_total, sizeof(ST_double));
+    if (!XZW || !XZWZX || !XZWZy || !XZWZX_L || !beta_temp) {
+        free(ZOmegaZ); free(ZOmegaZ_inv);
+        free(XZW); free(XZWZX); free(XZWZy); free(XZWZX_L); free(beta_temp);
+        return 920;
+    }
+
+    /* 1. XZW = (Z'X)' * W (K_total x K_iv) */
     for (i = 0; i < K_total; i++) {
@@ -96,3 +104,2 @@
-    /* 3. XZWZX = XZW * Z'X (K_total x K_total) */
-    ST_double *XZWZX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
+    /* 2. XZWZX = XZW * Z'X (K_total x K_total) */
     for (i = 0; i < K_total; i++) {
@@ -108,3 +115,2 @@
-    /* 4. XZWZy = XZW * Z'y (K_total x 1) */
-    ST_double *XZWZy = (ST_double *)calloc(K_total, sizeof(ST_double));
+    /* 3. XZWZy = XZW * Z'y (K_total x 1) */
     for (i = 0; i < K_total; i++) {
@@ -119,5 +125,3 @@
-    /* 5. Solve XZWZX * beta = XZWZy */
-    ST_double *XZWZX_L = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
-    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));
+    /* 4. Solve XZWZX * beta = XZWZy */
     memcpy(XZWZX_L, XZWZX, K_total * K_total * sizeof(ST_double));
@@ -126,3 +130,3 @@
         free(ZOmegaZ); free(ZOmegaZ_inv);
-        free(WZtX); free(XZW); free(XZWZX); free(XZWZy);
+        free(XZW); free(XZWZX); free(XZWZy);
         free(XZWZX_L); free(beta_temp);
@@ -160,22 +164,13 @@
-    /* Output the GMM Hessian inverse (X'ZWZ'X)^-1 if requested */
-    if (XZWZX_inv_out) {
-        /* XZWZX_L contains the Cholesky factor, compute the inverse */
-        ST_double *XZWZX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
-        if (XZWZX_inv) {
-            /* Re-compute Cholesky and inverse from XZWZX */
-            memcpy(XZWZX_inv, XZWZX, K_total * K_total * sizeof(ST_double));
-            if (ctools_cholesky(XZWZX_inv, K_total) == 0) {
-                ctools_invert_from_cholesky(XZWZX_inv, K_total, XZWZX_inv);
-                memcpy(XZWZX_inv_out, XZWZX_inv, K_total * K_total * sizeof(ST_double));
-            }
-            free(XZWZX_inv);
-        }
-    }
+    /* Output the GMM Hessian inverse (X'ZWZ'X)^-1: XZWZX_L holds its
+       Cholesky factor (the same factor the old re-factorization produced) */
+    ST_retcode out_rc = STATA_OK;
+    if (XZWZX_inv_out && ctools_invert_from_cholesky(XZWZX_L, K_total, XZWZX_inv_out) != 0)
+        out_rc = 920;
 
     free(ZOmegaZ); free(ZOmegaZ_inv);
-    free(WZtX); free(XZW); free(XZWZX); free(XZWZy);
+    free(XZW); free(XZWZX); free(XZWZy);
     free(XZWZX_L); free(beta_temp);
 
-    return STATA_OK;
+    return out_rc;
 }
@@ -1014,4 +1009,9 @@
     if (kclass != 1.0) {
         XtX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
         Xty = (ST_double *)calloc(K_total, sizeof(ST_double));
+        if (!XtX || !Xty) {
+            SF_error("civreghdfe: Memory allocation failed for X'X/X'y\n");
+            est_rc = 920;
+            goto est_cleanup;
+        }
 
@@ -1165,2 +1165,7 @@
         gmm_hessian_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
+        if (!gmm_hessian_inv) {
+            SF_error("civreghdfe: Memory allocation failed for the GMM Hessian\n");
+            est_rc = 920;
+            goto est_cleanup;
+        }
@@ -1187,4 +1192,5 @@
         if (gmm_rc != STATA_OK) {
-            if (verbose) SF_display("civreghdfe: GMM2S re-estimation failed, using 2SLS\n");
-            if (gmm_hessian_inv) { free(gmm_hessian_inv); gmm_hessian_inv = NULL; }
+            if (gmm_rc == 920) SF_error("civreghdfe: Memory allocation failed in GMM2S\n");
+            est_rc = gmm_rc;
+            goto est_cleanup;
         } else {
@@ -1219,2 +1225,7 @@
         gmm_hessian_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
+        if (!gmm_hessian_inv) {
+            SF_error("civreghdfe: Memory allocation failed for the CUE Hessian\n");
+            est_rc = 920;
+            goto est_cleanup;
+        }
@@ -1251,4 +1262,5 @@
         if (cue_rc != STATA_OK) {
-            if (verbose) SF_display("civreghdfe: CUE computation failed, using 2SLS\n");
-            if (gmm_hessian_inv) { free(gmm_hessian_inv); gmm_hessian_inv = NULL; }
+            /* beta may already hold a Newton iterate: do not post it as 2SLS */
+            SF_error("civreghdfe: CUE estimation failed (memory allocation)\n");
+            est_rc = cue_rc;
+            goto est_cleanup;
         } else {
```

Auditor's original suggestion: In ivest_compute_gmm2s: - Delete WZtX. - After allocating XZW, XZWZX, XZWZy, XZWZX_L and beta_temp, add `if (!XZW \|\| !XZWZX \|\| !XZWZy \|\| !XZWZX_L \|\| !beta_temp) { free(...all...); return 920; }` before first use. In ivest_compute_2sls: - Add `if (!XtX \|\| !Xty) { SF_error("civreghdfe: Memory allocation failed\n"); est_rc = 920; goto est_cleanup; }` after line 1016. - Add `if (!gmm_hessian_inv) { est_rc = 920; goto est_cleanup; }` after lines 1165 and 1219.

## F050 [P2, confirmed] Diagnostic test functions are still void: allocation or solve failures post 0 statistics (p = 1) or fall back to the iid Anderson LM labelled as Kleibergen-Paap, with no warning

Location: src/civreghdfe/civreghdfe_tests.c:599-606 and 616-664 (civreghdfe_compute_underid_test, Y/Zx and canonical-correlation failures leave 0), tests.c:677-678 (civ_kp_rk failure keeps the iid Anderson LM), 718-721 (civreghdfe_compute_hansen_j), 908-946 (civreghdfe_compute_cstat), 1001-1002 (civreghdfe_compute_endogtest_subset), 1124-1125 (civreghdfe_compute_redundant). The results are posted by estimate.c:1575-1676 and civreghdfe.ado:1399-1700.

```diff
--- a/src/civreghdfe/civreghdfe_tests.c
+++ b/src/civreghdfe/civreghdfe_tests.c
@@ -599,7 +599,11 @@ void civreghdfe_compute_underid_test(
     if (!Y || !Zx ||
         civ_partial_cols(X_endog, K_endog, Z, K_exog, w, N, Y) != 0 ||
         civ_partial_cols(Z + (size_t)K_exog * N, L, Z, K_exog, w, N, Zx) != 0) {
         free(Y); free(Zx);
+        /* not computed: missing, never 0 (p = 1) */
+        *underid_stat = SV_missval;
+        *cd_f = SV_missval;
+        *kp_f = SV_missval;
         return;
     }
@@ -658,7 +662,10 @@
         if (ok) {
             ST_double min_eval = ev[0];
             *underid_stat = (ST_double)N_eff * min_eval;
             if (min_eval > 0.0 && min_eval < 1.0)
                 *cd_f = min_eval / (1.0 - min_eval) * ((ST_double)df_resid / (ST_double)L);
+        } else {
+            *underid_stat = SV_missval;
+            *cd_f = SV_missval;
         }
@@ -676,5 +683,7 @@
         ST_double lm = 0.0, wald = 0.0;
-        if (civ_kp_rk(Y, K_endog, Zx, L, w, N, rk_om, 0, &lm) == 0) *underid_stat = lm;
+        /* never report the iid Anderson LM under the Kleibergen-Paap label */
+        *underid_stat = (civ_kp_rk(Y, K_endog, Zx, L, w, N, rk_om, 0, &lm) == 0) ? lm : SV_missval;
+        *kp_f = SV_missval;
         if (civ_kp_rk(Y, K_endog, Zx, L, w, N, rk_om, 1, &wald) == 0) {
@@ -720,2 +729,2 @@ void civreghdfe_compute_hansen_j(
                       : civ_gmm_eff_j(y, X_all, K_total, Z, K_iv, w, N, S, &j);
-    if (rc == 0) *J = j;
+    *J = (rc == 0) ? j : SV_missval;
 }
@@ -903,6 +912,8 @@ void civreghdfe_compute_cstat(
     if (n_orthog <= 0 || K_rest < K_total) return;
+    *cstat = SV_missval;   /* overwritten on success */
+    if (!(J_full < SV_missval)) return;
 
     ST_int *keep = (ST_int *)malloc((size_t)K_iv * sizeof(ST_int));
@@ -984,4 +995,5 @@ void civreghdfe_compute_endogtest_subset(
     if (n_endogtest <= 0 || n_endogtest > K_endog) return;
+    *endogtest_stat = SV_missval;   /* overwritten on success */
 
@@ -1116,4 +1128,5 @@ void civreghdfe_compute_redundant(
     if (n_redund <= 0 || K_endog <= 0 || K_rest < K_exog) return;
+    *redund_stat = SV_missval;   /* overwritten on success */
```

Auditor's original suggestion: Make each civreghdfe_compute_* return ST_retcode (0, 920, or -1 for a numerical failure), and set the outputs to SV_missval rather than 0 before any early return. In civreghdfe_estimate.c, propagate 920 as an error, or post SV_missval with an SF_display warning for numerical failures. In civreghdfe_compute_underid_test, set *underid_stat = SV_missval and *kp_f = SV_missval when civ_kp_rk fails on the non-iid path, rather than keeping the iid LM.

## F056 [P2, confirmed] cqreg_compute_quantile swallows an allocation failure and returns 0, so e(q_v), e(sum_rdev) and e(r2_p) are silently wrong with rc 0

Location: src/cqreg/cqreg_linalg.c:181-188 cqreg_compute_quantile; caller src/cqreg/cqreg_regress.c:825-826 cqreg_full_regression

```diff
--- a/src/cqreg/cqreg_linalg.h
+++ b/src/cqreg/cqreg_linalg.h
@@ -92,6 +92,7 @@
  * Compute sample quantile by selection (expected O(N)).
  * Creates a copy of y internally - does not modify input.
  * tau: quantile (0 < tau < 1), e.g., 0.5 for median
- * Returns: the tau-th sample quantile
+ * Stores the tau-th sample quantile in *q. Returns 0, or -1 when the
+ * working copy cannot be allocated (*q is then 0 and must not be used)
  */
-ST_double cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau);
+ST_int cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau, ST_double *q);
--- a/src/cqreg/cqreg_linalg.c
+++ b/src/cqreg/cqreg_linalg.c
@@ -181,8 +181,9 @@
-ST_double cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau)
+ST_int cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau, ST_double *q)
 {
-    if (N <= 0) return 0.0;
+    *q = 0.0;
+    if (N <= 0) return 0;
     /* Select from a copy so the response stays in observation order. */
     ST_double *y_work = (ST_double *)malloc((size_t)N * sizeof(ST_double));
     if (y_work == NULL) {
-        return 0.0;  /* Error: return 0 */
+        return -1;
     }
@@ -208,6 +209,6 @@
-    ST_double q = cqreg_select(y_work, N, k - 1);  /* Convert to 0-indexed */
+    *q = cqreg_select(y_work, N, k - 1);  /* Convert to 0-indexed */
 
     free(y_work);
-    return q;
+    return 0;
 }
--- a/src/cqreg/cqreg_regress.c
+++ b/src/cqreg/cqreg_regress.c
@@ -824,3 +824,9 @@
     /* Compute sample quantile and raw sum of deviations for pseudo R^2 */
-    state->q_v = cqreg_compute_quantile(state->y, N, quantile);
+    if (cqreg_compute_quantile(state->y, N, quantile, &state->q_v) != 0) {
+        ctools_error_alloc("cqreg");
+        cqreg_state_free(state);
+        ctools_filtered_data_free(&filtered);
+        main_debug_close();
+        return 920;
+    }
     state->sum_rdev = cqreg_sum_raw_deviations(state->y, N, state->q_v, quantile);
--- a/validation/cqreg_native.c
+++ b/validation/cqreg_native.c
@@ -69 +69,2 @@
-                assert(cqreg_compute_quantile(y,n,quantiles[q]) == sorted[k-1]);
+                double qv;
+                assert(cqreg_compute_quantile(y,n,quantiles[q],&qv) == 0 && qv == sorted[k-1]);
@@ -135 +136,2 @@
-    double q=cqreg_compute_quantile(y,n,.5);
+    double q;
+    if (cqreg_compute_quantile(y,n,.5,&q) != 0) q = 0.0;
```

Auditor's original suggestion: Report the failure instead of returning a substitute value. Change the function to `ST_int cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau, ST_double *out)`, returning -1 when malloc fails. In cqreg_full_regression, move the q_v/sum_rdev computation below `ctools_filtered_data_free(&filtered)` (line 841) to lower the peak, and on failure do `ctools_error_alloc("cqreg"); rc = 920; goto cleanup;`.

## F057 [P2, confirmed] cqreg_sandwich_product zeroes V on allocation failure and returns void, so cluster and robust-kernel VCEs post zero standard errors with rc 0

Location: src/cqreg/cqreg_vce.c:159-206 cqreg_sandwich_product; callers :393 (cqreg_vce_robust_kernel), :473 (cqreg_vce_cluster_kernel), :805 (cqreg_vce_cluster)

```diff
--- a/src/cqreg/cqreg_vce.c
+++ b/src/cqreg/cqreg_vce.c
@@ -159,21 +159,19 @@
-static void cqreg_sandwich_product(ST_double *V,
-                                   const ST_double *A,
-                                   const ST_double *B,
-                                   ST_int K)
+static ST_int cqreg_sandwich_product(ST_double *V,
+                                     const ST_double *A,
+                                     const ST_double *B,
+                                     ST_int K)
 {
     ST_int i, j, k;
     ST_double *AB = NULL;
 
     /* Compute size with overflow check */
     size_t kk_size;
     if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
-        memset(V, 0, (size_t)K * (size_t)K * sizeof(ST_double));
-        return;
+        return -1;
     }
 
     /* Allocate temporary for A * B */
     AB = (ST_double *)ctools_cacheline_alloc(kk_size);
     if (AB == NULL) {
-        memset(V, 0, kk_size);
-        return;
+        return -1;
     }
@@ -205,2 +203,3 @@
     ctools_aligned_free(AB);
+    return 0;
 }
@@ -393 +392 @@
-    cqreg_sandwich_product(V, H, XtX, K);
+    if (cqreg_sandwich_product(V, H, XtX, K) != 0) goto cleanup;  /* rc is still -1 */
@@ -473 +472 @@
-    cqreg_sandwich_product(V, H, M, K);
+    if (cqreg_sandwich_product(V, H, M, K) != 0) goto cleanup;  /* rc is still -1 */
@@ -805 +804,4 @@
-    cqreg_sandwich_product(V, XtX_inv, M, K);
+    if (cqreg_sandwich_product(V, XtX_inv, M, K) != 0) {
+        rc = -1;
+        goto cleanup;
+    }
```

Auditor's original suggestion: Make it `static ST_int cqreg_sandwich_product(...)` returning -1 on overflow or allocation failure and 0 on success. In each caller, write `if (cqreg_sandwich_product(V, H, XtX, K) != 0) goto cleanup;` with rc still -1; in cqreg_vce_cluster set `rc = -1; goto cleanup;`. cqreg_full_regression then reports the VCE failure instead of posting zeros.

## F060 [P2, confirmed] discrete + method(binsreg) + absorb(): dense N x U bin-indicator matrix and (U+K)^2 normal equations are built before the >500-bin capacity check, giving multi-GB allocations and effectively unbounded runtime on misused/continuous x

Location: src/cbinscatter/cbinscatter_binsreg.c:277-389 (adjust_bins_binsreg_hdfe), reached from src/cbinscatter/cbinscatter_impl.c:785-807 (do_compute_bins). The capacity is set at package/ado/cbinscatter.ado:246-251 and enforced only at cbinscatter_impl.c:474 (store_results).

```diff
--- a/src/cbinscatter/cbinscatter_impl.c
+++ b/src/cbinscatter/cbinscatter_impl.c
@@ -706,6 +706,11 @@
         if (rc != CBINSCATTER_OK) goto cleanup;
     }
 
+    /* The ado sized __cbinscatter_bins; refuse a group whose bins cannot be
+     * stored before the binsreg adjustment does dense N x bins work */
+    const ST_int bin_capacity = SF_row("__cbinscatter_bins");
+    ST_int bins_used = 0;
+
     /* Process each by-group */
     for (ST_int g = 0; g < num_groups; g++) {
         ByGroupResult *group = &results.groups[g];
@@ -792,6 +797,15 @@
                                            group, bin_ids);
         }
 
+        if (rc == CBINSCATTER_OK) {
+            if (group->num_bins > bin_capacity - bins_used) {
+                SF_error("cbinscatter: x has too many distinct values for option discrete\n");
+                rc = 134;
+            } else {
+                bins_used += group->num_bins;
+            }
+        }
+
         /* Binsreg method: regression-adjust the Y-means of the bins assigned
          * above, discrete or quantile (renumbered to index group->bins) */
         if (rc == CBINSCATTER_OK && bin_ids && group->bins != NULL && group->num_bins > 0) {
```

Auditor's original suggestion: Fail fast in C before any adjustment. In do_compute_bins, read the matrix capacity once on the main thread (ST_int bin_cap = SF_row("__cbinscatter_bins");) and keep ST_int rows_used = 0. Right after compute_bins_discrete/compute_bins_single_group returns OK, check: if (rows_used + group->num_bins > bin_cap) { SF_error("cbinscatter: x has more distinct values than discrete supports\n"); rc = CBINSCATTER_ERR_SYNTAX; free group buffers; goto cleanup; } rows_used += group->num_bins;. This runs before densify_bin_ids/adjust_bins_binsreg_hdfe.

## F061 [P2, confirmed] Block-count arithmetic overflows int for N within one block of 2^31: validity flags left uninitialized, or residualization silently skipped

Location: src/cbinscatter/cbinscatter_impl.c:350-351 (load_data validity-flag blocks) and src/cbinscatter/cbinscatter_resid.c:301-302 (ols_residualize residual blocks)

```diff
--- a/src/cbinscatter/cbinscatter_impl.c
+++ b/src/cbinscatter/cbinscatter_impl.c
@@ -349,4 +349,4 @@
         const ST_double missval = SV_missval;
         const ST_int block = 16384;
-        const ST_int num_blocks = (N + block - 1) / block;
+        const ST_int num_blocks = N / block + (N % block != 0);
         ST_int blk;
--- a/src/cbinscatter/cbinscatter_resid.c
+++ b/src/cbinscatter/cbinscatter_resid.c
@@ -300,4 +300,4 @@
     {
         enum { BLOCK = 2048 };
-        const ST_int num_blocks = (N + BLOCK - 1) / BLOCK;
+        const ST_int num_blocks = N / BLOCK + (N % BLOCK != 0);
         ST_int blk;
```

Auditor's original suggestion: Compute the block count without overflow at both sites, e.g. `const ST_int num_blocks = N / block + (N % block != 0);` (impl.c:351) and `const ST_int num_blocks = N / BLOCK + (N % BLOCK != 0);` (resid.c:302), or use `(ST_int)(((int64_t)N + BLOCK - 1) / BLOCK)`. i0 = blk * block stays below N, so nothing else changes.

## F062 [P2, confirmed] memset length K*K computed in int: overflows for >= 46,341 control columns and can memset a huge length

Location: src/cbinscatter/cbinscatter_resid.c:130 (compute_xtx_weighted), called from ols_residualize (line 264) with K = num_controls + 1

```diff
--- a/src/cbinscatter/cbinscatter_resid.c
+++ b/src/cbinscatter/cbinscatter_resid.c
@@ -129,3 +129,3 @@
     /* Initialize to zero */
-    memset(XtX, 0, K * K * sizeof(ST_double));
+    memset(XtX, 0, (size_t)K * (size_t)K * sizeof(ST_double));
```

Auditor's original suggestion: Widen before multiplying: `memset(XtX, 0, (size_t)K * (size_t)K * sizeof(ST_double));`. Or delete the memset, since lines 135-155 assign every element of XtX.

## F066 [P2, confirmed] XLSX import: a repeated <sst> element reallocates the shared-string table without resetting its count (heap overflow)

Location: src/cimport/cimport_xlsx.c:436-453 (shared_strings_callback, <sst> preallocation); the out-of-bounds writes are at 557-559 and the pointer fix-up at 547-548

```diff
--- a/src/cimport/cimport_xlsx.c
+++ b/src/cimport/cimport_xlsx.c
@@ -435,9 +435,12 @@ static bool shared_strings_callback(const xlsx_xml_event *event, void *user_data)
     if (event->type == XLSX_XML_START_ELEMENT) {
         if (strcmp(event->tag_name, "sst") == 0) {
-            /* Pre-allocate from uniqueCount if available */
+            /* Pre-allocate from uniqueCount once, before any string is stored;
+               a repeated <sst> must never replace or shrink the live table. */
             const char *unique_count = xlsx_xml_get_attr(event, "uniqueCount");
-            if (unique_count) {
+            if (unique_count && ctx->num_shared_strings == 0 && !ctx->shared_strings &&
+                !ctx->shared_string_lengths && !ctx->shared_strings_pool) {
                 uint32_t count = (uint32_t)strtoul(unique_count, NULL, 10);
                 if (count > 0 && count <= XLSX_MAX_SHARED_STRINGS) {
@@ -506,6 +509,8 @@ static bool shared_strings_callback(const xlsx_xml_event *event, void *user_data)
             if (ctx->num_shared_strings >= ctx->shared_strings_capacity) {
                 size_t new_cap = ctx->shared_strings_capacity == 0 ?
                                  1024 : ctx->shared_strings_capacity * 2;
+                if (new_cap <= ctx->num_shared_strings)
+                    new_cap = (size_t)ctx->num_shared_strings + 1024;
                 char **new_strings = (char **)realloc(ctx->shared_strings,
                                                        new_cap * sizeof(char *));
```

Auditor's original suggestion: Preallocate only once, and only while empty: `if (unique_count && !ctx->shared_strings && ctx->num_shared_strings == 0)`. Set capacity only after both allocations succeed. In the growth path compute `new_cap = max(capacity * 2, num_shared_strings + 1, 1024)`, so capacity can never fall below the count.

## F067 [P2, confirmed] An error after `cimport scan` leaves the whole parsed CSV context (file mmap plus per-row parse arenas) resident until the next ctools call

Location: package/ado/cimport.ado:379-382 (exit 103), 398-407 (frame copy/preserve), 411-429 (nvar==0 exit 0), 617-624 (`if `rc'` restore/exit); lifecycle src/ctools_runtime.c:210 and 237-238; context g_cimport_ctx (src/cimport/cimport_impl.c:62); sibling package/ado/_cio_import.ado:221-231 (frame copy/preserve after `cio scan` at :169)

```diff
--- a/package/ado/cimport.ado
+++ b/package/ado/cimport.ado
@@ -379,4 +379,5 @@
     if `: word count `explicit_names'' > `nvar' {
         di as error "too many variables specified"
+        capture plugin call ctools_plugin, "cimport clear"
         exit 103
     }
@@ -398,11 +399,17 @@
     if "`fast'" == "" {
         if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
             tempname __cimport_frame
-            frame copy `c(frame)' `__cimport_frame'
+            capture noisily frame copy `c(frame)' `__cimport_frame'
         }
         else {
-            preserve
+            capture noisily preserve
             local __cimport_snap 1
         }
+        if _rc {
+            local rc = _rc
+            capture plugin call ctools_plugin, "cimport clear"
+            exit `rc'
+        }
     }
@@ -426,4 +433,5 @@
         return local delimiters `"`macval(returned_delimiters)'"'
+        capture plugin call ctools_plugin, "cimport clear"
         if "`__cimport_frame'" != "" frame drop `__cimport_frame'
@@ -618,4 +626,6 @@
     local rc = _rc
     if `rc' {
+        * Free the parsed file before restoring, which may need the memory
+        capture plugin call ctools_plugin, "cimport clear"
         nobreak {
--- a/package/ado/_cio_import.ado
+++ b/package/ado/_cio_import.ado
@@ -223,9 +223,15 @@
         if "`snapshot'" == "" & c(stata_version) >= 16 & (_N > 0 | c(k) > 0) {
             tempname snapframe
-            frame copy `c(frame)' `snapframe'
+            capture noisily frame copy `c(frame)' `snapframe'
         }
         else {
-            preserve
+            capture noisily preserve
             local snap 1
         }
+        if _rc {
+            local rc = _rc
+            capture plugin call ctools_plugin, "cio clear"
+            macro drop CTOOLS_CIO_*
+            exit `rc'
+        }
     }
```

Auditor's original suggestion: In cimport.ado, call `capture plugin call ctools_plugin, "cimport clear"` before `exit 103` and inside the `if `rc'` block before `exit`. Put the frame copy/preserve under `capture noisily`, with a clear on failure, or move the snapshot before the scan. Do the same in _cio_import.ado around the snapframe/preserve block and the describe loop (`capture plugin call ctools_plugin, "cio clear"`).

## F068 [P2, confirmed] crangejoin keeps its full using-data index when _crangejoin_using fails after its plugin call: the clear guard global is set too late

Location: package/ado/crangejoin.ado:133-135 (_crangejoin_impl frame path), with the guard in the wrapper at crangejoin.ado:28 and :44

```diff
--- a/package/ado/crangejoin.ado
+++ b/package/ado/crangejoin.ado
@@ -130,8 +130,10 @@
         * and only their empty schema is appended to the master.
         tempname uframe
         frame create `uframe'
+        * Set before the plugin indexes the using data, so a failure after
+        * that call still has the wrapper release them
+        global CTOOLS_crangejoin_plugin "`__ctools_plugin'"
         frame `uframe': _crangejoin_using `"`using'"' `"`schema'"', `usingopts'
         frame drop `uframe'
-        global CTOOLS_crangejoin_plugin "`__ctools_plugin'"
         _crangejoin_master `"`schema'"', mastervars(`mastervars') lo(`lo') hi(`hi') ///
             renamed(`renamed') nmaster(`nmaster') byindices(`byindices') ///
```

Auditor's original suggestion: Set `global CTOOLS_crangejoin_plugin "`__ctools_plugin'"` before line 133. A rollback with state.ready==1 frees the state without writing, so this is safe. Alternatively, wrap the post-call steps of _crangejoin_using in `capture noisily` and call `crangejoin clear` on error.

## F003 [P3, latent] Division by zero in the store scheduler when storing a stata_data with zero variables

Location: src/ctools_data_store.c, ctools_io_data_store_impl: nvars is taken from data->nvars when var_indices == NULL (line 274) with no zero check, and the call ctools_io_chunks_per_column(nvars, nobs, MIN_OBS_PER_THREAD) is at line 433. The division is in src/ctools_transport.c:311, ctools_io_chunks_per_column: `(4 * nthreads + nvars - 1) / nvars`. The callers that can pass a zero-variable dataset are ctools_data_store (ctools_data_store.c:505), ctools_data_store_ex with var_indices == NULL (line 493) and ctools_data_store_sorted (line 498).

In ctools_io_data_store_impl, after resolving nvars, add: if (nvars == 0) return STATA_OK; (nothing to write). Also make ctools_io_chunks_per_column defensive: if (nvars == 0 || min_chunk_rows == 0) return 0;

## F007 [P3, latent] ctools_int64_to_str negates INT64_MIN (signed-overflow UB)

Location: src/ctools_types.c, ctools_int64_to_str (line 152 in the current tree): `return 1 + ctools_uint64_to_str((uint64_t)(-val), buf + 1);`

Negate in unsigned arithmetic: `return 1 + ctools_uint64_to_str(0 - (uint64_t)val, buf + 1);`

## F008 [P3, latent] Pairs radix sort orders only the low 32 bits of the size_t order_key, and its tie order follows input order rather than orig_idx (latent; SORT-12 unfixed)

Location: src/ctools_sort_pairs.c, ctools_sort_order_pairs: the fixed `for (int pass = 0; pass < 4; pass++) { int shift = pass * 8; ... }` radix loop, about lines 131-133. The contract it breaks is in src/ctools_sort_pairs.h at the ctools_order_pair_t typedef (`size_t order_key`) and in the doc comment on ctools_sort_order_pairs ("pairs with equal order_key values are ordered by orig_idx"), about lines 22-37. The only caller is src/cmerge/cmerge_impl.c, the preserve_order block of the plan phase (pairs fill loop and ctools_sort_order_pairs call, about lines 842-849).

Compute the maximum key once, then run `passes = max_key ? (bit_width(max_key) + 7) / 8 : 1` passes, capped at 8, instead of a fixed 4. Alternatively, narrow order_key to uint32_t and document both preconditions in ctools_sort_pairs.h: keys < 2^32, and orig_idx ascending in input order.

## F013 [P3, latent] Flat string buffer in streaming csort trusts the ado width hint: SF_sdata writes into the slot with no guard tail and no length check

Location: src/csort/csort_stream.c, stream_permute_string_var_flat (currently lines 211-237): `stride = (size_t)max_width + 1` (219), `calloc(nobs, stride)` (222), unguarded `SF_sdata(stata_var, obs_map[i], buf + inv_perm[i] * stride)` (227), SF_sdata/SF_sstore return codes ignored (227, 232). Hint source: src/csort/csort_impl.c parse_string_widths (about 217-233, `widths[j] = (int)val;` with no range check), passed through csort_main's stream branch (about 405-409) and csort_stream_sort's `nonkey_str_widths[nk] = str_widths[all_var_indices[i] - 1]` (about 951) into csort_stream_apply_permutation -> stream_permute_string_var (316).

Read each value into a local `char tmp[2046]`, reject `strlen(tmp) > max_width` with an error, then memcpy into the slot. Alternatively, allocate `nobs + (2046 + stride - 1) / stride` slots as ctools_data_io.c:54-55 does and check strlen after the direct read. Also clamp parsed widths to 0..2045 in parse_string_widths.

## F014 [P3, latent] Streaming csort's 16-wide read loop reads observation obs_map[i] + bb instead of obs_map[i + bb]

Location: src/csort/csort_stream.c, csort_stream_apply_permutation: the OpenMP Phase 1 scatter read at line 495 and the non-OpenMP read at line 624. Both are `SF_vdata(stata_var, (ST_int)(obs_map[i] + bb), &batch_vals[bb])` inside `for (i = 0; i < nobs_batch16; i += 16)`. The only caller is csort_stream_sort (around line 963), reached from csort_impl.c around line 408.

Use `(ST_int)obs_map[i + bb]` at lines 495 and 624.

## F015 [P3, latent] cmerge plugin argument parsers do not cross-check counts: negative counts, uninitialized/stale index arrays, and an unchecked keepusing count that indexes the using cache

Location: src/cmerge/cmerge_impl.c: cmerge_execute, the plan-consistency check (~1142-1153) and the static parse arrays (~1057-1066); cmerge_load_using's count parsing (~194-249). The faulting access is src/cmerge/cmerge_io.c: cmerge_write_keepusing_var_thread, `stata_variable *src = &g_using_cache.keepusing.data.vars[a->keepusing_idx];` followed by reads of src->type, src->nobs and src->data. The unchecked read of the destination index is the `SF_vdata(dest_idx, ...)` call in the shared-variable update/replace branch.

In load_using, reject nkeys < 0 and n_keepusing < 0, and after the token loop return 198 unless arg_idx > nkeys and keepusing_idx == n_keepusing. In execute, memset the three static arrays at entry; after parsing, require keepusing_idx == n_keepusing, the number of shared flags parsed == n_keepusing, n_keepusing == g_using_cache.n_keepusing, and 1 <= index <= SF_nvars() for every placeholder and for merge_var_idx (when > 0). Fail with 198 before any write.

## F018 [P3, confirmed] cdestring total converted/failed counters are int and overflow (signed UB) when nvars x nobs > INT_MAX

Location: src/cdestring/cdestring_impl.c, cdestring_main: the `int total_converted` / `int total_failed` declarations at the top of the function (currently lines 293-294), the "Aggregate results" loop that adds `var_converted[v]` / `var_failed[v]` into them (lines 466-469), and the `SF_scal_save("_cdestring_n_converted"/"_cdestring_n_failed", (double)total_...)` calls (lines 479-480). The values are used in package/ado/cdestring.ado at lines 275-276, 322-323 (the force note) and 370-373 (the verbose summary).

Declare the totals as `double total_converted = 0, total_failed = 0;` (or uint64_t), accumulate with `total_converted += (double)var_converted[v];`, and save them without the int cast.

## F035 [P3, confirmed] cio_days_civil signed overflow for Excel date-formatted serials near INT_MAX

Location: src/io/cio_context.c:148-150 (cio_days_civil) reached via cio_excel_date_text:160-180 from cio_xlsx.c:187 and cio_xls.c:312

Tighten the guard in cio_excel_date_text to `days < -715875.0 || days >= 2936550.0` (0000-01-01..9999-12-31) returning ".", or compute z/era in long long inside cio_days_civil.

Auditor's original suggestion: Tighten the guard to the domain cio_days_civil supports: `days < (double)INT_MIN + 146097 \|\| days > (double)INT_MAX - 715815` gives ".". Alternatively, compute `z` and `era` in `long long` inside cio_days_civil.

## F037 [P3, confirmed] cexport shp converts an arbitrary _ID double to uint32_t (UB for missing, negative or huge IDs)

Location: src/io/cio_shp.c:395-408 (cio_shp_export)

Before the conversion: `if (SF_is_missing(id) || id < 0 || id > (double)UINT32_MAX || id != floor(id)) { free(blob); rc = 459; goto done; }` (or 198).

Auditor's original suggestion: Validate before converting: `if (SF_is_missing(id) \|\| id < 0 \|\| id > (double)UINT32_MAX \|\| id != floor(id)) { free(blob); rc = 198; goto done; }`.

## F038 [P3, confirmed] cexport excel .xls: cell() column letters accumulate without a bound (signed int overflow); row parsed with atoi

Location: src/io/cexport_xls.c:308-320 (cexport_xls_main, cell= parsing); check at :326

Replace the hand-rolled parse with `int c, r; if (!ctools_excel_address(t + 5, 256, 65536, false, &c, &r, NULL)) { rc = 198; goto done; } col0 = c - 1; row0 = r - 1;` (optionally tighten the ado regex to ^[A-Z]{1,3}[0-9]{1,7}$).

Auditor's original suggestion: Bail out in the loop as soon as `col0 > 256` (rc 198). Parse the row with `strtol` with errno/range checks (1..65536). Optionally restrict the ado regex to `^[A-Z]{1,3}[0-9]{1,7}$`.

## F040 [P3, latent] cio_main leaves `kind`/`index` uninitialized when the argument string has only an operation

Location: src/io/cio.c:6-8 (cio_main)

Initialize `char op[24] = "", kind[24] = "";` and return 198 from workbook/scan/column/label/textnumbers/export when `!*kind`.

Auditor's original suggestion: `char op[24] = "", kind[24] = "";`, and return 198 for ops that require a kind when `!*kind`.

## F044 [P3, latent] creghdfe inverts the Cholesky factor in place, aliasing L and inv against ctools_invert_from_cholesky's documented contract

Location: src/creghdfe/creghdfe_regress.c:996 do_full_regression; contract in src/ctools_ols.h:37-45, implementation src/ctools_ols.c:65-98

Amend ctools_ols.h so the contract allows inv == L, since L is fully consumed into L_inv before inv is written. Optionally drop creghdfe's memcpy by factoring xtx_keep, which is dead after line 988, and inverting into inv_xx_x.

Auditor's original suggestion: Factor xtx_keep itself, which is dead after this point, and invert into the separate buffer. Drop the memcpy: `if (ctools_cholesky(xtx_keep, K_keep) != 0) { ... } stata_retcode inverse_status_0 = ctools_invert_from_cholesky(xtx_keep, K_keep, inv_xx_x);`. This adds no allocation and leaves inv_xx_x with the same contents as today.

## F048 [P3, confirmed] Signed overflow in the missing-value block count when more than INT_MAX-16383 observations are selected

Location: src/civreghdfe/civreghdfe_impl.c:520 (do_iv_regression, num_blocks = (N_total + block - 1) / block in ST_int)

```diff
--- a/src/civreghdfe/civreghdfe_impl.c
+++ b/src/civreghdfe/civreghdfe_impl.c
@@ -519,3 +519,3 @@
         const ST_int block = 16384;
-        const ST_int num_blocks = (N_total + block - 1) / block;
+        const ST_int num_blocks = N_total / block + (N_total % block != 0);
         ST_int blk;
```

Auditor's original suggestion: Replace line 520 with: const ST_int num_blocks = N_total / block + (N_total % block != 0);

## F049 [P3, confirmed] CUE copies beta onto itself with memcpy (overlapping memcpy is UB)

Location: src/civreghdfe/civreghdfe_estimate.c:287 (ivest_compute_cue: memcpy(beta, initial_beta, ...)), called at estimate.c:1246-1249 with initial_beta == beta

```diff
--- a/src/civreghdfe/civreghdfe_estimate.c
+++ b/src/civreghdfe/civreghdfe_estimate.c
@@ -286,2 +286,3 @@
     /* Initialize beta from input (2SLS estimate) */
-    memcpy(beta, initial_beta, K_total * sizeof(ST_double));
+    if (beta != initial_beta)
+        memcpy(beta, initial_beta, K_total * sizeof(ST_double));
```

Auditor's original suggestion: Use `if (beta != initial_beta) memcpy(beta, initial_beta, K_total * sizeof(ST_double));` (or memmove).

## F051 [P3, latent] civ_omega_with_kernel can turn a two-way-with-kernel builder into a kernel-free one whose intersection clusters were never built (cl12 = NULL, ncl12 = 0)

Location: src/civreghdfe/civreghdfe_omega.c:127-138 (civ_omega_with_kernel) with civ_omega_init:91-101 (cl12 built only for TWOWAY without a kernel) and civ_omega_build:362-363

In civ_omega_with_kernel, refuse to drop the kernel of a TWOWAY original that has no intersection ids: `if (c.kind == CIV_OMEGA_TWOWAY && om->kernel > 0 && kernel <= 0) kernel = CIVREGHDFE_KERNEL_BARTLETT;` (or add `if (!om->cl12) return 198;` in the TWOWAY kernel-free branch of civ_omega_build).

Auditor's original suggestion: In civ_omega_with_kernel, keep the original kernel for TWOWAY when the copy would need intersection ids that were not built: `if (c.kind == CIV_OMEGA_TWOWAY && om->kernel > 0 && kernel <= 0) kernel = CIVREGHDFE_KERNEL_BARTLETT;` (or add `if (!om->cl12) return 198;` in the TWOWAY kernel-free branch of civ_omega_build) so the state can never be used inconsistently.

## F052 [P3, latent] int products D*D and K1*K3 in the Kleibergen-Paap kernel overflow for K_endog*L > 46340

Location: src/civreghdfe/civreghdfe_tests.c:375/436/468-469 (civ_kp_rk: D = K*L and `i < D * D` in int), src/ctools_matrix.c:126 (ctools_matmul_ab memset(C, 0, K1 * K3 * sizeof) with an int product)

Add an up-front guard in civreghdfe (e.g. impl.c before ivest_compute_2sls): refuse K_iv > 46340 or K_endog*(K_iv-K_exog) > 46340 with r(198). Optionally use size_t counters in civ_kp_rk (`for (size_t i = 0; i < (size_t)D * D; i++)`), use `(size_t)K1 * K3` in ctools_matmul_ab's memset, and allocate B/T only when K > 1.

Auditor's original suggestion: Use size_t loop counters (`for (size_t i = 0; i < (size_t)D * D; i++)`) and `memset(C, 0, (size_t)K1 * K3 * sizeof(ST_double))` in ctools_matmul_ab. Allocate B and T only on the K > 1 path. Optionally refuse D > 46340 up front with rc = -1.

## F053 [P3, latent] A failed write of the d() variable is reported as success because ppml_save_effects overwrites output_rc

Location: src/cpplmhdfe/cpplmhdfe_irls.c:1562-1575 do_ppml_regression (d() store, then output_rc = ppml_save_effects at :1574); also :676-678 tag/certificate stores

Before line 1574 add `if (output_rc) goto ppml_cleanup;` (or use `if (!output_rc) output_rc = ppml_save_effects(...)`), and guard the certificate store at :678 with `!output_rc &&`.

Auditor's original suggestion: Keep the first error. Before line 1574 add `if (output_rc) goto ppml_cleanup;`, or write `if (!output_rc) output_rc = ppml_save_effects(...);`. At lines 677-678, guard the second store: `if(!output_rc && options->certificate_idx) output_rc = SF_vstore(...)`.

## F054 [P3, latent] K_x*K_x is computed in int for matrix buffer sizes and row-major indices (signed overflow for K_x >= 46341)

Location: src/cpplmhdfe/cpplmhdfe_irls.c:548,551 (first K_x*K_x allocations), also :770,807,810,885,906,1132,1183-1186,1199,1264,1412,1485; src/cpplmhdfe/cpplmhdfe_solvers.c:19

Reject K_x > 46340 with r(198) at the top of do_ppml_regression (or widen the products with `(size_t)K_x * K_x`), matching the civreghdfe cap proposed in F042.

Auditor's original suggestion: Widen before multiplying: `malloc((size_t)K_x * (size_t)K_x * sizeof(ST_double))` (likewise for K_keep, K_with_cons, and `(size_t)K * K` in ppml_qrsolve). Index with `(size_t)a * K + b`. Alternatively, reject K_x > 46340 with r(198) at the start of do_ppml_regression.

## F055 [P3, latent] ctools_invert_from_cholesky is called with L aliasing inv, against its documented contract

Location: src/cpplmhdfe/cpplmhdfe_irls.c:889, :908, :1202 do_ppml_regression

Amend the ctools_ols.h contract to state that inv may alias L (the implementation already guarantees this); do not factor xtx_keep in place at :1202, because the iterative refinement reads it.

Auditor's original suggestion: Factor into a separate buffer and invert into the destination. The existing workspaces already allow this: run Cholesky on xtx_k (or xtx) itself, which is not reused after the solve, then call `ctools_invert_from_cholesky(xtx_k, K_irls, inv_xx_k)`. Alternatively, change the header to document and guarantee in-place support.

## F058 [P3, confirmed] Uninitialized u/v arrays are read after a non-converged main solve

Location: src/cqreg/cqreg_regress.c:1008-1016 cqreg_full_regression (reads at src/cqreg/cqreg_types.c:174,186); u/v written only at src/cqreg/cqreg_fn.c:1612-1619, early return at :1496-1500

Move the `if (!state->converged) { ctools_error(...); rc = 430; goto cleanup; }` block from :1029-1033 to directly after `state->converged = ...` at :1008, before the get_objective and get_residuals calls.

Auditor's original suggestion: Move the `if (!state->converged) { ...; rc = 430; goto cleanup; }` check (lines 1029-1033) directly after line 1008, before the cqreg_ipm_get_objective and cqreg_ipm_get_residuals calls.

## F059 [P3, latent] C side never checks N > K; with N < K the crossover reads past the N-element idx array

Location: src/cqreg/cqreg_regress.c:756 (only guard N <= 0); OOB read in src/cqreg/cqreg_fn.c:1568 -> cqreg_basis_solve reading basis[bi] for bi < K (idx has N entries, :1544)

In cqreg_full_regression after the N <= 0 check add `if (N <= K) { ctools_filtered_data_free(&filtered); ctools_error("cqreg", "not enough observations"); return 2001; }`.

Auditor's original suggestion: In cqreg_full_regression, after `ST_int N = (ST_int)filtered.data.nobs;`, add `if (N <= K) { ctools_filtered_data_free(&filtered); ctools_error("cqreg", "not enough observations"); return 2001; }`, matching the ado. In cqreg_fn_solve_core, also return -1 at entry when `N < K`.

## F063 [P3, latent] By-group count derived as max-min+1 from a raw (ST_int) cast of the by value: UB/overflow and dropped groups if the dense-id contract from the ado is ever broken

Location: src/cbinscatter/cbinscatter_impl.c:316 (load_data by() cast) and 653-667 (do_compute_bins num_groups = max-min+1)

Hardening only: in load_data map any by value that is missing, < 1, > INT_MAX or non-integer to -1, and in do_compute_bins use num_groups = max_group with the min_group shift removed.

Auditor's original suggestion: Enforce the contract in C. In load_data map any by value that is not an integer in [1, INT_MAX] to -1 (`by_groups[i] = (SF_is_missing(val) \|\| val < 1 \|\| val > INT_MAX \|\| val != floor(val)) ? -1 : (ST_int)val;`). In do_compute_bins set num_groups = max_group, since ids are 1-based, and drop the min_group shift, so an emptied low group no longer hides the top group.

## F064 [P3, latent] compute_bin_statistics places the SE double accumulators after the int counts array, so they are misaligned whenever nquantiles is odd

Location: src/cbinscatter/cbinscatter_bins.c:174-181 (compute_bin_statistics block partition)

Place the ST_int counts array after sum_y2/sum_x2 (sum_y2 = sum_w + nq; sum_x2 = sum_y2 + nq; counts = (ST_int *)(sum_x2 + nq) when compute_se, else counts = (ST_int *)(sum_w + nq)). The allocation size is unchanged.

Auditor's original suggestion: Put the int array last: sum_y = block; sum_x = sum_y + nq; sum_w = sum_x + nq; if (compute_se) { sum_y2 = sum_w + nq; sum_x2 = sum_y2 + nq; counts = (ST_int *)(sum_x2 + nq); } else counts = (ST_int *)(sum_w + nq);. The total size is unchanged.

## F065 [P3, latent] Public store helpers call SF_error from OpenMP worker threads (cwinsor/cdestring parallel column stores), reachable with Break

Location: src/ctools_data_store.c:636-642 (ctools_store_filtered -> ctools_io_store_status -> SF_error at :3-7); reached from omp workers at src/cwinsor/cwinsor_impl.c:712-716 and :233-237 (cwinsor_undo_restore), and src/cdestring/cdestring_impl.c:452-457

Add a quiet variant (e.g. ctools_store_filtered_noreport) that returns ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE) without calling SF_error. Use it in the cwinsor/cdestring omp loops, and call ctools_io_store_status(store_rc) once on the calling thread after each region.

Auditor's original suggestion: Split out a quiet variant (e.g. ctools_store_filtered_noreport) that returns the rc without calling SF_error. Use it inside the parallel loops, and call ctools_io_store_status once on the calling thread after the region.

## F069 [P3, latent] cimport blob dereferences a NULL string cache for string columns that were not chosen as strL

Location: src/cimport/cimport_impl.c:2001-2007 cimport_blob (accepts any CIMPORT_COL_STRING column, then indexes col_cache[col-1].string_data, which cimport_build_column_cache:1516-1524 leaves NULL unless use_strl)

In cimport_blob's argument check, add `|| !ctx->columns[col-1].use_strl`. After cimport_build_column_cache, also return 198 when `!ctx->col_cache || !ctx->col_cache[col-1].string_data || !ctx->col_cache[col-1].string_data[row-1]`.

Auditor's original suggestion: Add `\|\| !ctx->columns[col-1].use_strl` to the argument check. After the cache build, also return 198 if `!ctx->col_cache \|\| !ctx->col_cache[col-1].string_data \|\| !ctx->col_cache[col-1].string_data[row-1]`.
