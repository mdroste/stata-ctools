/* Fixed-width string tokenization. The scan stages data/widths; write commits
 * only to ado-owned temporary variables. No Stata calls inside parsing loops. */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "csplit_impl.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_parse.h"
#include "ctools_runtime.h"

#define MAXSTR 2045
static struct {
    ctools_filtered_data input;
    char **delims;
    uint32_t *tokens;
    int lengths[MAXSTR + 1], next[MAXSTR + 1], head[256];
    char firstbytes[257];
    int ndelim, limit, trim, fields, verbose, ready;
    size_t nobs;
} state;

void csplit_cleanup_cache(void)
{
    ctools_filtered_data_free(&state.input);
    for (int i = 0; i < state.ndelim; i++) free(state.delims ? state.delims[i] : NULL);
    free(state.delims);
    free(state.tokens);
    memset(&state, 0, sizeof(state));
}

/* Match split's byte-string semantics: trim the initial string unless notrim;
 * default whitespace parsing also trims each remainder. Ties between parse()
 * strings go to the last listed delimiter, as in split. */
static int tokenize(const char *s, const char **starts, int *lengths)
{
    const char *end = s + strlen(s);
    if (state.trim) {
        while (s < end && *s == ' ') s++;
        while (end > s && end[-1] == ' ') end--;
    }
    int k = 0;
    do {
        const char *best = end;
        size_t delimiter_length = 0;
        if (state.ndelim == 1) {
            const char *p = strstr(s, state.delims[0]);
            if (p && p + state.lengths[0] <= end) {
                best = p; delimiter_length = state.lengths[0];
            }
        } else {
            /* Scan the string once. Only test delimiters whose first byte
             * matches, newest first to retain split's last-listed tie rule. */
            const char *p = s;
            while (p < end) {
                p += strcspn(p, state.firstbytes);
                if (p >= end) break;
                for (int d = state.head[(unsigned char)*p]; d >= 0; d = state.next[d]) {
                    size_t len = state.lengths[d];
                    if ((size_t)(end-p) >= len && !memcmp(p, state.delims[d], len)) {
                        best = p; delimiter_length = len; break;
                    }
                }
                if (delimiter_length) break;
                p++;
            }
        }
        starts[k] = s;
        lengths[k++] = (int)(best - s);
        if (!delimiter_length) break;
        s = best + delimiter_length;
        if (state.trim == 2) {
            while (s < end && *s == ' ') s++;
            while (end > s && end[-1] == ' ') end--;
        }
    } while (s < end && k < state.limit);
    return k;
}

ST_retcode csplit_main(const char *args)
{
    int rc = 0;
    if (!strcmp(args, "clear")) { csplit_cleanup_cache(); return 0; }
    if (!strncmp(args, "scan ", 5)) {
        csplit_cleanup_cache();
        args += 5;
        if (ctools_parse_next_int(&args, &state.ndelim) || state.ndelim < 1 || state.ndelim > 2045 ||
            ctools_parse_next_int(&args, &state.limit) || state.limit < 1 ||
            ctools_parse_next_int(&args, &state.trim) || state.trim < 0 || state.trim > 2 ||
            ctools_parse_next_int(&args, &state.verbose) || SF_nvars() != 1 || !SF_var_is_string(1)) {
            memset(&state, 0, sizeof(state)); return 198;
        }
        if (state.limit > MAXSTR + 1) state.limit = MAXSTR + 1;
        double start = ctools_timer_seconds();
        state.delims = calloc(state.ndelim, sizeof(*state.delims));
        if (!state.delims) { rc = 920; goto fail; }
        for (int b = 0; b < 256; b++) state.head[b] = -1;
        int nfirst = 0;
        for (int d = 0; d < state.ndelim; d++) {
            char name[64];
            snprintf(name, sizeof(name), "_csplit_delim%d", d+1);
            state.delims[d] = malloc(MAXSTR + 2);
            if (!state.delims[d]) { rc = 920; goto fail; }
            rc = SF_macro_use(name, state.delims[d], MAXSTR + 2);
            if (rc || !state.delims[d][0]) { rc = 198; goto fail; }
            state.lengths[d] = (int)strlen(state.delims[d]);
            unsigned char first = (unsigned char)state.delims[d][0];
            if (state.head[first] == -1) state.firstbytes[nfirst++] = (char)first;
            state.next[d] = state.head[first]; state.head[first] = d;
        }
        int var = 1;
        if (ctools_data_load(&state.input, &var, 1, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER)) { rc = 920; goto fail; }
        state.nobs = SF_nobs();
        state.tokens = malloc(state.input.data.nobs * sizeof(*state.tokens));
        if (!state.tokens) { rc = 920; goto fail; }
        double loaded = ctools_timer_seconds();
        int widths[MAXSTR + 1] = {0};
        #pragma omp parallel num_threads(ctools_get_max_threads()) if(state.input.data.nobs > 10000)
        {
            int localwidth[MAXSTR + 1] = {0}, maxfields = 0;
            const char *starts[MAXSTR + 1];
            int lengths[MAXSTR + 1];
            #pragma omp for schedule(static)
            for (size_t i = 0; i < state.input.data.nobs; i++) {
                int k = tokenize(state.input.data.vars[0].data.str[i], starts, lengths);
                if (k > maxfields) maxfields = k;
                char *row = state.input.data.vars[0].data.str[i], *out = row;
                size_t offset = state.ndelim == 1 ? (size_t)(starts[0]-row) : 0;
                state.tokens[i] = ((uint32_t)k << 16) | (uint32_t)offset;
                for (int j = 0; j < k; j++) {
                    if (lengths[j] > localwidth[j]) localwidth[j] = lengths[j];
                    if (state.ndelim == 1) {
                        /* Fixed delimiter width lets write walk the tokens in
                         * place. Avoid copying long fields just to remove it. */
                        ((char *)starts[j])[lengths[j]] = '\0';
                    } else {
                        /* Tokens plus terminators fit in the allocation: each
                         * consumed delimiter supplies at least one byte. */
                        memmove(out, starts[j], (size_t)lengths[j]);
                        out[lengths[j]] = '\0'; out += lengths[j] + 1;
                    }
                }
            }
            #pragma omp critical
            {
                if (maxfields > state.fields) state.fields = maxfields;
                for (int j = 0; j < maxfields; j++) if (localwidth[j] > widths[j]) widths[j] = localwidth[j];
            }
        }
        char sizes[(MAXSTR + 1) * 6 + 1];
        size_t used = 0;
        for (int j = 0; j < state.fields; j++)
            used += snprintf(sizes + used, sizeof(sizes) - used, "%d ", widths[j] ? widths[j] : 1);
        sizes[used] = '\0';
        rc = SF_macro_save("_csplit_widths", sizes);
        if (rc) goto fail;
        state.ready = 1;
        ctools_verbose("csplit", state.verbose, "load %.4fs; tokenize/size %.4fs; %d output fields",
            loaded-start, ctools_timer_seconds()-loaded, state.fields);
        return 0;
    }
    if (!strcmp(args, "write")) {
        if (!state.ready || SF_nvars() != state.fields || (size_t)SF_nobs() != state.nobs) { rc = 198; goto fail; }
        for (int j = 1; j <= state.fields; j++)
            if (!SF_var_is_string(j) || SF_var_is_strl(j)) { rc = 109; goto fail; }
        double start = ctools_timer_seconds();
        #pragma omp parallel num_threads(ctools_get_max_threads()) if(state.input.data.nobs > 10000)
        {

            int localrc = 0;
            #pragma omp for schedule(static)
            for (size_t i = 0; i < state.input.data.nobs; i++) {
                int count = (int)(state.tokens[i] >> 16);
                char *value = state.input.data.vars[0].data.str[i] + (state.tokens[i] & 65535);
                for (int j = 0; j < count; j++) {
                    int err = SF_sstore(j+1, state.input.obs_map[i], value);
                    if (err) localrc = err;
                    if (j+1 < count) {
                        value += strlen(value) + (state.ndelim == 1 ? state.lengths[0] : 1);
                        if (state.ndelim == 1 && state.trim == 2) while (*value == ' ') value++;
                    }
                }
            }
            #pragma omp critical
            { if (localrc) rc = localrc; }
        }
        ctools_verbose("csplit", state.verbose, "write %.4fs", ctools_timer_seconds()-start);
        csplit_cleanup_cache();
        return rc;
    }
    rc = 198;
fail:
    csplit_cleanup_cache();
    return rc;
}
