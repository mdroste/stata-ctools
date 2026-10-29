/*
    cdecode_impl.c
    Value-label decoding: write each observation's label text to a string variable.

    cdecode.ado reads the label through Mata (literal text, no label-save
    parsing), applies maxlength(), creates the empty destination variable and
    hands the mapping over in a local macro. The plugin only maps values to
    text. Anything this path does not cover (strL results, labels that need
    multibyte truncation, mappings too large for a macro) stays with Stata's
    decode, which the ado also falls back to if this entry point reports an
    error.

    Wire protocol, version 2:
        plugin call ctools_plugin <numeric source> <string destination> [marker] [if] [in],
            "cdecode v2 maplen=# width=# [ifcond [touse]] [fillall]"
    Local __cd_map holds maplen bytes: "<count> " followed by one entry per
    mapping, "<value> <bytes> <text> ", where <value> is an integer or .a-.z and
    <text> is exactly <bytes> bytes (at most width), and a closing "|". ifcond
    says the call carries an if condition to evaluate; with touse a numeric
    sample marker follows the destination in its place. Without fillall the
    destination must be all "" on entry and only observations with nonempty
    text are written; with fillall (a destination created without initial
    values, and no if or in) every observation is written.

    Earlier wrappers sent label-save text, which cannot represent every label;
    their calls are rejected rather than decoded.
*/
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_threads.h"
#include "ctools_parse.h"
#include "cdecode_impl.h"

/* Values above this span are found by binary search instead of a table. */
#define CDECODE_TABLE_MAX_SPAN ((int64_t)1 << 22)
/* Stata's largest label count; a larger header is not a label mapping. */
#define CDECODE_MAX_LABELS 65536
/* Bit patterns of . and .z; the 26 extended codes are 2^40 apart. */
#define CDECODE_SYSMISS_BITS UINT64_C(0x7FE0000000000000)
#define CDECODE_MISSING_STEP 40

typedef struct
{
    char *map;               /* owned copy of the macro; entries are terminated in place */
    size_t count;            /* integer-valued mappings */
    int64_t *values;         /* owned [count], ascending */
    const char **text;       /* owned [count], pointers into map ("" when the label is empty) */
    const char *missing[27]; /* [1..26]: text of .a-.z, NULL when unlabeled */
    int64_t table_min;       /* direct table covers table_min .. table_min + table_span - 1 */
    int64_t table_span;      /* 0 when lookups use binary search */
    int32_t *table;          /* owned [table_span]: index into text, -1 when unlabeled */
} cdecode_map;

static void cdecode_map_free(cdecode_map *m)
{
    free(m->map);
    free(m->values);
    free((void *)m->text);
    free(m->table);
    memset(m, 0, sizeof(*m));
}

typedef struct
{
    int64_t value;
    const char *text;
} cdecode_entry;

static int cdecode_entry_compare(const void *a, const void *b)
{
    int64_t x = ((const cdecode_entry *)a)->value, y = ((const cdecode_entry *)b)->value;
    return (x > y) - (x < y);
}

/* Parse the macro text in place. Returns 0, or a Stata return code. */
static int cdecode_map_parse(cdecode_map *m, size_t maplen, size_t width)
{
    char *p = m->map, *end = m->map + maplen, *stop;
    long declared = strtol(p, &stop, 10);
    if (stop == p || declared < 0 || declared > CDECODE_MAX_LABELS || stop >= end || *stop != ' ')
        return 198;
    p = stop + 1;

    cdecode_entry *entries = (cdecode_entry *)calloc((size_t)declared + 1, sizeof(*entries));
    if (!entries)
        return 920;
    size_t count = 0;
    for (long e = 0; e < declared; e++)
    {
        int missing_code = 0;
        int64_t value = 0;
        if (p + 2 < end && p[0] == '.' && p[1] >= 'a' && p[1] <= 'z' && p[2] == ' ')
        {
            missing_code = p[1] - 'a' + 1;
            p += 2;
        }
        else
        {
            long long parsed = strtoll(p, &stop, 10);
            if (stop == p || stop >= end || *stop != ' ')
                goto malformed;
            value = (int64_t)parsed;
            p = stop;
        }
        p++; /* the space after the value */
        long bytes = strtol(p, &stop, 10);
        if (stop == p || bytes < 0 || (size_t)bytes > width || stop >= end || *stop != ' ')
            goto malformed;
        p = stop + 1;
        /* The text is followed by the entry's separating space. */
        if ((size_t)(end - p) < (size_t)bytes + 1 || p[bytes] != ' ')
            goto malformed;
        p[bytes] = '\0';
        if (missing_code)
            m->missing[missing_code] = p;
        else
        {
            entries[count].value = value;
            entries[count].text = p;
            count++;
        }
        p += bytes + 1;
    }

    qsort(entries, count, sizeof(*entries), cdecode_entry_compare);
    for (size_t i = 1; i < count; i++)
        if (entries[i].value == entries[i - 1].value)
            goto malformed;
    m->values = (int64_t *)malloc((count + 1) * sizeof(int64_t));
    m->text = (const char **)malloc((count + 1) * sizeof(char *));
    if (!m->values || !m->text)
    {
        free(entries);
        return 920;
    }
    for (size_t i = 0; i < count; i++)
    {
        m->values[i] = entries[i].value;
        m->text[i] = entries[i].text;
    }
    m->count = count;
    free(entries);

    if (count > 0)
    {
        /* Label values are 32-bit integers, so the span cannot overflow. */
        int64_t span = m->values[count - 1] - m->values[0] + 1;
        if (span > 0 && span <= CDECODE_TABLE_MAX_SPAN)
        {
            m->table = (int32_t *)malloc((size_t)span * sizeof(int32_t));
            if (!m->table)
                return 920;
            memset(m->table, 0xFF, (size_t)span * sizeof(int32_t));
            for (size_t i = 0; i < count; i++)
                m->table[m->values[i] - m->values[0]] = (int32_t)i;
            m->table_min = m->values[0];
            m->table_span = span;
        }
    }
    return 0;

malformed:
    free(entries);
    return 198;
}

/* Label text for a stored value: NULL when it has no label. */
static inline const char *cdecode_lookup(const cdecode_map *m, double value)
{
    if (value >= SV_missval)
    {
        uint64_t bits;
        memcpy(&bits, &value, sizeof(bits));
        uint64_t offset = bits - CDECODE_SYSMISS_BITS;
        uint64_t code = offset >> CDECODE_MISSING_STEP;
        if (code < 1 || code > 26 || (offset & ((UINT64_C(1) << CDECODE_MISSING_STEP) - 1)))
            return NULL;
        return m->missing[code];
    }
    /* Labels attach to integers only; this range also keeps the cast defined. */
    if (!(value >= -2147483648.0 && value <= 2147483647.0))
        return NULL;
    int64_t key = (int64_t)value;
    if ((double)key != value || m->count == 0)
        return NULL;
    if (m->table)
    {
        int64_t slot = key - m->table_min;
        if (slot < 0 || slot >= m->table_span || m->table[slot] < 0)
            return NULL;
        return m->text[m->table[slot]];
    }
    size_t lo = 0, hi = m->count;
    while (lo < hi)
    {
        size_t mid = lo + (hi - lo) / 2;
        if (m->values[mid] < key)
            lo = mid + 1;
        else
            hi = mid;
    }
    return (lo < m->count && m->values[lo] == key) ? m->text[lo] : NULL;
}

static ST_retcode cdecode_reject_legacy(void)
{
    SF_error("cdecode: legacy plugin protocol retired; install the matching cdecode.ado\n");
    return 198;
}

ST_retcode cdecode_main(const char *args)
{
    double t_start = ctools_timer_seconds();
    while (args && (*args == ' ' || *args == '\t'))
        args++;
    if (!args || strncmp(args, "v2", 2) != 0 || (args[2] != ' ' && args[2] != '\0'))
        return cdecode_reject_legacy();

    int maplen = 0, width = 0;
    if (ctools_parse_int_option(args, "maplen", &maplen) <= 0 || maplen < 2 ||
        ctools_parse_int_option(args, "width", &width) <= 0 || width < 1 || width > 2045)
    {
        SF_error("cdecode: invalid arguments\n");
        return 198;
    }
    int has_if = ctools_parse_bool_option(args, "ifcond");
    int touse = has_if && ctools_parse_bool_option(args, "touse");
    int fill_all = ctools_parse_bool_option(args, "fillall");

    if (SF_nvars() != 2 + touse || SF_var_is_string(1) || !SF_var_is_string(2) ||
        SF_var_is_strl(2) || (touse && SF_var_is_string(3)))
    {
        SF_error("cdecode: expected a numeric source and a fixed-width string destination\n");
        return 198;
    }

    cdecode_map map;
    memset(&map, 0, sizeof(map));
    map.map = (char *)calloc((size_t)maplen + 2, 1);
    if (!map.map)
        return 920;
    /* Local __cd_map of the calling program (SPI names locals with a leading _). */
    if (SF_macro_use("___cd_map", map.map, maplen + 1) != 0 || strlen(map.map) != (size_t)maplen)
    {
        cdecode_map_free(&map);
        SF_error("cdecode: value-label mapping not available\n");
        return 198;
    }
    int rc = cdecode_map_parse(&map, (size_t)maplen, (size_t)width);
    if (rc)
    {
        cdecode_map_free(&map);
        if (rc == 198)
            SF_error("cdecode: malformed value-label mapping\n");
        return rc;
    }

    ST_int obs1 = SF_in1(), obs2 = SF_in2();
    size_t nobs = obs2 >= obs1 ? (size_t)(obs2 - obs1) + 1 : 0;

    /* The if condition is evaluated on the calling thread, once per
     * observation; the shared loader reads a sample marker in parallel. */
    uint8_t *selected = NULL;
    if (has_if && nobs > 0)
    {
        selected = (uint8_t *)calloc(nobs, 1);
        if (!selected)
        {
            cdecode_map_free(&map);
            return 920;
        }
        if (touse)
        {
            ctools_filtered_data marked;
            stata_retcode select_rc = ctools_data_select(&marked, 0, 0, CTOOLS_LOAD_MARKER);
            for (size_t k = 0; select_rc == STATA_OK && k < marked.data.nobs; k++)
                selected[marked.obs_map[k] - (perm_idx_t)obs1] = 1;
            ctools_filtered_data_free(&marked);
            if (select_rc != STATA_OK)
            {
                free(selected);
                cdecode_map_free(&map);
                return ctools_stata_rc(select_rc);
            }
        }
        else
        {
            for (size_t i = 0; i < nobs; i++)
                selected[i] = SF_ifobs(obs1 + (ST_int)i) != 0;
        }
    }

    int error = 0;
    size_t n_written = 0;
    const ST_IIIDp read_fn = (_stata_)->vdata;
    const ST_IIIS store_fn = (_stata_)->sstore;
#pragma omp parallel for schedule(static)                                                          \
    reduction(+ : n_written) if (nobs >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < nobs; i++)
    {
        int failed;
#pragma omp atomic read
        failed = error;
        if (failed || (selected && !selected[i]))
            continue;
        double value;
        int io_rc = read_fn(1, obs1 + (ST_int)i, &value);
        const char *text = io_rc ? NULL : cdecode_lookup(&map, value);
        if (!io_rc && text && text[0] != '\0')
        {
            io_rc = store_fn(2, obs1 + (ST_int)i, (char *)text);
            n_written++;
        }
        else if (!io_rc && fill_all)
        {
            io_rc = store_fn(2, obs1 + (ST_int)i, (char *)"");
        }
        if (io_rc)
        {
#pragma omp atomic write
            error = io_rc;
        }
    }

    free(selected);
    cdecode_map_free(&map);
    if (error)
    {
        if (error != 1)
            SF_error("cdecode: failed to read or write Stata data\n");
        return error == 1 ? 1 : 459;
    }

    SF_scal_save("_cdecode_n_decoded", (double)n_written);
    SF_scal_save("_cdecode_time_total", ctools_timer_seconds() - t_start);
    CTOOLS_SAVE_THREAD_INFO("_cdecode");
    return 0;
}

ST_retcode cdecode_scan_main(const char *args)
{
    (void)args;
    return cdecode_reject_legacy();
}
