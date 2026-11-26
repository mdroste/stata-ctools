/*
 * creshape_schema.c - the description of the reshaped dataset.
 *
 * The SPI moves values only. Storage types, display formats, variable and
 * value labels, characteristics and notes travel in a .dta file without
 * observations: the ado file saves one of the dataset as it is, this file
 * reads it, builds the description of the result and writes it for the ado
 * file to load (formats/ctools_dta_schema.h).
 *
 * What each variable of the result inherits follows reshape (measured
 * against StataNow 19.5, the default favor(speed) path):
 *
 * long   A stub's variable takes the widest storage type of its wide
 *        variables (byte < int < long < float < double; the longest str#),
 *        the format and value label of the last one in order of j, no
 *        variable label, every characteristic (a later variable's wins) and
 *        the notes of all of them in order. j gets the type compress would
 *        choose. The labels reshape wide set aside in _dta[__JValLab...],
 *        _dta[__JVarLab] and _dta[__XijVarLab...] are put back.
 * wide   Each new variable is a copy of the long variable's description,
 *        labelled "<j> <name>". The labels of j and of the long variables
 *        are set aside in the characteristics above.
 *
 * Both directions record _dta[ReS_Xij_n], _dta[ReS_Xij_wide#] and
 * _dta[ReS_Xij_long#] as reshape's summary does, and mark the result sorted
 * by i (and j): the store writes the observations in that order.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "creshape_internal.h"

/* Characteristics of the input grouped by owner: head[col] and next[] chain
 * the positions in in->chars, in file order; column nvars is _dta. */
typedef struct
{
    long *head, *next;
} schema_chains;

static int schema_chain(const creshape_state *st, const ctools_dta_schema *in, schema_chains *ch)
{
    ch->head = (long *)malloc((st->nvars + 1) * sizeof(long));
    ch->next = (long *)malloc((in->nchars ? in->nchars : 1) * sizeof(long));
    long *tail = (long *)malloc((st->nvars + 1) * sizeof(long));
    if (!ch->head || !ch->next || !tail)
    {
        free(tail);
        return 920;
    }
    for (size_t c = 0; c <= st->nvars; c++)
        ch->head[c] = tail[c] = -1;
    for (size_t i = 0; i < in->nchars; i++)
    {
        long col = !strcmp(in->chars[i].owner, "_dta")
                       ? (long)st->nvars
                       : creshape_names_find(st, in->chars[i].owner);
        ch->next[i] = -1;
        if (col < 0)
            continue;
        if (tail[col] < 0)
            ch->head[col] = (long)i;
        else
            ch->next[tail[col]] = (long)i;
        tail[col] = (long)i;
    }
    free(tail);
    return 0;
}

/* note1, note2, ...: the number, or 0 for another name (note0 is the count) */
static long schema_note_number(const char *name)
{
    if (strncmp(name, "note", 4) || name[4] < '1' || name[4] > '9')
        return 0;
    const char *p = name + 4;
    while (*p >= '0' && *p <= '9')
        p++;
    return *p ? 0 : atol(name + 4);
}

/* Copy the characteristics of input column col to owner. */
static int schema_copy_chars(const ctools_dta_schema *in, const schema_chains *ch, size_t col,
                             ctools_dta_schema *out, const char *owner)
{
    for (long i = ch->head[col]; i >= 0; i = ch->next[i])
        if (ctools_dta_schema_char_add(out, owner, in->chars[i].name, in->chars[i].text))
            return 920;
    return 0;
}

typedef struct
{
    long number;
    const char *text;
} schema_note;

static int schema_cmp_note(const void *a, const void *b)
{
    long x = ((const schema_note *)a)->number, y = ((const schema_note *)b)->number;
    return (x > y) - (x < y);
}

/* Long: the characteristics of a stub's wide variables, merged. A later
 * variable's characteristic replaces an earlier one of the same name. Notes
 * are gathered in order of the variables and of their numbers, each text
 * once, and numbered again from 1 (measured against reshape). */
static int schema_merge_chars(const creshape_state *st, const ctools_dta_schema *in,
                              const schema_chains *ch, const creshape_stub *sb,
                              ctools_dta_schema *out)
{
    /* collected in a list of its own, so a later variable's value replaces
     * an earlier one without searching every characteristic of the result */
    ctools_dta_schema merged;
    ctools_dta_schema_init(&merged);
    schema_note *notes = NULL; /* the notes of one variable */
    const char **kept = NULL;  /* the distinct note texts so far */
    size_t nnotes = 0, cap = 0, nkept = 0, kept_cap = 0;
    int rc = 0;
    for (size_t l = 0; l < st->nlevels && !rc; l++)
    {
        const int col = sb->cols[l];
        if (col < 0)
            continue;
        nnotes = 0;
        for (long i = ch->head[col]; i >= 0 && !rc; i = ch->next[i])
        {
            const ctools_dta_schema_char *c = &in->chars[i];
            const long number = schema_note_number(c->name);
            if (!strcmp(c->name, "note0"))
                continue;
            if (!number)
            {
                rc = ctools_dta_schema_char_set(&merged, sb->long_name, c->name, c->text);
                continue;
            }
            if (nnotes == cap)
            {
                size_t grown_cap = cap ? cap * 2 : 8;
                schema_note *grown = (schema_note *)realloc(notes, grown_cap * sizeof(*grown));
                if (!grown)
                {
                    rc = 920;
                    break;
                }
                notes = grown;
                cap = grown_cap;
            }
            notes[nnotes++] = (schema_note){number, c->text};
        }
        if (nnotes)
            qsort(notes, nnotes, sizeof(*notes), schema_cmp_note);
        for (size_t k = 0; k < nnotes && !rc; k++)
        {
            int seen = 0;
            for (size_t m = 0; m < nkept && !seen; m++)
                seen = !strcmp(kept[m], notes[k].text);
            if (seen)
                continue;
            if (nkept == kept_cap)
            {
                size_t grown_cap = kept_cap ? kept_cap * 2 : 8;
                const char **grown = (const char **)realloc(kept, grown_cap * sizeof(*grown));
                if (!grown)
                {
                    rc = 920;
                    break;
                }
                kept = grown;
                kept_cap = grown_cap;
            }
            kept[nkept++] = notes[k].text;
        }
    }
    for (size_t k = 0; k < nkept && !rc; k++)
    {
        char name[32];
        snprintf(name, sizeof(name), "note%zu", k + 1);
        rc = ctools_dta_schema_char_add(&merged, sb->long_name, name, kept[k]);
    }
    if (!rc && nkept)
    {
        char count[32];
        snprintf(count, sizeof(count), "%zu", nkept);
        rc = ctools_dta_schema_char_add(&merged, sb->long_name, "note0", count);
    }
    for (size_t i = 0; i < merged.nchars && !rc; i++)
        rc = ctools_dta_schema_char_add(out, sb->long_name, merged.chars[i].name,
                                        merged.chars[i].text);
    free(notes);
    free(kept);
    ctools_dta_schema_free(&merged);
    return rc;
}

static void schema_set_text(char *dst, size_t cap, const char *text)
{
    snprintf(dst, cap, "%s", text);
}

/* The wider of two storage types: numeric codes fall as the type widens,
 * str# codes are the width. */
static uint16_t schema_promote(uint16_t a, uint16_t b)
{
    const int a_str = a <= CTOOLS_DTA_CODE_STR_MAX, b_str = b <= CTOOLS_DTA_CODE_STR_MAX;
    if (a_str != b_str)
        return a;
    if (a_str)
        return a > b ? a : b;
    return a < b ? a : b;
}

/* ------------------------------------------------------------------------
 * The labels reshape carries across a round trip, in _dta characteristics
 * ------------------------------------------------------------------------ */

/* Wide: remember the labels of j and of the long variables. */
static int schema_save_labels(const creshape_state *st, const ctools_dta_schema *in,
                              const ctools_dta_schema *labels, ctools_dta_schema *dta)
{
    const ctools_dta_schema_var *j = &in->vars[st->jcol];
    int rc = 0;
    if (j->value_label[0])
    {
        rc = ctools_dta_schema_char_set(dta, "_dta", "__JValLabName", j->value_label);
        const long at = ctools_dta_schema_label_find(labels, j->value_label);
        size_t n = 0;
        if (!rc && at >= 0 && !ctools_dta_schema_label_count(labels, (size_t)at, &n) && n)
        {
            /* value `"text"' pairs, as reshape writes them; the file keeps
             * the entries by value, the missing codes last */
            size_t cap = 64, len = 0;
            char *text = (char *)malloc(cap);
            if (!text)
                return 920;
            text[0] = '\0';
            for (size_t k = 0; k < n; k++)
            {
                int32_t value;
                const char *label;
                if (ctools_dta_schema_label_entry(labels, (size_t)at, k, &value, &label))
                    break;
                if (value > 2147483620 || !label[0])
                    continue;
                const size_t need = strlen(label) + 32;
                if (len + need > cap)
                {
                    while (len + need > cap)
                        cap *= 2;
                    char *grown = (char *)realloc(text, cap);
                    if (!grown)
                    {
                        free(text);
                        return 920;
                    }
                    text = grown;
                }
                len += (size_t)snprintf(text + len, cap - len, " %d `\"%s\"' ", (int)value, label);
            }
            rc = ctools_dta_schema_char_set(dta, "_dta", "__JValLab", text);
            free(text);
        }
    }
    if (!rc && j->label[0])
        rc = ctools_dta_schema_char_set(dta, "_dta", "__JVarLab", j->label);
    size_t numbered = 0;
    for (size_t k = 0; k < st->nstubs && !rc; k++)
    {
        const creshape_stub *sb = &st->stubs[k];
        const char *label = in->vars[sb->long_col].label;
        if (!label[0])
            continue;
        char name[CTOOLS_DTA_NAME_BYTES + 16];
        if (strlen(sb->long_name) < 21)
        {
            snprintf(name, sizeof(name), "__XijVarLab%s", sb->long_name);
            rc = ctools_dta_schema_char_set(dta, "_dta", name, label);
        }
        else
        {
            /* a name too long for a characteristic's name: numbered, with
             * the name in front of the label */
            char count[32];
            char *text = creshape_join(sb->long_name, " ", label);
            if (!text)
                return 920;
            snprintf(name, sizeof(name), "__XijVarLab%zu", ++numbered);
            snprintf(count, sizeof(count), "%zu", numbered);
            rc = ctools_dta_schema_char_set(dta, "_dta", name, text);
            if (!rc)
                rc = ctools_dta_schema_char_set(dta, "_dta", "__XijVarLabTotal", count);
            free(text);
        }
    }
    return rc;
}

/* Read ` 1 `"one"'  2 `"two"' ` into values and texts (texts point into a
 * copy the caller frees through *buffer). Returns the number of pairs, or -1
 * when the text is not such a list. */
static long schema_parse_pairs(const char *text, int32_t **values_out, const char ***texts_out,
                               char **buffer)
{
    char *copy = creshape_strdup(text);
    size_t n = 0, cap = 16;
    int32_t *values = (int32_t *)malloc(cap * sizeof(*values));
    const char **texts = (const char **)malloc(cap * sizeof(*texts));
    long result = -1;
    if (!copy || !values || !texts)
        goto done;
    char *p = copy;
    for (;;)
    {
        while (*p == ' ')
            p++;
        if (!*p)
            break;
        char *end;
        long value = strtol(p, &end, 10);
        if (end == p)
            goto done;
        p = end;
        while (*p == ' ')
            p++;
        if (p[0] != '`' || p[1] != '"')
            goto done;
        p += 2;
        char *start = p;
        int depth = 1;
        while (*p)
        {
            if (p[0] == '`' && p[1] == '"')
            {
                depth++;
                p += 2;
            }
            else if (p[0] == '"' && p[1] == '\'')
            {
                if (--depth == 0)
                    break;
                p += 2;
            }
            else
                p++;
        }
        if (!*p)
            goto done;
        *p = '\0';
        p += 2;
        if (n == cap)
        {
            cap *= 2;
            int32_t *v = (int32_t *)realloc(values, cap * sizeof(*values));
            if (v)
                values = v;
            const char **t = (const char **)realloc(texts, cap * sizeof(*texts));
            if (t)
                texts = t;
            if (!v || !t)
                goto done;
        }
        values[n] = (int32_t)value;
        texts[n++] = start;
    }
    result = (long)n;
done:
    if (result < 0)
    {
        free(copy);
        free(values);
        free(texts);
        return -1;
    }
    *values_out = values;
    *texts_out = texts;
    *buffer = copy;
    return result;
}

/* Long: put the remembered labels back and forget them. */
static int schema_apply_labels(const creshape_state *st, ctools_dta_schema *dta,
                               ctools_dta_schema *out)
{
    ctools_dta_schema_var *j = &out->vars[st->ni];
    const char *name = ctools_dta_schema_char_get(dta, "_dta", "__JValLabName");
    int rc = 0;
    if (name && !st->string_j)
    {
        const char *pairs = ctools_dta_schema_char_get(dta, "_dta", "__JValLab");
        if (pairs && ctools_dta_schema_label_find(out, name) < 0)
        {
            int32_t *values = NULL;
            const char **texts = NULL;
            char *buffer = NULL;
            long n = schema_parse_pairs(pairs, &values, &texts, &buffer);
            if (n > 0)
                rc = ctools_dta_schema_label_define(out, name, (size_t)n, values, texts);
            if (n >= 0)
            {
                free(values);
                free(texts);
                free(buffer);
            }
            /* a list that cannot be defined leaves the name attached, as a
             * failed label define does in reshape */
            if (rc == 198 || rc == 110)
                rc = 0;
        }
        schema_set_text(j->value_label, sizeof(j->value_label), name);
        /* label values widens a display format to the longest label (in
         * bytes), and reshape attaches the label with it */
        const long at = ctools_dta_schema_label_find(out, name);
        size_t count = 0, widest = 0;
        if (at >= 0 && !ctools_dta_schema_label_count(out, (size_t)at, &count))
            for (size_t k = 0; k < count; k++)
            {
                int32_t value;
                const char *text;
                if (!ctools_dta_schema_label_entry(out, (size_t)at, k, &value, &text) &&
                    strlen(text) > widest)
                    widest = strlen(text);
            }
        if (widest > (size_t)atoi(j->format + 1) && widest <= CTOOLS_DTA_CODE_STR_MAX)
            snprintf(j->format, sizeof(j->format), "%%%zu.0g", widest);
        ctools_dta_schema_char_remove(dta, "_dta", "__JValLab");
        ctools_dta_schema_char_remove(dta, "_dta", "__JValLabName");
    }
    const char *label = ctools_dta_schema_char_get(dta, "_dta", "__JVarLab");
    if (!rc && label)
    {
        schema_set_text(j->label, sizeof(j->label), label);
        ctools_dta_schema_char_remove(dta, "_dta", "__JVarLab");
    }
    const char *total_text = ctools_dta_schema_char_get(dta, "_dta", "__XijVarLabTotal");
    const long total = total_text ? atol(total_text) : -1;
    for (size_t k = 0; k < st->nstubs && !rc; k++)
    {
        const creshape_stub *sb = &st->stubs[k];
        ctools_dta_schema_var *v = &out->vars[st->ni + 1 + k];
        char cname[CTOOLS_DTA_NAME_BYTES + 16];
        if (strlen(sb->long_name) < 21)
        {
            snprintf(cname, sizeof(cname), "__XijVarLab%s", sb->long_name);
            label = ctools_dta_schema_char_get(dta, "_dta", cname);
            if (label)
            {
                schema_set_text(v->label, sizeof(v->label), label);
                ctools_dta_schema_char_remove(dta, "_dta", cname);
            }
            continue;
        }
        const size_t len = strlen(sb->long_name);
        for (long i = 1; i <= total; i++)
        {
            snprintf(cname, sizeof(cname), "__XijVarLab%ld", i);
            label = ctools_dta_schema_char_get(dta, "_dta", cname);
            if (label && !strncmp(label, sb->long_name, len) && label[len] == ' ')
            {
                schema_set_text(v->label, sizeof(v->label), label + len + 1);
                ctools_dta_schema_char_remove(dta, "_dta", cname);
                break;
            }
        }
    }
    return rc;
}

/* _dta[ReS_Xij_n], _dta[ReS_Xij_wide#], _dta[ReS_Xij_long#] */
static int schema_report_chars(const creshape_state *st, ctools_dta_schema *dta)
{
    int rc = 0;
    for (size_t k = 0; k < st->nstubs && !rc; k++)
    {
        const creshape_stub *sb = &st->stubs[k];
        size_t need = 1;
        for (size_t r = 0; r < st->nreport; r++)
            need += strlen(sb->names[st->report_levels[r]]) + 1;
        char *list = (char *)malloc(need);
        if (!list)
            return 920;
        size_t at = 0;
        list[0] = '\0';
        for (size_t r = 0; r < st->nreport; r++)
            at += (size_t)snprintf(list + at, need - at, "%s%s", r ? " " : "",
                                   sb->names[st->report_levels[r]]);
        char name[48];
        snprintf(name, sizeof(name), "ReS_Xij_wide%zu", k + 1);
        rc = ctools_dta_schema_char_add(dta, "_dta", name, list);
        free(list);
        snprintf(name, sizeof(name), "ReS_Xij_long%zu", k + 1);
        if (!rc)
            rc = ctools_dta_schema_char_add(dta, "_dta", name, sb->long_name);
    }
    char count[32];
    snprintf(count, sizeof(count), "%zu", st->nstubs);
    if (!rc)
        rc = ctools_dta_schema_char_add(dta, "_dta", "ReS_Xij_n", count);
    return rc;
}

/* ------------------------------------------------------------------------
 * The result
 * ------------------------------------------------------------------------ */

int creshape_schema_write(creshape_state *st, const char *template_path, const char *out_path,
                          char *msg, size_t cap)
{
    ctools_dta_schema in, out, dta;
    schema_chains ch = {NULL, NULL};
    ctools_dta_schema_init(&out);
    ctools_dta_schema_init(&dta);
    int rc = ctools_dta_schema_read(template_path, &in, msg, cap);
    if (rc)
        return rc;

    /* the template must describe the dataset the plan saw */
    if (in.nvars != st->nvars)
        rc = 459;
    for (size_t c = 0; c < st->nvars && !rc; c++)
        if (strcmp(in.vars[c].name, st->names[c]) ||
            (in.vars[c].type <= CTOOLS_DTA_CODE_STR_MAX) != (st->is_str[c] != 0))
            rc = 459;
    if (rc)
    {
        snprintf(msg, cap, "the dataset changed while it was being reshaped");
        goto done;
    }
    if ((rc = schema_chain(st, &in, &ch)))
        goto oom;

    out.vars =
        (ctools_dta_schema_var *)calloc(st->out_nvars ? st->out_nvars : 1, sizeof(*out.vars));
    out.sort = (uint32_t *)calloc(st->ni + 1, sizeof(*out.sort));
    if (!out.vars || !out.sort)
    {
        rc = 920;
        goto oom;
    }
    out.nvars = st->out_nvars;
    memcpy(out.data_label, in.data_label, sizeof(out.data_label));

    for (size_t v = 0; v < st->out_nvars && !rc; v++)
    {
        const creshape_out *o = &st->out[v];
        ctools_dta_schema_var *ov = &out.vars[v];
        if (o->kind == CRESHAPE_OUT_COPY)
        {
            *ov = in.vars[o->col];
            rc = schema_copy_chars(&in, &ch, (size_t)o->col, &out, ov->name);
        }
        else if (o->kind == CRESHAPE_OUT_J)
        {
            ov->type = st->j_type;
            schema_set_text(ov->name, sizeof(ov->name), o->name);
            ctools_dta_schema_default_format(ov->type, ov->format);
        }
        else if (st->dir == CRESHAPE_LONG)
        {
            const creshape_stub *sb = &st->stubs[o->col];
            int have = 0;
            for (size_t l = 0; l < st->nlevels; l++)
            {
                if (sb->cols[l] < 0)
                    continue;
                const ctools_dta_schema_var *iv = &in.vars[sb->cols[l]];
                ov->type = have ? schema_promote(ov->type, iv->type) : iv->type;
                memcpy(ov->format, iv->format, sizeof(ov->format));
                memcpy(ov->value_label, iv->value_label, sizeof(ov->value_label));
                have = 1;
            }
            if (!have)
            {
                /* none of the stub's variables exists: all missing, a float */
                ov->type = CTOOLS_DTA_CODE_FLOAT;
                ctools_dta_schema_default_format(ov->type, ov->format);
            }
            schema_set_text(ov->name, sizeof(ov->name), o->name);
            rc = schema_merge_chars(st, &in, &ch, sb, &out);
        }
        else
        {
            const creshape_stub *sb = &st->stubs[o->col];
            *ov = in.vars[sb->long_col];
            schema_set_text(ov->name, sizeof(ov->name), o->name);
            snprintf(ov->label, sizeof(ov->label), "%.160s %.128s", st->levels[o->level].text,
                     sb->long_name);
            rc = schema_copy_chars(&in, &ch, (size_t)sb->long_col, &out, ov->name);
        }
    }
    if (rc)
        goto oom;

    /* the dataset's own characteristics, without the previous summary */
    for (long i = ch.head[st->nvars]; i >= 0 && !rc; i = ch.next[i])
    {
        const char *name = in.chars[i].name;
        if (!strcmp(name, "ReS_Xij_n") || !strncmp(name, "ReS_Xij_wide", 12) ||
            !strncmp(name, "ReS_Xij_long", 12))
            continue;
        rc = ctools_dta_schema_char_add(&dta, "_dta", name, in.chars[i].text);
    }
    /* value labels move over whole */
    out.labels = in.labels;
    out.nlabels = in.nlabels;
    out.labels_cap = in.labels_cap;
    in.labels = NULL;
    in.nlabels = in.labels_cap = 0;

    if (!rc)
        rc = st->dir == CRESHAPE_WIDE ? schema_save_labels(st, &in, &out, &dta)
                                      : schema_apply_labels(st, &dta, &out);
    if (!rc)
        rc = schema_report_chars(st, &dta);
    for (size_t i = 0; i < dta.nchars && !rc; i++)
        rc = ctools_dta_schema_char_add(&out, "_dta", dta.chars[i].name, dta.chars[i].text);
    if (rc)
        goto oom;

    /* sorted by i, and by j in long form: the order the store writes */
    for (size_t k = 0; k < st->ni; k++)
        out.sort[out.nsort++] = (uint32_t)k;
    if (st->dir == CRESHAPE_LONG)
        out.sort[out.nsort++] = (uint32_t)st->ni;

    rc = ctools_dta_schema_write(out_path, &out, msg, cap);
    goto done;

oom:
    if (rc == 920)
        snprintf(msg, cap, "out of memory");
    else if (!msg[0])
        snprintf(msg, cap, "cannot describe the result");
done:
    free(ch.head);
    free(ch.next);
    ctools_dta_schema_free(&in);
    ctools_dta_schema_free(&out);
    ctools_dta_schema_free(&dta);
    return rc;
}
