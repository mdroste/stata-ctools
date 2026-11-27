/*
    cfuzzyjoin_internal.h
    Shared types of the cfuzzyjoin modules

    cfuzzyjoin_grams.c   matchit's similarity functions: string -> grams
    cfuzzyjoin_vocab.c   gram vocabulary and per-string gram profiles
    cfuzzyjoin_lev.c     bounded Levenshtein distance (strgroup)
    cfuzzyjoin_match.c   match: similarity join of two datasets or columns
    cfuzzyjoin_link.c    link: record linkage with reclink's scores
    cfuzzyjoin_group.c   group: strgroup's grouping of similar strings
    cfuzzyjoin_impl.c    plugin calls
*/

#ifndef CFUZZYJOIN_INTERNAL_H
#define CFUZZYJOIN_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "stplugin.h"
#include "ctools_data.h"
#include "ctools_data_io.h"

#define CFJ_MODULE "cfuzzyjoin"

/* ----------------------------------------------------------------------------
 * Similarity functions (matchit's simf_*)
 * -------------------------------------------------------------------------- */

typedef enum
{
    CFJ_SIM_BIGRAM = 0,
    CFJ_SIM_NGRAM,
    CFJ_SIM_NGRAM_CIRC,
    CFJ_SIM_TOKEN,
    CFJ_SIM_COTOKEN,
    CFJ_SIM_SCOTOKEN,
    CFJ_SIM_FIRSTGRAM,
    CFJ_SIM_TKNGRAM,
    CFJ_SIM_SOUNDEX,
    CFJ_SIM_SOUNDEX_NARA,
    CFJ_SIM_SOUNDEX_FK,
    CFJ_SIM_SOUNDEX_EXT,
    CFJ_SIM_NYSIIS_FK,
    CFJ_SIM_TOKEN_SOUNDEX,
    CFJ_SIM_TOKENWRAP
} cfj_sim_kind;

typedef struct
{
    cfj_sim_kind kind;
    int n;               /* gram size of ngram, ngram_circ, firstgram, tkngram */
    int unit;            /* matchit's unitflag == 1: a gram counts once */
    cfj_sim_kind wrapped; /* tokenwrap: the function applied to each token */
} cfj_method;

/* Parse "name[,arg[,arg]]" as matchit's similmethod() passes it. 0 on
 * success; on failure *what names the problem. */
int cfj_method_parse(const char *text, cfj_method *m, const char **what);

/* The grams of one string, with repetitions, in the order matchit makes
 * them. Gram j is buf[off[j] .. off[j] + len[j]). */
typedef struct
{
    char *buf;
    size_t used, cap;
    uint32_t *off, *len;
    size_t n, ncap;
    char *work; /* scratch for phonetic keys and joined tokens */
    size_t wcap;
    /* token boundaries of the current string */
    uint32_t *tstart, *tlen;
    size_t ntok, tcap;
} cfj_grams;

void cfj_grams_init(cfj_grams *g);
void cfj_grams_free(cfj_grams *g);
/* Replace g's contents with the grams of s. 0, or -1 when out of memory. */
int cfj_extract(const cfj_method *m, const char *s, size_t len, cfj_grams *g);

/* Stata's soundex() and soundex_nara(): out holds at least 5 bytes; returns
 * the key length (0 when s has no ASCII letter). */
size_t cfj_soundex(const char *s, size_t len, int nara, char *out);

/* ----------------------------------------------------------------------------
 * Vocabulary and profiles
 * -------------------------------------------------------------------------- */

/* Distinct grams, numbered in order of first appearance over the string sets
 * in the order given, which does not depend on the thread count. */
typedef struct
{
    size_t n;
    char *bytes;    /* gram text, back to back */
    uint64_t *off;  /* n + 1 offsets into bytes */
} cfj_vocab;

/* Each string's distinct grams, ascending by id, and their counts */
typedef struct
{
    size_t n;
    uint64_t *off; /* n + 1 */
    uint32_t *id;
    uint32_t *cnt;
} cfj_profile;

typedef struct
{
    const char *const *s; /* strings (NULL counts as "") */
    size_t n;
} cfj_strset;

void cfj_vocab_free(cfj_vocab *v);
void cfj_profile_free(cfj_profile *p);
/* Build the vocabulary of nsets string sets and one profile per set. */
stata_retcode cfj_profiles_build(const cfj_method *m, const cfj_strset *sets, int nsets,
                                 cfj_vocab *vocab, cfj_profile *profiles, int threads);
/* Look a gram up in a vocabulary: its id, or -1. Builds an index on first
 * use, which cfj_vocab_free releases. */
typedef struct
{
    size_t mask;
    uint32_t *slot; /* id + 1, 0 empty */
} cfj_vocab_index;
stata_retcode cfj_vocab_index_build(const cfj_vocab *v, cfj_vocab_index *ix);
long long cfj_vocab_find(const cfj_vocab *v, const cfj_vocab_index *ix, const char *s, size_t len);
void cfj_vocab_index_free(cfj_vocab_index *ix);

uint64_t cfj_hash_bytes(const char *s, size_t len);

/* ----------------------------------------------------------------------------
 * Levenshtein distance
 * -------------------------------------------------------------------------- */

/* The Levenshtein distance of a and b (unit costs, bytes) when it is at most
 * k; otherwise some value above k. work holds 2 * (max length + 2) uint32. */
uint32_t cfj_lev_bounded(const unsigned char *a, uint32_t la, const unsigned char *b, uint32_t lb,
                         uint32_t k, uint32_t *work);
/* The exact distance (for tests) */
uint32_t cfj_lev_exact(const unsigned char *a, uint32_t la, const unsigned char *b, uint32_t lb,
                       uint32_t *work);

/* ----------------------------------------------------------------------------
 * Sub-commands
 * -------------------------------------------------------------------------- */

ST_retcode cfj_match_master(const char *args);
ST_retcode cfj_match_weights(const char *args, int after_master);
ST_retcode cfj_match_using(const char *args);
ST_retcode cfj_match_store(const char *args);
ST_retcode cfj_match_columns(const char *args);
void cfj_match_cleanup(void);

ST_retcode cfj_link_master(const char *args);
ST_retcode cfj_link_exclude(const char *args);
ST_retcode cfj_link_using(const char *args);
ST_retcode cfj_link_layout(const char *args);
ST_retcode cfj_link_store(const char *args);
void cfj_link_cleanup(void);

ST_retcode cfj_group(const char *args);

/* ----------------------------------------------------------------------------
 * Helpers (cfuzzyjoin_impl.c)
 * -------------------------------------------------------------------------- */

/* Load plugin variables first .. first + count - 1 (1-based) of every
 * observation in range, or of those the marker selects (CTOOLS_LOAD_MARKER:
 * the last plugin variable). */
stata_retcode cfj_load(ctools_filtered_data *fd, int first, int count, int flags);
/* One display line (SF_display must not see a newline before the end) */
void cfj_say(const char *fmt, ...);
/* Text of option name(...) in args, or "" */
int cfj_opt_text(const char *args, const char *name, char *out, size_t size);
int cfj_opt_double(const char *args, const char *name, double *out);
int cfj_opt_int(const char *args, const char *name, int *out);
int cfj_flag(const char *args, const char *word);
/* Numeric stores and reads through Stata's checked callbacks: the unchecked
 * ones SF_vstore and SF_vdata map to (SD_FASTMODE) assume double storage and
 * miss or corrupt other numeric types (see ctools_transport_internal.h) */
#define CFJ_VSTORE(i, j, v) ((_stata_)->safestore((i), (j), (v)))
#define CFJ_VDATA(i, j, d) ((_stata_)->safevdata((i), (j), (d)))

/* Stata's missing-value test for a loaded double */
static inline int cfj_missing(double x)
{
    return !(x < 8.988465674311579e307) || x != x;
}

#endif /* CFUZZYJOIN_INTERNAL_H */
