/* cfuzzyjoin's string kernels against references.
 *
 * Levenshtein: the bit-parallel and banded distances against the plain
 * dynamic program on random strings, exact and with every bound.
 * Soundex: Stata's soundex() and soundex_nara() on strings whose codes were
 * read from StataNow 19.5. Grams: Mata's tokens() rules and matchit's gram
 * functions on examples from matchit's help file and from Mata. Profiles:
 * gram ids in order of first appearance, counts, the unit flag, and the same
 * vocabulary for any task count.
 */
#undef NDEBUG
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"
#include "cfuzzyjoin/cfuzzyjoin_internal.h"
#include "cfuzzyjoin/cfuzzyjoin_lev.h"

#include "cfuzzyjoin/cfuzzyjoin_lev.c"
#include "cfuzzyjoin/cfuzzyjoin_grams.c"
#include "cfuzzyjoin/cfuzzyjoin_vocab.c"

static uint64_t rng = 0x9E3779B97F4A7C15ull;
static uint32_t next_u32(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return (uint32_t)(rng >> 16);
}

static uint32_t plain_lev(const unsigned char *a, uint32_t la, const unsigned char *b, uint32_t lb)
{
    uint32_t *prev = malloc((lb + 1) * sizeof(uint32_t)), *cur = malloc((lb + 1) * sizeof(uint32_t));
    assert(prev && cur);
    for (uint32_t j = 0; j <= lb; j++)
        prev[j] = j;
    for (uint32_t i = 1; i <= la; i++)
    {
        cur[0] = i;
        for (uint32_t j = 1; j <= lb; j++)
        {
            uint32_t v = prev[j - 1] + (a[i - 1] != b[j - 1]);
            if (prev[j] + 1 < v)
                v = prev[j] + 1;
            if (cur[j - 1] + 1 < v)
                v = cur[j - 1] + 1;
            cur[j] = v;
        }
        uint32_t *t = prev;
        prev = cur;
        cur = t;
    }
    uint32_t d = prev[lb];
    free(prev);
    free(cur);
    return d;
}

static void random_string(unsigned char *s, uint32_t n, uint32_t alphabet)
{
    for (uint32_t i = 0; i < n; i++)
        s[i] = (unsigned char)('a' + next_u32() % alphabet);
}

static void test_levenshtein(void)
{
    unsigned char a[300], b[300];
    uint32_t work[2 * 302];
    for (int trial = 0; trial < 20000; trial++)
    {
        uint32_t la = next_u32() % (trial < 10000 ? 70 : 260);
        uint32_t lb = next_u32() % (trial < 10000 ? 70 : 260);
        uint32_t alphabet = 2 + next_u32() % 6;
        random_string(a, la, alphabet);
        /* b: an edited copy of a, or a fresh string */
        if (next_u32() % 2)
        {
            lb = 0;
            for (uint32_t i = 0; i < la && lb < 299; i++)
            {
                uint32_t r = next_u32() % 10;
                if (r == 0)
                    continue;
                if (r == 1)
                    b[lb++] = (unsigned char)('a' + next_u32() % alphabet);
                b[lb++] = r == 2 ? (unsigned char)('a' + next_u32() % alphabet) : a[i];
            }
        }
        else
            random_string(b, lb, alphabet);
        uint32_t d = plain_lev(a, la, b, lb);
        assert(cfj_lev_exact(a, la, b, lb, work) == d);
        for (uint32_t k = 0; k <= d + 2; k++)
        {
            uint32_t got = cfj_lev_bounded(a, la, b, lb, k, work);
            if (d <= k)
                assert(got == d);
            else
                assert(got > k);
        }
        if (la <= 64)
        {
            cfj_myers p;
            cfj_myers_init(&p, a, la);
            assert(cfj_myers_dist(&p, b, lb, la > lb ? la : lb) == d);
            for (uint32_t k = 0; k <= d + 1; k++)
            {
                uint32_t got = cfj_myers_dist(&p, b, lb, k);
                assert(d <= k ? got == d : got > k);
            }
        }
    }
    printf("levenshtein ok\n");
}

static void sx(const char *s, const char *code, const char *nara)
{
    char out[8];
    cfj_soundex(s, strlen(s), 0, out);
    if (strcmp(out, code))
    {
        printf("soundex(%s) = %s, expected %s\n", s, out, code);
        assert(0);
    }
    cfj_soundex(s, strlen(s), 1, out);
    if (strcmp(out, nara))
    {
        printf("soundex_nara(%s) = %s, expected %s\n", s, out, nara);
        assert(0);
    }
}

static void test_soundex(void)
{
    /* read from Mata's soundex() and soundex_nara() */
    sx("", "", "");
    sx(" ", "", "");
    sx("123", "", "");
    sx("a", "A000", "A000");
    sx("A1b", "A100", "A100");
    sx("Ashcraft", "A226", "A261");
    sx("Tymczak", "T522", "T522");
    sx("Pfister", "P236", "P236");
    sx("Honeyman", "H555", "H555");
    sx("Lloyd", "L300", "L300");
    sx("  smith", "S530", "S530");
    sx("Smith, John", "S532", "S532");
    sx("1abc", "A120", "A120");
    sx("bhb", "B100", "B000");
    sx("bwb", "B100", "B000");
    sx("b1b", "B000", "B000");
    sx("b b", "B000", "B000");
    sx("wbwb", "W110", "W100");
    sx("Bbhb", "B100", "B000");
    sx("Bhab", "B100", "B100");
    sx("Abhb", "A110", "A100");
    sx("Ahxh", "A200", "A200");
    sx("bdbdbdbd", "B313", "B313");
    sx("\xc3\xa9lan", "L500", "L500");
    sx("Gutierrez", "G362", "G362");
    sx("Jackson", "J250", "J250");
    printf("soundex ok\n");
}

/* The grams of s under method m, joined with '|' */
static void grams(const char *spec, const char *s, const char *expected)
{
    cfj_method m;
    const char *what;
    assert(cfj_method_parse(spec, &m, &what) == 0);
    cfj_grams g;
    cfj_grams_init(&g);
    assert(cfj_extract(&m, s, strlen(s), &g) == 0);
    char out[512];
    size_t n = 0;
    for (size_t j = 0; j < g.n; j++)
    {
        if (j)
            out[n++] = '|';
        memcpy(out + n, g.buf + g.off[j], g.len[j]);
        n += g.len[j];
    }
    out[n] = '\0';
    if (strcmp(out, expected))
    {
        printf("%s(%s) = [%s], expected [%s]\n", spec, s, out, expected);
        assert(0);
    }
    cfj_grams_free(&g);
}

static void test_grams(void)
{
    /* Mata's tokens(), as observed */
    grams("token", "  john   smith  ", "john|smith");
    grams("token", "john\tsmith", "john\tsmith");
    grams("token", "john \"big jim\" smith", "john|big jim|smith");
    grams("token", "a\"b c\"d", "a\"b|c\"d");
    grams("token", "ab \"cd\"ef gh", "ab|cd|ef|gh");
    grams("token", "say `\"hi there\"' now", "say|hi there|now");
    grams("token", "unterminated \"quote here", "unterminated|quote here");
    grams("token", "\"", "");
    grams("token", " \"\" ", "");
    grams("token", "", "");
    /* matchit's help file: "John Smith" */
    grams("bigram", "John Smith", "Jo|oh|hn|n | S|Sm|mi|it|th");
    grams("ngram, 3", "John Smith", "Joh|ohn|hn |n S| Sm|Smi|mit|ith");
    grams("ngram_circ, 3", "John Smith", "Joh|ohn|hn |n S| Sm|Smi|mit|ith|th |h J| Jo");
    grams("token_soundex", "John Smith", "J500|S530");
    grams("soundex", "John Smith", "J525");
    grams("soundex_fk", "John Smith", "J5253");
    grams("soundex_ext", "John Smith", "J5|J52|J525|J5253");
    grams("nysiis_fk", "John Smith", "JANSNAT");
    grams("tokenwrap, \"nysiis_fk\"", "John Smith", "JAN|SNAT");
    /* matchit's help file shows S532255; its code (and matchit, in the Stata
     * suite) gives S53255: M5 T3 S2 N5 N5 */
    grams("tokenwrap, \"soundex_fk\"", "Johnatan Smithsonian", "J535|S53255");
    grams("firstgram, 3", "John Smith", "Joh|Smi");
    /* short strings are one gram; too few windows give the whole string */
    grams("bigram", "", "");
    grams("bigram", "a", "a");
    grams("bigram", "ab", "ab");
    grams("ngram, 3", "ab", "ab");
    grams("tkngram, 3", "ab", "ab");
    grams("cotoken", "b a c", "a b|a c");
    grams("scotoken", "b a c", "b|a|a b|a|c|a c");
    grams("scotoken", "solo", "solo");
    /* invalid specifications */
    cfj_method m;
    const char *what;
    assert(cfj_method_parse("nosuch", &m, &what) != 0);
    assert(cfj_method_parse("ngram", &m, &what) != 0);
    assert(cfj_method_parse("ngram, 0", &m, &what) != 0);
    assert(cfj_method_parse("tokenwrap, \"nosuch\"", &m, &what) != 0);
    assert(cfj_method_parse("bigram, 1", &m, &what) == 0 && m.unit == 1);
    assert(cfj_method_parse("ngram, 3, 1", &m, &what) == 0 && m.n == 3 && m.unit == 1);
    printf("grams ok\n");
}

static void test_profiles(void)
{
    const char *a[] = {"abab", "ba", "", "abc"};
    const char *b[] = {"cab", "abab"};
    cfj_strset sets[2] = {{a, 4}, {b, 2}};
    cfj_method m;
    const char *what;
    assert(cfj_method_parse("bigram", &m, &what) == 0);
    cfj_vocab v1, v2;
    cfj_profile p1[2], p2[2];
    assert(cfj_profiles_build(&m, sets, 2, &v1, p1, 1) == STATA_OK);
    assert(cfj_profiles_build(&m, sets, 2, &v2, p2, 7) == STATA_OK);
    /* first appearance: ab ba "" bc ca */
    const char *order[] = {"ab", "ba", "", "bc", "ca"};
    assert(v1.n == 5 && v2.n == 5);
    for (size_t g = 0; g < 5; g++)
    {
        size_t l = (size_t)(v1.off[g + 1] - v1.off[g]);
        assert(l == strlen(order[g]) && !memcmp(v1.bytes + v1.off[g], order[g], l));
        assert(v2.off[g + 1] == v1.off[g + 1] && !memcmp(v2.bytes, v1.bytes, (size_t)v1.off[5]));
    }
    /* "abab": ab x2, ba x1 */
    assert(p1[0].off[1] == 2 && p1[0].id[0] == 0 && p1[0].cnt[0] == 2 && p1[0].id[1] == 1 &&
           p1[0].cnt[1] == 1);
    for (int s = 0; s < 2; s++)
    {
        assert(p1[s].off[p1[s].n] == p2[s].off[p2[s].n]);
        size_t np = (size_t)p1[s].off[p1[s].n];
        assert(!memcmp(p1[s].id, p2[s].id, np * sizeof(uint32_t)));
        assert(!memcmp(p1[s].cnt, p2[s].cnt, np * sizeof(uint32_t)));
    }
    cfj_vocab_index ix;
    assert(cfj_vocab_index_build(&v1, &ix) == STATA_OK);
    assert(cfj_vocab_find(&v1, &ix, "ca", 2) == 4 && cfj_vocab_find(&v1, &ix, "zz", 2) == -1);
    cfj_vocab_index_free(&ix);
    /* unit flag: counts become 1 */
    assert(cfj_method_parse("bigram, 1", &m, &what) == 0);
    cfj_vocab v3;
    cfj_profile p3[2];
    assert(cfj_profiles_build(&m, sets, 2, &v3, p3, 3) == STATA_OK);
    assert(p3[0].cnt[0] == 1);
    cfj_vocab_free(&v1);
    cfj_vocab_free(&v2);
    cfj_vocab_free(&v3);
    for (int s = 0; s < 2; s++)
    {
        cfj_profile_free(&p1[s]);
        cfj_profile_free(&p2[s]);
        cfj_profile_free(&p3[s]);
    }
    printf("profiles ok\n");
}

int main(void)
{
    test_levenshtein();
    test_soundex();
    test_grams();
    test_profiles();
    return 0;
}
