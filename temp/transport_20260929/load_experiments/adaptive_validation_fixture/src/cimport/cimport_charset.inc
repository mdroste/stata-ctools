/* Portable C adaptation of Unicode ICU 67.1 charset recognition.
   Copyright (C) 1996-2014 IBM and others; (C) 2016-2020 Unicode, Inc.
   See src/io/vendor/icu_charset/LICENSE and README.md for provenance.
   Like ICU's InputStream detector, recognition samples at most 8,000 bytes. */
typedef struct {
  const char *name, *c1name;
  const unsigned char *map;
  const uint32_t (*grams)[64];
  size_t languages;
} charset_sbcs_model;
#include "cimport_charset_tables.inc"

static int charset_has(const uint32_t *table, size_t length, uint32_t value) {
  size_t lo = 0, hi = length;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (table[mid] < value)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo < length && table[lo] == value;
}
static int charset_sbcs_score(const unsigned char *data, size_t n,
                              const charset_sbcs_model *model) {
  int best = 0;
  for (size_t lang = 0; lang < model->languages; lang++) {
    uint32_t gram = 0;
    size_t count = 0, hits = 0;
    int space = 0;
    for (size_t i = 0; i <= n; i++) {
      unsigned char b = i == n ? 32 : model->map[data[i]];
      if (!b)
        continue;
      if (i == n || !(b == 32 && space)) {
        gram = ((gram << 8) | b) & 0xffffff;
        count++;
        hits += charset_has(model->grams[lang], 64, gram);
      }
      space = b == 32;
    }
    double fraction = (double)hits / (double)count;
    int score = fraction > 0.33 ? 98 : (int)(fraction * 300);
    if (score > best)
      best = score;
  }
  return best;
}
static int charset_utf8_score(const unsigned char *s, size_t n) {
  int valid = 0, invalid = 0;
  int bom = n >= 3 && s[0] == 0xef && s[1] == 0xbb && s[2] == 0xbf;
  for (size_t i = 0; i < n; i++) {
    unsigned char b = s[i];
    if (b < 128)
      continue;
    int trail = (b & 0xe0) == 0xc0   ? 1
                : (b & 0xf0) == 0xe0 ? 2
                : (b & 0xf8) == 0xf0 ? 3
                                     : 0;
    if (!trail) {
      invalid++;
      continue;
    }
    while (++i < n) {
      if ((s[i] & 0xc0) != 0x80) {
        invalid++;
        break;
      }
      if (!--trail) {
        valid++;
        break;
      }
    }
  }
  if (bom && !invalid)
    return 100;
  if (bom && valid > invalid * 10)
    return 80;
  if (valid > 3 && !invalid)
    return 100;
  if (valid && !invalid)
    return 80;
  if (!valid && !invalid)
    return 15;
  return valid > invalid * 10 ? 25 : 0;
}
static int charset_utf16_score(const unsigned char *s, size_t n, int little) {
  int score = 10;
  /* ICU InputStream supplies an 8,000-byte, zero-padded buffer. */
  for (size_t i = 0; i < 30; i += 2) {
    unsigned a = i < n ? s[i] : 0, b = i + 1 < n ? s[i + 1] : 0;
    unsigned cp = little ? a | b << 8 : a << 8 | b;
    if (!i && cp == 0xfeff)
      return 100;
    if (!cp)
      score -= 10;
    else if ((cp >= 32 && cp <= 255) || cp == 10)
      score += 10;
    if (score <= 0)
      return 0;
    if (score >= 100)
      return 100;
  }
  return score;
}
static int charset_utf32_score(const unsigned char *s, size_t n, int little) {
  int valid = 0, invalid = 0, bom = 0;
  for (size_t i = 0; i + 3 < n; i += 4) {
    uint32_t cp = 0;
    for (int j = 0; j < 4; j++)
      cp = cp << 8 | s[i + (little ? 3 - j : j)];
    if (!i && cp == 0xfeff)
      bom = 1;
    if (cp >= 0x10ffff || (cp >= 0xd800 && cp <= 0xdfff))
      invalid++;
    else
      valid++;
  }
  if (!valid && !invalid)
    return 0;
  if (bom && !invalid)
    return 100;
  if (bom && valid > invalid * 10)
    return 80;
  if (valid > 3 && !invalid)
    return 100;
  if (valid && !invalid)
    return 80;
  return valid > invalid * 10 ? 25 : 0;
}
/* Each sequence contains a byte count followed by its escape bytes. */
static const unsigned char charset_2022jp[][5] = {
    {4, 27, 36, 40, 67}, {4, 27, 36, 40, 68}, {3, 27, 36, 64}, {3, 27, 36, 65},
    {3, 27, 36, 66},     {3, 27, 38, 64},     {3, 27, 40, 66}, {3, 27, 40, 72},
    {3, 27, 40, 73},     {3, 27, 40, 74},     {3, 27, 46, 65}, {3, 27, 46, 70}};
static const unsigned char charset_2022cn[][5] = {
    {4, 27, 36, 41, 65}, {4, 27, 36, 41, 71}, {4, 27, 36, 42, 72},
    {4, 27, 36, 41, 69}, {4, 27, 36, 43, 73}, {4, 27, 36, 43, 74},
    {4, 27, 36, 43, 75}, {4, 27, 36, 43, 76}, {4, 27, 36, 43, 77},
    {2, 27, 78},         {2, 27, 79}};
static const unsigned char charset_2022kr[][5] = {{4, 27, 36, 41, 67}};
static int charset_2022_score(const unsigned char *s, size_t n,
                              const unsigned char (*sequences)[5],
                              size_t count) {
  int hits = 0, misses = 0, shifts = 0;
  for (size_t i = 0; i < n; i++) {
    if (s[i] == 27) {
      size_t j = 0;
      for (; j < count; j++) {
        size_t length = sequences[j][0];
        if (length <= n - i && !memcmp(s + i, sequences[j] + 1, length)) {
          hits++;
          i += length - 1;
          break;
        }
      }
      if (j < count)
        continue;
      misses++;
    }
    if (s[i] == 14 || s[i] == 15)
      shifts++;
  }
  if (!hits)
    return 0;
  int score = (100 * hits - 100 * misses) / (hits + misses);
  if (hits + shifts < 5)
    score -= (5 - hits - shifts) * 10;
  return score > 0 ? score : 0;
}
static int charset_mbcs_score(const unsigned char *s, size_t n, int kind,
                              const uint32_t *common, size_t common_length) {
  size_t i = 0;
  int total = 0, doubles = 0, bad = 0, hits = 0;
  while (i < n) {
    unsigned first = s[i++];
    uint32_t ch = first;
    int error = 0;
    int single = kind == 0   ? first <= 0x7f || (first >= 0xa1 && first <= 0xdf)
                 : kind == 1 ? first <= 0x80
                 : kind == 4 ? first <= 0x7f || first == 0xff
                             : first <= 0x8d;
    if (!single) {
      if (i == n)
        break;
      unsigned second = s[i++];
      ch = ch << 8 | second;
      if (!kind)
        error = second < 0x40;
      else if (kind == 4)
        error = second < 0x40 || second == 0x7f || second == 0xff;
      else if (kind == 1) {
        if (first >= 0x81 && first <= 0xfe &&
            !((second >= 0x40 && second <= 0x7e) ||
              (second >= 80 && second <= 0xfe))) {
          error = 1;
          if (second >= 0x30 && second <= 0x39) {
            if (i == n)
              break;
            unsigned third = s[i++];
            if (third >= 0x81 && third <= 0xfe) {
              if (i == n)
                break;
              unsigned fourth = s[i++];
              if (fourth >= 0x30 && fourth <= 0x39) {
                ch = ch << 16 | third << 8 | fourth;
                error = 0;
              }
            }
          }
        }
      } else if ((first >= 0xa1 && first <= 0xfe) || first == 0x8e)
        error = second < 0xa1;
      else if (first == 0x8f) {
        if (i == n)
          break;
        unsigned third = s[i++];
        ch = ch << 8 | third;
        error = third < 0xa1;
      }
    }
    total++;
    if (error)
      bad++;
    else if (ch > 255) {
      doubles++;
      hits += charset_has(common, common_length, ch);
    }
    if (bad >= 2 && bad * 5 >= doubles)
      return 0;
  }
  if (doubles <= 10 && !bad)
    return !doubles && total < 10 ? 0 : 10;
  if (doubles < 20 * bad)
    return 0;
  double scale = 90 / log((double)((float)doubles / 4));
  int score = (int)(log(hits + 1) * scale + 10);
  return score > 100 ? 100 : score;
}
const char *cimport_detect_charset(const char *data, size_t length) {
  const unsigned char *s = (const unsigned char *)data;
  size_t n = length < 8000 ? length : 8000;
  const char *name = "UTF-8";
  int best = 0, score;
  /* Stable ascending sort followed by reversal favors the later recognizer
     when confidences tie. Disabled EBCDIC recognizers are omitted. */
#define CHARSET_CANDIDATE(label, expression)                                   \
  do {                                                                         \
    score = (expression);                                                      \
    if (score && score >= best) {                                              \
      best = score;                                                            \
      name = (label);                                                          \
    }                                                                          \
  } while (0)
  CHARSET_CANDIDATE("UTF-8", charset_utf8_score(s, n));
  CHARSET_CANDIDATE("UTF-16BE", charset_utf16_score(s, n, 0));
  CHARSET_CANDIDATE("UTF-16LE", charset_utf16_score(s, n, 1));
  CHARSET_CANDIDATE("UTF-32BE", charset_utf32_score(s, n, 0));
  CHARSET_CANDIDATE("UTF-32LE", charset_utf32_score(s, n, 1));
#define CHARSET_MBCS(label, idx)                                               \
  CHARSET_CANDIDATE(label, charset_mbcs_score(s, n, idx, charset_common_##idx, \
                                              sizeof(charset_common_##idx) /   \
                                                  sizeof(uint32_t)))
  CHARSET_MBCS("Shift_JIS", 0);
  CHARSET_CANDIDATE("ISO-2022-JP",
                    charset_2022_score(s, n, charset_2022jp, 12));
  CHARSET_CANDIDATE("ISO-2022-CN",
                    charset_2022_score(s, n, charset_2022cn, 11));
  CHARSET_CANDIDATE("ISO-2022-KR", charset_2022_score(s, n, charset_2022kr, 1));
  CHARSET_MBCS("GB18030", 1);
  CHARSET_MBCS("EUC-JP", 2);
  CHARSET_MBCS("EUC-KR", 3);
  CHARSET_MBCS("Big5", 4);
  int c1 = 0;
  for (size_t i = 0; i < n; i++)
    if (s[i] >= 0x80 && s[i] <= 0x9f) {
      c1 = 1;
      break;
    }
  for (size_t j = 0;
       j < sizeof(charset_sbcs_models) / sizeof(charset_sbcs_models[0]); j++) {
    const charset_sbcs_model *m = &charset_sbcs_models[j];
    CHARSET_CANDIDATE(c1 && *m->c1name ? m->c1name : m->name,
                      charset_sbcs_score(s, n, m));
  }
#undef CHARSET_MBCS
#undef CHARSET_CANDIDATE
  /* The native Java decoder does not support this ICU logical-Hebrew name;
     its automatic-encoding path consequently falls back to UTF-8. */
  return !strcmp(name, "ISO-8859-8-I") ? "UTF-8" : name;
}
