/* Locale numeric parsing is independent of string storage and CSV quoting. */
typedef struct {
  const char *tag;
  uint32_t decimal, group;
  const char *negative, *exponent;
} cimport_locale_profile;
#include "cimport_locale_profiles.inc"
static cimport_locale_profile cimport_active_locale;
static bool cimport_locale_enabled;
static int cimport_locale_compare(const void *key, const void *entry) {
  return strcmp((const char *)key,
                ((const cimport_locale_profile *)entry)->tag);
}
int cimport_configure_locale(const char *name, const char *decimal,
                             const char *group) {
  cimport_locale_enabled = false;
  if (!*name && !*decimal && !*group)
    return 0;
  char tag[128];
  size_t n = strlen(name);
  if (n >= sizeof(tag))
    return 198;
  if (!n)
    strcpy(tag, "en-us");
  else {
    for (size_t i = 0; i <= n; i++)
      tag[i] = name[i] == '_' ? '-' : (char)tolower((unsigned char)name[i]);
  }
  const cimport_locale_profile *profile = bsearch(
      tag, cimport_locale_profiles,
      sizeof(cimport_locale_profiles) / sizeof(cimport_locale_profiles[0]),
      sizeof(cimport_locale_profiles[0]), cimport_locale_compare);
  if (!profile)
    return 198;
  cimport_active_locale = *profile;
  if (*decimal)
    cimport_active_locale.decimal = (unsigned char)*decimal;
  if (*group)
    cimport_active_locale.group = (unsigned char)*group;
  cimport_locale_enabled = true;
  return 0;
}
static int cimport_locale_codepoint(const char *src, size_t n, uint32_t *cp) {
  if (!n)
    return 0;
  unsigned char b = (unsigned char)src[0];
  if (b < 128) {
    *cp = b;
    return 1;
  }
  int width = b >= 0xc2 && b <= 0xdf   ? 2
              : b >= 0xe0 && b <= 0xef ? 3
              : b >= 0xf0 && b <= 0xf4 ? 4
                                       : 0;
  if (!width || (size_t)width > n)
    return 0;
  uint32_t value = b & ((1U << (7 - width)) - 1);
  for (int i = 1; i < width; i++) {
    unsigned char next = (unsigned char)src[i];
    if ((next & 0xc0) != 0x80)
      return 0;
    value = (value << 6) | (next & 63);
  }
  if (value > 0x10ffff || (value >= 0xd800 && value <= 0xdfff))
    return 0;
  *cp = value;
  return width;
}
static int cimport_locale_digit(uint32_t cp) {
  if (cp > 0xffff)
    return -1; /* DecimalFormat uses BMP Character.digit(char). */
  for (size_t i = 0;
       i < sizeof(cimport_digit_starts) / sizeof(cimport_digit_starts[0]); i++)
    if (cp >= cimport_digit_starts[i] && cp - cimport_digit_starts[i] < 10)
      return (int)(cp - cimport_digit_starts[i]);
  return -1;
}
static bool cimport_locale_number(const char *src, int len, double *value,
                                  double missing) {
  char stack[256], *out = stack;
  if (len < 0)
    return false;
  if ((size_t)len + 1 > sizeof(stack)) {
    out = malloc((size_t)len + 1);
    if (!out)
      return false;
  }
  size_t p = 0, used = 0, n = (size_t)len;
  bool valid = true, point = false, exponent = false, group = false;
  const char *prefix = cimport_active_locale.negative,
             *exp = cimport_active_locale.exponent;
  size_t prefixlen = strlen(prefix), explen = strlen(exp);
  if (prefixlen <= n && !memcmp(src, prefix, prefixlen)) {
    out[used++] = '-';
    p = prefixlen;
  }
  if (n - p == 3 && !memcmp(src + p, "\xe2\x88\x9e", 3)) {
    *value = missing;
    if (out != stack)
      free(out);
    return true;
  }
  while (p < n) {
    if (explen && explen <= n - p) {
      bool match = true;
      for (size_t i = 0; i < explen; i++) {
        unsigned char ch = (unsigned char)src[p + i];
        if (ch < 128)
          ch = (unsigned char)toupper(ch);
        if (ch != (unsigned char)exp[i]) {
          match = false;
          break;
        }
      }
      if (match) {
        if (exponent || group) {
          valid = false;
          break;
        }
        out[used++] = 'e';
        p += explen;
        exponent = true;
        /* The exponent sign is the locale's minus character. */
        const char *sign = prefix + prefixlen;
        while (sign > prefix && ((unsigned char)sign[-1] & 0xc0) == 0x80)
          sign--;
        if (sign > prefix)
          sign--;
        size_t signlen = (size_t)(prefix + prefixlen - sign);
        if (signlen && signlen <= n - p && !memcmp(src + p, sign, signlen)) {
          out[used++] = '-';
          p += signlen;
        }
        continue;
      }
    }
    uint32_t cp;
    int width = cimport_locale_codepoint(src + p, n - p, &cp);
    if (!width) {
      valid = false;
      break;
    }
    int digit = cimport_locale_digit(cp);
    if (digit >= 0) {
      out[used++] = (char)('0' + digit);
      group = false;
    } else if (cp == cimport_active_locale.decimal && !point && !exponent) {
      out[used++] = '.';
      point = true;
      group = false;
    } else if (cp == cimport_active_locale.group && !point && !exponent)
      group = true;
    else {
      valid = false;
      break;
    }
    p += (size_t)width;
  }
  if (group)
    valid = false;
  out[used] = 0;
  if (valid)
    valid = ctools_parse_double_with_separators(out, (int)used, value, missing,
                                                '.', 0);
  if (out != stack)
    free(out);
  return valid;
}
