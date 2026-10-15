#include "excel_address.h"
#include <stddef.h>
bool ctools_excel_address(const char *text, int max_col, int max_row, bool absolute, int *col,
                          int *row, const char **end)
{
    if (!text || !col || !row || max_col < 1 || max_row < 1)
        return false;
    const char *p = text;
    int c = 0, r = 0;
    if (absolute && *p == '$')
        p++;
    while ((*p >= 'A' && *p <= 'Z') || (*p >= 'a' && *p <= 'z'))
    {
        int digit = (*p >= 'a' ? *p - 'a' : *p - 'A') + 1;
        if (digit > max_col || c > (max_col - digit) / 26)
            return false;
        c = c * 26 + digit;
        p++;
    }
    if (!c)
        return false;
    if (absolute && *p == '$')
        p++;
    while (*p >= '0' && *p <= '9')
    {
        int digit = *p++ - '0';
        if (digit > max_row || r > (max_row - digit) / 10)
            return false;
        r = r * 10 + digit;
    }
    if (!r || (!end && *p))
        return false;
    *col = c;
    *row = r;
    if (end)
        *end = p;
    return true;
}
