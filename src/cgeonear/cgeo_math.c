/* Runtime math functions that match Stata's; see cgeo_math.h. */
#include "cgeo_math.h"

#ifdef _WIN32
#include <windows.h>
#endif

static double cgeo_libm_sin(double x)
{
    return sin(x);
}
static double cgeo_libm_cos(double x)
{
    return cos(x);
}
static double cgeo_libm_atan2(double y, double x)
{
    return atan2(y, x);
}

cgeo_fn1 cgeo_sin_fn = cgeo_libm_sin;
cgeo_fn1 cgeo_cos_fn = cgeo_libm_cos;
cgeo_fn2 cgeo_atan2_fn = cgeo_libm_atan2;

int cgeo_math_init(void)
{
    static int done = 0, found = 0;
    if (done)
        return found;
    done = 1;
#ifdef _WIN32
    /* Stata for Windows calls the Universal C Runtime, which is part of the
     * system on Windows 10 and later; the module is already loaded in
     * Stata's process. Without it the C library's functions stay. */
    HMODULE ucrt = GetModuleHandleA("ucrtbase.dll");
    if (!ucrt)
        ucrt = LoadLibraryA("ucrtbase.dll");
    if (ucrt)
    {
        FARPROC s = GetProcAddress(ucrt, "sin");
        FARPROC c = GetProcAddress(ucrt, "cos");
        FARPROC a = GetProcAddress(ucrt, "atan2");
        if (s && c && a)
        {
            cgeo_sin_fn = (cgeo_fn1)(void (*)(void))s;
            cgeo_cos_fn = (cgeo_fn1)(void (*)(void))c;
            cgeo_atan2_fn = (cgeo_fn2)(void (*)(void))a;
            found = 1;
        }
    }
#else
    /* Elsewhere Stata and the plugin share the system's math library. */
    found = 1;
#endif
    return found;
}
