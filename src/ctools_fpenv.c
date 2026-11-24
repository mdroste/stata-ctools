/*
 * ctools_fpenv.c - every thread computes with Stata's x87 precision (Windows GCC)
 *
 * MinGW-w64's exp(), log() and related functions evaluate on the x87 unit,
 * so the last bit of their result depends on the x87 precision control.
 * Stata's thread runs with the Windows default control word, 0x027F
 * (53-bit precision, round to nearest, exceptions masked), as does every
 * thread Windows creates. winpthreads, under libgomp and the ctools thread
 * pool, calls _fpreset() as each thread starts. MinGW's _fpreset() is a bare
 * fninit, which selects 64-bit precision (0x037F). A worker's exp() then
 * differed in the last bit from the calling thread's for about one argument
 * in six, so a result depended on which thread computed it. The bootstrap
 * hands replicates to threads dynamically, so its results changed from run
 * to run and with threads().
 *
 * The linker takes this object's definition before MinGW's archive member,
 * so a thread starts in the state Stata's thread runs in. A thread Windows
 * creates directly (without winpthreads) already starts in that state.
 */
#if defined(_WIN32) && defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))

void _fpreset(void);

__attribute__((used)) void _fpreset(void)
{
    const unsigned short control = 0x027F;
    __asm__ __volatile__("fninit\n\tfldcw %0" : : "m"(control));
}

#else
typedef int ctools_fpenv_unused; /* ISO C wants a declaration in every file */
#endif
