#include "ctools_file.h"
#include <assert.h>
#include <unistd.h>
static const char *fail_destination;
static int force_copy;
#ifdef __APPLE__
static int test_rename_excl(const char *from, const char *to, unsigned flags)
{
    if (fail_destination && !strcmp(to, fail_destination)) {
        fail_destination = NULL; errno = EIO; return -1;
    }
    if (force_copy) { errno = ENOTSUP; return -1; }
    return renamex_np(from, to, flags);
}
#define renamex_np test_rename_excl
#elif defined(__linux__)
#include <sys/syscall.h>
#include <stdarg.h>
/* Match syscall's variadic declaration when the implementation is included. */
static long test_syscall(long number, ...)
{
    assert(number == SYS_renameat2);
    va_list args;
    va_start(args, number);
    int dir1 = va_arg(args, int);
    const char *from = va_arg(args, const char *);
    int dir2 = va_arg(args, int);
    const char *to = va_arg(args, const char *);
    int flags = va_arg(args, int);
    va_end(args);
    if (fail_destination && !strcmp(to, fail_destination)) {
        fail_destination = NULL; errno = EIO; return -1;
    }
    if (force_copy) { errno = ENOTSUP; return -1; }
    return syscall(number, dir1, from, dir2, to, flags);
}
#define syscall test_syscall
#endif
static int test_link(const char *from, const char *to)
{
    if (force_copy) { errno = ENOTSUP; return -1; }
    return link(from, to);
}
#define link test_link
#include "formats/ctools_output.c"

static void write_file(const char *path, const char *text)
{
    FILE *fp = fopen(path, "wb"); assert(fp);
    assert(fwrite(text, 1, strlen(text), fp) == strlen(text)); assert(!fclose(fp));
}
static void contents(const char *path, const char *text)
{
    char buf[100] = {0}; FILE *fp = fopen(path, "rb"); assert(fp);
    fread(buf, 1, sizeof(buf)-1, fp); fclose(fp); assert(!strcmp(buf, text));
}
static void prepare(ctools_output *out, const char *path, int replace)
{
    assert(!ctools_output_prepare(out, path, replace)); write_file(out->temporary, "new");
}
int main(int argc, char **argv)
{
    assert(argc == 2); assert(!chdir(argv[1])); umask(0027);
    ctools_output a, b; ctools_output *pair[] = {&a, &b};
    prepare(&a, "mode", 0); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    struct stat st; assert(!stat("mode", &st)); assert((st.st_mode & 0777) == 0640);
    assert(!chmod("mode", 0664)); prepare(&a, "mode", 1);
    assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    assert(!stat("mode", &st)); assert((st.st_mode & 0777) == 0664);
    char name[251]; memset(name, 'x', 250); name[250] = 0;
    prepare(&a, name, 0); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a); contents(name, "new");
    force_copy = 1; prepare(&a, "fat", 0); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    contents("fat", "new"); prepare(&a, "fat", 0); assert(ctools_output_commit(&a) == 602); ctools_output_cleanup(&a);
    force_copy = 0;
    write_file("b", "old-b"); prepare(&a, "a", 0); prepare(&b, "b", 0);
    assert(ctools_output_commit_group(pair, 2) == 602); assert(access("a", F_OK)); contents("b", "old-b");
    ctools_output_cleanup(&a); ctools_output_cleanup(&b);
    for (int existing = 0; existing < 2; existing++) {
        if (existing) write_file("a", "old-a");
        prepare(&a, "a", existing); prepare(&b, "b", 1);
        fail_destination = "b";
        assert(ctools_output_commit_group(pair, 2) == 603);
        if (existing) contents("a", "old-a"); else assert(access("a", F_OK));
        contents("b", "old-b"); ctools_output_cleanup(&a); ctools_output_cleanup(&b);
    }
    prepare(&a, "a", 1); prepare(&b, "b", 1); assert(!ctools_output_commit_group(pair, 2));
    contents("a", "new"); contents("b", "new"); ctools_output_cleanup(&a); ctools_output_cleanup(&b);
#ifdef __APPLE__
    write_file("acl", "old");
    assert(!system("chmod +a 'everyone allow read' acl"));
    acl_t acl = acl_get_file("acl", ACL_TYPE_EXTENDED); assert(acl);
    char *before = acl_to_text(acl, NULL); assert(before); acl_free(acl);
    prepare(&a, "acl", 1); assert(!ctools_output_commit(&a)); ctools_output_cleanup(&a);
    acl = acl_get_file("acl", ACL_TYPE_EXTENDED); assert(acl);
    char *after = acl_to_text(acl, NULL); assert(after); assert(!strcmp(before, after));
    acl_free(before); acl_free(after); acl_free(acl);
#endif
    write_file("alias", "original"); prepare(&a, "alias", 1); prepare(&b, "./alias", 1);
    assert(ctools_output_commit_group(pair, 2) == 198); contents("alias", "original");
    ctools_output_cleanup(&a); ctools_output_cleanup(&b);
    puts("PASS output rollback, permissions/ACLs, long names, and no-hardlink fallback");
    return 0;
}
