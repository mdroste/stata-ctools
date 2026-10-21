#include "ctools_output.h"
#include "ctools_file.h"
#ifdef __APPLE__
#include <sys/acl.h>
#elif defined(__linux__)
#include <sys/syscall.h>
/* Some native harnesses request only POSIX declarations. */
extern long syscall(long number, ...);
#include <sys/xattr.h>
#endif

/* A short sibling name leaves the destination's component limit untouched. */
ST_retcode ctools_output_prepare(ctools_output *out, const char *filename, bool replace)
{
    memset(out, 0, sizeof(*out));
    out->replace = replace;
    out->destination = strdup(filename);
    const char *slash = strrchr(filename, '/');
#ifdef _WIN32
    const char *backslash = strrchr(filename, '\\');
    if (backslash && (!slash || backslash > slash)) slash = backslash;
#endif
    size_t parent = slash ? (size_t)(slash - filename) + 1 : 0;
#ifdef _WIN32
    if (!slash && strlen(filename) >= 2 && filename[1] == ':') parent = 2;
#endif
    out->temporary = malloc(parent + 32);
    if (!out->destination || !out->temporary) {
        ctools_output_cleanup(out); return 920;
    }
    memcpy(out->temporary, filename, parent);
    strcpy(out->temporary + parent, ".ctools.XXXXXX");
#ifdef _WIN32
    wchar_t *wide = ctools_wpath(out->temporary);
    int fd = -1;
    if (wide && !_wmktemp_s(wide, wcslen(wide) + 1)) {
        fd = _wopen(wide, _O_CREAT | _O_EXCL | _O_BINARY | _O_WRONLY, _S_IREAD | _S_IWRITE);
        if (fd >= 0) {
            size_t n = strlen(out->temporary);
            size_t w = wcslen(wide);
            for (int i = 0; i < 6; i++) out->temporary[n - 6 + i] = (char)wide[w - 6 + i];
            _close(fd);
        }
    }
    free(wide);
#else
    int fd = mkstemp(out->temporary);
    if (fd >= 0) {
        /* Recreate the empty reservation exclusively to capture the caller's
         * umask/default ACL without changing the process-wide umask. No data
         * is written until permissions have been restricted to 0600. */
        close(fd);
        if (unlink(out->temporary)) { out->owns_temporary = true; fd = -1; }
        else {
            fd = open(out->temporary, O_CREAT | O_EXCL | O_WRONLY, 0666);
            if (fd >= 0) {
                out->owns_temporary = true;
                struct stat st;
                if (fstat(fd, &st) || fchmod(fd, 0600)) { close(fd); fd = -1; }
                else {
                    out->mode = st.st_mode & 0777;
#ifdef __APPLE__
                    out->creation_acl = acl_get_fd_np(fd, ACL_TYPE_EXTENDED);
                    int acl_error = errno;
                    /* macOS reports ENOENT for a file with no extended ACL. */
                    if (!out->creation_acl && acl_error == ENOENT)
                        out->creation_acl = acl_init(0);
                    if (out->creation_acl) {
                        acl_t empty = acl_init(0);
                        if (!empty || acl_set_fd_np(fd, empty, ACL_TYPE_EXTENDED)) {
                            close(fd); fd = -1;
                        }
                        if (empty) acl_free(empty);
                    } else if (acl_error != ENOTSUP) { close(fd); fd = -1; }
#endif
                    if (fd >= 0) close(fd);
                }
            }
        }
    }
#endif
    if (fd < 0) { ctools_output_cleanup(out); return 603; }
    out->owns_temporary = true;
    return 0;
}

#ifndef _WIN32
static int output_copy_acl(const char *from, int fd)
{
    int rc = 0;
#ifdef __APPLE__
        acl_t acl = acl_get_file(from, ACL_TYPE_EXTENDED);
        if (!acl && errno == ENOENT) acl = acl_init(0);
        if (acl) { if (acl_set_fd_np(fd, acl, ACL_TYPE_EXTENDED)) rc = -1; acl_free(acl); }
        else if (errno != ENOTSUP) rc = -1;
#elif defined(__linux__)
        ssize_t n = getxattr(from, "system.posix_acl_access", NULL, 0);
        if (n > 0) {
            void *acl = malloc((size_t)n);
            if (!acl || getxattr(from, "system.posix_acl_access", acl, (size_t)n) != n ||
                fsetxattr(fd, "system.posix_acl_access", acl, (size_t)n, 0)) rc = -1;
            free(acl);
        } else if (n < 0 && errno != ENODATA && errno != ENOTSUP) rc = -1;
        else if (fremovexattr(fd, "system.posix_acl_access") && errno != ENODATA && errno != ENOTSUP) rc = -1;
#endif
    return rc;
}

static int output_metadata(ctools_output *out)
{
    struct stat st;
    int exists = lstat(out->destination, &st) == 0;
    if (!exists && errno != ENOENT) return -1;
    if (exists && !out->replace) { errno = EEXIST; return -1; }
    if (exists && !S_ISREG(st.st_mode)) { errno = EINVAL; return -1; }
    int fd = open(out->temporary, O_WRONLY);
    if (fd < 0) return -1;
    int rc = 0;
    if (exists) {
        struct stat temp;
        if (fstat(fd, &temp) || (temp.st_gid != st.st_gid && fchown(fd, -1, st.st_gid))) rc = -1;
        if (output_copy_acl(out->destination, fd)) rc = -1;
    }
#ifdef __APPLE__
    if (!exists && out->creation_acl && acl_set_fd_np(fd, out->creation_acl, ACL_TYPE_EXTENDED)) rc = -1;
#endif
    if (!rc && fchmod(fd, exists ? (st.st_mode & 0777) : out->mode)) rc = -1;
    close(fd);
    return rc;
}

/* Last-resort no-overwrite copy for filesystems with neither exclusive rename
 * nor hard links. The destination is private until complete; unlike rename,
 * its existence is visible during the copy. Remove it on any ordinary error. */
static int output_copy_exclusive(const char *from, const char *to)
{
    int in = open(from, O_RDONLY);
    if (in < 0) return -1;
    int out = open(to, O_WRONLY | O_CREAT | O_EXCL, 0600);
    if (out < 0) { close(in); return -1; }
    char buffer[65536];
    int rc = 0;
#ifdef __APPLE__
    acl_t empty = acl_init(0);
    if (!empty || (acl_set_fd_np(out, empty, ACL_TYPE_EXTENDED) && errno != ENOTSUP)) rc = -1;
    if (empty) acl_free(empty);
#endif
    ssize_t n;
    while (!rc && (n = read(in, buffer, sizeof(buffer))) != 0) {
        if (n < 0) { if (errno == EINTR) continue; rc = -1; break; }
        for (ssize_t at = 0; at < n;) {
            ssize_t written = write(out, buffer + at, (size_t)(n - at));
            if (written < 0 && errno == EINTR) continue;
            if (written <= 0) { rc = -1; break; }
            at += written;
        }
        if (rc) break;
    }
    struct stat st;
    if (!rc && (fstat(in, &st) || output_copy_acl(from, out) || fchmod(out, st.st_mode & 0777) || fsync(out))) rc = -1;
    if (close(out)) rc = -1;
    close(in);
    if (rc) unlink(to);
    else unlink(from);
    return rc;
}
#endif

static int output_move(const char *from, const char *to, bool replace)
{
#ifdef _WIN32
    return ctools_move_file(from, to, replace);
#else
    if (replace) return rename(from, to);
#ifdef __APPLE__
    if (!renamex_np(from, to, RENAME_EXCL)) return 0;
    if (errno != ENOTSUP && errno != EINVAL && errno != ENOSYS) return -1;
#elif defined(__linux__) && defined(SYS_renameat2)
    if (!syscall(SYS_renameat2, AT_FDCWD, from, AT_FDCWD, to, 1 /* RENAME_NOREPLACE */)) return 0;
    if (errno != ENOTSUP && errno != EINVAL && errno != ENOSYS) return -1;
#endif
    if (!link(from, to)) { unlink(from); return 0; }
    if (errno != ENOTSUP && errno != EOPNOTSUPP && errno != EPERM && errno != ENOSYS) return -1;
    return output_copy_exclusive(from, to);
#endif
}

ST_retcode ctools_output_commit(ctools_output *out)
{
#ifndef _WIN32
    if (output_metadata(out)) return errno == EEXIST ? 602 : 603;
#endif
    if (output_move(out->temporary, out->destination, out->replace))
        return errno == EEXIST ? 602 : 603;
    out->owns_temporary = false;
    return 0;
}

/* Compare directory-resolved destinations before a multi-file publication.
 * Windows also treats different letter casing as the same destination. */
#ifndef _WIN32
static int output_parent_stat(const char *path, struct stat *st)
{
    const char *slash = strrchr(path, '/');
    char *parent = slash ? strndup(path, slash == path ? 1 : (size_t)(slash - path)) : strdup(".");
    if (!parent) return -1;
    int rc = stat(parent, st);
    free(parent);
    return rc;
}
#endif
static int output_same_destination(const char *a, const char *b)
{
#ifdef _WIN32
    wchar_t *left = ctools_wpath(a), *right = ctools_wpath(b);
    int same = left && right ? !_wcsicmp(left, right) : -1;
#else
    const char *left = strrchr(a, '/'), *right = strrchr(b, '/');
    left = left ? left + 1 : a;
    right = right ? right + 1 : b;
    struct stat sa, sb;
    if (stat(a, &sa) == 0 && stat(b, &sb) == 0 && sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino)
        return 1;
    if (strcmp(left, right)) return 0;
    if (output_parent_stat(a, &sa) || output_parent_stat(b, &sb)) return -1;
    int same = sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
#endif
#ifdef _WIN32
    free(left); free(right);
#endif
    return same;
}

/* Rollback ordinary errors across multiple output files. Backups remain on
 * disk (with a diagnostic) if rollback itself fails. There is no filesystem-
 * independent atomic commit across names, nor a crash/concurrent-writer lock. */
ST_retcode ctools_output_commit_group(ctools_output **outputs, size_t count)
{
    ctools_output *backup = calloc(count, sizeof(*backup));
    bool *saved = calloc(count, sizeof(*saved)), *published = calloc(count, sizeof(*published));
    ST_retcode rc = 0;
    if (!backup || !saved || !published) { rc = 920; goto done; }
    for (size_t i = 0; i < count; i++) {
        ctools_output *out = outputs[i];
        for (size_t j = 0; j < i; j++) {
            int same = output_same_destination(out->destination, outputs[j]->destination);
            if (same) { rc = same < 0 ? 603 : 198; goto rollback; }
        }
#ifndef _WIN32
        if (output_metadata(out)) { rc = errno == EEXIST ? 602 : 603; goto rollback; }
#else
        wchar_t *wide = ctools_wpath(out->destination);
        if (!wide) { rc = 603; goto rollback; }
        DWORD attributes = GetFileAttributesW(wide), error = GetLastError();
        free(wide);
        if (attributes != INVALID_FILE_ATTRIBUTES && (!out->replace || (attributes & FILE_ATTRIBUTE_DIRECTORY))) {
            rc = out->replace ? 603 : 602; goto rollback;
        }
        if (attributes == INVALID_FILE_ATTRIBUTES && error != ERROR_FILE_NOT_FOUND && error != ERROR_PATH_NOT_FOUND) {
            rc = 603; goto rollback;
        }
#endif
        if (out->replace && (rc = ctools_output_prepare(&backup[i], out->destination, true))) goto rollback;
    }
    for (size_t i = 0; i < count; i++) {
        ctools_output *out = outputs[i];
        if (out->replace) {
            if (!output_move(out->destination, backup[i].temporary, true)) saved[i] = true;
            else {
#ifdef _WIN32
                if (GetLastError() != ERROR_FILE_NOT_FOUND)
#else
                if (errno != ENOENT)
#endif
                { rc = 603; goto rollback; }
            }
        }
        /* Metadata was applied before moving the original away. */
        if (output_move(out->temporary, out->destination, false)) {
            rc = errno == EEXIST ? 602 : 603; goto rollback;
        }
        out->owns_temporary = false;
        published[i] = true;
    }
    goto done;
rollback:
    for (size_t i = count; i-- > 0;) {
        if (saved[i]) {
            if (output_move(backup[i].temporary, outputs[i]->destination, true)) {
                backup[i].owns_temporary = false;
                fprintf(stderr, "ctools: rollback failed; original retained at %s\n", backup[i].temporary);
                rc = 603;
            }
        } else if (published[i] && ctools_remove(outputs[i]->destination)) rc = 603;
    }
done:
    if (backup) for (size_t i = 0; i < count; i++) ctools_output_cleanup(&backup[i]);
    free(backup); free(saved); free(published);
    return rc;
}

void ctools_output_cleanup(ctools_output *out)
{
    if (out->owns_temporary) ctools_remove(out->temporary);
#ifdef __APPLE__
    if (out->creation_acl) acl_free(out->creation_acl);
#endif
    free(out->temporary);
    free(out->destination);
    memset(out, 0, sizeof(*out));
}
