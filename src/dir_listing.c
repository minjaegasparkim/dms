/*
 * Directory listing for DMS: names of the entries of one directory (POSIX opendir/readdir;
 * Linux and macOS). Called from Fortran through ISO_C_BINDING (module tdep_search) to find the
 * files of temperature-dependent optical data, <base>_<T>K or <base>_<T>K.<extension>.
 * One directory is open at a time; the routines are called during the set-up, never from an
 * OpenMP region.
 */
#include <dirent.h>
#include <string.h>

static DIR *dms_dir = NULL;

/* open a directory; 0 on success, -1 if it cannot be opened */
int dms_dir_open(const char *path)
{
    if (dms_dir != NULL) {
        closedir(dms_dir);
        dms_dir = NULL;
    }
    dms_dir = opendir(path);
    return dms_dir != NULL ? 0 : -1;
}

/* next entry: copies at most maxlen-1 characters of the name to name (NUL-terminated) and returns
   the full length of the name, or -1 at the end of the directory */
int dms_dir_next(char *name, int maxlen)
{
    struct dirent *e;
    size_t n;

    if (dms_dir == NULL || maxlen < 1) return -1;
    e = readdir(dms_dir);
    if (e == NULL) return -1;
    n = strlen(e->d_name);
    if (n > (size_t)(maxlen - 1)) {
        memcpy(name, e->d_name, (size_t)(maxlen - 1));
        name[maxlen - 1] = '\0';
    } else {
        memcpy(name, e->d_name, n + 1);
    }
    return (int)n;
}

void dms_dir_close(void)
{
    if (dms_dir != NULL) {
        closedir(dms_dir);
        dms_dir = NULL;
    }
}
