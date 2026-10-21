"""Remove compiled artifacts identically under POSIX and MSYS2 make."""
from pathlib import Path
import shutil
import sys


def clean(build, source):
    build, source = Path(build), Path(source)
    objects = build / 'obj'
    if objects.exists():
        shutil.rmtree(objects)
    for directory, patterns in ((build, ('ctools*.plugin', '*.plugin.build.json', '*.o')),
                                (source, ('*.o',))):
        for pattern in patterns:
            for path in directory.glob(pattern):
                path.unlink()


if __name__ == '__main__':
    clean(*sys.argv[1:])
