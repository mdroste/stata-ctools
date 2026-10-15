"""ASan/UBSan coverage for binary codecs, corrupt files, and cache cleanup."""
from native.support.compiler import compile_native
from native.support.runner import harness
from pathlib import Path
import os
import platform
import subprocess
import tempfile
from run_io_formats import ROOT, shapefiles, codec_fixtures
from io_excel_format_fixtures import fixtures

SOURCE = harness("io_formats_native_source.c")


def main():
    vendor = ROOT / 'src/vendor'
    with tempfile.TemporaryDirectory(prefix='ctools-binary-asan-') as tmp:
        directory = Path(tmp)
        shapefiles(directory)
        codec_fixtures(directory)
        fixtures(directory)
        cfile = directory / 'test.c'
        cfile.write_text(SOURCE)
        sources = [*sorted((vendor / 'readstat').rglob('*.c')),
                   *sorted((vendor / 'zlib').glob('*.c')),
                   *sorted((vendor / 'libxls/src').glob('*.c')),
                   *sorted((ROOT/'src/vendor/miniz').glob('*.c')),
                   *sorted((ROOT/'src/vendor/libdeflate').rglob('*.c')),
                   ROOT / 'src/cexport/cexport_parse.c',
                   *[ROOT/'src'/name for name in ['cimport/cimport_xlsx.c','formats/xlsx_xml.c','formats/xlsx_zip.c','ctools_types.c','ctools_arena.c','ctools_threads.c']]]
        flags = ['-std=c11', '-O1', '-g', '-DHAVE_ZLIB=1', '-DSD_FASTMODE',
                 '-D_GNU_SOURCE', '-D_DARWIN_C_SOURCE', '-pthread', '-fno-fast-math', '-ffp-contract=off',
                 '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-I', str(ROOT / 'src'),
                 '-I', str(vendor / 'readstat'), '-I', str(vendor / 'zlib'),
                 '-I', str(vendor / 'libxls/include'), '-I', str(vendor / 'libxls/include/libxls')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip', '-liconv']
        else:
            flags += ['-DSYSTEM=STUNIX', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', '-lm']
        binary = directory / 'test'
        compile_native([os.environ.get('CC', 'cc'), *flags, str(cfile), *map(str, sources), '-o', str(binary)], check=True)
        subprocess.run([str(binary), str(directory)], check=True, timeout=60,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})
        print('PASS binary codecs, truncation, and cache cleanup (ASan/UBSan)')


if __name__ == '__main__':
    main()
