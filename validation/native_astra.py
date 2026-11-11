"""Small native regressions for the November 6 audit (no Stata license needed)."""
import os
import platform
import shlex
import subprocess
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import unit, harness


def test_astra_native():
    compiler = shlex.split(os.environ.get('CTOOLS_AUDIT_CC', os.environ.get('CC', 'cc')))
    with tempfile.TemporaryDirectory(prefix='ctools-astra-') as tmp:
        folder = Path(tmp)
        flags = ['-std=gnu11', '-O1', '-UNDEBUG', '-fno-sanitize-recover=all', '-fno-fast-math', '-ffp-contract=off',
                 '-ffunction-sections', '-fdata-sections', '-I' + str(ROOT / 'src')]
        if platform.system() == 'Windows':
            flags += ['-DSYSTEM=STWIN32', '-Wl,--gc-sections']
        elif platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip', '-fsanitize=undefined']
        else:
            flags += ['-DSYSTEM=STUNIX', '-D_GNU_SOURCE', '-Wl,--gc-sections',
                      '-fsanitize=' + os.environ.get('CTOOLS_SANITIZERS', 'address,undefined,float-cast-overflow')]
        for name in ('astra_parquet', 'astra_merge', 'astra_stream', 'astra_encode'):
            binary = folder / (name + ('.exe' if os.name == 'nt' else ''))
            subprocess.run([*compiler, *flags, str(unit(folder, name)), *shlex.split(os.environ.get('CTOOLS_AUDIT_LDFLAGS', '')), '-lm', '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True, timeout=60)
            print('PASS ' + name, flush=True)
            if name == 'astra_stream' and platform.system() == 'Linux':
                subprocess.run([*compiler, *flags, '-fopenmp', str(unit(folder, name)), '-lm', '-o', str(binary)], check=True)
                subprocess.run([str(binary)], check=True, timeout=60)
                print('PASS astra_stream OpenMP', flush=True)
        # These two fixtures use only ordinary allocation. The shared POSIX
        # fault allocator's aligned entry point must never run on Windows.
        preamble = ''
        if os.name == 'nt':
            preamble = '#include <stdlib.h>\n#include <malloc.h>\nstatic int posix_memalign(void **p,size_t a,size_t n){abort();}\n'
        for bootstrap in (False, True):
            source = folder / 'publication.c'
            source.write_text(preamble + ('#define BOOTSTRAP\n' if bootstrap else '') + harness('memory_sampling'))
            binary = folder / ('publication.exe' if os.name == 'nt' else 'publication')
            subprocess.run([*compiler, *flags, str(source), '-lm', '-o', str(binary)], check=True)
            for failure in range(0, 12 if bootstrap else 11):
                subprocess.run([str(binary), '0', '0', str(failure)], check=True, timeout=60, stdout=subprocess.DEVNULL)
            print('PASS sampler publications bootstrap=' + str(bootstrap), flush=True)
        source.write_text(preamble + (ROOT / 'validation/native/alloc.h').read_text() + '\n' + harness('cpsmatch_paths_oom'))
        subprocess.run([*compiler, *flags, str(source), '-lm', '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60)
        print('PASS matching allocation and publication failures', flush=True)
        source.write_text(harness('newcommands_native_source'))
        subprocess.run([*compiler, *flags, str(source), str(ROOT / 'src/ctools_parse.c'), '-lm', '-o', str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60)
