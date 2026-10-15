"""UBSan checks for CSV formatting and parsing boundaries, independent of Stata."""
from native.support.runner import harness
from pathlib import Path
import tempfile
from native.support.runner import run_test

SOURCE = harness("io_parity_native_source.c")

if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='ctools-io-parity-') as directory:
        run_test(Path(directory), 'io_parity', SOURCE, modules=('ctools_types.c',))
