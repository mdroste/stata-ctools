"""Complete native Parquet parser/decoder corpus and failure injection."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test
from validation.support.parquet_native_fixtures import corpus

MODULES=('vendor/libdeflate/lib/deflate_decompress.c','vendor/libdeflate/lib/crc32.c',
         'vendor/libdeflate/lib/utils.c','vendor/libdeflate/lib/x86/cpu_features.c',
         'vendor/libdeflate/lib/arm/cpu_features.c')


def source(folder):
    body=(ROOT/'validation/native/parquet_paths.c').read_text(encoding='utf-8')
    body=body.replace('#include "parquet_alloc.h"',
        '#include "'+(ROOT/'validation/native/parquet_alloc.h').as_posix()+'"')
    return ('#define PARQUET_FIXTURES "'+folder.as_posix()+'"\n'+
            '#define PARQUET_CASES "'+(folder/'cases.h').as_posix()+'"\n'+body)


def test_parquet_paths_native():
    with tempfile.TemporaryDirectory(prefix='ctools-parquet-paths-') as tmp:
        folder=corpus(Path(tmp)/'corpus')
        run_test(Path(tmp),'parquet_paths',source(folder),modules=MODULES)
