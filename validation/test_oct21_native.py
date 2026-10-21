"""Regression tests for the October 21 platform/data-integrity audit."""
from pathlib import Path
import os
import platform
import struct
import subprocess
import tempfile
import zipfile
from native.support.compiler import compile_native

ROOT = Path(__file__).resolve().parents[1]


def xlsx_fixtures(folder):
    for name, method in [('valid.xlsx', zipfile.ZIP_DEFLATED), ('stored.xlsx', zipfile.ZIP_STORED)]:
        with zipfile.ZipFile(folder / name, 'w', compression=method, compresslevel=0) as archive:
            archive.writestr('sheet.xml', b'<v>123</v>')
    original = (folder / 'valid.xlsx').read_bytes()
    data = bytearray(original)
    data[data.index(b'<v>123') + 3] ^= 8
    (folder / 'bitflip.xlsx').write_bytes(data)
    for name, field in [('crc.xlsx', 16), ('length.xlsx', 24)]:
        data = bytearray(original)
        offset = data.index(b'PK\x01\x02') + field
        struct.pack_into('<I', data, offset, struct.unpack_from('<I', data, offset)[0] + 1)
        (folder / name).write_bytes(data)


def main():
    with tempfile.TemporaryDirectory(prefix='ctools-oct21-') as tmp:
        directory = Path(tmp)
        xlsx_fixtures(directory)
        flags = ['-std=c11', '-O1', '-g', '-DSD_FASTMODE', '-D_GNU_SOURCE', '-D_DARWIN_C_SOURCE',
                 '-pthread', '-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-I', str(ROOT/'src')]
        if platform.system() == 'Darwin':
            flags += ['-DSYSTEM=APPLEMAC', '-Wl,-dead_strip']
        else:
            flags += ['-DSYSTEM=STUNIX', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', '-lm']
        sources = [ROOT/'validation/native/harnesses/oct21_xlsx.c',
                   *[ROOT/'src'/p for p in ['formats/xlsx_xml.c', 'ctools_types.c', 'ctools_arena.c', 'ctools_threads.c']],
                   *sorted((ROOT/'src/vendor/miniz').glob('*.c')),
                   *sorted((ROOT/'src/vendor/libdeflate').rglob('*.c'))]
        binary = directory/'xlsx'
        compile_native([os.environ.get('CC', 'cc'), *flags, *map(str, sources), '-o', str(binary)], check=True)
        subprocess.run([str(binary), str(directory)], check=True,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})

        output = directory/'output'
        compile_native([os.environ.get('CC', 'cc'), *flags,
                        str(ROOT/'validation/native/harnesses/oct21_output.c'), '-o', str(output)], check=True)
        subprocess.run([str(output), str(directory)], check=True,
                       env={**os.environ, 'UBSAN_OPTIONS': 'halt_on_error=1'})


if __name__ == '__main__':
    main()
