"""Offline differential tests for binary formats; invoke Stata through the stata shell alias."""
import argparse
import os
import platform
import tempfile
from pathlib import Path
import re
import shlex
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def shapefiles(directory):
    """Small ESRI records, built independently of either reader/writer."""
    def write(name, kind, bodies):
        records = b''.join(struct.pack('>II', i + 1, len(body) // 2) + body
                           for i, body in enumerate(bodies))
        header = struct.pack('>7I', 9994, 0, 0, 0, 0, 0, (100 + len(records)) // 2)
        header += struct.pack('<2I8d', 1000, kind, 1, 2, 3, 4, 7, 8, 9, 10)
        (directory / (name + '.shp')).write_bytes(header + records)
    xy = struct.pack('<4d', 1, 2, 3, 4)
    for kind in (0, 1, 11, 21):
        body = struct.pack('<I', kind)
        if kind:
            body += struct.pack('<2d', 1, 2)
            if kind == 11:
                body += struct.pack('<d', 7)
            if kind in (11, 21):
                body += struct.pack('<d', 9)
        write('shape' + str(kind), kind, [body])
    for kind in (3, 5, 8, 13, 15, 18, 23, 25, 28, 31):
        multi = kind in (8, 18, 28)
        body = struct.pack('<I4d', kind, 1, 2, 3, 4)
        body += struct.pack('<I', 2) if multi else struct.pack('<3I', 1, 2, 0)
        if kind == 31:
            body += struct.pack('<I', 2)
        body += xy
        if kind in (13, 15, 18, 31):
            body += struct.pack('<4d', 7, 8, 7, 8)
        if kind in (13, 15, 18, 23, 25, 28, 31):
            body += struct.pack('<4d', 9, 10, 9, 10)
        write('shape' + str(kind), kind, [body])
    body = struct.pack('<I4d4I', 3, 1, 2, 3, 4, 2, 3, 0, 1)
    body += struct.pack('<6d', 1, 2, 3, 4, 2, 3)
    write('multipart', 3, [body])
    write('nullrecord', 1, [struct.pack('<I', 0), struct.pack('<I2d', 1, 1, 2)])
    (directory / 'corrupt.sav').write_bytes(b'$FL2' + b'\0' * 30)
    (directory / 'corrupt.shp').write_bytes(b'\0' * 100)
    (directory / 'corrupt.dbf').write_bytes(b'\x03\0')
    (directory / 'corrupt.xls').write_bytes(b'\xd0\xcf\x11\xe0')
    (directory / 'memory.csv').write_text('id,text\n1,first\n2,second\n')


def codec_fixtures(directory):
    vendor = ROOT / 'src/io/vendor'
    sources = [*sorted((vendor / 'readstat').rglob('*.c')),
               *sorted((vendor / 'zlib').glob('*.c')),
               # vendor/zlib/zlib.h maps ReadStat's zlib calls onto miniz.
               *sorted((ROOT / 'src/cimport/miniz').glob('*.c'))]
    with tempfile.TemporaryDirectory(prefix='ctools-format-fixtures-') as tmp:
        binary = Path(tmp) / 'fixtures'
        command = [os.environ.get('CC', 'cc'), '-std=c11', '-O2', '-DHAVE_ZLIB=1',
                   '-D_GNU_SOURCE', '-D_DARWIN_C_SOURCE', '-I', str(vendor / 'readstat'),
                   '-I', str(vendor / 'zlib'), str(ROOT / 'validation/io_format_fixtures.c'),
                   *map(str, sources), '-lm']
        if platform.system() == 'Darwin':
            command += ['-liconv']
        subprocess.run([*command, '-o', str(binary)], check=True)
        subprocess.run([str(binary), str(directory)], check=True)

def failure_fixtures(directory):
    headers = {
        'sas.sas7bdat': (directory/'long.sas7bdat').read_bytes()[:16],
        'spss.sav': b'$FL2'+bytes(12),
        'sasxport5.xpt': b'HEADER RECORD***',
        'sasxport8.v8xpt': b'HEADER RECORD***',
        'dbase.dbf': b'\x03'+bytes(15),
        'xls.xls': bytes.fromhex('d0cf11e0a1b11ae1')+bytes(8),
        'xlsx.xlsx': b'PK\x03\x04'+bytes(12),
        'shp.shp': struct.pack('>I',9994)+bytes(12),
    }
    for name, header in headers.items():
        kind, ext = name.split('.')
        (directory/f'{kind}_invalid.{ext}').write_bytes(b'invalid file header\n')
        (directory/f'{kind}_truncated.{ext}').write_bytes(header)
        (directory/f'{kind}_longinvalid.{ext}').write_bytes(b'not valid\0'*200)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plugin-path', default='')
    args = parser.parse_args()
    directory = ROOT / 'temp/io_allformats'
    directory.mkdir(parents=True, exist_ok=True)
    for sub in ('native', 'result'):
        (directory / sub).mkdir(exist_ok=True)
    shapefiles(directory)
    codec_fixtures(directory)
    failure_fixtures(directory)
    log = directory / 'regressions.log'
    log.unlink(missing_ok=True)
    driver = directory / 'driver.do'
    driver.write_text(f'''clear all
capture log close _all
log using "{log}", text replace
capture noisily do "{ROOT}/validation/validate_io_formats.do" "{args.plugin_path}"
local rc = _rc
di "IO_FORMATS_DRIVER_RC=`rc'"
log close
exit, clear
''')
    subprocess.run(['/bin/zsh', '-lic', 'stata -q -b do ' + shlex.quote(str(driver))],
                   cwd=ROOT, check=True)
    output = log.read_text()
    print('\n'.join(line for line in output.splitlines()
                    if re.match(r'^(PASS:|FAIL:|KNOWN GAP:|DETAIL:|REFERENCE FAILED:|IO_FORMATS_)', line)))
    if not re.search(r'^IO_FORMATS_DRIVER_RC=0$', output, re.M):
        raise SystemExit(f'Binary format regressions failed: {log}')


if __name__ == '__main__':
    main()
