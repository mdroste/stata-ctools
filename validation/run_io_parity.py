"""Exact Stata differential regressions; Stata execution uses the stata shell alias."""
from stata_runner import launch_driver
from pathlib import Path
import argparse
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def multichunk():
    """Over 2 MB, so cimport parses it in parallel chunks: quoted delimiters,
    padded and leading-zero numbers, exponents, ragged rows, and a column
    that becomes text only near the end of the file."""
    rows = ['id,num,padded,lead0,sci,text,late,ragged']
    for i in range(1, 60001):
        text = '"has, comma"' if i % 7 == 0 else '"q ""x"" y"' if i % 11 == 0 else f'plain{i}'
        # id stays in int range: native's %12.0g-for-long heuristic is a separate gap.
        row = [str(i % 30000 + 1), f'{(i * 7919 % 100003) / 97:.6f}', f' {i % 500} ', '000' + str(i % 1000),
               f'{i % 97 + 1}.5e-{i % 30 + 1}', text, 'word' if i == 59995 else str(i % 1000), str(i % 3)]
        if i % 1000 == 0:
            row = row[:5]
        elif i % 1500 == 0:
            row.append('extra')
        rows.append(','.join(row))
    return '\n'.join(rows) + '\n'

def fixtures():
    p = ROOT / 'temp/io_parity'
    p.mkdir(parents=True, exist_ok=True)
    cases = {
        'headers': 'First Name,1,,a-b,a b,a_b,UPPER,if,_all,über\n1,2,3,4,5,6,7,8,9,10\n11,12,13,14,15,16,17,18,19,20\n',
        'noheaders': '1,2,abc\n3,4,def\n',
        'numbers': 'small,decimal,precise,missing,na,empty\n1,1.1,1.12345678901234,.a,NA,\n2,1.2,2.12345678901234,.z,1,\n',
        'quoted': 'value,text\n"1.25","hello"\n"2.5","a""b"\n',
        'exponents': 'value\n1e-2\n2e-3\n',
        'outlier': 'id,value\n' + ''.join(f'{i},{1000000000 if i == 10001 else 1}\n' for i in range(1, 20002)),
        # Whitespace-only fields are text; padded numbers and missings are numeric
        # (padding turns .a into system missing).
        'whitespace': 'sp1,lead,trail,both,sp2,tab,qsp,qempty,spa,atr,xa,qa,trailtab\n'
                      '1,1,1,1,1,1,1,1,1,1,1,1,1\n'
                      ' , 2,2 , 2 ,  ,\t," ",""," .a",.a ,.a,".a",2\t\n'
                      '3,3,3,3,3,3,3,3,3,3,3,3,3\n',
        'digits': 'lead,sci\n0000000000000000007.9770,2.59462e-18\n00000000000000.0851289740626,7.7651552348e-14\n'
                  '1167.411680745998,-4.356865436322338e-21\n0.00000000000000000012,3.474157e-30\n',
        'multichunk': multichunk(),
    }
    for name, data in cases.items():
        (p / (name + '.csv')).write_text(data)
    return p

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plugin-path', default='')
    args = parser.parse_args()
    p = fixtures()
    log = p / 'regressions.log'
    log.unlink(missing_ok=True)
    driver = p / 'driver.do'
    driver.write_text(f'''clear all
capture log close _all
log using "{log}", text replace
capture noisily do "{ROOT}/validation/validate_io_parity.do" "{args.plugin_path}"
local rc = _rc
di "IO_PARITY_DRIVER_RC=`rc'"
log close
exit, clear
''')
    launch_driver(driver, log)
    text = log.read_text()
    print('\n'.join(line for line in text.splitlines() if re.match(r'^(PASS:|FAIL:|METADATA:|REFERENCE FAILED:|IO_PARITY_)', line)))
    if not re.search(r'^IO_PARITY_DRIVER_RC=0$', text, re.M):
        raise SystemExit(f'Parity regressions failed: {log}')

if __name__ == '__main__':
    main()
