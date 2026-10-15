"""Sample a warmed command for five seconds inside a bounded ten-second loop."""
import shlex
import subprocess
import sys
import time
from pathlib import Path

root = Path(__file__).resolve().parents[2]
out = Path(__file__).resolve().parent
build, command = sys.argv[1:]
folder = out / build
name = f"profile_{command}_{build}"
log = out / f"{name}.log"
driver = out / f"{name}.do"
if command == "encode":
    setup = '\n'.join(['set obs 200000', 'gen str12 x="value_"+string(mod(_n-1,10000), "%05.0f")'])
    operation = 'cencode x, gen(y) threads(8)'
    cleanup = 'drop y\nlabel drop y'
else:
    setup = '\n'.join(['set obs 1000000', 'set seed 928143', 'gen double x=rnormal()', 'gen long g=floor((_n-1)/1000)', 'gen double shuffle=runiform()', 'sort shuffle', 'drop shuffle'])
    operation = 'cwinsor x, by(g) suffix(_w) threads(8)'
    cleanup = 'drop x_w'
driver.write_text(f'''clear all
set more off
set linesize 255
log using "{log}", text replace
program define run_profile
adopath ++ "{folder}"
program perf_clock, plugin using("{out}/encode_1m/clock.plugin")
{setup}
quietly {operation}
{cleanup}
di "PROFILE_READY"
plugin call perf_clock, "start"
scalar profile_elapsed = 0
while scalar(profile_elapsed) < 10 {{
    quietly {operation}
    {cleanup}
    plugin call perf_clock, "stop"
    scalar profile_elapsed = scalar(__ctools_perf_seconds)
}}
end
capture noisily run_profile
local rc=_rc
di "PROFILE_COMPLETE_RC=`rc'"
log close
exit, clear
''')
log.unlink(missing_ok=True)
p = subprocess.Popen(['/bin/zsh','-lic', 'stata -q -b do ' + shlex.quote(str(driver))], cwd=root)
deadline = time.monotonic()+30
while p.poll() is None and time.monotonic() < deadline:
    if log.exists() and '\nPROFILE_READY\n' in log.read_text(): break
    time.sleep(.05)
ps = subprocess.check_output(['ps','-axo','pid=,args='], text=True)
matches = [int(line.split()[0]) for line in ps.splitlines() if '/stata-mp ' in line and str(driver) in line]
if len(matches) != 1: raise RuntimeError('Cannot identify profile process')
subprocess.run(['sample',str(matches[0]),'5','1','-file',str(out/f'{name}.txt')], check=True)
p.wait(timeout=30)
if '\nPROFILE_COMPLETE_RC=0\n' not in log.read_text(): raise RuntimeError('Incomplete profile')
print('PROFILE_COMPLETE', name)
