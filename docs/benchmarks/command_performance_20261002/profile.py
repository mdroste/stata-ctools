from pathlib import Path
import subprocess,sys,time
root=Path.cwd();out=root/'temp/command_performance_20261002';tag=sys.argv[1]
driver=out/f'profile_{tag}.do';log=out/f'profile_{tag}.log'
driver.write_text(f'''clear all
set more off
set linesize 255
adopath ++ "{out}/{tag}"
log using "{log}", text replace
set obs 1000000
set seed 928143
gen long g=floor((_n-1)/1000)
gen double key=mod(_n-1,1000)
gen double x=100+5*rnormal()
gen double shuffle=runiform()
sort shuffle
drop shuffle
crangestat (mean) result=x, interval(key -5 5) by(g) threads(8)
drop result
di "PROFILE_READY"
forvalues i=1/400 {{
    quietly crangestat (mean) result=x, interval(key -5 5) by(g) threads(8)
    quietly drop result
}}
di "PROFILE_COMPLETE"
log close
exit, clear
''')
log.unlink(missing_ok=True)
p=subprocess.Popen(['/bin/zsh','-lic',f'stata -q -b do "{driver}"'])
while p.poll() is None:
 if log.exists() and '\nPROFILE_READY\n' in log.read_text():break
 time.sleep(.1)
ps=subprocess.check_output(['ps','-axo','pid=,args='],text=True)
matches=[int(line.split()[0]) for line in ps.splitlines() if '/stata-mp ' in line and str(driver) in line]
if len(matches)!=1:raise RuntimeError('Cannot identify profile process: '+repr(matches))
subprocess.run(['sample',str(matches[0]),'5','1','-file',str(out/f'profile_{tag}.txt')],check=True)
p.wait(timeout=90)
assert '\nPROFILE_COMPLETE\n' in log.read_text()
print('PROFILE_COMPLETE',tag)
