"""Profile only a newly started benchmark through oldstata; never a user's process."""
from pathlib import Path
import subprocess,time,shlex,sys,re,json
root=Path(sys.argv[1]).resolve();driver=root/'run.do';log=root/'transport.log'
log.unlink(missing_ok=True)
start=time.monotonic()
with (root/'wrapper_resources.txt').open('w') as output:
 process=subprocess.Popen(['/usr/bin/time','-l','/bin/zsh','-lic','oldstata -q -b do '+shlex.quote(str(driver))],cwd=root,stdout=output,stderr=subprocess.STDOUT)
 pid=None
 while process.poll() is None and time.monotonic()-start<90:
  text=log.read_text() if log.exists() else ''
  if '\nTRANSPORT,' in text:
   for line in subprocess.check_output(['ps','-axo','pid,args'],text=True).splitlines():
    if '/stata-mp -q -b do '+str(driver) in line:
     pid=int(line.split()[0]);break
   if pid:break
  time.sleep(.2)
 if pid:
  sample=subprocess.run(['sample',str(pid),'5','1','-file',str(root/'cpu_profile.txt')],capture_output=True,text=True)
  (root/'sample_status.txt').write_text(sample.stdout+sample.stderr)
  print('SAMPLED',pid,'status',sample.returncode,flush=True)
 else: print('PROFILE_NOT_STARTED',flush=True)
 result=process.wait()
text=log.read_text() if log.exists() else ''
if result or '\nTRANSPORT_COMPLETE RC=0\n' not in text:raise SystemExit('profile workload failed')
(root/'profile_metadata.json').write_text(json.dumps({'pid':pid,'elapsed_monotonic':time.monotonic()-start,'exit_code':result,'sample_duration_seconds':5,'profiled_timings_excluded_from_benchmark_evidence':True},indent=2)+'\n')
print('PROFILE_COMPLETE',str(root),flush=True)
