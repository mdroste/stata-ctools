"""Sequential runner. Every Stata invocation uses the user's oldstata wrapper."""
from pathlib import Path
import json,subprocess,sys,time,shlex,os
for directory in sys.argv[1:]:
 root=Path(directory).resolve(); m=json.loads((root/'manifest.json').read_text())
 for run in m['runs']:
  log=Path(run['log'])
  if log.exists() and '\nTRANSPORT_COMPLETE RC=0\n' in log.read_text():
   print('ALREADY_COMPLETE',run['name'],flush=True);continue
  print('START',run['name'],flush=True)
  started=time.monotonic()
  shell='oldstata -q -b do '+shlex.quote(run['driver'])
  with Path(run['resources']).open('w') as f:
   r=subprocess.run(['/usr/bin/time','-l','/bin/zsh','-lic',shell],cwd=Path(run['driver']).parent,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,**m.get('requested_runtime_environment',{})))
  text=log.read_text() if log.exists() else ''
  (Path(run['driver']).parent/'execution.json').write_text(json.dumps({'elapsed_monotonic_seconds':time.monotonic()-started,'exit_code':r.returncode,'requested_runtime_environment':m.get('requested_runtime_environment',{})},indent=2)+'\n')
  success=r.returncode==0 and '\nTRANSPORT_COMPLETE RC=0\n' in text
  print('COMPLETE' if success else 'FAILED',run['name'],round(time.monotonic()-started,2),flush=True)
  if not success:
   print(text[-3000:]);raise SystemExit(1)
 summary = [sys.executable,m['summarize_script']] if m.get('summarize_script') else [sys.executable,'validation/benchmark_transport_campaign.py','summarize',str(root)]
 subprocess.run(summary,check=True)
