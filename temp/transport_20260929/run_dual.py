from pathlib import Path
import json, subprocess, sys, shlex, time
root=Path(__file__).resolve().parent/'dual_impl_control'
for block in range(2):
 d=root/f'b{block}'; start=time.monotonic()
 shell='oldstata -q -b do '+shlex.quote(str(d/'run.do'))
 print('DUAL_START',block,flush=True)
 with (d/'resources.txt').open('w') as f:
  r=subprocess.run(['/usr/bin/time','-l','/bin/zsh','-lic',shell],cwd=d,stdout=f,stderr=subprocess.STDOUT)
 (d/'execution.json').write_text(json.dumps(dict(exit_code=r.returncode,elapsed_monotonic=time.monotonic()-start),indent=2)+'\n')
 if r.returncode: raise SystemExit(r.returncode)
 print('DUAL_DONE',block,flush=True)
subprocess.run([sys.executable,str(root.parent/'prepare_dual_impl.py'),'summarize',str(root)],check=True)
