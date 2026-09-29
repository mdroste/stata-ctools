from pathlib import Path
import hashlib,json,os,subprocess,time
ROOT=Path('/Users/Mike/Documents/GitHub/stata-ctools')
OUT=ROOT/'temp/transport_20260929/final_native_checks_final'
SOURCE=ROOT/'temp/transport_20260929/production_final/src'
CC=os.environ.get('FINAL_CHECK_CC','/usr/bin/clang')
LIBOMP=os.environ.get('FINAL_CHECK_LIBOMP','/private/tmp/ctools-memory-audit/libomp')
SAN=os.environ.get('FINAL_CHECK_SANITIZERS','undefined')
LABEL='asan_ubsan' if 'address' in SAN else 'ubsan'
SUITES=['test_transport_native','test_transport_scheduling','test_transport_adaptive_native','test_transport_store_native']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
manifest={str(p.relative_to(SOURCE)):sha(p) for p in SOURCE.rglob('*') if p.is_file()}
(OUT/'source_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
report={'source':str(SOURCE),'source_io_sha256':sha(SOURCE/'ctools_data_io.c'),'sanitizers':SAN,'compiler':CC,'compiler_version':subprocess.check_output([CC,'--version'],text=True),'libomp':LIBOMP,'libomp_archive_sha256':sha(Path(LIBOMP)/'lib/libomp.a'),'results':[]}
env=dict(os.environ,CC=CC,LIBOMP_PREFIX=LIBOMP,CTOOLS_TEST_SOURCE=str(SOURCE),CTOOLS_SANITIZERS=SAN,KMP_BLOCKTIME='0',OMP_WAIT_POLICY='PASSIVE')
# A local macOS ASan run, if available, is memory-error evidence only; leak
# detection is not claimed. Leak detection is separately configured in Linux CI.
if 'address' in SAN:env['ASAN_OPTIONS']='detect_leaks=0:halt_on_error=1'
for suite in SUITES:
 path=ROOT/'validation'/f'{suite}.py';log=OUT/f'{suite}_{LABEL}.log';started=time.monotonic()
 with log.open('w') as handle:
  try:
   done=subprocess.run(['python3',str(path)],cwd=ROOT,env=env,stdout=handle,stderr=subprocess.STDOUT,timeout=180)
   rc=done.returncode
  except subprocess.TimeoutExpired:
   rc='timeout';handle.write('\nSUITE TIMEOUT after180 seconds\n')
 result={'suite':suite,'exit_code':rc,'elapsed_seconds':round(time.monotonic()-started,3),'log':str(log),'test_sha256':sha(path)}
 report['results'].append(result);print(json.dumps(result),flush=True)
 (OUT/f'{LABEL}_results.json').write_text(json.dumps(report,indent=2)+'\n')
print('DONE '+LABEL,flush=True)
