"""Validate complete ABBA benchmark logs and retain every measured observation."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

CASES=('ip_sorted','ip_shuffled','ip_single','ip_skew','ip_string',
       'sp_short','sp_long','sp_multi','rj_sorted','rj_string','rj_dense','rj_wide')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    parser.add_argument('--csv',type=Path,required=True)
    parser.add_argument('--cases',nargs='+',default=CASES,choices=CASES)
    args=parser.parse_args()
    rows=[];timings=defaultdict(list)
    for label in ('b1','c1','c2','b2'):
        text=(args.directory/(label+'.log')).read_text()
        if f'\nPERFORMANCE_COMPLETE label={label}\n' not in text or '\nPERFORMANCE_DRIVER_RC=0\n' not in text:
            raise SystemExit(f'Incomplete or failed run: {label}')
        seen=set()
        for line in text.splitlines():
            if not line.startswith('PERF,'):continue
            _,run,case,threads,rep,seconds,n=line.split(',')
            threads,rep,n=int(threads),int(rep),int(n);seconds=float(seconds)
            if run!=label or case not in args.cases or threads not in (1,4,12) or rep not in range(1,7) or seconds<=0:
                raise SystemExit(f'Invalid benchmark record: {line}')
            key=(case,threads,rep)
            if key in seen:raise SystemExit(f'Duplicate benchmark record: {line}')
            seen.add(key)
            rows.append((run,case,threads,rep,seconds,n))
            timings[label[0],case,threads].append(seconds)
        if len(seen)!=len(args.cases)*3*6:raise SystemExit(f'Missing records: {label}: {len(seen)}')
    args.csv.parent.mkdir(parents=True,exist_ok=True)
    with args.csv.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(('run','case','threads','rep','seconds','output_rows'));writer.writerows(rows)
    print('| Case | Before, 12 threads | After, 12 threads | Speedup | After, 1 thread | After, 4 threads |')
    print('|---|---:|---:|---:|---:|---:|')
    for case in args.cases:
        b=median(timings['b',case,12]);c=median(timings['c',case,12])
        one=median(timings['c',case,1]);four=median(timings['c',case,4])
        print(f'| {case} | {b:.4f} | {c:.4f} | {b/c:.2f}× | {one:.4f} | {four:.4f} |')
    print(f'\nRetained all {len(rows)} observations; no timing outliers removed.')

if __name__=='__main__':
    main()
