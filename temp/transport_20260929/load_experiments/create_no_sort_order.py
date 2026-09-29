from pathlib import Path
import shutil,re,json,hashlib,difflib
root=Path(__file__).resolve().parent.parent;base=root/'baseline/src';out=root/'load_experiments';dest=out/'no_sort_order/src';shutil.copytree(base,dest,dirs_exist_ok=True)
original={}
def edit(rel,new):
 p=dest/rel;original[rel]=p.read_text();p.write_text(new)
p=dest/'ctools_types.h';s=p.read_text();s=s.replace('#define CTOOLS_LOAD_SKIP_IF     0x01  /* Skip SF_ifobs checks, assume all observations pass */','#define CTOOLS_LOAD_SKIP_IF     0x01  /* Skip SF_ifobs checks, assume all observations pass */\n#define CTOOLS_LOAD_NO_SORT_ORDER 0x02 /* Omit identity sort_order; caller will not sort */')
s=s.replace('    - If CTOOLS_LOAD_SKIP_IF flag is set: skips SF_ifobs() checks entirely','    - If CTOOLS_LOAD_SKIP_IF flag is set: skips SF_ifobs() checks entirely\n    - If CTOOLS_LOAD_NO_SORT_ORDER is set: leaves data.sort_order NULL.\n      Use only when no subsequent operation requires a sort/permutation order.\n      Without this flag the identity sort_order is always populated.')
s=s.replace('    @param flags        [in] CTOOLS_LOAD_CHECK_IF (0) or CTOOLS_LOAD_SKIP_IF','    @param flags        [in] Bitwise OR of CTOOLS_LOAD_* flags (0 checks SF_ifobs)')
edit('ctools_types.h',s)
p=dest/'ctools_data_io.c';s=p.read_text();s=s.replace('static stata_retcode init_data_structure(stata_data *data, size_t nvars, size_t nobs)','static stata_retcode init_data_structure(stata_data *data, size_t nvars, size_t nobs, int flags)')
s=s.replace('    /* Allocate sort order array (cache-line aligned, overflow-safe) */','    /* Read-only/column-transform callers can omit the unused row order. */\n    if (flags & CTOOLS_LOAD_NO_SORT_ORDER) return STATA_OK;\n\n    /* Allocate sort order array (cache-line aligned, overflow-safe) */',1)
s=s.replace('init_data_structure(&result->data, nvars, 0)', 'init_data_structure(&result->data, nvars, 0, flags)')
s=s.replace('init_data_structure(&result->data, nvars, n_filtered)', 'init_data_structure(&result->data, nvars, n_filtered, flags)')
edit('ctools_data_io.c',s)
files=['creghdfe/creghdfe_regress.c','civreghdfe/civreghdfe_impl.c','cpplmhdfe/cpplmhdfe_irls.c','cqreg/cqreg_regress.c','cbinscatter/cbinscatter_impl.c','cpsmatch/cpsmatch_impl.c','cencode/cencode_impl.c','cdestring/cdestring_impl.c','csplit/csplit_impl.c','cexport/cexport_impl.c','cexport/cexport_xlsx.c']
for rel in files:
 p=dest/rel;s=p.read_text();start=s.index('ctools_data_load(&');depth=1;end=start+len('ctools_data_load(')
 while depth:
  if s[end]=='(':depth+=1
  if s[end]==')':depth-=1
  end+=1
 call=s[start:end];at=call.rindex(',')+1;value=call[at:-1].strip();assert value in ('0','CTOOLS_LOAD_CHECK_IF'),(rel,value)
 replacement=call[:at]+f' {value} | CTOOLS_LOAD_NO_SORT_ORDER)'
 s=s[:start]+replacement+s[end:];edit(rel,s)
core=''.join(''.join(difflib.unified_diff(original[r].splitlines(keepends=True),(dest/r).read_text().splitlines(keepends=True),fromfile='a/src/'+r,tofile='b/src/'+r)) for r in ['ctools_types.h','ctools_data_io.c'])
callers=''.join(''.join(difflib.unified_diff(original[r].splitlines(keepends=True),(dest/r).read_text().splitlines(keepends=True),fromfile='a/src/'+r,tofile='b/src/'+r)) for r in files)
(out/'no_sort_order_core.patch').write_text(core);(out/'no_sort_order_callers.patch').write_text(callers)
p=out/'manifest.json';m=json.loads(p.read_text());m['no_sort_order']={'source':str(dest.resolve()),'sha256':hashlib.sha256((dest/'ctools_data_io.c').read_bytes()).hexdigest(),'callers':files};p.write_text(json.dumps(m,indent=2)+'\n')
print(dest)
