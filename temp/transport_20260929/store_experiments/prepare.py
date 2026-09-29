from pathlib import Path
import shutil, hashlib, json
root = Path(__file__).resolve().parents[3]
base = root/'src'
out = Path(__file__).resolve().parent
source = (base/'ctools_data_io.c').read_text()
variants = {}
anchor = '''    /* Wide rows: keep a bounded row range hot while writing each column.'''
serial = '''    /* Experiment: eliminate pool dispatch for small complete stores. */
    if (nvars && nobs <= 65536 / nvars) {
        for (size_t j = 0; j < nvars; j++) {
            int err = store_single_variable(&data->vars[j],
                var_indices ? var_indices[j] : (int)j + 1, obs1, nobs, order);
            if (err) return store_status(err);
        }
        return STATA_OK;
    }

'''
assert source.count(anchor)==1
variants['store_tiny_serial'] = source.replace(anchor,serial+anchor)
start = source.index('            double gathered[1024];')
end = source.index('        } else {', start)
replacement = '''            for (i = 0; i < nobs; i++) {
                #if defined(__GNUC__) || defined(__clang__)
                if (i + 16 < nobs)
                    __builtin_prefetch(values + order[i + 16], 0, 0);
                #endif
                if (store_fn(stata_var_idx, (ST_int)(i + obs1), values[order[i]]))
                    return STATA_ERR_STATA_WRITE;
            }
'''
variants['store_direct_numeric'] = source[:start]+replacement+source[end:]
old = '''            for (size_t j = 0; j < nvars; j++) {
                size_t share = nobs / per_col, extra = nobs % per_col;
                size_t dest = 0;
                for (size_t c = 0; c < per_col; c++) {
                    size_t count = share + (c < extra);'''
new = '''            size_t share = nobs / per_col, extra = nobs % per_col;
            for (size_t c = 0; c < per_col; c++) {
                size_t dest = c * share + (c < extra ? c : extra);
                size_t count = share + (c < extra);
                for (size_t j = 0; j < nvars; j++) {'''
assert source.count(old)==1
interleaved = source.replace(old,new)
old_tail = '''                    ck->error = &error;
                    dest += count;'''
assert interleaved.count(old_tail)==1
variants['store_interleaved_chunks'] = interleaved.replace(old_tail,'''                    ck->error = &error;''')
start = source.index('stata_retcode ctools_store_filtered_rowpar(')
end = source.index('\n/*',start)
new_safety = (out/'numeric_map_safety.c').read_text()
variants['store_map_safety'] = source[:start]+new_safety+source[end:]
old = 'if (few_columns || (store_tile_ok && string_bytes >= IO_MIN_TILED_STRING_BYTES)) {'
new = 'if (few_columns || (store_tile_ok && (string_bytes >= IO_MIN_TILED_STRING_BYTES || (nvars >= 64 && string_bytes == 0)))) {'
assert source.count(old) == 1
variants['store_numeric_tiles'] = source.replace(old, new)
manifest = {'baseline_sha256':hashlib.sha256(source.encode()).hexdigest(),'variants':{}}
for name,text in variants.items():
    target=out/name/'src'
    shutil.copytree(base,target,dirs_exist_ok=True)
    (target/'ctools_data_io.c').write_text(text)
    manifest['variants'][name]={'source':str(target),'sha256':hashlib.sha256(text.encode()).hexdigest()}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
