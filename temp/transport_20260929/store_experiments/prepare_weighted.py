from pathlib import Path
import shutil, hashlib, json
out=Path(__file__).resolve().parent
base=out.parent/'baseline'/'src'
source=(base/'ctools_data_io.c').read_text()
helper='''/* Estimate tiny transfer work before paying pool/task wakeup overhead.
 * One numeric callback costs one unit; strings include callback overhead
 * and a term proportional to their declared bytes. Unknown widths use the
 * largest supported fixed width. Saturate the sum rather than overflowing. */
static size_t io_column_work(int is_string, size_t width)
{
    if (!is_string) return 1;
    if (width < 1 || width > STATA_STR_MAXLEN) width = STATA_STR_MAXLEN;
    return 2 + (width + 7) / 8;
}

static size_t io_add_work(size_t total, size_t column)
{
    return total > SIZE_MAX - column ? SIZE_MAX : total + column;
}

static int io_serial_work(size_t nobs, size_t nvars, size_t row_work)
{
    if (nvars < 2 || row_work == 0) return 0;
    size_t budget = IO_SERIAL_BUDGET;
    return nobs <= budget / row_work;
}

'''
anchor='static int execute_io_parallel(ctools_var_io_args *args, size_t nvars,'
assert source.count(anchor)==1
s=source.replace(anchor,helper+anchor)
old='''    ctools_persistent_pool *pool = nvars >= 2 ? ctools_get_global_pool() : NULL;'''
new='''    size_t row_work = 0;
    for (size_t j = 0; j < nvars; j++) {
        int is_string = func == load_variable_thread
            ? args[j].is_string : args[j].var->type == STATA_TYPE_STRING;
        size_t width = func == load_variable_thread
            ? (args[j].str_width > 0 ? (size_t)args[j].str_width : 0)
            : args[j].var->str_maxlen;
        row_work = io_add_work(row_work, io_column_work(is_string, width));
    }
    int serial = nvars >= 2 && io_serial_work(args[0].nobs, nvars, row_work);
    ctools_persistent_pool *pool = nvars >= 2 && !serial ? ctools_get_global_pool() : NULL;'''
assert s.count(old)==1
s=s.replace(old,new)
old='''        ctools_persistent_pool *pool = parallel_count >= 2 ? ctools_get_global_pool() : NULL;'''
new='''        size_t row_work = 0;
        for (size_t j = 0; j < parallel_count; j++) {
            int width = thread_args[j].str_width;
            row_work = io_add_work(row_work, io_column_work(thread_args[j].is_string,
                                  width > 0 ? (size_t)width : 0));
        }
        int serial = io_serial_work(n_filtered, parallel_count, row_work);
        ctools_persistent_pool *pool = parallel_count >= 2 && !serial ? ctools_get_global_pool() : NULL;'''
assert s.count(old)==1
s=s.replace(old,new)
variants={
    'tiny_weighted':s.replace('IO_SERIAL_BUDGET','32768'),
    'tiny_adaptive':s.replace('IO_SERIAL_BUDGET','nvars > (SIZE_MAX - 16384) / 256 ? SIZE_MAX : 16384 + 256 * nvars'),
}
manifest={'baseline_sha256':hashlib.sha256(source.encode()).hexdigest(),'variants':{}}
for name,text in variants.items():
    target=out/name/'src'
    shutil.copytree(base,target,dirs_exist_ok=True)
    (target/'ctools_data_io.c').write_text(text)
    manifest['variants'][name]={'source':str(target.resolve()),'sha256':hashlib.sha256(text.encode()).hexdigest()}
(out/'weighted_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
