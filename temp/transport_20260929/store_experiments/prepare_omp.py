from pathlib import Path
import shutil,hashlib,json
out=Path(__file__).resolve().parent
base=out.parent/'baseline'/'src'
tile=out.parent/'load_experiments'/'numeric_tiles'/'src'

def store_omp(s):
 old='''    ctools_persistent_pool *chunk_pool = per_col >= 2 ? ctools_get_global_pool() : NULL;
    size_t store_tasks;
    if (chunk_pool && ctools_safe_mul_size(nvars, per_col, &store_tasks) == 0) {'''
 new='''    #ifdef _OPENMP
    int chunk_schedule_ok = per_col >= 2;
    #else
    ctools_persistent_pool *chunk_pool = per_col >= 2 ? ctools_get_global_pool() : NULL;
    int chunk_schedule_ok = chunk_pool != NULL;
    #endif
    size_t store_tasks;
    if (chunk_schedule_ok && ctools_safe_mul_size(nvars, per_col, &store_tasks) == 0) {'''
 assert s.count(old)==1;s=s.replace(old,new)
 old='''            if (ctools_persistent_pool_submit_batch(chunk_pool, store_chunk_thread,
                                    chunks, built, sizeof(ctools_store_chunk)) != 0)
                record_io_error(&error, STATA_ERR_MEMORY);
            else if (ctools_persistent_pool_wait(chunk_pool))
                record_io_error(&error, STATA_ERR_MEMORY);'''
 new='''            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            for (size_t t = 0; t < built; t++) store_chunk_thread(&chunks[t]);
            #else
'''+old+'''
            #endif'''
 assert s.count(old)==1;s=s.replace(old,new)
 old='''    rc = execute_io_parallel(thread_args, nvars, store_variable_thread);'''
 new='''    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    for (size_t j = 0; j < nvars; j++) store_variable_thread(&thread_args[j]);
    rc = STATA_OK;
    for (size_t j = 0; j < nvars; j++) {
        if (thread_args[j].error) { rc = thread_args[j].error; break; }
    }
    #else
'''+old+'''
    #endif'''
 assert s.count(old)==1;s=s.replace(old,new)
 return s

def all_omp(s):
 s=store_omp(s)
 start=s.index('    ctools_persistent_pool *pool = nvars >= 2 ? ctools_get_global_pool() : NULL;')
 end=s.index('    return STATA_OK;',start)
 old=s[start:end]
 new='''    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(nvars >= 2)
    for (size_t j = 0; j < nvars; j++) func(&args[j]);
    for (size_t j = 0; j < nvars; j++)
        if (args[j].error) return args[j].error;
    #else
'''+old+'''    #endif
'''
 s=s[:start]+new+s[end:]
 start=s.index('        ctools_persistent_pool *pool = parallel_count >= 2 ? ctools_get_global_pool() : NULL;')
 end=s.index('        for (size_t j = 0; j < parallel_count; j++) {\n            if (thread_args[j].error)',start)
 old=s[start:end]
 new='''        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(parallel_count >= 2)
        for (size_t j = 0; j < parallel_count; j++) load_filtered_variable_thread(&thread_args[j]);
        #else
'''+old+'''        #endif
'''
 s=s[:start]+new+s[end:]
 # Multi-million-row numeric chunks also stay on the same runtime.
 old='''    ctools_persistent_pool *pool = ctools_get_global_pool();
    if (!pool) return 0;'''
 new='''    #ifndef _OPENMP
'''+old+'''
    #endif'''
 assert s.count(old)==1;s=s.replace(old,new)
 old='''    if (ctools_persistent_pool_submit_batch(pool, load_chunk_thread,
                chunks, built, sizeof(ctools_io_chunk)) != 0)
        record_io_error(&error, STATA_ERR_MEMORY);
    else if (ctools_persistent_pool_wait(pool))
        record_io_error(&error, STATA_ERR_MEMORY);'''
 new='''    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    for (size_t t = 0; t < built; t++) load_chunk_thread(&chunks[t]);
    #else
'''+old+'''
    #endif'''
 assert s.count(old)==1;s=s.replace(old,new)
 return s

manifest={}
for name,parent,text in [
 ('numeric_tiles_ompstores',tile,store_omp((tile/'ctools_data_io.c').read_text())),
 ('baseline_ompcols',base,all_omp((base/'ctools_data_io.c').read_text())),
]:
 target=out/name/'src';shutil.copytree(parent,target,dirs_exist_ok=True)
 (target/'ctools_data_io.c').write_text(text)
 manifest[name]={'source':str(target.resolve()),'sha256':hashlib.sha256(text.encode()).hexdigest(),'parent':str(parent.resolve())}
(out/'omp_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
