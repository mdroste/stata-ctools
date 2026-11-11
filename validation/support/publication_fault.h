/* Test-only SPI publication interposition for a separate serial plugin.
 * Include after the command's normal headers. Never ship this header in a
 * plugin: reading injection scalars on worker threads is not supported. */
#ifdef _OPENMP
#error Publication fault-injection plugins must be built without OpenMP
#endif
static int publication_fault(char *name)
{
    double left;
    if (!SF_scal_use(name,&left) && left>0) {
        SF_scal_save(name,left-1);
        if(left<=1)return 459;
    }
    return 0;
}
static ST_retcode fault_macro_save(char *key,char *value)
{
    int rc=publication_fault("__ctools_test_macro_after");
    return rc?rc:SF_macro_save(key,value);
}
static ST_retcode fault_sstore(ST_int col,ST_int row,char *value)
{
    int rc=publication_fault("__ctools_test_store_after");
    return rc?rc:SF_sstore(col,row,value);
}
static ST_retcode fault_vstore(ST_int col,ST_int row,ST_double value)
{
    int rc=publication_fault("__ctools_test_store_after");
    return rc?rc:SF_vstore(col,row,value);
}
static ST_retcode fault_scal_save(char *key,ST_double value)
{
    int rc=publication_fault("__ctools_test_scalar_after");
    return rc?rc:SF_scal_save(key,value);
}
#undef SF_macro_save
#undef SF_sstore
#undef SF_vstore
#undef SF_scal_save
#define SF_macro_save(key,value) fault_macro_save((key),(value))
#define SF_sstore(col,row,value) fault_sstore((col),(row),(value))
#define SF_vstore(col,row,value) fault_vstore((col),(row),(value))
#define SF_scal_save(key,value) fault_scal_save((key),(value))
