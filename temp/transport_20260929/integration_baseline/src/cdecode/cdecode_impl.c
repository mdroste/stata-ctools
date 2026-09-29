/* Legacy wire protocol: textual label-save parsing cannot preserve Stata value
 * labels. Current cdecode.ado uses the native engine and stages all outputs.
 * Reject old wrappers instead of silently decoding corrupt/truncated text. */
#include "cdecode_impl.h"
ST_retcode cdecode_main(const char *args)
{
    (void)args;
    SF_error("cdecode: legacy plugin protocol retired; install the matching cdecode.ado\n");
    return 198;
}
ST_retcode cdecode_scan_main(const char *args)
{
    return cdecode_main(args);
}
