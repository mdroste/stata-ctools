
readstat_error_t readstat_convert(char *dst, size_t dst_len, const char *src, size_t src_len, iconv_t converter);

readstat_error_t readstat_convert_untrimmed(char *dst,size_t cap,const char *src,size_t len,iconv_t c);
