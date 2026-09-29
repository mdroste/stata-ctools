from pathlib import Path
base=Path(__file__).with_name('test_load_edges.py').read_text()
base=base.replace('    ctools_destroy_global_pool();\n    return 0;\n}',r'''    ctools_destroy_global_pool();
    shape=0;
    int wide_ns[]={31254,62503,4100,8195};
    int wide_ks[]={64,64,512,512};
    for(int c=0;c<4;c++) {
      nrows=wide_ns[c];ncols=wide_ks[c];filter=c%2;
      for(int fail=0;fail<=3;fail++) {
        fail_allocation=fail;
        assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_CHECK_IF)==STATA_OK);
        size_t out=0;
        for(int row=3;row<=nrows-2;row++) if(selected(row)) {
          for(int j=1;j<=ncols;j++)assert(fd.data.vars[j-1].data.dbl[out]==row+0.125*j);
          out++;
        }
        assert(out==fd.data.nobs);
        #ifdef _OPENMP
        assert(fail_allocation==0);
        #endif
        ctools_filtered_data_free(&fd);
      }
      bad_row=777;
      assert(ctools_data_load(&fd,NULL,0,3,nrows-2,CTOOLS_LOAD_CHECK_IF)==STATA_ERR_STATA_READ);
      ctools_filtered_data_free(&fd);bad_row=0;
    }
    ctools_destroy_global_pool();
    return 0;
}''',1)
exec(compile(base,str(Path(__file__).with_name('test_load_edges.py')),'exec'))
