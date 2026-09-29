#include <stdio.h>
#include <stdlib.h>
int main(void) { fputs("ASAN_MAIN\n", stderr); char *p=malloc(16); if (!p) return 1; p[0]=42; fprintf(stderr,"VALUE=%d\n",p[0]); free(p); fputs("ASAN_DONE\n", stderr); return 0; }
