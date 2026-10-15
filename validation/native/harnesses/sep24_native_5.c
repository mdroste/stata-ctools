
#include "ctools_matrix.c"
int main(void) {
 double X[]={1,2,3,4},D[]={1.0/30},e[]={1,-1,2,1},w[]={1,2,1,2},V[1]; int ids[]={0,0,1,1};
 ctools_vce_data d={.X_eff=X,.D=D,.resid=e,.weights=w,.N=4,.K=1,.weight_type=1,.normalize_weights=1};
 for(int cluster=0;cluster<2;cluster++) {
  calls=0;fail_at=0;V[0]=123;
  assert((cluster?ctools_vce_cluster(&d,ids,2,1,V):ctools_vce_robust(&d,1,V))==0);
  int allocations=calls;assert(V[0]>0 && V[0]!=123 && live_count==0);
  for(int f=1;f<=allocations;f++) {calls=0;fail_at=f;V[0]=123;
   assert((cluster?ctools_vce_cluster(&d,ids,2,1,V):ctools_vce_robust(&d,1,V))==920);
   assert(V[0]==123 && live_count==0);
  }
 }
 fail_at=0; memset(e,0,sizeof(e));assert(ctools_vce_robust(&d,1,V)==0 && V[0]==0);
 V[0]=123;assert(ctools_vce_cluster(&d,ids,1,1,V)==498 && V[0]==123);
 e[0]=1e308;V[0]=123;assert(ctools_vce_robust(&d,1,V)==498 && V[0]==123);
 assert(live_count==0);return 0;
}
