
#include "cbinscatter/cbinscatter_resid.c"
int main(void) {
 double gram[9]={0},rhs[3]={0},beta[3];
 for(int i=0;i<5;i++){double x[]={1,1e-9*i,2e-9*i},y=3+2*i;
  for(int j=0;j<3;j++){rhs[j]+=x[j]*y;for(int k=0;k<3;k++)gram[j*3+k]+=x[j]*x[k];}}
 assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==0);
 int allocations=calls;assert(live_count==0);
 for(int i=0;i<5;i++)assert(fabs(beta[0]+1e-9*i*beta[1]+2e-9*i*beta[2]-(3+2*i))<1e-10);
 for(int f=1;f<=allocations;f++) {calls=0;fail_at=f;beta[0]=beta[1]=beta[2]=123;
  assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==920);
  assert(beta[0]==123 && beta[1]==123 && beta[2]==123 && live_count==0);}
 fail_at=0;gram[0]=-1;
 assert(cbinscatter_solve_with_collinearity(gram,rhs,3,beta)==199);
 assert(live_count==0);return 0;
}
