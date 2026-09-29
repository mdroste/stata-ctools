/* Native homogeneous simplex separation on residualized suspect regressors.
 * Included by cpplmhdfe_irls.c to share its projection/precision helpers. */
static int ppml_simplex_tableau(double *X,int n,int k,int maxiter,int *separated)
{
    double *A=calloc((size_t)n*k,sizeof(double));
    double *column=malloc((size_t)n*sizeof(double));
    double *row=malloc((size_t)k*sizeof(double));
    double *cost=malloc((size_t)n*sizeof(double));
    int *basic=malloc((size_t)n*sizeof(int));
    int *nonbasic=malloc((size_t)k*sizeof(int));
    int *used=calloc((size_t)n,sizeof(int));
    int rc=0,rank=0,m=0;
    if(!A||!column||!row||!cost||!basic||!nonbasic||!used){rc=-920;goto done;}
    for(int j=0;j<k;j++) {
        int p=-1;
        for(int i=0;i<n;i++)if(fabs(X[(size_t)i*k+j])>1.4210854715202004e-14){p=i;break;}
        if(p<0)continue;
        double pivot=-1/X[(size_t)p*k+j];
        for(int l=0;l<k;l++)row[l]=X[(size_t)p*k+l]*pivot;
        for(int i=0;i<n;i++)column[i]=X[(size_t)i*k+j];
        A[(size_t)p*k+rank]=1;
        for(int i=0;i<n;i++)if(i!=p)for(int l=0;l<=rank;l++) {
            double z=A[(size_t)i*k+l]+column[i]*pivot*A[(size_t)p*k+l];
            A[(size_t)i*k+l]=fabs(z)<1.4210854715202004e-14?0:z;
        }
        memset(A+(size_t)p*k,0,k*sizeof(double));
        for(int i=0;i<n;i++)for(int l=0;l<k;l++) {
            double z=X[(size_t)i*k+l]+column[i]*row[l];
            X[(size_t)i*k+l]=fabs(z)<1.4210854715202004e-14?0:z;
        }
        nonbasic[rank++]=p;used[p]=1;
    }
    for(int i=0;i<n;i++){cost[i]=1;if(!used[i])basic[m++]=i;}
    if(!rank)goto done;
    if(!m){for(int i=0;i<n;i++)separated[i]=1;rc=n;goto done;}
    for(int i=0;i<m;i++)for(int j=0;j<rank;j++)X[(size_t)i*rank+j]=A[(size_t)basic[i]*k+j];
    for(int iter=0;iter<maxiter;iter++) {
        int enter=-1;double gain=0;
        for(int j=0;j<rank;j++) {
            double reduced=cost[nonbasic[j]];
            for(int i=0;i<m;i++)reduced-=cost[basic[i]]*X[(size_t)i*rank+j];
            if(reduced>gain+1.4210854715202004e-14){gain=reduced;enter=j;}
        }
        if(enter<0){for(int i=0;i<n;i++){separated[i]=cost[i]==0;rc+=separated[i];}goto done;}
        int leave=-1;double maximum=0;
        for(int i=0;i<m;i++)if(X[(size_t)i*rank+enter]>maximum){maximum=X[(size_t)i*rank+enter];leave=i;}
        if(leave<0) {
            cost[nonbasic[enter]]=0;
            for(int i=0;i<m;i++)if(X[(size_t)i*rank+enter]<0)cost[basic[i]]=0;
            continue;
        }
        double pivot=X[(size_t)leave*rank+enter];
        for(int j=0;j<rank;j++)row[j]=X[(size_t)leave*rank+j]/pivot;
        for(int i=0;i<m;i++)column[i]=X[(size_t)i*rank+enter];
        for(int i=0;i<m;i++)for(int j=0;j<rank;j++) {
            double z=i==leave?row[j]:X[(size_t)i*rank+j]-column[i]*row[j];
            X[(size_t)i*rank+j]=fabs(z)<1.7763568394002505e-15?0:z;
        }
        for(int i=0;i<m;i++)X[(size_t)i*rank+enter]=i==leave?1/pivot:-column[i]/pivot;
        int temp=basic[leave];basic[leave]=nonbasic[enter];nonbasic[enter]=temp;
    }
    rc=-430;
done:
    free(A);free(column);free(row);free(cost);free(basic);free(nonbasic);free(used);
    return rc;
}

static int ppml_simplex(HDFE_State *S,const double *y,const double *X,const double *w_user,
    int weight_type,int N,int K,int *separated)
{
    double *data=malloc((size_t)N*K*sizeof(double));
    double *full=malloc((size_t)N*K*sizeof(double));
    int *suspect=malloc((size_t)K*sizeof(int));
    int *rows=malloc((size_t)N*sizeof(int));
    int *mask=calloc((size_t)N,sizeof(int));
    double *tableau=NULL;
    int rc=0,ns=0,nz=0;
    double saved_tol=S->tolerance;
    memset(separated,0,(size_t)N*sizeof(int));
    if(!data||!full||!suspect||!rows||!mask){rc=-920;goto done;}
    for(int i=0;i<N;i++){S->weights[i]=w_user?w_user[i]:1;if(y[i]==0)rows[nz++]=i;}
    if(!nz)goto done;
    S->tolerance=1e-9;
    ppml_update_fe_weights(S,S->weights,N);
    memcpy(full,X,(size_t)N*K*sizeof(double));
    HDFE_SolveResult pr=ppml_reference_partial(S,full,N,K,1);
    if(pr.status){rc=-pr.status;goto done;}
    for(int i=0;i<N;i++)S->weights[i]=y[i]>0?(weight_type==2&&w_user?w_user[i]:1):0;
    ppml_update_fe_weights(S,S->weights,N);
    memcpy(data,X,(size_t)N*K*sizeof(double));
    /* The positive-sample null space must be resolved more accurately than
     * the sign threshold; otherwise projection noise can flag extra zeros. */
    S->tolerance=fmin(1e-12,ppml_options.simplex_tol*.1);
    pr=partial_out_columns(S,data,N,K,1);
    if(pr.status){rc=-pr.status;goto done;}
    double largest=0;
    for(int j=0;j<K;j++)largest=fmax(largest,ppml_quad_dot(data+(size_t)j*N,data+(size_t)j*N,S->weights,N));
    for(int j=0;j<K;j++) {
        if(ppml_quad_dot(full+(size_t)j*N,full+(size_t)j*N,NULL,N)==0)continue;
        double norm=ppml_quad_dot(data+(size_t)j*N,data+(size_t)j*N,S->weights,N);
        if(norm<=fmax(DBL_MIN,largest*1e-10)){suspect[ns++]=j;continue;}
        norm=sqrt(norm);
        for(int i=0;i<N;i++)data[(size_t)j*N+i]/=norm;
        for(int k=j+1;k<K;k++)for(int pass=0;pass<2;pass++) {
            double dot=ppml_quad_dot(data+(size_t)j*N,data+(size_t)k*N,S->weights,N);
            for(int i=0;i<N;i++)data[(size_t)k*N+i]-=dot*data[(size_t)j*N+i];
        }
    }
    if(!ns)goto done;
    tableau=malloc((size_t)nz*ns*sizeof(double));
    if(!tableau){rc=-920;goto done;}
    for(int i=0;i<nz;i++)for(int j=0;j<ns;j++) {
        double z=data[(size_t)suspect[j]*N+rows[i]];
        tableau[(size_t)i*ns+j]=fabs(z)<=ppml_options.simplex_tol?0:z;
    }
    rc=ppml_simplex_tableau(tableau,nz,ns,ppml_options.simplex_maxiter,mask);
    if(rc>=0)for(int i=0;i<nz;i++)separated[rows[i]]=mask[i];
done:
    S->tolerance=saved_tol;
    free(data);free(full);free(suspect);free(rows);free(mask);free(tableau);
    return rc;
}
