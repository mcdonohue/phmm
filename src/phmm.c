/* phmm.c */

/* Input data should repeat the col's that corresp to random effects */
/* Ties are split */
/* To compile, use
gcc -lm reff2.c reffaux.c arms.c -lm
or
gcc -lm reff.c reffvar.c remnewt.c  /home/rxu/recipes/lubksb.c  /home/rxu/recipes/ludcmp.c  /home/rxu/recipes/nrutil.c arms.c  /home/rxu/recipes/mylubksb.c  /home/rxu/recipes/myludcmp.c
*/
/* Format of data: id, x=event time, delta=event indicator, cluster,
z=covariates (fixed effect), w= r.eff. covariates */
/*   x, delta, z, w are sorted by x=time */
/*   xx, ddelta, zz, ww, cluster, rank are sorted by cluster */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "phmm.h"
#include "time.h"

void phmm(double *x, double *zv, double *wv, int *delta,
  int *cluster,
  int *svarcov,
  double *sSigma, double *sinvSigma, double *sdetSigma,  
  int *sNINIT, 
  int *sMAXSTEP,
  double *smaxtime,
  int *sCONVERG,
  int *semstep,
  int *sGbs,
  int *sGbsvar,
  int *snobs,
  int *sncov,
  int *snreff,
  double *bhat2,
  double *sdbhat2,
  double *steps2,
  double *svar,
  int *sverbose,
  double *bridgeC, 
  double *laplacetot,
  double *llaplace, 
  double *limport, 
  double *lbridge,
  double *lambda,
  double *Lambda)
{
  int i, j, k, nclust;

  int varcov = *svarcov;
  int NINIT = *sNINIT;
  int MAXSTEP = *sMAXSTEP;
  double maxtime = *smaxtime;
  int CONVERG = *sCONVERG;
  int emstep = *semstep;
  int Gbs = *sGbs;
  int Gbsvar = *sGbsvar;
  int nreff = *snreff;
  int nobs = *snobs;
  int ncov = *sncov;
  int verbose = *sverbose;

  double *laplace;
	
  double *xx, *zz[ncov+1], *ww[nreff+1], *B[nreff+1], *a[nreff+1];
  double *z[ncov+1], *w[nreff+1], *omega, *Lambexp;
  int *ddelta, *rank, *clust_start, *clust_sizes;
  double **bhat, **sdbhat, **steps, **var, **Sigma, **invSigma;
  /* to be returned by R wrapper */
	
  /*   EM variables */
  double detSigma = *sdetSigma;
  double *sum0, *sum1[ncov+1], *sum2[ncov+1][ncov+1];
  double *eb[nreff+1], *v[nreff+1], *condvar[nreff+1][nreff+1];
  double *invcondv[nreff+1][nreff+1], *detcondv;
  double *sumb[nreff+1], *sumbb[nreff+1], *sumv[nreff+1][nreff+1];
  float betahat[ncov+1];
	
  xx= (double *)R_alloc((nobs+1), sizeof(double));
  Lambexp= (double *)R_alloc((nobs+1), sizeof(double));
  omega= (double *)R_alloc((nobs+1), sizeof(double));
  for (i=1;i<=ncov;i++) {
    z[i]= (double *)R_alloc((nobs+1), sizeof(double));
    zz[i]= (double *)R_alloc((nobs+1), sizeof(double));
  }
  for (i=1;i<=nreff;i++) {
    w[i]= (double *)R_alloc((nobs+1), sizeof(double));
    ww[i]= (double *)R_alloc((nobs+1), sizeof(double));
  }
	
  ddelta= (int *)R_alloc((nobs+1), sizeof(int));
  rank= (int *)R_alloc((nobs+1), sizeof(int));
  clust_sizes= (int *)R_alloc((nobs+1), sizeof(int));
	
  k = 1;
  for (i=1;i<=ncov;i++) {
    for (j=1;j<=nobs;j++) {
    	z[i][j]=zv[k];
    	k++;
    }
    k++;
  }
	
  k = 1;
  for (i=1;i<=nreff;i++) {
    for (j=1;j<=nobs;j++) {
    	w[i][j]=wv[k];
    	k++;
    }
    k++;
  }
	
	/*for (i=1; i<=nobs; i++)  Rprintf("%.5f %d %d %.1f %.1f %.1f %.1f %.1f %.1f\n", 
  x[i], delta[i], cluster[i],
  z[1][i], z[2][i], z[3][i], z[4][i], z[5][i], 
  w[1][i]);*/
	
  sum0= (double *)R_alloc((nobs+2), sizeof(double));
  for (i=1;i<=ncov;i++) {
    sum1[i]= (double *)R_alloc((nobs+2), sizeof(double));
    for (j=1;j<=ncov;j++)
    	sum2[i][j]= (double *)R_alloc((nobs+2), sizeof(double));
  }
	
  sum0[nobs+1]=0;   /* set all sums to 0 */
  for (i=1;i<=ncov;i++) {
    sum1[i][nobs+1]=0;
    for (j=1;j<=ncov;j++)
    	sum2[i][j][nobs+1]=0;
  }
  Lambda[0]=0;
	
  /* Sort data: */
	Sort( nobs, x, delta, cluster, ncov, z, nreff, w );
  Clust( nobs, rank, xx, x, ddelta, delta, ncov, zz, z, nreff, ww,
    w, cluster, clust_sizes, &nclust, verbose);
	/* sort data by cluster and create
  cluster_sizes[] */
	
  clust_start= (int *)R_alloc((nclust+2),sizeof(int));
  for (i=1;i<=nreff;i++) {
    B[i]= (double *)R_alloc((nclust+1),sizeof(double));
    a[i]= (double *)R_alloc((nclust+1),sizeof(double));
  }
  clust_start[1]=1;
  for (i=1; i<=nclust; i++) {
    for (j=1; j<=nreff; j++) {
      B[j][i]=a[j][i]=0;
      for (k=clust_start[i]; k<clust_start[i]+clust_sizes[i]; k++)
        a[j][i]+=ddelta[k]*ww[j][k];
      }
    clust_start[i+1] = clust_start[i] + clust_sizes[i];
  }
	
  bhat = dmatrix2(bhat2, nreff, nclust);
  sdbhat = dmatrix2(sdbhat2, nreff, nclust);
  steps = dmatrix2(steps2, ncov+nreff, MAXSTEP+1);
  Sigma = dmatrix2(sSigma, nreff+1, nreff+1);
  invSigma = dmatrix2(sinvSigma, nreff+1, nreff+1);
  var = dmatrix2(svar, ncov+nreff+nobs+1, ncov+nreff+nobs+1);
	
	detcondv= (double *)R_alloc((nclust+1),sizeof(double));
  for (i=1;i<=nreff;i++) {
    sumb[i]= (double *)R_alloc((nclust+1),sizeof(double));
    sumbb[i]= (double *)R_alloc((nclust+1),sizeof(double));
    eb[i]= (double *)R_alloc((nclust+1),sizeof(double));
    v[i]= (double *)R_alloc((nclust+1),sizeof(double));
    for (j=1;j<=nreff;j++) {
    	sumv[i][j]= (double *)R_alloc((nclust+1),sizeof(double));
    	condvar[i][j]= (double *)R_alloc((nclust+1),sizeof(double));
    	invcondv[i][j]= (double *)R_alloc((nclust+1),sizeof(double));
    }
  }
	laplace= (double *)R_alloc((nclust+1),sizeof(double));
	
  /* EM iterative algorithm: */
  EM( ncov, nreff, Sigma, nobs, omega,
    betahat, sum0, sum1, z, delta,
    sum2, lambda, Lambda, Lambexp,
    steps, NINIT, emstep, MAXSTEP, maxtime, nclust,
    sumb, sumbb, B, clust_start,
    ww, Gbs, Gbsvar, a, rank, CONVERG, bhat,
    sdbhat, sumv, invSigma, eb, v, condvar, detSigma, 
    verbose, varcov);
	
	Var(Gbsvar, nobs, nreff, ncov, NINIT, betahat,
    nclust, sumb, sumbb, B, clust_start,
    ww, Lambexp, lambda, a, 
    Sigma, invSigma, rank,
    zz, delta, ddelta, omega, z, var);
	
	Likelihood(ncov, nreff, Sigma, nobs, omega,
	  betahat, sum0, sum1, z, delta, ddelta,
	  sum2, lambda, Lambda, Lambexp,
	  steps, NINIT, emstep, MAXSTEP, nclust,
	  sumb, sumbb, B, clust_start,
	  ww, zz, Gbs, Gbsvar, a, rank, CONVERG, bhat,
	  sdbhat, sumv, invSigma, eb, v, condvar, detSigma, 
	  verbose, varcov, llaplace, limport, lbridge, 
	  laplacetot, bridgeC, invcondv, detcondv, laplace);
	
	return;
}
/* end phmm() */


/* Function definitions:
=====================*/


void Sort(int nobs, double *x, int *delta, int *cluster,
  int ncov, double **z, int nreff, double **w ) /* x, delta, z, w are sorted by time */
{
	double temp1;
	int temp2, i, j, d;
	
	
  for ( i=1; i<=nobs-1; i++ ) {
    for ( j=i+1; j<=nobs; j++ ) {
      if ((x[i]>x[j]) || ((x[i]==x[j])&&(delta[i]<delta[j]))) {
      temp1=x[i]; x[i]=x[j]; x[j]=temp1;
      temp2=delta[i]; delta[i]=delta[j]; delta[j]=temp2;
      temp2=cluster[i]; cluster[i]=cluster[j]; cluster[j]=temp2;
      for (d=1;d<=ncov;d++) {
      temp1=z[d][i]; z[d][i]=z[d][j]; z[d][j]=temp1;
      }
      for (d=1;d<=nreff;d++) {
      temp1=w[d][i]; w[d][i]=w[d][j]; w[d][j]=temp1;
      }
      }
    }
  }
}   /* end Sort() */


void Clust(int nobs, int *rank, double *xx, double *x, int *ddelta, int *delta,
  int ncov, double **zz, double **z, int nreff, double **ww,
  double **w, int *cluster, int *clust_sizes, int *pnclust, int verbose)
  /* xx, ddelta, cluster, rank, zz, ww are sorted by clusters */
{
	double temp1;
	int temp2, i, j, d, k, flag;
	int nclust = *pnclust;
	
	for (i=1; i<=nobs; i++) {
  rank[i]=i;
  xx[i]=x[i]; ddelta[i]=delta[i];
  for (d=1; d<=ncov; d++)
  	zz[d][i]=z[d][i];
  for (d=1; d<=nreff; d++)
  	ww[d][i]=w[d][i];
	}
	
	for ( i=1; i<=nobs-1; i++ ) {
    for ( j=i+1; j<=nobs; j++ ) {
      	if (cluster[i]>cluster[j]) {
        temp1=xx[i]; xx[i]=xx[j]; xx[j]=temp1;
        temp2=ddelta[i]; ddelta[i]=ddelta[j]; ddelta[j]=temp2;
        temp2=cluster[i]; cluster[i]=cluster[j]; cluster[j]=temp2;
        temp2=rank[i]; rank[i]=rank[j]; rank[j]=temp2;
        for (d=1;d<=ncov;d++) {
          temp1=zz[d][i]; zz[d][i]=zz[d][j]; zz[d][j]=temp1;
        }
        for (d=1;d<=nreff;d++) {
          temp1=ww[d][i]; ww[d][i]=ww[d][j]; ww[d][j]=temp1;
        }
      }
    }
	}
	
	/* for (i=1; i<=nobs; i++)  Rprintf("%d %.1f %.1f\n", cluster[i], zz[1][i], zz[2][i]);
	Rprintf("\n"); */
	
	/* Create the vector clust_sizes: */
	j=k=1;
	for (i=1; i<=nobs-1; i++) {
  if (cluster[i]==cluster[i+1]) {
  	k++; flag=1;   /* flag=1 means the cluster contains more than one obs*/
  	continue;
  }
  clust_sizes[j]=k;
  j++; k=1; flag=0;
	}
	/* next two lines takes care of the last cluster */
	if (flag) { clust_sizes[j]=k; nclust=j; }
	else { clust_sizes[j]=1; nclust=j; }
	if(verbose){
  Rprintf("\nNo. of clusters: %d\n  Cluster sizes: ", nclust);
  for (i=1; i<=nclust; i++) Rprintf("%d ", clust_sizes[i]);
  Rprintf("\n");
	}
	
	*pnclust = nclust;
	
	/*   for (i=1; i<=nobs; i++)  Rprintf("%.5f %d %d %.1f %.1f %.1f %.1f %.1f %.1f\n", 
  xx[i], ddelta[i], cluster[i],
  zz[1][i], zz[2][i], zz[3][i], zz[4][i], zz[5][i], 
  ww[1][i]);
Rprintf("\n"); */
	
}   /* end Clust() */


void EM(int ncov, int nreff, double *Sigma[nreff+1],
  int nobs, double *omega,
  float *betahat, double *sum0, double *sum1[ncov+1], double **z, int *delta,
  double *sum2[ncov+1][ncov+1], double *lambda, double *Lambda, double *Lambexp,
  double **steps, int NINIT, int emstep, int MAXSTEP, double maxtime, int nclust,
  double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, int Gbs, int Gbsvar, double **a, int *rank, int CONVERG, double **bhat,
  double **sdbhat, double *sumv[nreff+1][nreff+1], double *invSigma[nreff+1], 
  double *eb[nreff+1], double *v[nreff+1], double *condvar[nreff+1][nreff+1], 
  double detSigma, int verbose, int varcov )
{
	int i, d;
	double dif;
	time_t start,end;
	
	time (&start);

	for ( i=1; i<=nobs; i++ )
  omega[i]=1;  /* omega is expect. of exp(b*w) */
	for (d=1; d<=ncov; d++)
  betahat[d]=0;
	
	Betahat(betahat, ncov, nobs, omega, sum0,
  	sum1, z, delta, sum2,
  	lambda, Lambda, Lambexp);   /* starting value for betahat and lambda (Cox) */
	
	if(verbose){
  Rprintf("\nstep  ");
  for (d=1; d<=ncov; d++)   Rprintf(" beta%d ", d);
  for ( d=1; d<=nreff; d++ )   Rprintf(" var%d ", d);
  
  /*if(nreff>1)
  	for ( d=1; d<=nreff; d++ )
  	for ( dd=1; dd<=nreff; dd++ )
  	if (d<dd) Rprintf("   cov%d%d ", d, dd);*/
  Rprintf(":\n   0 ");
	}
	for (d=1; d<=ncov; d++){
  if(verbose) Rprintf("%.4f ", betahat[d]);
  steps[emstep][d-1]=betahat[d];
	}
	for ( d=1; d<=nreff; d++ ){
  if(verbose) Rprintf("%.4f ", Sigma[d][d]);
  steps[emstep][ncov+d-1]=Sigma[d][d];
	}
	/*if(nreff>1)
	for ( d=1; d<=nreff; d++ )
	for ( dd=1; dd<=nreff; dd++ )
	if (d<dd){
  if(verbose) Rprintf("%.4f ", Sigma[d][dd]);
  steps[emstep][ncov+d-1]=Sigma[d][dd];
	}*/
	if(verbose) Rprintf("\n");
  	
	while (1) {
  /* Test for EM convergence: */
  emstep++;
  if (emstep > MAXSTEP) {
  	/*if(verbose) Rprintf("\n");
  	if(verbose) Rprintf("\n E(b) and sd:\n");*/
  	for (i=1; i<=nclust; i++) {
  for (d=1; d<=nreff; d++){
  /*if(verbose) Rprintf("%f %f  ", sumbb[d][i]/Gbs, sqrt(sumb[d][i]/Gbs - (sumbb[d][i]/Gbs)*(sumbb[d][i]/Gbs)) );*/
  bhat[i-1][d-1]=sumbb[d][i]/Gbs;
  sdbhat[i-1][d-1]=sqrt(sumb[d][i]/Gbs - (sumbb[d][i]/Gbs)*(sumbb[d][i]/Gbs));
  }
  /*if(verbose) Rprintf("\n");*/
  	}
  	return;
  }
  /* Gibbs E-step: */
  
  Estep(ncov, nreff, Sigma, nobs, NINIT, nclust,
  	  sumb, sumbb, sumv, B, ww, Gbs, rank,
  	  CONVERG, clust_start, omega, Lambexp,
  	  a, invSigma, eb, v, condvar, detSigma, varcov );   /* also updates Sigma (M-step) */
  	  
  /* M-step for beta */
  Betahat(betahat, ncov, nobs, omega, sum0,
  sum1, z, delta, sum2, lambda, Lambda, Lambexp);   /* updates betahat and lambda */
  
  /* Print some results: */
  if(verbose) Rprintf("%4d ", emstep);
  for (d=1; d<=ncov; d++){
  	if(verbose) Rprintf("%.4f ", betahat[d]);
  	steps[emstep][d-1]=betahat[d];
  }
  for ( d=1; d<=nreff; d++ ){
  	if(verbose) Rprintf("%.4f ", Sigma[d][d]);
  	steps[emstep][ncov+d-1]=Sigma[d][d];
  }
  /*if(nreff>1)
   for ( d=1; d<=nreff; d++ )
   for ( dd=1; dd<=nreff; dd++ )
   if (d<dd){
   if(verbose) Rprintf("%.4f ", Sigma[d][dd]);
   steps[emstep][ncov+d-1]=Sigma[d][dd];
   }*/
  if(verbose) Rprintf("\n");
  
  time (&end);
  dif = difftime (end,start);
  if(dif > maxtime)
  	emstep = CONVERG++;
  
  /* Update Gbs */
  if (emstep > CONVERG)
  	Gbs = Gbsvar;
	}
}   /* end EM() */


double BetaZ( int i, float *bbeta, int ncov, double **z )
{   /* used by betahat() */
  int d;
  double a=0;

  for (d=1;d<=ncov;d++){
    a+=z[d][i]*bbeta[d];
  }
  return a;
}

void Betahat(float *betahat, int ncov, int nobs, double *omega, double *sum0,
  double *sum1[ncov+1], double **z, int *delta, double *sum2[ncov+1][ncov+1],
  double *lambda, double *Lambda, double *Lambexp)  /* program for Cox regression */
{
	int i;

	
	/** SHOULDN'T THE OMEGAS BE INPUT HERE RATHER THAN IN Estep() ? **/
	mnewt(ncov, 50, betahat, 0.0001*ncov, 0.0001*ncov,
  nobs, omega, sum0, sum1, z, delta, sum2);
	
	/*   for (d=1; d<=ncov; d++)
  Rprintf("betahat%d = %.4f \n", d, betahat[d]);
	*/
	
	for (i=1; i<=nobs; i++) {
    lambda[i] = delta[i]/sum0[i];
    Lambda[i] = Lambda[i-1]+lambda[i];
    Lambexp[i] = Lambda[i]*exp(BetaZ(i, betahat, ncov, z));
	}
}
  	 
double logdens(double bb, void *para)
/* log-density for */
{   /* the Gibbs E-step */
  struct dens_para *dpar;
  double y, sum=0, *alpha, *b;
  int i, d, l, nreff;
  int *rank, *clust_start;
  double *Lambexp, **ww, **a, **Sigma, **invSigma;

  dpar = para;
  i = dpar->i;   /* i = cluster */
  d = dpar->d;   /* d = random effect */
  alpha = dpar->alpha;
  b = dpar->b;
  clust_start = dpar->clust_start;
  Lambexp = dpar->Lambexp;
  rank = dpar->rank;
  ww = dpar->ww;
  a = dpar->a;
  Sigma = dpar->Sigma;
  nreff = dpar->nreff;
  invSigma = dpar->invSigma;

  /*Rprintf("bb=%.4f\n", bb);*/

  for (l=clust_start[i]; l<clust_start[i+1]; l++)
  sum += Lambexp[rank[l]]*exp(alpha[l]+bb*ww[d][l]-b[d]*ww[d][l]);

  y = a[d][i]*bb - sum - bb*bb/(2*Sigma[d][d]);
  return y;
}

void Estep(int ncov, int nreff, double *Sigma[nreff+1], int nobs, int NINIT, int nclust,
  double **sumb, double **sumbb, double *sumv[nreff+1][nreff+1], double **B, double **ww,
  int Gbs, int *rank, int CONVERG, int *clust_start, double *omega, double *Lambexp,
  double **a, double *invSigma[nreff+1], double *eb[nreff+1],
  double *v[nreff+1], double *condvar[nreff+1][nreff+1],
  double detSigma, int varcov)
{
	int i, l, d, in, dd;
	int err, ninit = 4, dometrop = 0,  npoint = 100, ncent = 0,
  neval, nsamp=1;
	double xl = -100.0, xr = 100.0, xprev = 0.0, xsamp,
  xcent, qcent, convex = 1., temp;
	double b[nreff+1], sum[nreff+1], alpha[nobs+1], oomega[nobs+1];
	struct dens_para para;
	double myxinit[nreff+1][NINIT];
	double *vtemp[nreff+1];
	
	for (d=1;d<=nreff;d++)
  vtemp[d]= (double *)R_alloc((nreff+1),sizeof(double));
	
	/* Initializations: */
	for (l=1; l<=nobs; l++)   oomega[l] = 0;
	for (d=1; d<=nreff; d++)   sum[d]=0;
	for (d=1; d<=nreff; d++)
  for (i=1; i<=nclust; i++){
  	sumb[d][i] = sumbb[d][i] = 0;
  	/* sumb for E(b^2) and sumbb for E(b) */
  	for (dd=1; dd<=nreff; dd++)  sumv[d][dd][i] = 0;
  }
  	
  	/* Gibbs sampler: */  	
  	{
  register int i, d, g;
  for (i=1; i<=nclust; i++) {
  for (d=1; d<=nreff; d++)
  b[d] = B[d][i];
  
  for (l=clust_start[i]; l<clust_start[i+1]; l++) {
  alpha[l] = 0;
  for (d=1; d<=nreff; d++)
  alpha[l] += b[d]*ww[d][l];
  }
  for (g=1; g<=Gbs; g++) {  /* MCMC loop  */
  for (d=1; d<=nreff; d++) {
  /* initialize parameters for arms(): */
  para.i = i;
  para.d = d;
  para.nreff = nreff;
  para.alpha = alpha;
  para.b = b;
  para.clust_start = clust_start;
  para.Lambexp = Lambexp;
  para.ww = ww;
  para.a = a;
  para.Sigma = Sigma;
  para.rank = rank;
  para.invSigma = invSigma;
  
  xr = 10 * sqrt(Sigma[d][d]);
  xl = -xr;
  /* THIS IS A SENSITIVE SPOT! */
  for(in = 0; in < ninit; in++)
  myxinit[d][in] = sqrt(Sigma[d][d]) *
  (in +.5 - ninit*.5);
  
  /* sample current r.eff element using ARS: */
  err = arms(myxinit[d], ninit, &xl, &xr, logdens, &para,
  &convex, npoint, dometrop, &xprev, &xsamp, nsamp,
  &qcent, &xcent, ncent, &neval);
  
  if (err>0) {
  error("ARMS error code = %d\n",err);
  }
  
  /** Rprintf("emstep = %d, cluster = %d, iteration = %d, b = %.4f\n",
  emstep, i, g, xsamp); **/
  
  for (l=clust_start[i]; l<clust_start[i+1]; l++){ 
  alpha[l] += (xsamp*ww[d][l] - b[d]*ww[d][l]);
  /* if(i==1&g==1) Rprintf("\nalpha[l]: %d", alpha[l]); */
  }
  b[d] = xsamp;
  sumbb[d][i] += b[d];
  sumb[d][i] += b[d]*b[d];
  sum[d] += b[d]*b[d];
  }  /* for d */
  
   	  for (d=1; d<=nreff; d++)
  for (dd=1; dd<=nreff; dd++)  sumv[d][dd][i]+=b[d]*b[dd];
  
  for (l=clust_start[i]; l<clust_start[i+1]; l++)
  oomega[l] += exp(alpha[l]);
  }  /* for g */

  	for (d=1; d<=nreff; d++)
  B[d][i] = b[d]; /* B keeps current set of r.effs */

  	for (d=1; d<=nreff; d++) {
  eb[d][i] = sumbb[d][i]/Gbs;
  
  v[d][i] = sumb[d][i]/Gbs-(sumbb[d][i]/Gbs)*(sumbb[d][i]/Gbs);
  for (dd=1; dd<=nreff; dd++) {
  condvar[d][dd][i] = 
  sumv[d][dd][i]/Gbs-(sumbb[d][i]/Gbs)*(sumbb[dd][i]/Gbs);
  }
  
  	}
  }  /* for i */
  	}   /* register block */

   for (i=1; i<=nobs; i++)   omega[rank[i]] = oomega[i]/Gbs;
   for (d=1; d<=nreff; d++) {
	   temp = 0;
	   for (i=1; i<=nclust; i++)   temp += sumv[d][d][i];
	   Sigma[d][d] = temp/(nclust*Gbs);
	   for (dd=1; dd<=nreff; dd++)   
  if (d!=dd)   Sigma[d][dd]=0;
   }

  
   /* Invert Sigma */

   for (d=1; d<=nreff; d++)
	   for (dd=1; dd<=nreff; dd++) vtemp[d][dd]=Sigma[d][dd];
   detSigma = Invmatrix( vtemp, nreff, nreff, invSigma );
   
/*   Rprintf("\nGbs: %d", Gbs);
   Rprintf("\nnclust: %d", nclust);
   for(i=1; i<=20; i++) Rprintf("\nww[1][i]: %d", ww[1][i]);
   Rprintf("\nb[1]: %d", b[1]);
   Rprintf("\nxsamp: %d", xsamp);
   Rprintf("\nalpha[2]: %d", alpha[2]);
   Rprintf("\noomega[2]: %d", oomega[2]);
   Rprintf("\nomega[rank[2]]: %d", omega[rank[2]]);
   Rprintf("\nrank: %d", rank[1]);
   Rprintf("\nsum[1]: %d", sum[1]);
   Rprintf("\nSigma[1][1]: %.4f", Sigma[1][1]);   */

}   /* end Estep() */

void Likelihood(int ncov, int nreff, double *Sigma[nreff+1],
  int nobs, double *omega,
  float *betahat, double *sum0, double *sum1[ncov+1], double **z, int *delta, int *ddelta,
  double *sum2[ncov+1][ncov+1], double *lambda, double *Lambda, double *Lambexp,
  double **steps, int NINIT, int emstep, int MAXSTEP, int nclust,
  double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, double **zz, int Gbs, int Gbsvar, double **a, int *rank, int CONVERG, double **bhat,
  double **sdbhat, double *sumv[nreff+1][nreff+1], double *invSigma[nreff+1], 
  double *eb[nreff+1], double *v[nreff+1], double *condvar[nreff+1][nreff+1], 
  double detSigma, int verbose, int varcov, double *pllaplace, double *plimport, 
  double *plbridge, double *laplacetot, double *bridgeC, double *invcondv[nreff+1][nreff+1], 
  double *detcondv, double *laplace)
{
	double incrB, logratio, uu[3], temp;
	int i, g, d, dd;
	long idum[1]={-10};
	double lbridge = *plbridge;
	double llaplace = *pllaplace;
	double limport = *plimport;
	
	llaplace = Laplace(nreff, ncov, nclust, condvar, detcondv, 
  	  invcondv, laplace, eb, Sigma,
  	  clust_start, invSigma, rank, ww, zz, ddelta, 
  	  betahat, z, lambda, Lambexp);
	
	*laplacetot = llaplace;
	
	limport = Importance(ncov, nreff, Sigma,
   nobs, omega,
   betahat, sum0, sum1, z, delta, ddelta,
   sum2, lambda, Lambda, Lambexp,
   steps, NINIT, emstep, MAXSTEP, nclust,
   sumb, sumbb, B, clust_start,
   ww, zz, Gbs, Gbsvar, a, rank, CONVERG, bhat,
   sdbhat, sumv, invSigma, 
   eb, v, condvar, 
   detSigma, verbose, varcov, 
   llaplace, bridgeC, invcondv, 
   detcondv, laplace);
   
	if(nreff <= 2){   
	
	lbridge = *laplacetot - *bridgeC;
	
	for (i=1; i<=nclust; i++) {
  incrB = 0;
  for (g=1; g<=Gbsvar; g++) {	
  	if (nreff==1) {
  uu[1] = gasdev(idum)*sqrt(v[1][i]) + eb[1][i]; 
  logratio = -(uu[1]-eb[1][i])*(uu[1]-eb[1][i])/2/v[1][i] 
  -log(2*PI*v[1][i])/2 -logJointprob(uu, i, ncov, nreff, Sigma,
  clust_start, invSigma, rank, ww, zz,
  ddelta, betahat, z, lambda, Lambexp);   /*nreff=1 only*/
  	} /*end if*/
  	if (nreff==2) {
  uu[1] = gasdev(idum)*sqrt(condvar[1][1][i]) + eb[1][i]; 
  uu[2] = gasdev(idum)*sqrt(condvar[2][2][i] - 
  condvar[1][2][i]*condvar[1][2][i]/condvar[1][1][i]) + eb[2][i] 
  +(uu[1]-eb[1][i])*condvar[1][2][i]/condvar[1][1][i];
  temp = 0;
  for (d=1; d<=nreff; d++) {
  temp += (uu[d]-eb[d][i])*(uu[d]-eb[d][i])*invcondv[d][d][i];
  for (dd=1; dd<d; dd++)   
  temp += 
  2*(uu[d]-eb[d][i])*(uu[dd]-eb[dd][i])*invcondv[d][dd][i];
  }
  logratio = -temp/2 -nreff*log(2*PI)/2 - log(detcondv[i])/2 
  -logJointprob(uu, i, ncov, nreff, Sigma,
  clust_start, invSigma, rank, ww, zz,
  ddelta, betahat, z, lambda, Lambexp);	
  	} /*end if*/
  	incrB += 1/( 1+exp( laplace[i] + logratio ) );
  }  /*end g*/
  lbridge += log(incrB/Gbsvar);
	}  /*end i*/
   
   *plbridge = lbridge;
   } /* end if nreff <= 2 */
   
   *pllaplace = llaplace;
   *plimport = limport;
}

double logJointprob(double *bb, int i, int ncov, int nreff, double *Sigma[nreff+1],
  int *clust_start, double *invSigma[nreff+1], int *rank, double **ww, double **zz,
  int *ddelta, float *betahat, double **z, double *lambda, double *Lambexp)
/* called in Laplace() and Importance() */
{
	double loglik, temp, temp1;
	int l, d, dd;
	
	loglik = 0;
	for (l=clust_start[i]; l<clust_start[i+1]; l++) {
  temp1 = temp = 0;
  for (d=1; d<=nreff; d++)   temp += bb[d]*ww[d][l];
  if (ddelta[l]) {	
  	for (d=1; d<=ncov; d++)   temp1 += betahat[d]*zz[d][l];
  	loglik += log(lambda[rank[l]])+ temp1+ temp;
  }
  loglik -= Lambexp[rank[l]]*exp(temp);
	} /*end l*/
   temp = 0; temp1=1;
   for (d=1; d<=nreff; d++) {
	   temp += bb[d]*bb[d]*invSigma[d][d] /*Sigma[d][d]*/;
	   temp1 *= Sigma[d][d];
	   for (dd=1; dd<d; dd++)   temp += 2*bb[d]*bb[dd]*invSigma[d][dd];
   }
   loglik -= nreff*log(2*PI)/2 +log(/*detSigma*/temp1)/2 + temp/2;
   return loglik;
}



double Laplace(int nreff, int ncov, int nclust, double *condvar[nreff+1][nreff+1], double *detcondv, 
  double *invcondv[nreff+1][nreff+1], double *laplace, double *eb[nreff+1], double *Sigma[nreff+1], 
  int *clust_start, double *invSigma[nreff+1], int *rank, double **ww, double **zz, int *ddelta, 
  float *betahat, double **z, double *lambda, double *Lambexp)
{
	int i, d, dd;
	double *vtemp[nreff+1], *vtemp1[nreff+1], temp1[nreff+1];
	
	for (d=1;d<=nreff;d++) {
  vtemp[d]= (double *)R_alloc((nreff+1),sizeof(double));
  vtemp1[d]= (double *)R_alloc((nreff+1),sizeof(double));
	}
	
	double laplacetot;
	laplacetot = 0;
	for (i=1; i<=nclust; i++) {
  for (d=1; d<=nreff; d++)
  	for (dd=1; dd<=nreff; dd++)  vtemp[d][dd]=condvar[d][dd][i]; 
  detcondv[i] = Invmatrix( vtemp, nreff, nreff, vtemp1 );
  for (d=1; d<=nreff; d++)
  	for (dd=1; dd<=nreff; dd++)  invcondv[d][dd][i]=vtemp1[d][dd];
  for (d=1; d<=nreff; d++)   temp1[d] = eb[d][i];
  laplace[i] = nreff*log(2*PI)/2 +log(detcondv[i])/2 +logJointprob(temp1, 
     i, ncov, nreff, Sigma,
     clust_start, invSigma, rank, ww, zz,
     ddelta, betahat, z, lambda, Lambexp);
  /* see twin/hyptest?? */
  laplacetot += laplace[i];
	} /*end i*/	
   return laplacetot;
}


void Nulllik(int nreff, int ncov, int nclust, double *condvar[nreff+1][nreff+1], double *detcondv, 
  	 double *invcondv[nreff+1][nreff+1], double *laplace, double *eb[nreff+1], double *Sigma[nreff+1], 
  	 int *clust_start, double *invSigma[nreff+1], int *rank, double **ww, double **zz, int *ddelta, 
  	 float *betahat, double **z, double *lambda, double *Lambexp, double loglik0)
/* log partial lik - log full lik = #events if no ties */
/* called in EM() step 0 */
{
	double temp1;
	int i, d, l;
	
	loglik0 = 0;
	for (i=1; i<=nclust; i++) 
  for (l=clust_start[i]; l<clust_start[i+1]; l++) {
  	temp1 = 0;
  	for (d=1; d<=ncov; d++)   temp1 += betahat[d]*zz[d][l];
  	if (ddelta[l])   loglik0 += log(lambda[rank[l]])+ temp1;
  	loglik0 -= Lambexp[rank[l]];
  } /*end l*/
}

double Importance(int ncov, int nreff, double *Sigma[nreff+1],
  int nobs, double *omega,
  float *betahat, double *sum0, double *sum1[ncov+1], double **z, int *delta, int *ddelta,
  double *sum2[ncov+1][ncov+1], double *lambda, double *Lambda, double *Lambexp,
  double **steps, int NINIT, int emstep, int MAXSTEP, int nclust,
  double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, double **zz, int Gbs, int Gbsvar, double **a, int *rank, int CONVERG, double **bhat,
  double **sdbhat, double *sumv[nreff+1][nreff+1], double *invSigma[nreff+1], 
  double *eb[nreff+1], double *v[nreff+1], double *condvar[nreff+1][nreff+1], 
  double detSigma, int verbose, int varcov, 
  double import, double *pbridgeC, double *invcondv[nreff+1][nreff+1], 
  double *detcondv, double *laplace)
{
	int g, i, d, dd, l, in;
	double incrA, logratio, incrC, temp;
	int err, ninit = 4, dometrop = 0,  npoint = 100, ncent = 0, neval, nsamp=1;
	double xl, xr, xprev = 0.0, xsamp, xcent, qcent, convex = 1.;
	double b[nreff+1], alpha[nobs+1];
	struct dens_para para;
	double myxinit[nreff+1][NINIT];
	
	double bridgeC = *pbridgeC;
	
	/*Gbs = MSAMP;*/
	bridgeC = 0;
	
	for (i=1; i<=nclust; i++) {
  
  incrA = incrC = 0;
  for (d=1; d<=nreff; d++) 
	  b[d] = B[d][i];
  
  for (l=clust_start[i]; l<clust_start[i+1]; l++) {
  	alpha[l] = 0;
  	for (d=1; d<=nreff; d++)   
  alpha[l] += b[d]*ww[d][l];
  }
  
  /* Gibbs sampler: */
  for (g=1; g<=Gbsvar; g++) {	/* MCMC loop	*/
  	
   for (d=1; d<=nreff; d++) {	
  /* initialize parameters for arms(): */
  /* initialize parameters for arms(): */
  para.i = i;
  para.d = d;
  para.nreff = nreff;
  para.alpha = alpha;
  para.b = b;
  para.clust_start = clust_start;
  para.Lambexp = Lambexp;
  para.ww = ww;
  para.a = a;
  para.Sigma = Sigma;
  para.rank = rank;
  para.invSigma = invSigma;
  
  xr = 10 * sqrt(Sigma[d][d]);
  xl = -xr;
  /* THIS IS A SENSITIVE SPOT! */
  for(in = 0; in < ninit; in++) 
  myxinit[d][in] = sqrt(Sigma[d][d]) * 
  (in +.5 - ninit*.5);
  
  /* sample current r.eff element using ARS: */
  err = arms(myxinit[d], ninit, &xl, &xr, logdens, &para,
  &convex, npoint, dometrop, &xprev, &xsamp, nsamp,
  &qcent, &xcent, ncent, &neval);
  
  if (err>0) {
  error("ARMS error code = %d\n",err);
  }
  
  for (l=clust_start[i]; l<clust_start[i+1]; l++) 
   alpha[l] += (xsamp*ww[d][l] - b[d]*ww[d][l]);
  B[d][i] = b[d] = xsamp;
  /*sumbb[d][i] += b[d];
  sumb[d][i] += b[d]*b[d];
  sum[d] += b[d]*b[d];*/
   }  /* for d */

   temp = 0;
   for (d=1; d<=nreff; d++) {
  	temp += (b[d]-eb[d][i])*(b[d]-eb[d][i])*invcondv[d][d][i];
  	for (dd=1; dd<d; dd++)   
  temp +=   
  2*(b[d]-eb[d][i])*(b[dd]-eb[dd][i])*invcondv[d][dd][i];
  }
  logratio = -temp/2 -nreff*log(2*PI)/2 - log(detcondv[i])/2 
  -logJointprob(b, i, ncov, nreff, Sigma,
   clust_start, invSigma, rank, ww, zz,
   ddelta, betahat, z, lambda, Lambexp);	
  temp = exp( laplace[i] + logratio );
/*Rprintf("%f ", temp);*/
  incrA += temp;
  incrC += 1/( 1+1/temp );
/*Rprintf("%f ", incrA);*/
  }	/* for g */

	import -= log( incrA/Gbsvar );
   	bridgeC += log(incrC/Gbsvar);
	}  /* for i */
   
   *pbridgeC = bridgeC;
  
   return import;
}

double Invmatrix ( double **a, int N, int M, double **y )
/* N is size of the matrix, M is # of col's from the inverse matrix */
{
	double d, col[N+1];
	int i, j, indx[N+1];
  
	myludcmp(a, N, indx, &d); 
	for (j=1; j<=N; j++)   d*=a[j][j];
	
  for (j=1; j<=M; j++) {
    for (i=1; i<=N; i++) col[i]=0.0;
    col[j]=1.0;
    mylubksb(a, N, indx, col);
  
    for (i=1; i<=N; i++) y[i][j]=col[i];
  }
	
	return d;
}

double **dmatrix2(double *array, int ncol, int nrow)
{
  register int i;
  register double **pointer;
	
  pointer = (double **)R_alloc(nrow, sizeof(double *));
  for (i=0; i<nrow; i++) {
    pointer[i] = array;
    array += ncol;
  }
  return(pointer);
}