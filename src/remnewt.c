#include <R.h>
#include <Rinternals.h>
#define NRANSI
#include "phmm.h"

void mnewt(int n, int ntrial, float *x, float tolx, float tolf,
   int nobs, double *omega, double *sum0, double **sum1, double **z,
   int *delta, double *sum2[n+1][n+1])
{
    void lubksb(float **a, int n, int *indx, float b[]);
    void ludcmp(float **a, int n, int *indx, float *d);
    int k,i,*indx;
    float errx,errf,d,*fvec,**fjac,*p;

    indx=(int *)R_alloc(n, sizeof(int));
    p=(float *)R_alloc(n, sizeof(float));
    fvec=(float *)R_alloc(n, sizeof(float));
    fjac=matrix(1,n,1,n);
    errx=tolx+1;
    for (k=1;k<=ntrial;k++) {            /*Rprintf("trial no.1");*/
            usrfun(x,n,fvec,fjac,nobs,omega,sum0,sum1,z,delta,sum2);
        if (errx <= tolx) return;
        errf=0.0;
        for (i=1;i<=n;i++) errf += fabs(fvec[i]);
        if (errf <= tolf) return;
        for (i=1;i<=n;i++) p[i] = -fvec[i];
        ludcmp(fjac,n,indx,&d);
        lubksb(fjac,n,indx,p);
        errx=0.0;
        for (i=1;i<=n;i++) {
            errx += fabs(p[i]);
            x[i] += p[i];
        }
    }
    return;
}

/* BETTER MOVED IN THE SAME FILE WITH   mnewt() */
void usrfun(float *beta, int n, float *score, float **fjac,
   int nobs, double *omega, double *sum0, double *sum1[n], double **z,
   int *delta, double *sum2[n+1][n+1])
{           /* used by mnewt(); updated for frailties */

   int i, j, k;
   double temp;

/** Rprintf("usrfun "); **/
   for (j=1; j<=n; j++)   {
      score[j]=0;
      for (k=1; k<=n; k++)   fjac[j][k]=0;
   }

   for ( i=nobs; i>=1; i-- )  {
      temp=exp(BetaZ(i, beta, n, z))*omega[i];
      sum0[i] = sum0[i+1]+temp;
      for (j=1;j<=n;j++) {
         sum1[j][i] = sum1[j][i+1]+temp*z[j][i];
         if (delta[i])   score[j]+= (z[j][i]-sum1[j][i]/sum0[i]);
      }
      for (j=1;j<=n;j++)
         for (k=1;k<=n;k++)  {
            sum2[j][k][i] = sum2[j][k][i+1]+temp*z[j][i]*z[k][i];
            if (delta[i])
              fjac[j][k]+= ((sum1[j][i]/sum0[i])*(sum1[k][i]/sum0[i]) -
            sum2[j][k][i]/sum0[i]);
         }
   }
}   /* end usrfun() */

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **)R_alloc((nrow+1),sizeof(float*));
	if (!m) error("allocation failure 1 in matrix()");
	m += 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *)R_alloc((nrow*ncol+1),sizeof(float));
	if (!m[nrl]) error("allocation failure 2 in matrix()");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
