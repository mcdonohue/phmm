/* header file for phmm.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "arms.h"

/*#define ALLOC(a,b) S_alloc(a,b)*/

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
  double *Lambda);

void phmm2(int *sii,
  int *sNINIT, 
  int *sGbs,
  int *snreff,
  int *snobs,
  int *snclust,
  double *salpha,
  double *sb,
  double *sLambexp,
  double *sclust_start,
  double *wwv,
  double *av,
  double *sSigma,
  double *srank,
  double *sinvSigma,
  double *ssum,
  double *ssumb,
  double *ssumbb,
  double *ssumv,
  double *oomega);

struct dens_para {  /* used in logdens() */
  int i, d, *rank, *clust_start, nreff;

  double *alpha, *b, *Lambexp, **ww, **a, **Sigma, **invSigma;

};

void mnewt(int n, int ntrial, float *x, float tolx, float tolf,
  int nobs, double *omega, double *sum0, double **sum1, double **z,
  int *delta, double *sum2[n+1][n+1]);

void Sort(int nobs, double *x, int *delta, int *cluster,
  int ncov, double **z, int nreff, double **w );  /* for sorting data by x=time */
void Clust(int nobs, int *rank, double *xx, double *x, int *ddelta, int *delta,
  int ncov, double **zz, double **z, int nreff, double **ww,
  double **w, int *cluster, int *clust_sizes, int *nclust, int verbose);

  /* for sorting data by x=time */
void EM(int ncov, int nreff, double *Sigma[nreff+1],
  int nobs, double *omega,
  float *betahat, double *sum0, double *sum1[ncov+1], double **z, int *delta,
  double *sum2[ncov+1][ncov+1], double *lambda, double *Lambda, double *Lambexp,
  double **steps, int NINIT, int emstep, int MAXSTEP, double maxtime, int nclust,
  double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, int Gbs, int Gbsvar, double **a, int *rank, int CONVERG, double **bhat,
  double **sdbhat, double *sumv[nreff+1][nreff+1], double *invSigma[nreff+1], 
  double *eb[nreff+1], double *v[nreff+1], double *condvar[nreff+1][nreff+1], 
  double detSigma, int verbose, int varcov);

void Estep(int ncov, int nreff, double *Sigma[nreff+1], int nobs, int NINIT, int nclust,
  double **sumb, double **sumbb, double *sumv[nreff+1][nreff+1], double **B, double **ww,
  int Gbs, int *rank, int CONVERG, int *clust_start, double *omega, double *Lambexp,
  double **a, double *invSigma[nreff+1], double *eb[nreff+1],
  double *v[nreff+1], double *condvar[nreff+1][nreff+1],
  double detSigma, int varcov);

void Betahat(float *betahat, int ncov, int nobs, double *omega, double *sum0,
  double *sum1[ncov+1], double **z, int *delta, double *sum2[ncov+1][ncov+1],
  double *lambda, double *Lambda, double *Lambexp);

void Var(int Gbs, int Nobs, int nreff, int ncov, int NINIT, float *betahat,
  int nclust, double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, double *Lambexp, double *lambda, double **a, 
  double *Sigma[nreff+1], double *invSigma[nreff+1], int *rank,
  double **zz, int *delta, int *ddelta, double *omega, double **z, double **var);

void sVar(int *Gbs, int *nobs, int *nreff, int *ncov, int *NINIT, float *betahat,
  int *nclust, double *sumbv, double *sumbbv, double *sB, int *clust_start,
  double *sww, double *Lambexp, double *lambda, double *sa, 
  double *sSigma, double *sinvSigma, int *rank,
  double *szz, int *delta, int *ddelta, double *omega, double *sz,
  double *svar);

double BetaZ( int, float *, int, double ** );  /* beta*z in the linear predictor */
void usrfun(float *beta, int n, float *score, float **fjac,
  int nobs, double *omega, double *sum0, double *sum1[n], double **z,
  int *delta, double *sum2[n+1][n+1]);

double logdens(double, void *); /* log-density for Gibbs in arms() */
double Invmatrix( double **a, int N, int M, double **y );

void mylubksb(double **, int n, int *indx, double *);

void myludcmp(double **, int n, int *indx, double *);

double **dmatrix2(double *array, int ncol, int nrow);

float **matrix(long nrl, long nrh, long ncl, long nch);

extern float ran2( long * );

extern float gasdev( long * );

extern float ran1( long * );

extern double qgaus( double (*func)(double,int), double, double, int);

double Laplace(int nreff, int ncov, int nclust, double *condvar[nreff+1][nreff+1], double *detcondv, 
  double *invcondv[nreff+1][nreff+1], double *laplace, double *eb[nreff+1], double *Sigma[nreff+1], 
  int *clust_start, double *invSigma[nreff+1], int *rank, double **ww, double **zz, int *ddelta, 
  float *betahat, double **z, double *lambda, double *Lambexp);

void Nulllik(int nreff, int ncov, int nclust, double *condvar[nreff+1][nreff+1], double *detcondv, 
  double *invcondv[nreff+1][nreff+1], double *laplace, double *eb[nreff+1], double *Sigma[nreff+1], 
  int *clust_start, double *invSigma[nreff+1], int *rank, double **ww, double **zz, int *ddelta, 
  float *betahat, double **z, double *lambda, double *Lambexp, double loglik0);

double Jointprob( double, int);

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
  double laplacetot, double *bridgeC, double *invcondv[nreff+1][nreff+1], 
  double *detcondv, double *laplace);

double logJointprob(double *bb, int i, int ncov, int nreff, double *Sigma[nreff+1],
  int *clust_start, double *invSigma[nreff+1], int *rank, double **ww, double **zz,
  int *ddelta, float *betahat, double **z, double *lambda, double *Lambexp);

void Likelihood(int ncov, int nreff, double *Sigma[nreff+1],
  int nobs, double *omega,
  float *betahat, double *sum0, double *sum1[ncov+1], double **z, int *delta, int *ddelta,
  double *sum2[ncov+1][ncov+1], double *lambda, double *Lambda, double *Lambexp,
  double **steps, int NINIT, int emstep, int MAXSTEP, int nclust,
  double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, double **zz, int Gbs, int Gbsvar, double **a, int *rank, int CONVERG, double **bhat,
  double **sdbhat, double *sumv[nreff+1][nreff+1], double *invSigma[nreff+1], 
  double *eb[nreff+1], double *v[nreff+1], double *condvar[nreff+1][nreff+1], 
  double detSigma, int verbose, int varcov, double *llaplace, double *limport, 
  double *lbridge, double *laplacetot, double *bridgeC, double *invcondv[nreff+1][nreff+1], 
  double *detcondv, double *laplace);