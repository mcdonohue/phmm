/* reffvar.c */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "phmm.h" 

void sVar(int *sGbs, int *snobs, int *snreff, int *sncov, int *sNINIT, float *sbetahat,
  int *snclust, double *ssumb, double *ssumbb, double *Bv, int *sclust_start,
  double *wwv, double *sLambexp, double *slambda, double *av, 
  double *sSigma, double *sinvSigma, int *srank,
  double *szz, int *sdelta, int *sddelta, double *somega, double *sz,
  double *svar)
{
  int i, j, k, d, dd;
  
  /*variables from R*/
  int NINIT = *sNINIT;
  int Gbs = *sGbs;
  int nreff = *snreff;
  int nobs = *snobs;
  int ncov = *sncov;
  int nclust = *snclust;
  double *alpha;
  double *b;
  double *Lambexp;
  double *lambda;
  int *delta;
  int *ddelta;
  double *omega;
  int *clust_start;
  double *ww[ncov+1];
  double *zz[ncov+1];
  double *z[ncov+1];
  double *a[nreff+1];
  double *B[nreff+1];
  double *Sigma[nreff+1];
  int *rank;
  double **invSigma;
  double *sumb[nreff+1];
  double *sumbb[nreff+1];
  double myxinit[nreff+1][NINIT];
  float betahat[ncov+1];
  double **var;
  
  b = (double *)R_alloc((nreff+1), sizeof(double));
  Lambexp = (double *)R_alloc((nobs+1), sizeof(double));
  lambda = (double *)R_alloc((nobs+1), sizeof(double));
  delta = (int *)R_alloc((nobs+1), sizeof(double));
  ddelta = (int *)R_alloc((nobs+1), sizeof(double));
  omega = (double *)R_alloc((nobs+1), sizeof(double));
  clust_start = (int *)R_alloc((nclust+2), sizeof(int));
  rank = (int *)R_alloc((nobs+1), sizeof(int));
  invSigma= dmatrix2(sinvSigma, nreff+1, nreff+1);  
  var = dmatrix2(svar, ncov+nreff+nobs+1, ncov+nreff+nobs+1);
    
  for (i=1;i<=nreff;i++) {
    B[i]= (double *)R_alloc((nclust+1),sizeof(double));
    a[i]= (double *)R_alloc((nclust+1),sizeof(double));
    Sigma[i] = (double*)R_alloc((nreff+1), sizeof(double));
    ww[i] = (double*)R_alloc((nobs+1), sizeof(double));
    sumb[i]= (double *)R_alloc((nclust+1),sizeof(double));
    sumbb[i]= (double *)R_alloc((nclust+1),sizeof(double));
  }
    
  for (i=1;i<=ncov;i++) {
    z[i]= (double *)R_alloc((nobs+1), sizeof(double));
    zz[i]= (double *)R_alloc((nobs+1), sizeof(double));
  }
    
  /*get data from R*/
    
  // for (i=1; i<=nobs; i++){ Rprintf("alpha=%.4f ", alpha[i]);}
  
  // for (i=1; i<=nreff; i++){ Rprintf("b=%.4f ", b[i]);}
  
  k=0;
  for (i=1; i<=nobs; i++){
    Lambexp[i] = sLambexp[k];
    lambda[i] = slambda[k];
    delta[i] = sdelta[k];
    ddelta[i] = sddelta[k];
    omega[i] = somega[k];
    k++;
  }

  k=0;
  for (i=1; i<=nclust+1; i++) {
    clust_start[i] = sclust_start[k];
    k++;
  }

  k=0;
  for (i=1; i<=nreff; i++){
    for (j=1; j<=nobs; j++){
      ww[i][j] = wwv[k];
      k++;
    }
  }

  k=0;
  for (i=1; i<=ncov; i++){
    for (j=1; j<=nobs; j++){
      zz[i][j] = szz[k];
      z[i][j] = sz[k];
      k++;
    }
  }

  k=0;
  for (i=1; i<=nclust; i++){
    for (j=1; j<=nreff; j++){
      a[j][i] = av[k];
      B[j][i] = Bv[k];
      k++; 
    }
  }

  k=0;
  for (i=1; i<=nobs; i++){
    rank[i] = srank[k];
    k++;
  }

  k=0;
  for (i=1; i<=nclust; i++){
    for (j=1; j<=nreff; j++){
      sumb[j][i] = ssumb[k];
      sumbb[j][i] = ssumbb[k];
      k++; 
    }
  }

  Rprintf("\n\n betahat:\n");
  for (i=0; i<=ncov; i++)  Rprintf("%.4f ", sbetahat[i]);

  k=0;
  for (i=1; i<=ncov; i++){
    betahat[i] = sbetahat[k];
    k++;
  }

  k=0;  
  for (dd=1; dd<=nreff; dd++)
    for (d=1; d<=nreff; d++){  
	 	  Sigma[dd][d] = sSigma[k];
	 	  invSigma[dd][d] = sinvSigma[k];
      k++;
	  }

  Var(Gbs, nobs, nreff, ncov, NINIT, betahat,
    nclust, sumb, sumbb, B, clust_start,
    ww, Lambexp, lambda, a, 
    Sigma, invSigma, rank,
    zz, delta, ddelta, omega, z,
    var);
};

/* Fisher's Observed Information and its Inverse */
void Var(int Gbs, int nobs, int nreff, int ncov, int NINIT, float *betahat,
  int nclust, double **sumb, double **sumbb, double **B, int *clust_start,
  double **ww, double *Lambexp, double *lambda, double **a, 
  double *Sigma[nreff+1], double *invSigma[nreff+1], int *rank,
  double **zz, int *delta, int *ddelta, double *omega, double **z,
  double **var) 
{
  int g, i, d, j, l, ll, nfail, dd, fail[nobs+1], in;
  double ooomega[nobs+1], ebetaz[nobs+1], sum[nreff+1];
  double sbeta[ncov+1], ssigma2[nreff+1], slambda[nobs+1], ssum0, ssum1;
  double *Iobs[ncov+nreff+nobs+1];
  double l2beta[ncov+1][ncov+1], l2sigma2[nreff+1], l2lambda[nobs+1], l2blamb[ncov+1][nobs+1];
  int err, ninit = 4, dometrop = 0,  npoint = 100, ncent = 0, neval, nsamp=1;
  double xinit[NINIT], xl, xr, xprev = 0.0, xsamp, xcent, qcent, convex = 1., temp;
  double b[nreff+1], alpha[nobs+1], oomega[nobs+1];
  struct dens_para para;
  double myxinit[nreff+1][NINIT];
  double msbeta[ncov+1], mssigma2[nreff+1], mslambda[nobs+1];
  double varLamb[nobs+1], sumcov;
  
  for (i=1;i<=ncov+nreff+nobs;i++) { 
    Iobs[i]= (double *)R_alloc((ncov+nreff+nobs+1),sizeof(double));
  }

  for (j=1; j<=nobs; j++){
    ebetaz[j] = exp(BetaZ(j, betahat, ncov, z)); 
  }   

  varLamb[0] = 0;

  /* Initializations: */
  for (l=1; l<=nobs; l++){
    oomega[l] = 0; 
    alpha[l] = 0; 
  }

  for (d=1; d<=nreff; d++) 
    for (i=1; i<=nclust; i++)   sumb[d][i] = sumbb[d][i] = 0;   
  /* sumb for E(b^2) and sumbb for E(b) */

  for (d=1; d <= nreff; d++)
    mssigma2[d] = l2sigma2[d] = 0;
  for (d=1; d<=ncov; d++) {
    msbeta[d] = 0;
    for (dd=1; dd<=ncov; dd++)   l2beta[d][dd] = 0;
  }
  for (i=1; i<=nobs; i++)   mslambda[i] = 0;

  for (l=1; l<=ncov+nreff+nobs; l++)  
  for (j=1; j<=ncov+nreff+nobs; j++) 
    Iobs[l][j] = 0;

  // x, delta, z, w are sorted by time   
  // xx, ddelta, cluster, rank, zz, ww are sorted by clusters     
  
  // for (i=1; i<=nobs; i++)  Rprintf("\n%d: %d", i, ddelta[i]);
  // for (d=1; d<=nreff; d++){
  //   Rprintf("\n%d: ", d);    
  //   for (l=clust_start[i]; l<clust_start[i+1]; l++) {
  //     Rprintf(" %d", ww[d][l]);
  //   }
  // }
  
  // for (d=1; d<=nreff; d++){
  //   Rprintf("\n%d: ", d);
  //   for (dd=1; dd<=nreff; dd++)
  //   Rprintf(" %f", Sigma[d][dd]);
  // }
  // 
  // for (d=1; d<=nreff; d++){
  //   Rprintf("\n%d: ", d);
  //   for (dd=1; dd<=nreff; dd++)
  //   Rprintf(" %f", invSigma[d][dd]);
  // }
  
  // for (i=1; i<=nclust; i++) {
  //   Rprintf("\n%d: ", i);
  //   for (d=1; d<=nreff; d++)
  //   Rprintf(" %f", B[d][i]);    
  //   for (l=clust_start[i]; l<clust_start[i+1]; l++) {
  //     alpha[l] = 0;
  //     for (d=1; d<=nreff; d++)   
  //       alpha[l] += b[d]*ww[d][l];
  //     Rprintf("\n alpha=%f ", alpha[l]);
  //   }   
  // }

  Rprintf("\n\n alpha:\n");
  for (i=1; i<=nobs; i++){ Rprintf("alpha=%.4f ", alpha[i]);}
  Rprintf("\n\n b:\n");
  for (i=1; i<=nreff; i++){ Rprintf("b=%.4f ", b[i]);}
  Rprintf("\n\n Lambexp:\n");
  for (i=1; i<=nobs; i++){ Rprintf("Lambexp=%.4f ", Lambexp[i]);} 
  Rprintf("\n\n ww:\n");
  for (i=1; i<=nobs; i++)  Rprintf("w1=%.4f w2=%.4f \n", ww[1][i], ww[2][i]);
  Rprintf("\n\n omega:\n");
  for (i=1; i<=nobs; i++)  Rprintf("omega=%.4f \n", omega[i]);
  Rprintf("\n\n a:\n");
  for (i=1; i<=nclust; i++)  Rprintf("%.4f %.4f \n", a[1][i], a[2][i]);
  Rprintf("\n\n betahat:\n");
  for (i=1; i<=ncov; i++)  Rprintf("%.4f ", betahat[i]);

  /* Gibbs sampler: */
  for (g=1; g<=Gbs; g++) {  /* MCMC loop  */
    for (d=1; d<=ncov; d++)   sbeta[d] = 0;
    for (d=1; d<=nreff; d++)    sum[d] = 0;   
    /* sum[] for b^2 within each Gibbs iteration */

    for (i=1; i<=nclust; i++) {
      for (d=1; d<=nreff; d++)
      b[d] = B[d][i];    
      for (l=clust_start[i]; l<clust_start[i+1]; l++) {
        alpha[l] = 0;
        for (d=1; d<=nreff; d++)   
          alpha[l] += b[d]*ww[d][l];
      }

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

        // if((xinit[0] <= xl) || (xinit[ninit-1] >= xr)){
        //   /* initial points do not satisfy bounds */
        //   return 1003;
        // }
        
        // Rprintf("\n%f %f %f %f", myxinit[d][0], xl, xr, myxinit[3]);

        /* sample current r.eff element using ARS: */
        err = arms(myxinit[d], ninit, &xl, &xr, logdens, &para, 
          &convex, npoint, dometrop, &xprev, &xsamp, nsamp,  
          &qcent, &xcent, ncent, &neval);
        // Rprintf(" here");

        if (err>0) {
          error("ARMS error code = %d in reffvar.c\n", err);
        }

        for (l=clust_start[i]; l<clust_start[i+1]; l++) 
          alpha[l] += (xsamp*ww[d][l] - b[d]*ww[d][l]);

        B[d][i] = b[d] = xsamp;
        sumbb[d][i] += b[d];
        sumb[d][i] += b[d]*b[d];
        sum[d] += b[d]*b[d];
      }  /* for d */

      for (l=clust_start[i]; l<clust_start[i+1]; l++) {
        temp = exp(alpha[l]);
        ooomega[rank[l]] = temp;   /* ooomega is current sampled exp(b*w), ordered by time */
        oomega[l] += temp;   /* oomega will be E[exp(b*w)]*Gbs */ 
        for (d=1; d<=ncov; d++) {
          sbeta[d] += zz[d][l]*(ddelta[l]-Lambexp[rank[l]]*temp);
        }
      }  /* for l */
    }  /* for i */

    /* compute ssigma2 */
    for (d=1; d<=nreff; d++) 
    ssigma2[d] = -0.5*(nclust-sum[d]/Sigma[d][d])/Sigma[d][d];

    /* compute slambda */
    ssum0 = l = 0;
    for (j=nobs; j>=1; j--) {
      ssum0 += ebetaz[j]*ooomega[j];
      if (delta[j]) {
        l++;
        slambda[l] = 1/lambda[j] - ssum0;
        fail[l]=j;  /* fail is the link between slambda and time */
      }
    }   /* for j */
    nfail = l;

    /* update Imis */
    for (l=1; l<=ncov; l++) {
      msbeta[l] += sbeta[l];
      for (d=1; d<=ncov; d++)   
        Iobs[l][d] += sbeta[l]*sbeta[d];
      for (d=1; d<=nreff; d++)   
        Iobs[l][ncov+d] += sbeta[l]*ssigma2[d];
      for (j=1; j<=nfail; j++)   
        Iobs[l][ncov+nreff+j] += sbeta[l]*slambda[j];
    }
    for (l=1; l<=nreff; l++) { 
      mssigma2[l] += ssigma2[l];
      for (d=1; d<=ncov; d++)   
      Iobs[ncov+l][d] += ssigma2[l]*sbeta[d];
      for (d=1; d<=nreff; d++)   
      Iobs[ncov+l][ncov+d] += ssigma2[l]*ssigma2[d];
      for (j=1; j<=nfail; j++)   
      Iobs[ncov+l][ncov+nreff+j] += ssigma2[l]*slambda[j];
    }
    for (l=1; l<=nfail; l++) { 
      mslambda[l] += slambda[l];
      for (d=1; d<=ncov; d++)   
      Iobs[ncov+nreff+l][d] += slambda[l]*sbeta[d];
      for (d=1; d<=nreff; d++)   
      Iobs[ncov+nreff+l][ncov+d] += slambda[l]*ssigma2[d];
      for (j=1; j<=nfail; j++)   
      Iobs[ncov+nreff+l][ncov+nreff+j] += slambda[l]*slambda[j];
    }
    /* compute l2sigma2 (another version) */
    /* for (d=1; d <= nreff; d++) {
    l2sigma2[d] += -0.5*(-nclust + 2.*sum[d]/sigma2[d])/
    (sigma2[d]*sigma2[d])/Gbs;
    } */
  }  /* for g */

  /* Print mean of the scores */
  /*
   Rprintf("\n mean of the scores:\n");
   for (d=1; d<=ncov; d++)   Rprintf("%.4f ", msbeta[d]/Gbs);
   Rprintf("\n");
   for (d=1; d <= nreff; d++)   Rprintf("%.4f ", mssigma2[d]/Gbs);
   Rprintf("\n");
   */
  /*
   for (i=1; i<=nfail; i++)   Rprintf("%.4f ", mslambda[i]/Gbs);
   Rprintf("\n");
   */

  /* print E(b) */
  /*
   Rprintf("\n E(b) and sd:\n");
   for (i=1; i<=nclust; i++) {
  for (d=1; d<=nreff; d++)   Rprintf("%f %f    ", sumbb[d][i]/Gbs, sqrt(sumb[d][i]/Gbs - (sumbb[d][i]/Gbs)*(sumbb[d][i]/Gbs)) );
  Rprintf("\n");
   }
   */


  /* Print Imis for beta, sigma2 */
  /*
   Rprintf("\n\n Imis:\n");
   for (l=1; l <= ncov+nreff; l++) {
  for( ll=1; ll <= ncov+nreff; ll++)
  Rprintf("%9.2f  ", Iobs[l][ll]/Gbs);
  Rprintf("\n");
   }
   Rprintf("\n");
   */

  for (j=1; j<=nobs; j++)  
  omega[rank[j]] = oomega[j]/Gbs;

  /* compute 2nd derivatives */

  for (j=1; j<=nobs; j++)
  for (d=1; d<=ncov; d++) 
  for (dd=1; dd<=ncov; dd++)   l2beta[d][dd] -= z[d][j]*z[dd][j]*Lambexp[j]*omega[j];

  for (d=1; d<=nreff; d++) {
    sum[d] = 0;
    for (i=1; i<=nclust; i++)   sum[d] += sumb[d][i];
  }   

  for (d=1; d<=nreff; d++)   l2sigma2[d] = -0.5*(-nclust+ 2*(sum[d]/Gbs)/Sigma[d][d]) / (Sigma[d][d]*Sigma[d][d]);

  for (l=1; l<=nfail; l++)   l2lambda[l] = -1/(lambda[fail[l]]*lambda[fail[l]]);

  for (d=1; d<=ncov; d++) {
    ssum1 = l = 0;
    for (j=nobs; j>=1; j--) {
   ssum1 += z[d][j]*ebetaz[j]*omega[j];  
   if (delta[j]) {
   l++;
   l2blamb[d][l] = - ssum1;
   } 
    }  /* for j */ 
  }   /* for d */


  /* Print Iaug */
  /*
   Rprintf("Iaug:\n");
   for (d=1; d <= ncov; d++)  {
  for (dd=1; dd <= ncov; dd++) 
  Rprintf("%7.2f ", -l2beta[d][dd]);
  Rprintf("\n");
   }
   for  (d=1; d <= nreff; d++)
   Rprintf("%7.2f ", -l2sigma2[d]);
   Rprintf("\n");
   */

  /* finish computing Iobs */
  for (l=1; l<=ncov+nreff+nobs; l++)  
  for (j=1; j<=ncov+nreff+nobs; j++)   Iobs[l][j] = -Iobs[l][j]/Gbs;

  for (l=1; l<=ncov; l++) {
    for (j=l; j<=ncov; j++)   Iobs[l][j] -= l2beta[l][j];
    for (j=1; j<=nfail; j++)   Iobs[l][ncov+nreff+j] -= l2blamb[l][j];
  }
  for (l=1; l<=nreff; l++)   Iobs[ncov+l][ncov+l] -= l2sigma2[l];
  for (l=1; l<=nfail; l++)   Iobs[ncov+nreff+l][ncov+nreff+l] -= l2lambda[l];

  for (l=1; l<=ncov+nreff+nobs; l++)  
  for (j=l+1; j<=ncov+nreff+nobs; j++)   Iobs[j][l] = Iobs[l][j];

  /* Print out Iobs */
  /*
   Rprintf("\nIobs:\n");
   for (i=1; i<=ncov+nreff+20; i++) {
  for (j=1; j<= ncov+nreff; j++)   
  Rprintf("%7.2f ", Iobs[i][j]);
  Rprintf("\n");
   }
   */


  /* get the first ncov+nreff col's of Iobs^(-1) */
  Invmatrix( Iobs, ncov+nreff+nfail, ncov+nreff, var );

  /*Rprintf("\n Var-cov matrix:\n");
  for (i=1; i<=ncov+nreff; i++) {
    for (j=1; j<=ncov+nreff; j++)   Rprintf("%.4f ", var[i][j]);
    Rprintf("\n");
  }*/

  /* Calculate var of Lambda() */
  /*
   for (i=1; i<=nfail; i++) {
  sumcov = 0;
  for (j=ncov+nreff+1; j<ncov+nreff+i; j++)
  sumcov += var[ncov+nreff+i][j];
  varLamb[i] = varLamb[i-1] + var[ncov+nreff+i][ncov+nreff+i] + 2*sumcov;
   } 
   */

  /*Rprintf("\n standard error:\n");
  for (d=1; d<=ncov; d++)   Rprintf("%.4f ", sqrt(var[d][d]));
  Rprintf("\n");
  for (d=ncov+1; d <= ncov+nreff; d++)   Rprintf("%.4f ", sqrt(var[d][d]));
  Rprintf("\n");*/
  /*
   for (j=1; j<=nfail; j++) 
   Rprintf("%f ", sqrt(varLamb[j]));
   Rprintf("\n");
   */
}