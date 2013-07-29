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
  double *soomega)
{
  int i, j, k, g, in, dd, l, d;
  int NINIT = *sNINIT;
  int err, ninit = 4, dometrop = 0,  npoint = 100, ncent = 0, neval, nsamp=1;
  double xinit[NINIT], xl = -100.0, xr = 100.0, xprev = 0.0, xsamp, xcent, qcent, convex = 1., temp;
  struct dens_para para;

  /*variables from R*/
  int ii = *sii;   
  int Gbs = *sGbs;
  int nreff = *snreff;
  int nobs = *snobs;
  int nclust = *snclust;
  double *alpha;
  double *b;
  double *Lambexp;
  int *clust_start;
  double *ww[nreff+1];
  double *a[nreff+1];
  double *Sigma[nreff+1];
  int *rank;
  double **invSigma;
  double *sum;
  double *sumb[nreff+1];
  double *sumbb[nreff+1];
  double *sumv[nreff+1][nreff+1];
  double *oomega;
  double myxinit[nreff+1][NINIT];

  alpha = (double *)R_alloc((nobs+1), sizeof(double));
  b = (double *)R_alloc((nreff+1), sizeof(double));
  Lambexp = (double *)R_alloc((nobs+1), sizeof(double));
  clust_start = (int *)R_alloc((nclust+2), sizeof(int));
  
  for (i=1; i<=nreff; i++){
    ww[i] = (double*)R_alloc((nobs+1), sizeof(double));
  }
  
  for (i=1;i<=nreff;i++) {
    a[i]= (double *)R_alloc((nclust+1),sizeof(double));
  }
  for (i=1; i<=nreff; i++){
    Sigma[i] = (double*)R_alloc((nreff+1), sizeof(double));
  }
  /* Sigma = dmatrix2(sSigma, nreff+1, nreff+1); */
  
  rank = (int *)R_alloc((nobs+1), sizeof(int));
  invSigma= dmatrix2(sinvSigma, nreff+1, nreff+1);
  sum = (double *)R_alloc((nreff+1), sizeof(double));

  for (i=1;i<=nreff;i++){
    sumb[i]= (double *)R_alloc((nclust+1),sizeof(double));
    sumbb[i]= (double *)R_alloc((nclust+1),sizeof(double));
  }
  
  for (i=1;i<=nreff;i++)  {
    for (j=1;j<=nreff;j++) {
      sumv[i][j]= (double *)R_alloc((nclust+1),sizeof(double));
    }
  }
    
  oomega = (double *)R_alloc((nobs+1), sizeof(double));

  /*get data from R*/
  
  k=0;
  for (i=1; i<=nobs; i++){
    alpha[i] = salpha[k];
    k++;
  }
  
  // for (i=1; i<=nobs; i++){ Rprintf("alpha=%.4f ", alpha[i]);}

  k=0;
  for (i=1; i<=nreff; i++){
    b[i] = sb[k];
    k++;
  }
  
  // for (i=1; i<=nreff; i++){ Rprintf("b=%.4f ", b[i]);}
  
  k=0;
  for (i=1; i<=nobs; i++){
    Lambexp[i] = sLambexp[k];
    k++;
  }
  // Rprintf("\n");
  // for (i=1; i<=nobs; i++){ Rprintf("Lambexp=%.4f ", Lambexp[i]);} 
  
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

  // for (i=1; i<=nobs; i++)  Rprintf("w1=%.4f w2=%.4f \n", ww[1][i], ww[2][i]);

  k=0;
  for (i=1; i<=nclust; i++){
    for (j=1; j<=nreff; j++){
      a[j][i] = av[k];
      k++; 
    }
  }
  
  // Rprintf("\n\n a:\n");
  // for (i=1; i<=nclust; i++)  Rprintf("%.4f %.4f \n", a[1][i], a[2][i]);
  
  k=0;
  for (i=1; i<=nobs; i++){
    rank[i] = srank[k];
    k++;
  }

  k=0;
  for (i=1; i<=nreff; i++){
    sum[i] = ssum[k];
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
    
  k=0;
  for (i=1; i<=nclust; i++)  
  for (dd=1; dd<=nreff; dd++)
    for (d=1; d<=nreff; d++){  
	 	  sumv[d][dd][i] = ssumv[k];
      k++;
		}
  
  k=0;
  for (i=1; i<=nobs; i++){
    oomega[i] = soomega[k];
    k++; 
  }
  
  k=0;  
  for (dd=1; dd<=nreff; dd++)
    for (d=1; d<=nreff; d++){  
	 	  Sigma[dd][d] = sSigma[k];
      k++;
	  }

  /* for (i=1; i<=nreff; i++)  Rprintf("%.4f %.4f\n", Sigma[1][i], Sigma[2][i]);*/

  for (g=1; g<=Gbs; g++) {  /* MCMC loop  */
		for (d=1; d<=nreff; d++) {
			/* initialize parameters for arms(): */
			para.i = ii;
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
				myxinit[d][in] = sqrt(Sigma[d][d]) * (in +.5 - ninit*.5);
			
      /* sample current r.eff element using ARS: */
      /*   Rprintf("cluster = %d, reff = %d, nreff = %d\n", ii, d, nreff);
      for (i=1; i<=nobs; i++){ Rprintf("alpha=%.4f ", alpha[i]);}
      Rprintf("\n");
      for (i=1; i<=nreff; i++){ Rprintf("%.4f ", b[i]);}
      Rprintf("\n");
      for (i=1; i<=nobs; i++){ Rprintf("Lambexp=%.4f ", Lambexp[i]);}
      Rprintf("\n");
      for (i=1; i<=nobs; i++)  Rprintf("%.1f %.1f\n", ww[1][i], ww[2][i]);
      Rprintf("\n");
      for (i=1; i<=nclust; i++)  Rprintf("%.1f %.1f\n", a[1][i], a[2][i]);
      Rprintf("Sigma\n");
      for (i=1; i<=nreff; i++)  Rprintf("%.1f %.1f\n", Sigma[1][i], Sigma[2][i]);
      Rprintf("invSigma\n");
      for (i=1; i<=nreff; i++)  Rprintf("%.1f %.1f\n", invSigma[1][i], invSigma[2][i]);
      Rprintf("\n");
      for (i=1; i<=nobs; i++){ Rprintf("%.4f ", rank[i]);}
      Rprintf("\n");*/

			err = arms(myxinit[d], ninit, &xl, &xr, logdens, &para,
					   &convex, npoint, dometrop, &xprev, &xsamp, nsamp,
					   &qcent, &xcent, ncent, &neval);
			
			if (err>0) {
				error("ARMS error code = %d\n",err);
			}

			for (l=clust_start[ii]; l<clust_start[ii+1]; l++){ 
				alpha[l] += (xsamp*ww[d][l] - b[d]*ww[d][l]);
				/* if(i==1&g==1) Rprintf("\nalpha[l]: %d", alpha[l]); */
			}
			b[d] = xsamp;
			sumbb[d][ii] += b[d];
			sumb[d][ii] += b[d]*b[d];
			sum[d] += b[d]*b[d];
	  }  /* for d */

    for (d=1; d<=nreff; d++)
    for (dd=1; dd<=nreff; dd++)  sumv[d][dd][ii]+=b[d]*b[dd];

    for (l=clust_start[ii]; l<clust_start[ii+1]; l++)
      oomega[l] += exp(alpha[l]);
  }  /* for g */
  
  /*return the results to R*/
  
  for (i=1; i<=nobs; i++){
    salpha[i-1] = alpha[i];
  }
  
  for (i=1; i<=nreff; i++){
    sb[i-1] = b[i];
  }
  
  for (i=1; i<=nobs; i++){
    sLambexp[i-1] = Lambexp[i];
  }
  
  for (i=1; i<=nreff; i++){
    ssum[i-1] = sum[i];
  }
  k=0;
  for (i=1; i<=nclust; i++){
    for (j=1; j<=nreff; j++){
      ssumb[k] = sumb[j][i];
      k++; 
    }
  }
  k=0;
  for (i=1; i<=nclust; i++){
    for (j=1; j<=nreff; j++){
      ssumbb[k] = sumbb[j][i];
      k++; 
    }
  }
  
  k=0;
  for (i=1; i<=nclust; i++)   
  for (dd=1; dd<=nreff; dd++)
    for (d=1; d<=nreff; d++){ 
			ssumv[k] = sumv[d][dd][i];
      k++;
		}
  
  for (i=1; i<=nobs; i++){
    soomega[i-1] = oomega[i];
  }
  // Rprintf("\n\n oomega:\n");
  // for (l=1; l <= nobs; l++){
  //    Rprintf("%.4f  ", oomega[l]);
  //    Rprintf("\n");
  // }
}