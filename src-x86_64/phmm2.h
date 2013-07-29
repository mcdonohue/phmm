/* header file for phmm.c */

/*#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "arms.h"

/*#define ALLOC(a,b) S_alloc(a,b)*/

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


