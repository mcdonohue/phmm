#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "phmm.h"

static R_NativePrimitiveArgType phmm_t[]={REALSXP, REALSXP, REALSXP, INTSXP,
	INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
	REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
	REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
	REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

static const R_CMethodDef CEntries[] ={
    CDEF(phmm),
    {NULL,NULL,0}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_phmm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL,NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}