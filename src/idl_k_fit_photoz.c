#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to fit photometric redshifts to a set of
   maggies */

/********************************************************************/
IDL_LONG idl_k_fit_photoz
  (int      argc,
   void *   argv[])
{
	IDL_LONG nk,nv,nz,maxiter,*niter,ngalaxy,verbose,nprior;
	float *coeffs,*photoz,*rmatrix,*zvals,*maggies,*maggies_ivar,*chi2;
	float *lprior, *zprior, tolerance;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	photoz=(float *)argv[i]; i++;
	coeffs=(float *)argv[i]; i++;
	rmatrix=(float *)argv[i]; i++;
	nk=*((IDL_LONG *)argv[i]); i++;
	nv=*((IDL_LONG *)argv[i]); i++;
	lprior=(float *)argv[i]; i++;
	zprior=(float *)argv[i]; i++;
	nprior=*((IDL_LONG *)argv[i]); i++;
	zvals=(float *)argv[i]; i++;
	nz=*((IDL_LONG *)argv[i]); i++;
	maggies=(float *)argv[i]; i++;
	maggies_ivar=(float *)argv[i]; i++;
	ngalaxy=*((IDL_LONG *)argv[i]); i++;
	tolerance=*((float *)argv[i]); i++;
	maxiter=*((IDL_LONG *)argv[i]); i++;
	niter=(IDL_LONG *)argv[i]; i++;
	chi2=(float *)argv[i]; i++;
	verbose=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) k_fit_photoz(photoz,coeffs,rmatrix,nk,nv,lprior,zprior, 
                                 nprior,zvals,nz,maggies,maggies_ivar,ngalaxy,
                                 tolerance,maxiter,niter,chi2,verbose);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

