#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "kcorrect.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code which converts to redshift from time;
   uses the utitilities in ztransform.c */

/********************************************************************/
IDL_LONG idl_t2z(int      argc,
                  void *   argv[])
{
  IDL_LONG ngals;
	float *z,omega0,omegal0,*t;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	t=((float *)argv[i]); i++;
	omega0=*((float *)argv[i]); i++;
	omegal0=*((float *)argv[i]); i++;
	z=((float *)argv[i]); i++;
	ngals=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
  for(i=0;i<ngals;i++) 
    z[i]=t2z(t[i],omega0,omegal0);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

