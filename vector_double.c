/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  vector_double.c : Vector operations for data type double
 *                                                                     
 *  Copyright (C) 1995-2000 by Mark Austin, Xiaoguang Chen, and Wane-Jang Lin
 *  Institute for Systems Research,                                           
 *  University of Maryland, College Park, MD 20742                                   
 *                                                                     
 *  This software is provided "as is" without express or implied warranty.
 *  Permission is granted to use this software on any computer system
 *  and to redistribute it freely, subject to the following restrictions:
 * 
 *  1. The authors are not responsible for the consequences of use of
 *     this software, even if they arise from defects in the software.
 *  2. The origin of this software must not be misrepresented, either
 *     by explicit claim or by omission.
 *  3. Altered versions must be plainly marked as such and must not
 *     be misrepresented as being the original software.
 *  4. This software may not be sold or included in commercial software
 *     products without a license. 
 *  5. This notice is to remain intact.
 *                                                                    
 *  Written by: Mark Austin                                  July 1993 - May 1997
 *  ============================================================================= 
 */

/* #define DEBUG */

#include <stdio.h>
#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"

/*
 *  =================================================================
 *  Print Vector [length], where length = no items in Vector
 *  =================================================================
 */

#ifdef __STDC__
void PrintVectorDouble(VECTOR *v)
#else
void PrintVectorDouble(v)
VECTOR *v;
#endif
{
int i;


    if(v->cpVectorName != NULL) 
       printf ("\nVECTOR \"%s\" \n\n", v->cpVectorName);
    else 
       printf("\nVECTOR : \"UNTITLED\" \n\n");

    for(i = 1; i <= v->iLength; i++) {
       printf(" %3d          ",i);
       printf(" %16.9e\n", v->uVector.da[i-1]);
    }
    printf("\n");
}


/*
 *  =================================================
 *  Allocate and Free Vector data structure
 *  =================================================
 */

#ifdef __STDC__
double *dVectorAlloc(int length)
#else
double *dVectorAlloc( length)
int length;
#endif
{
double *array;

      array = (double *) MyCalloc( length, sizeof(double));

      return (array);
}

#ifdef __STDC__
void VectorFreeDouble(VECTOR *v)
#else
void VectorFreeDouble(v)
VECTOR *v;
#endif
{
   free ((char *) v->uVector.da);
   free ((char *) v->cpVectorName);
   free ((char *) v);
   v = (VECTOR *)NULL;
}


/*
 *  =================
 *  Vector Operations
 *  =================
 */

#ifdef __STDC__
VECTOR *VectorAddDouble(VECTOR *v1, VECTOR *v2)
#else
VECTOR *VectorAddDouble(v1, v2)
VECTOR *v1;
VECTOR *v2;
#endif
{
VECTOR *v3;
int i;

    /* [a] : Check Dimensions of Vector */

       if(v1->iLength != v2->iLength) {
          printf("FATAL ERROR >> Execution halted in VectorAdd()\n");
          printf("FATAL ERROR >> Problem : v1->iLength = %4d\n", v1->iLength);
          printf("FATAL ERROR >> Problem : v2->iLength = %4d\n", v2->iLength);
          FatalError("Inconsistent Dimensions",(char *)NULL);
       }

    /* [b] : Add Matrices */

       v3 = VectorAlloc((char *) NULL, DOUBLE_ARRAY, v2->iLength);
       for(i = 1; i <= v2->iLength; i++)
           v3->uVector.da[i-1] = v1->uVector.da[i-1] + v2->uVector.da[i-1];

       return(v3);
}

#ifdef __STDC__
VECTOR *VectorSubDouble(VECTOR *v1, VECTOR *v2)
#else
VECTOR *VectorSubDouble(v1, v2)
VECTOR *v1;
VECTOR *v2;
#endif
{
VECTOR *v3;
int i;

    /* [a] : Check Dimensions of Vector */

       if(v1->iLength != v2->iLength) {
          printf("FATAL ERROR >> Execution halted in VectorSub()\n");
          printf("FATAL ERROR >> Problem : v1->iLength = %4d\n", v1->iLength);
          printf("FATAL ERROR >> Problem : v2->iLength = %4d\n", v2->iLength);
          FatalError("Inconsistent Dimensions",(char *)NULL);
       }

    /* [b] : Add Matrices */

       v3 = VectorAlloc((char *) NULL, DOUBLE_ARRAY, v2->iLength);
       for(i = 1; i <= v2->iLength; i++)
           v3->uVector.da[i-1] = v1->uVector.da[i-1] - v2->uVector.da[i-1];

       return(v3);
}
