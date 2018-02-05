/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  vector.c : Functions for Vector Operations and Printing
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
 *  Written by: Mark Austin                                              May 1997
 *  ============================================================================= 
 */

#include <stdio.h>
#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"

/* #define DEBUG */

/*
 *  =================================================================
 *  Print a Vector [length]
 *  
 *  where --  length = number of items in vector.
 *  =================================================================
 */

#ifdef __STDC__
void VectorPrint(VECTOR *v)
#else
void VectorPrint(v)
VECTOR *v;
#endif
{

       switch((int) v->eType) {
              case INTEGER_ARRAY:
                   (void) PrintVectorInteger(v);
                   break;
              case DOUBLE_ARRAY:
                   (void) PrintVectorDouble(v);
                   break;
              default:
                   FatalError("In VectorPrint() : Undefined v->eType",(char *)NULL);
                   break;
       }
}


/*
 *  =================================================
 *  Allocate Vector data structure, and array[length]
 *  =================================================
 */

#ifdef __STDC__
VECTOR *VectorAlloc(char *name, DATA_TYPE type, int length)
#else
VECTOR *VectorAlloc(name, type, length)
char       *name;
DATA_TYPE   type;
int       length;
#endif
{
VECTOR *v;
int     i;

   /* Allocate matrix parent data structure, and name */

      v = (VECTOR *) MyMalloc(sizeof(VECTOR));
      if(name != NULL)
         v->cpVectorName = SaveString(name);
      else 
         v->cpVectorName = (char *) NULL;

   /* Set parameters and allocate memory for matrix array */

      v->eType       = type;
      v->iLength    = length;

      switch((int) v->eType) {
          case INTEGER_ARRAY:
               v->uVector.ia = iVectorAlloc(length);
               break;
          case DOUBLE_ARRAY:
               v->uVector.da = dVectorAlloc(length);
               break;
          default:
               break;
      }

      return (v);
}

#ifdef __STDC__
void VectorFree(VECTOR *v)
#else
void VectorFree(v)
VECTOR *v;
#endif
{
      switch((int) v->eType) {
          case INTEGER_ARRAY:
               VectorFreeInteger(v);
               break;
          case DOUBLE_ARRAY:
               VectorFreeDouble(v);
               break;
          default:
               FatalError("In VectorFree() : Undefined v->eType",(char *)NULL);
               break;
      }
}


/*
 *  =================
 *  Vector Operations
 *  =================
 */

#ifdef __STDC__
VECTOR *VectorAdd(VECTOR *v1, VECTOR *v2)
#else
VECTOR *VectorAdd(v1, v2)
VECTOR *v1;
VECTOR *v2;
#endif
{
int i,j;

    switch((int) v1->eType) {
        case INTEGER_ARRAY:
             v1 = VectorAddInteger(v1, v2);
             break;
        case DOUBLE_ARRAY:
             v1 = VectorAddDouble(v1, v2);
             break;
        default:
             FatalError("In VectorAdd() : Undefined v->eType",(char *)NULL);
             break;
    }

    return (v1);
}

#ifdef __STDC__
VECTOR *VectorSub(VECTOR *v1, VECTOR *v2)
#else
VECTOR *VectorSub(v1, v2)
VECTOR *v1;
VECTOR *v2;
#endif
{

    switch((int) v1->eType) {
        case INTEGER_ARRAY:
             v1 = VectorSubInteger(v1, v2);
             break;
        case DOUBLE_ARRAY:
             v1 = VectorSubDouble(v1, v2);
             break;
        default:
             FatalError("In VectorSub() : Undefined v->eType",(char *)NULL);
             break;
    }

    return (v1);
}
