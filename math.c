/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  math.c : Math Functions
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
 *  Written by: Mark Austin and Wane-Jang Lin                            May 1997
 *  ============================================================================= 
 */

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stddef.h>
#include "defs.h"
#include "miscellaneous.h"
#include "units.h"

/*
#define DEBUG
*/

extern int errno;

double	errcheck();

QUANTITY *QuantityUnitsLess(q)
QUANTITY *q;
{
    if(CheckUnits() == OFF) {
       printf("***** You should set units on to use this function \n");
       printf(" **** Check QDimenLess() in input file \n");
       printf(" **** Fail to set units on in QuantityUnitsLess() ");
       exit(1);
    }
    ZeroUnits(q->dimen);
    return (q);
}

/* --------------------------------------------------------- */
/* Random() : Generate Uniformly Distributed Random Numbers. */
/* --------------------------------------------------------- */

static     int  flag = 0;
QUANTITY *Random()
{
static     long int  seed;
long int   a = 16807;
long int   m = 2147483647;
long int   n = 127773;
long int   r = 2836;
long int   tmp_seed;
QUANTITY   *q;

        if( flag == 0 )      seed = time(NULL);
	flag = 1;

	tmp_seed = a*( seed%n ) - r*( seed/n );
	if( tmp_seed >= 0 )
	    seed = tmp_seed;
	else
	    seed = tmp_seed + m;

	q = (QUANTITY *) MyMalloc(sizeof(QUANTITY));
	q->value = errcheck( ((double) seed)/m, "random");
	if(CheckUnits() == ON) {
	    q->dimen = (DIMENSIONS *) MyMalloc(sizeof(DIMENSIONS));
            ZeroUnits(q->dimen);
        }
	else
	    q->dimen = (DIMENSIONS *)NULL;

        return(q);
}


QUANTITY *Log(q)
QUANTITY *q;
{
        q->value = errcheck(log(q->value), "log");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}

QUANTITY *Log10(q)
QUANTITY *q;
{
	q->value = errcheck(log10(q->value), "log10");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}

QUANTITY *Exp(q)
QUANTITY *q;
{
        q->value = errcheck(exp(q->value), "exp");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}

QUANTITY *Sin(q)
QUANTITY *q;
{
        q->value = errcheck(sin(q->value), "sin");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}
QUANTITY *Cos(q)
QUANTITY *q;
{
        q->value = errcheck(cos(q->value), "cos");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}
QUANTITY *Tan(q)
QUANTITY *q;
{
        q->value = errcheck(tan(q->value), "tan");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}
QUANTITY *Atan(q)
QUANTITY *q;
{
        q->value = errcheck(atan(q->value), "atan");
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}
QUANTITY *Integer(q)
QUANTITY *q;
{
        q->value = (double) (long) q->value;
        if(CheckUnits() == ON) 
           ZeroUnits(q->dimen);
        return(q);
}
QUANTITY *Fabs(q)
QUANTITY *q;
{
     q->value = fabs(q->value);
     return(q);
}

QUANTITY *Sqrt(q)
QUANTITY *q;
{

#ifdef DEBUG
    printf(" enter Sqrt()\n"); 
#endif 

     q->value = errcheck(sqrt(q->value), "sqrt");
   
    switch(CheckUnits()) {
      case ON:
          UnitsPowerRep( q->dimen, q->dimen, 0.5, YES );
      break;
      case OFF:
      break;
    }
     
#ifdef DEBUG
    printf(" Leaving Sqrt()\n"); 
#endif 
     
     return(q);
}

QUANTITY *Pow(q, value)
QUANTITY     *q;
double    value;
{

#ifdef DEBUG
       printf(" Enter Pow()\n"); 
       printf("x = %lf \n", q1->value);
       printf("y = %lf \n", value);
#endif

     if(value > 0) 
         q->value = errcheck(pow(q->value, value), "exponentiation");
     else {
       if(value == 0) q->value = 1.0;
       else
         q->value = 1.0/errcheck(pow(q->value, value), "exponentiation");
     }

    switch(CheckUnits()) {
      case ON:
         UnitsPowerRep( q->dimen, q->dimen, value, YES );
      break;
      case OFF:
      break;
    }
     
#ifdef DEBUG
    printf(" Leaving Pow()\n"); 
#endif 
     return(q);
}

double errcheck(d, s)
double d;
char   *s;
{
#ifdef DEBUG 
   printf(" enter errcheck() \n");
   printf(" d = %lf \n", d);
#endif
   if(errno == EDOM)  {
      errno = 0;
      EXECUTION_ERROR(s, "argument out of domain");
   }
   else if (errno == ERANGE)   {
   	errno = 0;
	EXECUTION_ERROR(s, "result out of range");
   }

#ifdef DEBUG 
   printf(" d = %lf \n", d);
   printf(" leaving  errcheck() \n");
#endif
   return d;
}
