/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  code.c : Functions for Stack Machine
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
 *  Written by: Mark Austin, Xiaoguang Chen, and Wane-Jang Lin         March 2000
 *  ============================================================================= 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "code.h"
#include "vector.h"
#include "fe_functions.h"
#include "y.tab.h"

/* Setup Program Stack */

#define  NSTACK   1500
#define  NPROG    5000
#define  NFRAME      5

/* Setup Stack Frame for Program Machine */

static DATUM  stack[NSTACK];
static DATUM        *stackp;

Inst     prog[NPROG];
Inst     *progp;
Inst     *pc;
Inst     *progbase = prog;
int      returning;

static int iBREAK_FLAG = 0;

typedef struct Frame {
        SYMBOL     *sp;
        Inst    *retpc;
        DATUM    *argn;
        int      nargs;
} Frame;

Frame frame_in_code[NFRAME];
Frame *fp;


/* Initialise the pointers stackp and progp to top of their stacks. */

int Init_Code()
{
   progp     = progbase;
   stackp    = stack;
   fp        = frame_in_code;
   returning = 0;
}

/* Push a data type Datum onto the interpreter and
   increment stackp to point to the address to which will
   be assigned the next element to be pushed onto stack.*/

#ifdef __STDC__
int Push( DATUM d )
#else
int Push(d)
DATUM d;
#endif
{
    if(stackp >= &stack[NSTACK])
       EXECUTION_ERROR("ERROR >> stack overflow", (char *) 0);
       *stackp++ = d;  /* *stackp=d; stackp++; */
}

DATUM Pop()
{
    if(stackp == stack)
       EXECUTION_ERROR("ERROR >> stack underflow ", (char *) 0);
       return *--stackp;  /* stackp--; return *stackp */
}

int Pop_Eval()
{
   Pop();
}

#ifdef __STDC__
int Execute( Inst *p )
#else
int Execute(p)
Inst *p;
#endif
{
    for(pc = p; *pc != STOP && !returning; ) {
        pc = pc+1;
        if(Check_Break())
           break;
        (*(*(pc-1)))();
    }
}

/* Push_constant() pushes a data type Datum onto the 
   interpreter stack, and increments the pointer to prog.*/

int Push_Constant()
{
DATUM d;

    d.q        = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
    if(CheckUnits() == ON)
        d.q->dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
    else
        d.q->dimen = (DIMENSIONS *)NULL;

    d.q->value = ((SYMBOL *) *pc)->u.value;
    pc = pc+1;
    Push(d);
}

int Push_Variable()
{
DATUM d;

    d.sym = (SYMBOL *)(*pc);
    pc = pc+1;
    Push(d);
}

int Push_String()
{
DATUM d;

    d.sym = (SYMBOL *) SaveString((char *)(*pc));
    pc = pc+1;
    Push(d);
}

int Push_Matrix()
{
DATUM d;

    d.sym = (SYMBOL *)(*pc);
    pc = pc+1;
    Push(d);
}


/*
 *  ---------------------------------------------------------
 *  Callback Functions for Looping and Conditional Constructs
 *  ---------------------------------------------------------
 */

int While_Code()
{
DATUM d;
Inst *savepc = pc;

     Execute(savepc + 2);                /* test condition */
     d = Pop();

     while(d.q->value) {
	 if(Check_Break())  break;
	 Execute(*((Inst **)(savepc)));  /* body of code */
	 if(returning || Check_Break())  break;
	 Execute(savepc + 2);            /* test condition */
         if(CheckUnits()==ON) {
             free((char *)d.q->dimen->units_name);
             free((char *)d.q->dimen);
         }
         free((char *)d.q);
	 d = Pop();
     }
     if(CheckUnits()==ON) {
         free((char *)d.q->dimen->units_name);
         free((char *)d.q->dimen);
     }
     free((char *)d.q);

     if(returning == 0)
        pc = *((Inst **)(savepc + 1));   /* next statement */

     After_Break();
}

int If_Code() {
DATUM d;
Inst *savepc = pc;

     Execute(savepc + 3);               /* condition part */
     d = Pop();
     if(d.q->value) {
       Execute(*((Inst **)(savepc)));
     }
     else if(*((Inst **)(savepc + 1))) {
	Execute(*((Inst **)(savepc+1)));
     }
     if(CheckUnits()==ON) {
         free((char *)d.q->dimen->units_name);
         free((char *)d.q->dimen);
     }
     free((char *)d.q);
	
     if(!returning)
	pc = *((Inst **)(savepc + 2));   /* next statement */
}

int For_Code()
{
DATUM d;
Inst *savepc = pc;

     Execute(*((Inst **) savepc)) ;     /* initiation */ 
     Execute(*((Inst **) (savepc+1))) ; /* cond */

     d = Pop();
     while(d.q->value) {
        if(Check_Break()) break;
        Execute(*((Inst **)(savepc+2))); /* body of loop */
        Execute(*((Inst **)(savepc+3))); /* increments */
        if(returning || Check_Break()) break;
        Execute(*((Inst **) (savepc+1))) ; /* cond */

        if(CheckUnits()==ON) {
           free((char *)d.q->dimen->units_name);
           free((char *)d.q->dimen);
        }

        free((char *)d.q);
        d = Pop();
     }
     if(CheckUnits()==ON) {
         free((char *)d.q->dimen->units_name);
         free((char *)d.q->dimen);
     }
     free((char *)d.q);

     if(!returning)
        pc = *((Inst **)(savepc + 4));   /* next statement */
     After_Break();
}

/* ========================================= */
/* Put Function or Procedure in Symbol Table */
/* ========================================= */

#ifdef __STDC__
int define( SYMBOL *sp )
#else
int define(sp)
SYMBOL *sp;
#endif
{
     sp->u.defn = (Inst) progbase; /* start of code         */
     progbase   = progp;           /* next code starts here */
}

/* Call a function */

int Call()
{
SYMBOL *sp = (SYMBOL *) pc[0];

     if(fp++ >= &frame_in_code[NFRAME - 1])
        EXECUTION_ERROR(sp->cpSymName, "ERROR >> Call nested too deeply", (char *) 0);

     fp->sp    = sp;
     fp->nargs = (int) pc[1];
     fp->retpc = pc + 2;
     fp->argn  = stackp - 1;

     Execute((Inst *) sp->u.defn);
     returning = 0;
}

int Ret()
{
int i;
DATUM d;

     for(i = 0; i < fp->nargs; i++) {
         d = Pop();
         if(CheckUnits()==ON) {
             free((char *)d.q->dimen->units_name);
             free((char *)d.q->dimen);
         }
         free((char *)d.q);
     }

     pc = (Inst *) fp->retpc;
     --fp;
     returning = 1;
}

/*  
 *  ----------------------------------------------------------------------
 *  The function code() builds up the machine stack prog, adding data
 *  types Inst to it. In practice the arguments to code() are pointers
 *  to the functions defined later on in code.c, and pointers to
 *  data types Symbol, which have been coerced by the use of the 
 *  cast Inst.
 *  ----------------------------------------------------------------------
 */ 

#ifdef __STDC__
Inst *Code( Inst f )
#else
Inst *Code(f)
Inst f;
#endif
{
Inst *oprogp = progp;

   if(progp >= &prog[NPROG])
      EXECUTION_ERROR("ERROR >> program too big", (char *) 0);

    *progp++ = f;
    return oprogp;
}


/* 
 *  ------------------------------------
 *  Functions for Engineering Quantities
 *  ------------------------------------
 */ 

int Quantity_Add()
{
DATUM d1, d2, d3;
int    UnitsType;

   d2 = Pop();  /* second item, d2.q */
   d1 = Pop();  /* first item,  d1.q */

   switch(CheckUnits()) {
     case ON:
        if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE) {
           UnitsType = CheckUnitsType();     
           if(d2.q->dimen->units_type == UnitsType ||
              d2.q->dimen->units_type == SI_US     ) {
              d2.q->value = d2.q->value + d1.q->value;
              free((char *) d1.q->dimen->units_name);
              free((char *) d1.q->dimen);
              free((char *) d1.q);
              Push(d2);
           }
           else {
              d1.q->value = d2.q->value + d1.q->value;
              free((char *) d2.q->dimen->units_name);
              free((char *) d2.q->dimen);
              free((char *) d2.q);
              Push(d1);
           }
        }
        else
           FatalError("In Add() : Inconsistent Dimensions",(char *)NULL);
        break;
     case OFF:
        d2.q->value = d2.q->value + d1.q->value;
        free((char *) d1.q);
        Push(d2);
        break;
     default:
        FatalError("In Add() : CheckUnits is not ON or OFF",(char *)NULL);
        break;
   }
}

int Quantity_Sub()
{
DATUM d1, d2;
int UnitsType;

    d2 = Pop();   /* 2nd quantity, d2.q */
    d1 = Pop();   /* 1st quantity, d1.q */

    switch( CheckUnits() )  {
      case ON:
        if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE) {
           UnitsType = CheckUnitsType();     
           if(d2.q->dimen->units_type == UnitsType ||
              d2.q->dimen->units_type == SI_US ) {
              d2.q->value = d1.q->value - d2.q->value;
              free((char *) d1.q->dimen->units_name);
              free((char *) d1.q->dimen);
              free((char *) d1.q);
              Push(d2);
           }
           else {
              d1.q->value = d1.q->value - d2.q->value;
              free((char *) d2.q->dimen->units_name);
              free((char *) d2.q->dimen);
              free((char *) d2.q);
              Push(d1);
           }
        }
        else
           FatalError("In Sub() : Inconsistent Dimensions",(char *)NULL);
        break;
      case OFF:
        d1.q->value = d1.q->value - d2.q->value;
        free((char *) d2.q);
        Push(d1);
        break;
     default:
        FatalError("In Sub() : CheckUnits is not ON or OFF",(char *)NULL);
        break;
    }
}

int Quantity_Mul()
{
DATUM  d1, d2, d3;
int        length;

    d2 = Pop();   /* 2nd quantity, d2.q */
    d1 = Pop();   /* 1st quantity, d1.q */
    d3.q         = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));

    switch( CheckUnits() ) {
      case ON:
         d3.q->value  = d1.q->value*d2.q->value;
         d3.q->dimen  = UnitsMult( d1.q->dimen, d2.q->dimen );
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         free((char *) d1.q->dimen);
         break;
      case OFF:
         d3.q->value  = d1.q->value*d2.q->value;
         break;
      default:
         FatalError("In Mul() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    free((char *) d1.q);
    Push(d3);
}

int Quantity_Div()
{
DATUM d1, d2, d;
int      length;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */
    d.q        = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));

    if(d2.q->value != 0)
       d.q->value = d1.q->value/d2.q->value;
    else
       FatalError("In Quantity_Div() : Attempt to divide by zero",(char *)NULL);

    switch(CheckUnits()) {
      case ON:
         d.q->dimen = UnitsDiv( d1.q->dimen, d2.q->dimen, NO );
         free((char *) d1.q->dimen->units_name);
         free((char *) d1.q->dimen);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         break;
      case OFF:
         break;
      default:          
         FatalError("In Div() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d1.q);
    free((char *) d2.q);
    Push(d);
}

int Quantity_Mod()
{
DATUM d1, d2, d;
int      length;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */
    d.q        = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));

    if(d2.q->value != 0) {
       if( ( (double)(int)d1.q->value != d1.q->value ) ||
           ( (double)(int)d2.q->value != d2.q->value ) )
          FatalError("In Quantity_Mod() : operands of % (mod) must be integers",(char *)NULL);
       else
          d.q->value = ((int) d1.q->value) % ((int) d2.q->value);
    }
    else
       FatalError("In Quantity_Mod() : Attempt to divide by zero",(char *)NULL);

    switch(CheckUnits()) {
      case ON:
         d.q->dimen = UnitsDiv( d1.q->dimen, d2.q->dimen, NO );
         free((char *) d1.q->dimen->units_name);
         free((char *) d1.q->dimen);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         break;
      case OFF:
         break;
      default:          
         FatalError("In Mod() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d1.q);
    free((char *) d2.q);
    Push(d);
}

int Quantity_Negate()
{
DATUM d1;

    d1 = Pop();
    d1.q->value = -d1.q->value;
    Push(d1);
}

#ifdef __STDC__
int dComp_Le( double x, double y )
#else
int dComp_Le(x, y)
double x, y;
#endif
{
static int FLAG; 

   FLAG = (x < y );
   if(dComp_Eq(x,y)) FLAG = 1;
   return (FLAG);
}

#ifdef __STDC__
int dComp_Ge( double x, double y )
#else
int dComp_Ge(x, y)
double x, y;
#endif
{
static int FLAG; 

   FLAG = (x > y );
   if(dComp_Eq(x,y)) FLAG = 1;
   return (FLAG);
}

#ifdef __STDC__
int dComp_Eq( double x, double y )
#else
int dComp_Eq(x,y)
double x,y;
#endif
{
static int FLAG; 
double ERROR;
       
   ERROR = x - y;
   if((ABS(ERROR) == 0.0) || (ABS(ERROR) != 0.0 && ABS(ERROR) < 1E-15))
       FLAG = 1;
   else
       FLAG = 0;
   return (FLAG);
}

#ifdef __STDC__
dComp_Ne( double x, double y )
#else
dComp_Ne(x,y)
double x,y;
#endif
{
static int FLAG; 
double ERROR;
       
   ERROR = x - y;
   if(ABS(ERROR) != 0.0 && ABS(ERROR) > 1E-15)
       FLAG = 1;
   else
       FLAG = 0;
   return (FLAG);
}

int Quantity_Gt()
{
DATUM d1, d2;

    d2     = Pop();  /* 2nd quantity, d2.q */
    d1     = Pop();  /* 1st quantity, d1.q */

    switch( CheckUnits() ) {
      case ON:
         if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE)
            d1.q->value = (double) (d1.q->value > d2.q->value);
         else
            FatalError("In Quantity_Gt() : Inconsistent Dimensions",(char *)NULL);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         d1.q->dimen->units_name = (char *)NULL;
         break;
      case OFF:
         d1.q->value = (double) (d1.q->value > d2.q->value);
         break;
      default:          
         FatalError("In Quantity_Gt() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Lt()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    switch( CheckUnits() ) {
      case ON:
         if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE) {
            d1.q->value = (double) (d1.q->value < d2.q->value);
         }
         else
            FatalError("In Quantity_Lt() : Inconsistent Dimensions",(char *)NULL);

         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         d1.q->dimen->units_name = (char *)NULL;
         break;
      case OFF:
         d1.q->value = (double) (d1.q->value < d2.q->value);
         break;
      default:          
         FatalError("In Quantity_Lt() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Ge()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    switch( CheckUnits() ) {
      case ON:
         if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE)
            d1.q->value = (double) dComp_Ge(d1.q->value, d2.q->value);
         else
            FatalError("In Quantity_Ge() : Inconsistent Dimensions",(char *)NULL);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         d1.q->dimen->units_name = (char *)NULL;
         break;
      case OFF:
         d1.q->value = (double) dComp_Ge(d1.q->value, d2.q->value);
         break;
      default:          
         FatalError("In Quantity_Ge() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Le()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    switch( CheckUnits() ) {
      case ON:
         if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE)
            d1.q->value = (double) dComp_Le(d1.q->value, d2.q->value);
         else
            FatalError("In Quantity_Le() : Inconsistent Dimensions",(char *)NULL);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         d1.q->dimen->units_name = (char *)NULL;
         break;
      case OFF:
         d1.q->value = (double) dComp_Le(d1.q->value, d2.q->value);
         break;
      default:          
         FatalError("In Quantity_Le() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Eq()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    switch( CheckUnits() ) {
      case ON:
         if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE)
            d1.q->value = (double) dComp_Eq(d1.q->value, d2.q->value);
         else
            FatalError("In Quantity_Eq() : Inconsistent Dimensions",(char *)NULL);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         d1.q->dimen->units_name = (char *)NULL;
         break;
      case OFF:
         d1.q->value = (double) dComp_Eq(d1.q->value, d2.q->value);
         break;
      default:          
         FatalError("In Quantity_Eq() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Ne()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    switch( CheckUnits() ) {
      case ON:
         if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE)
            d1.q->value = (double) dComp_Ne(d1.q->value, d2.q->value);
         else
            FatalError("In Quantity_Ne() : Inconsistent Dimensions",(char *)NULL);
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d1.q->dimen->units_name);
         d1.q->dimen->units_name = (char *)NULL;
         break;
      case OFF:
         d1.q->value = (double) dComp_Ne(d1.q->value, d2.q->value);
         break;
      default:          
         FatalError("In Quantity_Ne() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_And()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    d1.q->value = (double) (d1.q->value != 0 && d2.q->value != 0);
    if(CheckUnits()==ON) {
       free((char *) d2.q->dimen->units_name);
       free((char *) d2.q->dimen);
       free((char *) d1.q->dimen->units_name);
       d1.q->dimen->units_name = (char *)NULL;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Or()
{
DATUM d1, d2;

    d2 = Pop();  /* 2nd quantity, d2.q */
    d1 = Pop();  /* 1st quantity, d1.q */

    d1.q->value = (double) (d1.q->value != 0 || d2.q->value != 0);
    if(CheckUnits()==ON) {
       free((char *)d2.q->dimen->units_name);
       free((char *)d2.q->dimen);
       free((char *)d1.q->dimen->units_name);
       d1.q->dimen->units_name = (char *)NULL;
    }
    free((char *) d2.q);
    Push(d1);
}

int Quantity_Yes()
{
DATUM d1;

     d1.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
     if(CheckUnits() == ON){
        d1.q->dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
        d1.q->dimen->units_name = (char *)NULL;
     }
     else
        d1.q->dimen = (DIMENSIONS *)NULL;
     d1.q->value = (double) (TRUE);
     Push(d1);
}

int Quantity_Not()
{
DATUM d1;

     d1 = Pop();
     d1.q->value = (double) (d1.q->value == 0);

     if( CheckUnits() == ON ) {
        free((char *) d1.q->dimen->units_name);
        d1.q->dimen->units_name = (char *)NULL;
     }
     Push(d1);
}

int Quantity_Power()
{
DATUM     d1, d2, d;
int          length;

  /* [a] : Pop quantities from stack and allocate memory for result */

     d2  = Pop();  /* 2nd quantity, d2.q */
     d1  = Pop();  /* 1st quantity, d1.q */
     d.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));

  /* [a] : Pop quantities from stack and allocate memory for result */

     switch( CheckUnits() ) {
       case ON:
            if(((d1.q->dimen->length_expnt != 0) ||
                (d1.q->dimen->mass_expnt   != 0) ||
                (d1.q->dimen->time_expnt   != 0) ||
                (d1.q->dimen->temp_expnt   != 0)) &&
               ((d2.q->dimen->length_expnt == 0) &&
                (d2.q->dimen->mass_expnt   == 0) &&
                (d2.q->dimen->time_expnt   == 0) &&
                (d2.q->dimen->temp_expnt   == 0))) { 
                 d.q->value = pow(d1.q->value, d2.q->value);
                 d.q->dimen = UnitsPower( d1.q->dimen, d2.q->value, NO );
            } else if (((d1.q->dimen->length_expnt == 0) &&
                (d1.q->dimen->mass_expnt   == 0) &&
                (d1.q->dimen->time_expnt   == 0) &&
                (d1.q->dimen->temp_expnt   == 0)) &&
               ((d2.q->dimen->length_expnt == 0) ||
                (d2.q->dimen->mass_expnt   == 0) ||
                (d2.q->dimen->time_expnt   == 0) ||
                (d2.q->dimen->temp_expnt   == 0))) { 
                 d.q->value = pow(d1.q->value, d2.q->value);
                 d.q->dimen = UnitsPower( d2.q->dimen, 1.0, NO );
            } else if (((d1.q->dimen->length_expnt != 0) ||
                (d1.q->dimen->mass_expnt   != 0) ||
                (d1.q->dimen->time_expnt   != 0) ||
                (d1.q->dimen->temp_expnt   != 0)) &&
               ((d2.q->dimen->length_expnt != 0) ||
                (d2.q->dimen->mass_expnt   != 0) ||
                (d2.q->dimen->time_expnt   != 0) ||
                (d2.q->dimen->temp_expnt   != 0))) { 
                 FatalError("In Quantity_Power() : Illegal Power Operation",
                           (char *) NULL);
            }
            free((char *) d2.q->dimen->units_name);
            free((char *) d2.q->dimen);
            free((char *) d1.q->dimen->units_name);
            free((char *) d1.q->dimen);
            break;
       case OFF:
            d.q->value = pow(d1.q->value, d2.q->value);
            break;
       default:          
            FatalError("In Quantity_Power() : CheckUnits is not ON or OFF",
                        (char *)NULL);
            break;
     }
     free((char *) d2.q);
     free((char *) d1.q);
     Push(d);
}

int Quantity_Affirm()
{
    Push(Pop());
}

int Quantity_Extract()
{
DATUM d1, d2, d3, d;
int i, j, length;
int iMin, iMax;

     d3 = Pop();  /* input matrix, d3.m */
     d2 = Pop();  /* index j , d2.q */
     d1 = Pop();  /* index i , d1.q */
     i   = (int) d1.q->value;
     j   = (int) d2.q->value;
     d.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));

     /* Check arguments are defined */

     if((i < 1) || (j<1) || (i > d3.m->iNoRows) || (j > d3.m->iNoColumns)) {
        printf("ERROR >> In Quantity_Extract() : Array cpMatrixName = %s\n", d3.m->cpMatrixName);
        printf("ERROR >> In Quantity_Extract() : Array dimensions = [%4d][%4d]\n",
                                                       d3.m->iNoRows, d3.m->iNoColumns);
        printf("ERROR >> In Quantity_Extract() : Arguments [i][j] = [%4d][%4d]\n", i,j);
        FatalError(" Array arguments a[i][j] out of bounds",(char *)NULL);
     }

     if(d3.m->eType == DOUBLE_ARRAY) {
        switch( (int) d3.m->eRep ) {
            case INDIRECT:
                 d.q->value = d3.m->uMatrix.daa[i-1][j-1];
                 break;
            case SKYLINE:
                 iMin = MIN(i,j);
                 iMax = MAX(i,j);
                 if((iMax - iMin + 1) <= d3.m->uMatrix.daa[iMax-1][0])
                    d.q->value = d3.m->uMatrix.daa[iMax-1][iMax-iMin+1];
                 else
                    d.q->value = 0.0;
                 break;
            default:          
                 FatalError("In Quantity_Extract() : Cannot extract matrix type",
                            (char *) NULL);
        }

        switch( CheckUnits() ) {
          case ON:
             d.q->dimen =  UnitsMult( &(d3.m->spColUnits[j-1]), &(d3.m->spRowUnits[i-1]) );
             free((char *)d2.q->dimen->units_name);
             free((char *)d2.q->dimen);
             free((char *)d1.q->dimen->units_name);
             free((char *)d1.q->dimen);
             break;
          case OFF:
             break;
          default:          
             FatalError("In Quantity_Extract() : CheckUnits is not ON or OFF",(char *)NULL);
             break;
        }
     }
     free((char *)d1.q);
     free((char *)d2.q);
     MatrixFree( d3.m );
     Push(d);
}

/*
 *  ======================================================
 *  Evaluation of Variables and Strings
 *  
 *  Purpose : Transfer contents from Symbol Table to Stack
 *  ======================================================
 */

int String_Eval()
{
DATUM  d1;

    d1 = Pop();
    Push(d1);
}

int Variable_Eval()
{
DATUM  d1, d2;
int         i;
int   length1;

   /* variable name, d1.sym, stored in symbol table, cannot be freed */
    d1 = Pop();

    if(d1.sym->type != VAR && d1.sym->type != QUAN)
       FatalError("Attempt to evaluate non-variable\n", (char *)NULL);

    d2.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
    if( d1.sym->u.q != (QUANTITY *)NULL ) {
        switch( CheckUnits() ) {
          case ON:
             d2.q->value  = d1.sym->u.q->value;
             d2.q->dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             UnitsCopy( d2.q->dimen, d1.sym->u.q->dimen );
             break;
          case OFF:
             d2.q->value  = d1.sym->u.q->value;
             d2.q->dimen  = (DIMENSIONS *)NULL;
             break;
          default:          
             FatalError("In Variable_Eval() : CheckUnits is not ON or OFF",(char *)NULL);
             break;
        }
    }
    else
      if( CheckUnits()==ON )
          d2.q->dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));

    Push(d2);
}

int Matrix_Eval()
{
DATUM  d1, d2;
int iNoRows;
int iNoColumns;
int i, j;

   /* matrix, stored in symbol table, cannot be free */
    d1 = Pop();
    if(d1.sym->type != MATX) {
       printf("FATAL ERROR >> Execution halted in Matrix_Eval()\n");
       printf("FATAL ERROR >> Problem : %s->type != MATX\n", d1.sym->cpSymName);
       printf("FATAL ERROR >> Recommendation : Check input file !! \n");
       FatalError("Assignment to non-matrix",(char *)NULL);
    }
    d2.m = MatrixCopy(d1.sym->u.m);
    Push(d2);
}

int Matrix_Build()
{
DATUM    d1, d2, d3;
int iNoRows, row_length;
int             i,j,r,s;
int              length;
int        UNITS_SWITCH;

    d1 = Pop();  /* constant= no of rows = address of d1.sym */
    iNoRows = (int) d1.sym;

    d2.m = (MATRIX *) MyCalloc(1,sizeof(MATRIX));
    d2.m->cpMatrixName    = (char *) NULL;         /* matrix name */
    d2.m->eType    = DOUBLE_ARRAY;
    d2.m->eRep     = INDIRECT;
    d2.m->iNoRows  = iNoRows;
    d2.m->uMatrix.daa = (double **) MyCalloc(iNoRows, sizeof(double *));

    UNITS_SWITCH = CheckUnits();
    for(i = 1; i <= iNoRows; i++) {
        d3  = Pop();    /* constant= no of columns = address of d3.sym */
        row_length = (int) d3.sym;
        r = iNoRows - i;
        d2.m->uMatrix.daa[r]    = (double *) MyCalloc(row_length, sizeof(double));

        if(i == 1) {
          d2.m->iNoColumns    = row_length;
          if( UNITS_SWITCH == ON ) {
              d2.m->spColUnits = (DIMENSIONS *) MyCalloc(row_length, sizeof(DIMENSIONS));
              d2.m->spRowUnits = (DIMENSIONS *) MyCalloc(iNoRows, sizeof(DIMENSIONS));
          }
          else {
              d2.m->spColUnits = (DIMENSIONS *)NULL;
              d2.m->spRowUnits = (DIMENSIONS *)NULL;
          }
        }
        
        /* check for no of columns in matrix */

        if(i != 1 && d2.m->iNoColumns != row_length ) {
           printf(" In Matrix_Build(): \n");
           FatalError("==> no of columns in each row of a matrix must be same",(char *)NULL);
        }

        for(j = 1; j <= row_length; j++) {
            d3  = Pop();    /* matrix items, d3.q */
            s = row_length - j;
            d2.m->uMatrix.daa[r][s] = (double) d3.q->value;

        /* Units of matrix elements are stored in column_units_buf as an default */

           switch( UNITS_SWITCH ) {
             case ON:
                if( (row_length != 1) || ((iNoRows == 1) && (row_length == 1))  ) {
                   UnitsCopy( &(d2.m->spColUnits[s]), d3.q->dimen );
                }
                else {
                   ZeroUnits( &(d2.m->spColUnits[s]) );
                   d2.m->spColUnits[s].units_type = SI_US;
                }

                if((row_length == 1) && (iNoRows > 1) ) { /* a vector */
                   break;
                }
                else {
                   free((char *) d3.q->dimen->units_name);
                   free((char *) d3.q->dimen);
                   free((char *) d3.q);
                }
                break;
             case OFF:
                if((row_length == 1) && (iNoRows > 1) )   /* a vector */
                   break;
                else
                   free((char *) d3.q);
                break;
             default:          
                FatalError("In Matrix_Build() : CheckUnits is not ON or OFF",(char *)NULL);
                break;
           }
        }
        if( UNITS_SWITCH == ON ) {
            if((row_length == 1) && (iNoRows > 1) ) { /* a vector */
               UnitsCopy( &(d2.m->spRowUnits[r]), d3.q->dimen );
               free((char *) d3.q->dimen->units_name);
               free((char *) d3.q->dimen);
               free((char *) d3.q);
            } 
            else {
               ZeroUnits( &(d2.m->spRowUnits[r]) );
               d2.m->spRowUnits[r].units_type = SI_US;
            }
        }
   }
   Push(d2);
}

int Assign_Matrix()
{
DATUM        d1, d2, d3;
int iNoRows, iNoColumns;
int                i, j;
int              length;
int                 *ld;
int          iColHeight;
int        UNITS_SWITCH;

    d1     = Pop();  /* matrix name , d1.sym */
    d2     = Pop();  /* matrix contents , d2.m */

    iNoRows    = d2.m->iNoRows;
    iNoColumns = d2.m->iNoColumns;

    /* Check for assignment to "variable" or "matrix" */

    if(d1.sym->type != VAR && d1.sym->type != MATX)
       EXECUTION_ERROR("Assignment to non-variable\n", d1.sym->cpSymName);
 
    /* Check for change in "dimensions" and "type" of matrix */

    if(d1.sym->type == MATX && d1.sym->u.m != (MATRIX *) NULL)  {
       if(d1.sym->u.m->iNoRows    != iNoRows     ||
          d1.sym->u.m->iNoColumns != iNoColumns  ||
          d1.sym->u.m->eType      != d2.m->eType ||
          d1.sym->u.m->eRep       != d2.m->eRep  ) {
             MatrixFree(d1.sym->u.m);
             d1.sym->u.m = (MATRIX *) NULL;
          }
    }

    /* Allocate new matrix if needed */

    if( d1.sym->u.m == (MATRIX *) NULL || d1.sym->type != MATX) {

       d1.sym->type = MATX;
       switch((int) d2.m->eRep ) {
           case INDIRECT:
           {
               switch((int) d2.m->eType) {
                   case DOUBLE_ARRAY:
                        d1.sym->u.m 
                        = MatrixAllocIndirect((char *) d1.sym->cpSymName,
                          DOUBLE_ARRAY, iNoRows, iNoColumns);
                        break;
                   default:
                        FatalError("In Assign_Matrix() : Undefined Matrix Type",
                                  (char *)NULL);
                        break;
               }
               break;
           }
           case SKYLINE:
           {
               ld = (int *)iVectorAlloc( d2.m->iNoRows );
               for( i=0 ; i<d2.m->iNoRows ; i++ )
                   ld[i] = (int) d2.m->uMatrix.daa[i][0];
               switch((int) d2.m->eType) {
                   case DOUBLE_ARRAY:
                        d1.sym->u.m 
                        = MatrixAllocSkyline((char *) d1.sym->cpSymName,
                          DOUBLE_ARRAY, iNoRows, iNoColumns, ld);
                        free((char *) ld);
                        break;
                   default:
                        FatalError("In Assign_Matrix() : Undefined Matrix Type",
                        (char *)NULL);
                        break;
               }
               break;
           }
           default:
               FatalError("In Assign_Matrix() : Undefined Matrix eRep",
               (char *)NULL);
               break;
       } /* end of switch( d2.m->eRep ) */

    } /* end of if */

    /* Copy contents of Matrix d2.m to d1.sym->u.m with no reallocation of memory */

    UNITS_SWITCH = CheckUnits();
    switch((int) d2.m->eRep ) {
      case INDIRECT:
      {
         switch((int) d2.m->eType) {
           case DOUBLE_ARRAY:
           {
             for(i = 1; i <= iNoRows; i++) {
                if( UNITS_SWITCH == ON )
                    UnitsCopy( &(d1.sym->u.m->spRowUnits[i-1]), &(d2.m->spRowUnits[i-1]) ); 

                for(j = 1; j <= iNoColumns; j++)
                   d1.sym->u.m->uMatrix.daa[i-1][j-1]
                   = d2.m->uMatrix.daa[i-1][j-1];
             }
             if( UNITS_SWITCH == ON ) {
                 for(j = 1; j <= iNoColumns; j++)
                    UnitsCopy( &(d1.sym->u.m->spColUnits[j-1]), &(d2.m->spColUnits[j-1]) ); 
             }
             break;
           }
           default:
                FatalError("In Assign_Matrix() : Undefined Matrix Type",
               (char *)NULL);
                break;
         } /* end of case INDIRECT, switch( eType ) */
         break;
      } /* end of case INDIRECT */
      case SKYLINE:
      {
         switch((int) d2.m->eType) {
           case DOUBLE_ARRAY:
           {
             for(i = 1; i <= iNoRows; i++) {
                if( UNITS_SWITCH == ON )
                    UnitsCopy(&(d1.sym->u.m->spRowUnits[i-1]),&(d2.m->spRowUnits[i-1])); 

                if( d2.m->uMatrix.daa[i-1][0] != d1.sym->u.m->uMatrix.daa[i-1][0] ) {
                   free( d1.sym->u.m->uMatrix.daa[i-1] );
                   d1.sym->u.m->uMatrix.daa[i-1] 
                   = (double *)MyCalloc((d2.m->uMatrix.daa[i-1][0]+1), sizeof(double));
                   d1.sym->u.m->uMatrix.daa[i-1][0] = d2.m->uMatrix.daa[i-1][0];
                }
                for(j = 1; j <= d2.m->uMatrix.daa[i-1][0]; j++) 
                   d1.sym->u.m->uMatrix.daa[i-1][j]  = d2.m->uMatrix.daa[i-1][j];
             }
             if( UNITS_SWITCH == ON ) {
                 for(j = 1; j <= iNoColumns; j++)
                    UnitsCopy(&(d1.sym->u.m->spColUnits[j-1]), &(d2.m->spColUnits[j-1])); 
             }
             break;
           } /* end of case SKYLINE, DOUBLE_ARRAY */
           default:
                FatalError("In Assign_Matrix() : Undefined Matrix Type",
                          (char *)NULL);
                break;
         } /* end of case SKYLINE , switch( eType ) */
         break;
      } /* end of case SKYLINE */
      default:
           FatalError("In Assign_Matrix() : Undefined Matrix eRep",(char *)NULL);
           break;
    } /* end of switch(d2.m->eRep) */

    MatrixFree(d2.m);
    Push(d1);
}

int Assign_Matrix_Item()
{
DATUM  d1, d2, d3, d4, d5;
int row, column;
int i, j, length;
int iMin, iMax;
int UNITS_SWITCH;

DIMENSIONS *dm1;
DIMENSIONS *dm_ptr;

    d1 = Pop();         /* matrix name , d1.sym    */
    d2 = Pop();         /* matrix item , d2.q      */
    d3 = Pop();         /* matrix iNoColumns, d3.q */
    d4 = Pop();         /* matrix iNoRows   , d4.q */

    UNITS_SWITCH = CheckUnits();

    row    = (int) d4.q->value;
    column = (int) d3.q->value;
    if( UNITS_SWITCH == ON ) {
      free((char *) d3.q->dimen->units_name);
      free((char *) d3.q->dimen);
      free((char *) d4.q->dimen->units_name);
      free((char *) d4.q->dimen);
    }
    free((char *) d3.q);
    free((char *) d4.q);

    /* Check Variable type and Matrix dimensions */

    if(d1.sym->type != MATX)
       FatalError("In Assign_Matrix_Item : Assignment to non-matrix",(char *)NULL);

    if((row < 1) || (row > d1.sym->u.m->iNoRows))
       FatalError("In Assign_Matrix_Item : row out of array bounds",(char *)NULL);
     
    if((column < 1) || (column > d1.sym->u.m->iNoColumns))
       FatalError("In Assign_Matrix_Item : column out of array bounds",(char *)NULL);

    /* Update Matrix Item */

    switch( d1.sym->u.m->eRep ) {
       case INDIRECT :
          switch((int) d1.sym->u.m->eType) {
              case DOUBLE_ARRAY:

                switch( UNITS_SWITCH ) {
                  case  ON:
                  /* [a]  calculate the units of m[row][col] */


                     dm1 = UnitsMult( &(d1.sym->u.m->spRowUnits[row-1]),
                                      &(d1.sym->u.m->spColUnits[column-1]) );

                     if(dm1->units_name  != (char *)NULL  ||
                        dm1->length_expnt != 0     ||
                        dm1->mass_expnt != 0       ||
                        dm1->time_expnt != 0       ||
                        dm1->temp_expnt != 0 ) {
                     
                        /* check units of m[row][col]    */
                     
                        if(SameUnits(dm1, d2.q->dimen) != TRUE) {
                           printf("\n **** WARNING : in  Assign_Matrix_Item(): \n");
                           printf(" =====> %s[%d][%d] and %15.4e %s have inconsistent units \n\n",
                                  d1.sym->u.m->cpMatrixName, row, column, d2.q->value/d2.q->dimen->scale_factor,
                                  d2.q->dimen->units_name);
                           printf(" =====> This assignment will NOT overwrite the previous units\n");
                        }
                        else
                          /* Assign the value to matrix element */
                          d1.sym->u.m->uMatrix.daa[row-1][column-1] = d2.q->value;
                     }
                     else {

                       /* [x] m[row][column] has no unit                */
                       /* Save units in the col_units_buf as an default */

                       if(d2.q->dimen->units_name != NULL    ||
                          d2.q->dimen->length_expnt != 0     ||
                          d2.q->dimen->mass_expnt != 0       ||
                          d2.q->dimen->time_expnt != 0       ||
                          d2.q->dimen->temp_expnt != 0 ) {
/*
                          printf("\n ****** WARNING: Matrix %s element %s[%d][%d] \n",
                                 d1.sym->u.m->cpMatrixName,d1.sym->u.m->cpMatrixName, row, column);
                          printf("****WARNING: assign a units to a single element might affect the units of other elements\n\n");
*/
                          d1.sym->u.m->uMatrix.daa[row-1][column-1] = d2.q->value;
                          UnitsCopy( &(d1.sym->u.m->spColUnits[column-1]), d2.q->dimen );
                       }
                       else /* both matrix and d2.q are non-dimensional */
                          d1.sym->u.m->uMatrix.daa[row - 1][column - 1] = d2.q->value;
                     }

                     free((char *) d2.q->dimen->units_name);
                     free((char *) d2.q->dimen);
                     free((char *) dm1->units_name);
                     free((char *) dm1);

                     break;   /* end of case UNITS_SWITCH==ON */
                  case OFF:
                     d1.sym->u.m->uMatrix.daa[row-1][column-1] = d2.q->value;
                     break;
                  default:          
                     FatalError("In Assign_Matrix_Item() : CheckUnits is not ON or OFF",(char *)NULL);
                     break;
                }
                break;  /* end of case INDIRECT, DOUBLE_ARRAY */
              default:
                FatalError("In Assign_Matrix_Item() : Undefined m->eType",(char *)NULL);
                break;
          }
          break;  /* end of case INDIRECT */
       case SKYLINE:
          switch((int) d1.sym->u.m->eType) {
              case DOUBLE_ARRAY:

                switch( UNITS_SWITCH ) {
                   case ON:
                     /* [a]  calculate the units of m[row][col] */

                     dm1 = UnitsMult( &(d1.sym->u.m->spRowUnits[row-1]),
                                      &(d1.sym->u.m->spColUnits[column-1]) );

                     if(dm1->units_name  != NULL  ||
                        dm1->length_expnt != 0     ||
                        dm1->mass_expnt != 0       ||
                        dm1->time_expnt != 0       ||
                        dm1->temp_expnt != 0 ) {
                     
                        /* check units of m[row][col]    */
                     
                        if(SameUnits(dm1, d2.q->dimen) != TRUE) {
                           printf("\n **** WARNING : in  Assign_Matrix_Item(): \n");
                           printf(" =====> %s[%d][%d] and %15.4e %s have inconsistent units \n\n",
                                  d1.sym->u.m->cpMatrixName, row, column, d2.q->value/d2.q->dimen->scale_factor,
                                  d2.q->dimen->units_name);
                           printf(" =====> This assignment will NOT overwrite the previous units\n");
                        }
                        else {
                          /* Assign the value to matrix element */
                           iMin = MIN(row, column);
                           iMax = MAX(row, column);
                           if((iMax-iMin+1) <= d1.sym->u.m->uMatrix.daa[iMax-1][0])
                               d1.sym->u.m->uMatrix.daa[iMax-1][iMax-iMin+1] = d2.q->value;
                           else {
                               d1.sym->u.m->uMatrix.daa[iMax-1]
                               = (double *) realloc(d1.sym->u.m->uMatrix.daa[iMax-1],
                                 (iMax-iMin+1));
                               for( j=(int) d1.sym->u.m->uMatrix.daa[iMax-1][0]+1; j < iMax-iMin+1 ; j++ )
                                   d1.sym->u.m->uMatrix.daa[iMax-1][j] = 0.0;
                               d1.sym->u.m->uMatrix.daa[iMax-1][0] = iMax-iMin+1;
                               d1.sym->u.m->uMatrix.daa[iMax-1][iMax-iMin+1] = d2.q->value;
                           }
                        }
                     }
                     free((char *) d2.q->dimen->units_name);
                     free((char *) d2.q->dimen);
                     free((char *) dm1->units_name);
                     free((char *) dm1);
                     break;  /* end of case UNITS_SWITCH == ON */
                  case OFF:
                     iMin = MIN(row, column);
                     iMax = MAX(row, column);
                     if((iMax-iMin+1) <= d1.sym->u.m->uMatrix.daa[iMax-1][0])
                         d1.sym->u.m->uMatrix.daa[iMax-1][iMax-iMin+1] = d2.q->value;
                     else {
                         d1.sym->u.m->uMatrix.daa[iMax-1] 
                         = (double *) realloc(d1.sym->u.m->uMatrix.daa[iMax-1], (iMax-iMin+1));
                         for( j=(int) d1.sym->u.m->uMatrix.daa[iMax-1][0]+1; j < iMax-iMin+1 ; j++ )
                             d1.sym->u.m->uMatrix.daa[iMax-1][j] = 0.0;
                         d1.sym->u.m->uMatrix.daa[iMax-1][0] = iMax-iMin+1;
                         d1.sym->u.m->uMatrix.daa[iMax-1][iMax-iMin+1] = d2.q->value;
                     }
                     break;  /* end of case UNITS_SWITCH == OFF */
                  default:
                     FatalError("In Matrix_Assign_Item(): CheckUnits() is not ON or OFF",
                     (char *)NULL);
                     break;
                  }
                  break;  /* end of case SKYLINE, DOUBLE_ARRAY */
              default:
                  FatalError("In Assign_Matrix_Item() : Undefined m->eType",(char *)NULL);
                  break;
          }
          break;  /* end of case SKYLINE */
       default:
          FatalError("In Assign_Matrix_Item() : Undefined m->eRep",(char *)NULL);
          break;
    }  /* end of switch of m->eRep */

    free((char *) d2.q);

    Push(d1);
}

int Assign_Quantity()
{
DATUM     d1, d2;
int       length;
int UNITS_SWITCH;

    d1 = Pop();   /* variable name , d1.sym */
    d2 = Pop();   /* value, scale factor and dimensions, d2.q */

    if(d1.sym->type != VAR && d1.sym->type != QUAN) {
       FatalError("In Assign_Quantity() : Assignment to non-variable",(char *)NULL);
    }

    UNITS_SWITCH = CheckUnits();
    if(d1.sym->type != QUAN) {
       d1.sym->u.q   = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
       d1.sym->type  = QUAN;
       if( UNITS_SWITCH == ON )
          d1.sym->u.q->dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
    }

    switch( UNITS_SWITCH ) {
      case ON:
         d1.sym->u.q->value = d2.q->value;
         UnitsCopy( d1.sym->u.q->dimen,  d2.q->dimen );
         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         break;
      case OFF:
         d1.sym->u.q->value             = d2.q->value;
         break;
      default:          
         FatalError("In Assign_Quantity() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d2.q);
    Push(d1);
}


/*
 *  ------------------------------------------
 *  Functions to Print Quantities and Matrices
 *  ------------------------------------------
 */

int Print_Expr()
{
DATUM d;
int UnitsType;

    d = Pop();

    switch(CheckUnits()) {
      case ON:
         UnitsSimplify(d.q->dimen);
         if(d.q->dimen->units_name != NULL) {
           UnitsType = CheckUnitsType(); 
           switch(UnitsType) { 
             case SI:
                  if(!strcmp(d.q->dimen->units_name, "deg_F") )
                     d.q->value = ConvertTempUnits(d.q->dimen->units_name, d.q->value, US); 
                  break;
             case US:
                  if(!strcmp(d.q->dimen->units_name, "deg_C") ) 
                     d.q->value = ConvertTempUnits(d.q->dimen->units_name, d.q->value, SI); 
                  break;
           }
           printf("%10.4g ", d.q->value/d.q->dimen->scale_factor);
           printf("%s ", d.q->dimen->units_name);
           free((char *) d.q->dimen->units_name);
         } else {
           printf("%10.4g ", d.q->value);
         }
         fflush(stdout);
         free((char *) d.q->dimen);
         break;
      case OFF:
         printf("%10.4g ", d.q->value);
         break;
      default:          
         FatalError("In Print_Expr() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    free((char *) d.q);
}

/*
 *  =================================================
 *  Print_Dimen_Expr() : Print Dimensional Expression
 *  =================================================
 */

int Print_Dimen_Expr() {
DATUM d1, d2;
int UnitsType;

#ifdef DEBUG
       printf("*** Enter Print_Dimen_Expr()\n");
#endif

    d2 = Pop(); /* units    */
    d1 = Pop(); /* quantity */

#ifdef DEBUG
       printf("*** d1.q->value = %10.5f\n", d1.q->value );
#endif

    switch(CheckUnits()) {
      case ON:
           if(SameUnits(d1.q->dimen, d2.q->dimen) == TRUE) {
           if(d2.q->dimen->units_name != NULL) {
              UnitsType = CheckUnitsType(); 
              switch(UnitsType) {
                 case SI:
                      if(!strcmp(d2.q->dimen->units_name, "deg_F")) 
                         d1.q->value = ConvertTempUnits(d1.q->dimen->units_name, d1.q->value, US);
                      break;
                 case US:
                      if(!strcmp(d2.q->dimen->units_name, "deg_C")) 
                         d1.q->value = ConvertTempUnits(d1.q->dimen->units_name, d1.q->value, SI);
                      break;
              }
              printf("%10.4g ", d1.q->value/d2.q->dimen->scale_factor);
              printf("%s ", d2.q->dimen->units_name);
              free((char *) d2.q->dimen->units_name);
           } else {
              printf("%10.4g ", d1.q->value);
           }
           fflush(stdout);
           free((char *) d1.q->dimen->units_name);
           free((char *) d1.q->dimen);
           free((char *) d2.q->dimen);
           } else {
             printf(" Wanted units_name   = %s\n", d2.q->dimen->units_name);
             printf(" quantity units_name = %s\n", d1.q->dimen->units_name);
             FatalError(" In Print_Dimen_Expr(): ",
                        "Try to print quantity with inconsistent units",
                        (char *) NULL);
           }
           break;
      case OFF:
           printf("%10.4g ", d1.q->value);
           break;
      default:          
           FatalError("In Print_Dimen_Expr() : CheckUnits is not ON or OFF",(char *)NULL);
           break;
    }
    free((char *) d1.q);
    free((char *) d2.q);

#ifdef DEBUG
       printf("*** Leave Print_Dimen_Expr()\n");
#endif

}

int Print_String()
{
    printf("%s", (char *) *pc);
    pc = pc+1;
    fflush(stdout);
}


/* 
 *  ---------------------------------------------
 *  Call Built-in Expression and Matrix Functions 
 *  ---------------------------------------------
 */

int Bltin_Break()
{
     iBREAK_FLAG = (int) 1;
}

int Check_Break()
{
    return(iBREAK_FLAG);
}

int After_Break()
{
    iBREAK_FLAG = (int) 0;
}

int Bltin_Quantity()
{
DATUM  d;

  /* random number function, return a new quantity */
     d.q  = (*(QUANTITY * (*)()) *pc)();
     pc = pc+1;
     Push(d);
}

int Bltin1_Quantity()
{
DATUM  d;

  /* mathematic functions, return the old quantity */
     d = Pop();  /* parameter value, d.q */
     d.q  = (*(QUANTITY * (*)()) *pc)(d.q);
     pc = pc+1;
     Push(d);
}

int Bltin2_Quantity()
{
DATUM d2, d;

  /* mathematic function, Pow(quantity,value), return old quantity */
     d2 = Pop();  /* parameter 2 value, d2.q */
     d  = Pop();  /* parameter 1 value, d.q */
     d.q = (*(QUANTITY * (*)()) *pc)(d.q, (double) d2.q->value);

     if( CheckUnits() == ON ) {
       free((char *) d2.q->dimen->units_name);
       free((char *) d2.q->dimen);
     }
     free((char *) d2.q);
     pc = pc+1;
     Push(d);
}

int Bltin3_Quantity()
{
DATUM d1, d2;

     d1 = Pop();  /* parameter matrix, d1.m */
     d2.q = (*(QUANTITY * (*)()) *pc)(d1.m);
     MatrixFree(d1.m);
     pc = pc+1;
     Push(d2);
}

int Bltin_Matrix()
{
DATUM d;

     d.m = (*(MATRIX * (*)()) *pc)();
     pc = pc+1;
     Push(d);
}

int Bltin1_Matrix()
{
DATUM d1;
DATUM d2;

     d1 = Pop();  /* matrix, d.m */
     d2.m = (*(MATRIX * (*)()) *pc)(d1.m);
     MatrixFree(d1.m);
     pc = pc+1;
     Push(d2);
}

int Bltin2_Matrix()
{
DATUM d1, d2, d3;

     d2 = Pop();
     d1 = Pop();
     d3.m = (*(MATRIX * (*)()) *pc)(d1.m, d2.m);
     MatrixFree( d1.m );
     MatrixFree( d2.m );
     pc = pc+1;
     Push(d3);
}

int Bltin3_Matrix()
{
DATUM d1, d2, d3, d4;

     d3 = Pop();
     d2 = Pop();
     d1 = Pop();
     d4.m = (*(MATRIX * (*)()) *pc)(d1.m, d2.m, d3.m);
     MatrixFree( d1.m );
     MatrixFree( d2.m );
     MatrixFree( d3.m );
     pc = pc+1;
     Push(d4);
}

int BltinVar_Matrix()
{
DATUM d1, d[5], d3;
int no_args, i;

    d1 = Pop();  /* constant = address of d1.sym */
    no_args = (int) d1.sym;

    switch(no_args) {
      case 0:
        d3.m = (*(MATRIX * (*)()) *pc)((MATRIX *) NULL);
      break;
      case 1:
        d[0]  = Pop();
        d3.m = (*(MATRIX * (*)()) *pc)(d[0].m, (MATRIX *) NULL);
        MatrixFree(d[0].m);
      break;
      
      case 2:
        d[1]  = Pop();
        d[0]  = Pop();
        d3.m = (*(MATRIX * (*)()) *pc)(d[0].m, d[1].m, (MATRIX *) NULL);
        MatrixFree(d[1].m);
        MatrixFree(d[0].m);
      break;

      case 3:
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.m = (*(MATRIX * (*)()) *pc)(d[0].m, d[1].m, d[2].m, (MATRIX *) NULL);
        MatrixFree(d[2].m);
        MatrixFree(d[1].m);
        MatrixFree(d[0].m);
      break;

      case 4:
        d[3]  = Pop();
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.m = (*(MATRIX * (*)()) *pc)(d[0].m,d[1].m, d[2].m, d[3].m, (MATRIX *) NULL);
        MatrixFree(d[3].m);
        MatrixFree(d[2].m);
        MatrixFree(d[1].m);
        MatrixFree(d[0].m);
      break;

      case 5:
        d[4]  = Pop();
        d[3]  = Pop();
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.m = (*(MATRIX * (*)()) *pc)(d[0].m, d[1].m,d[2].m, d[3].m, d[4].m, (MATRIX *) NULL);
        MatrixFree(d[4].m);
        MatrixFree(d[3].m);
        MatrixFree(d[2].m);
        MatrixFree(d[1].m);
        MatrixFree(d[0].m);
      break;
    }

     pc = pc+1;
     Push(d3);
}

int BltinVar_Quantity()
{
DATUM d1, d[5], d3;
int no_args, i;
int UNITS_SWITCH;

    d1 = Pop();
    no_args = (int) d1.sym;

    UNITS_SWITCH = CheckUnits();
    switch(no_args) {
      case 1:
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)(d[0].q, (QUANTITY *) NULL);
        if( UNITS_SWITCH == ON ) {
            free((char *) d[0].q->dimen->units_name);
            free((char *) d[0].q->dimen);
        }
        free((char *) d[0].q);
      break;
      
      case 2:
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)(d[0].q, d[1].q, (QUANTITY *) NULL);
        if( UNITS_SWITCH == ON ) {
            free((char *) d[1].q->dimen->units_name);
            free((char *) d[0].q->dimen->units_name);
            free((char *) d[1].q->dimen);
            free((char *) d[0].q->dimen);
        }
        free((char *) d[1].q);
        free((char *) d[0].q);
      break;

      case 3:
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)(d[0].q, d[1].q, d[2].q, (QUANTITY *) NULL);
        if( UNITS_SWITCH == ON ) {
            free((char *) d[2].q->dimen->units_name);
            free((char *) d[1].q->dimen->units_name);
            free((char *) d[0].q->dimen->units_name);
            free((char *) d[2].q->dimen);
            free((char *) d[1].q->dimen);
            free((char *) d[0].q->dimen);
        }
        free((char *) d[2].q);
        free((char *) d[1].q);
        free((char *) d[0].q);
      break;

      case 4:
        d[3]  = Pop();
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)(d[0].q,d[1].q, d[2].q, d[3].q, (QUANTITY *) NULL);
        if( UNITS_SWITCH == ON ) {
            free((char *) d[3].q->dimen->units_name);
            free((char *) d[2].q->dimen->units_name);
            free((char *) d[1].q->dimen->units_name);
            free((char *) d[0].q->dimen->units_name);
            free((char *) d[3].q->dimen);
            free((char *) d[2].q->dimen);
            free((char *) d[1].q->dimen);
            free((char *) d[0].q->dimen);
        }
        free((char *) d[3].q);
        free((char *) d[2].q);
        free((char *) d[1].q);
        free((char *) d[0].q);
      break;

      case 5:
        d[4]  = Pop();
        d[3]  = Pop();
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)(d[0].q, d[1].q,d[2].q, d[3].q, d[4].q, (QUANTITY *) NULL);
        if( UNITS_SWITCH == ON ) {
            free((char *) d[4].q->dimen->units_name);
            free((char *) d[3].q->dimen->units_name);
            free((char *) d[2].q->dimen->units_name);
            free((char *) d[1].q->dimen->units_name);
            free((char *) d[0].q->dimen->units_name);
            free((char *) d[4].q->dimen);
            free((char *) d[3].q->dimen);
            free((char *) d[2].q->dimen);
            free((char *) d[1].q->dimen);
            free((char *) d[0].q->dimen);
        }
        free((char *) d[4].q);
        free((char *) d[3].q);
        free((char *) d[2].q);
        free((char *) d[1].q);
        free((char *) d[0].q);
      break;
    }

     pc = pc+1;
     Push(d3);
}

int BltinVar_String()
{
DATUM d1, d[5], d3;
int no_args, i;

    d1 = Pop();
    no_args = (int) d1.sym;

    switch(no_args) {
      case 1:
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)((char *) d[0].sym, (char *) NULL);
        free((char *) d[0].sym);
      break;
      
      case 2:
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)((char *) d[0].sym, (char *) d[1].sym, (char *) NULL);
        free((char *) d[1].sym);
        free((char *) d[0].sym);
      break;

      case 3:
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)((char *) d[0].sym, (char *) d[1].sym,(char *) d[2].sym, (char *) NULL);
        free((char *) d[2].sym);
        free((char *) d[1].sym);
        free((char *) d[0].sym);
      break;

      case 4:
        d[3]  = Pop();
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)((char *) d[0].sym, (char *) d[1].sym, (char *) d[2].sym,
                            (char *) d[3].sym, (char *) NULL);
        free((char *) d[3].sym);
        free((char *) d[2].sym);
        free((char *) d[1].sym);
        free((char *) d[0].sym);
      break;

      case 5:
        d[4]  = Pop();
        d[3]  = Pop();
        d[2]  = Pop();
        d[1]  = Pop();
        d[0]  = Pop();
        d3.q = (*(QUANTITY * (*)()) *pc)((char *) d[0].sym,(char *) d[1].sym, (char*) d[2].sym,
                            (char *) d[3].sym,(char *) d[4].sym, (char *) NULL);
        free((char *) d[4].sym);
        free((char *) d[3].sym);
        free((char *) d[2].sym);
        free((char *) d[1].sym);
        free((char *) d[0].sym);
      break;
    }
     pc = pc+1;
     Push(d3);
}


/* 
 * ------------------------------------
 * Functions for Engineering Dimensions
 * ------------------------------------
 */ 

/* Push_Dimension() pushes a data type Datum onto the 
   interpreter stack. Increment the pointer to prog. */

int Push_Dimension()
{
DATUM    d;
int length;

    d.q        = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
    switch(CheckUnits()) {
      case ON:
         d.q->dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
         UnitsCopy( d.q->dimen, ((SYMBOL *) *pc)->u.q->dimen );
         break;
      case OFF:
         d.q->dimen = (DIMENSIONS *)NULL; 
         break;
      default:          
         FatalError("In Push_Dimension() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
    pc = pc + 1;
    Push(d);
}

int Push_Dimensionless()
{
DATUM d;

    d.q        = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
    if(CheckUnits() == ON ) {
        d.q->dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
        d.q->dimen->units_name   = (char *)NULL;
        d.q->dimen->units_type   = SI_US; /* A NUMBER is default as SI_US for now */
        d.q->dimen->scale_factor = 1.0;
        d.q->dimen->length_expnt = 0.0;
        d.q->dimen->mass_expnt   = 0.0;
        d.q->dimen->time_expnt   = 0.0;
        d.q->dimen->temp_expnt   = 0.0;
        d.q->dimen->radian_expnt = 0.0;
    }
    else
        d.q->dimen = (DIMENSIONS *)NULL;

    Push(d);
}

int Dimension_Eval()
{
DATUM    d1, d2;
int      length;
int   UnitsType;

    d1 = Pop();   /* quantity number, d1.q */
    d2 = Pop();   /* dimensions     , d2.q */

    switch(CheckUnits()) {
      case ON:
         UnitsType = CheckUnitsType();
         if( (d2.q->dimen->units_type != UnitsType) &&
             (d2.q->dimen->units_type != SI_US    ) &&
             (!strcmp(d2.q->dimen->units_name, "deg_C") ||
              !strcmp(d2.q->dimen->units_name, "deg_F")) ) {

            d2.q->value = d1.q->value*d2.q->dimen->scale_factor;
            d2.q->value = ConvertTempUnits(d2.q->dimen->units_name, d2.q->value, UnitsType);
            d1.q->value = d2.q->value;

         } else
            d1.q->value = d1.q->value*d2.q->dimen->scale_factor;

         UnitsSimplify( d2.q->dimen );
         UnitsCopy( d1.q->dimen, d2.q->dimen );

         free((char *) d2.q->dimen->units_name);
         free((char *) d2.q->dimen);
         free((char *) d2.q);
         Push(d1);
         break;
      case OFF:
         free((char *) d2.q);
         Push(d1);
         break;
      default:          
         FatalError("In Dimension_Eval() : CheckUnits is not ON or OFF",(char *)NULL);
         break;
    }
}

int Dimension_Mult()
{
DATUM d1, d2, d3;
int       length;

   if(CheckUnits() == OFF)
      FatalError("You should set units on to use this function",
                 "Fail to set units on in Dimension_Mult()",(char *)NULL);

   d2 = Pop();
   d1 = Pop();
   d3.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
   d3.q->value = d1.q->value * d2.q->value;
   d3.q->dimen = UnitsMult( d1.q->dimen, d2.q->dimen );
   free((char *) d1.q->dimen->units_name);
   free((char *) d1.q->dimen);
   free((char *) d1.q);
   free((char *) d2.q->dimen->units_name);
   free((char *) d2.q->dimen);
   free((char *) d2.q);
   Push(d3);
}

int Dimension_Div()
{
DATUM d1, d2, d3;
int       length;

   if(CheckUnits() == OFF)
      FatalError("You should set units on to use this function",
                 "Fail to set units on Dimension_Div()",(char *)NULL);

   d2 = Pop();
   d1 = Pop();
   d3.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
   d3.q->value = d1.q->value/d2.q->value;
   d3.q->dimen = UnitsDiv( d1.q->dimen, d2.q->dimen, NO );
   free((char *) d1.q->dimen->units_name);
   free((char *) d1.q->dimen);
   free((char *) d1.q);
   free((char *) d2.q->dimen->units_name);
   free((char *) d2.q->dimen);
   free((char *) d2.q);
   Push(d3);
}

int Dimension_Power()
{
DATUM d1, d2, d3;
double      temp;
int       length;

   if(CheckUnits() == OFF)
      FatalError("You should set units on to use this function",
                 "Fail to set units on in Dimension_Power()",(char *)NULL);

   d2 = Pop();  /* power value */
   d1 = Pop();  /* dimension */
   d3.q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
   d3.q->value = pow(d1.q->value, d2.q->value);
   d3.q->dimen = UnitsPower( d1.q->dimen, d2.q->value, NO );
   free((char *) d1.q->dimen->units_name);
   free((char *) d1.q->dimen);
   free((char *) d1.q);
   free((char *) d2.q->dimen->units_name);
   free((char *) d2.q->dimen);
   free((char *) d2.q);
   Push(d3);
}

/* 
 *  ------------------------
 *  Builtin Matrix Functions 
 *  ------------------------
 */ 

int Bltin_Matrix_Add()
{
DATUM d1, d2, d3;

   d2 = Pop();
   d1 = Pop();
   d1.m = MatrixAddReplace(d1.m, d2.m);
   MatrixFree(d2.m);
   Push(d1);
}

int Bltin_Matrix_Sub()
{
DATUM d1, d2, d3;

   d2 = Pop();
   d1 = Pop();
   d1.m = MatrixSubReplace(d1.m, d2.m);
   MatrixFree(d2.m);
   Push(d1);
}

int Bltin_Matrix_Mult()
{
DATUM d1, d2, d3;

   d2 = Pop();
   d1 = Pop();
   d3.m = MatrixMult(d1.m, d2.m);
   MatrixFree( d1.m );
   MatrixFree( d2.m );
   Push(d3);
}

int Bltin_Matrix_Quan_Mult()
{
DATUM d1, d2, d3;

   d1 = Pop();      /* quantity */
   d2 = Pop();      /* matrix   */
   d3.m = MatrixQuanMult(d1.q,d2.m);
   MatrixFree(d2.m);
   if( CheckUnits() == ON ) {
     free((char *) d1.q->dimen->units_name);
     free((char *) d1.q->dimen);
   }
   free((char *) d1.q);
   Push(d3);
}

int Bltin_Quan_Matrix_Mult() 
{
DATUM d1, d2, d3;

   d2 = Pop();      /* matrix   */
   d1 = Pop();      /* quantity */
   d3.m = MatrixQuanMult(d1.q,d2.m);
   MatrixFree(d2.m);
   if( CheckUnits() == ON ) {
     free((char *) d1.q->dimen->units_name);
     free((char *) d1.q->dimen);
   }
   free((char *) d1.q);
   Push(d3);
}

int Bltin_Matrix_Quan_Div()
{
DATUM d1, d2, d3;

   d1 = Pop();      /* quantity */
   d2 = Pop();      /* matrix   */
   d3.m = MatrixQuanDiv(d2.m, d1.q);
   MatrixFree(d2.m);
   if( CheckUnits() == ON ) {
     free((char *) d1.q->dimen->units_name);
     free((char *) d1.q->dimen);
   }
   free((char *) d1.q);
   Push(d3);
}

int Bltin_Matrix_Power()
{
DATUM d1, d2, d3;

   d1 = Pop();      /* quantity */
   d2 = Pop();      /* matrix   */
   d3.m = MatrixPower(d2.m, d1.q);
   MatrixFree(d2.m);
   if( CheckUnits() == ON ) {
     free((char *) d1.q->dimen->units_name);
     free((char *) d1.q->dimen);
   }
   free((char *) d1.q);
   Push(d3);
}

int Bltin_Matrix_Trans()
{
DATUM d1, d2, d3;

   d2   = Pop();      /* matrix  */
   d3.m = MatrixTranspose(d2.m);
   MatrixFree(d2.m);
   Push(d3);
}

int Bltin_Matrix_Negate()
{
DATUM d1;

   d1 = Pop();
   d1.m = MatrixNegateReplace(d1.m);
   Push(d1);
}

int Bltin_Matrix_Affirm()
{
   Push(Pop());
}

/* 
 * --------------------------------
 * Builtin Finite Element Functions
 * --------------------------------
*/ 

/* func(quantity, matrix) , such as AddNode, FixNode, NodeLoad */
int Bltin_Node_Quant()
{
DATUM d1, d2;

     d1 = Pop();        /* matrix            , d1.m */ 
     d2 = Pop();        /* quantity : node No, d2.q */
     (*(int (*)())*pc)(d2.q->value, d1.m);

     if(CheckUnits()==ON) {
        free((char *)d2.q->dimen->units_name);
        free((char *)d2.q->dimen);
     }
     free((char *)d2.q);
     MatrixFree(d1.m);
     pc = pc + 1;
     Push(d2);
}

int Bltin_Link_Node() 
{
DATUM d1, d2, d3;

     d2 = Pop();  /* matrix of nodes which will be linked together */
     d1 = Pop();  /* matrix of restriction condition of each dof   */
     Link_Node(d1.m, d2.m);
     MatrixFree( d1.m );
     MatrixFree( d2.m );
     pc = pc+1;
     Push(d1);
}

int Bltin_Add_Elmt()
{
DATUM d1, d2, d3;

     d1 = Pop();       /* name string of elmt attribute, d1.sym */
     d2 = Pop();       /* connection matrix of elmt: nodes connected to matrix, d2.m */
     d3 = Pop();       /* elmt No., d3.q */
     Add_Elmt(d3.q->value, d2.m, (char *) d1.sym);

     if(CheckUnits()==ON) {
        free((char *)d3.q->dimen->units_name);
        free((char *)d3.q->dimen);
     }
     free((char *)d3.q);
     MatrixFree(d2.m);
     free((char *) d1.sym);
     pc = pc+1;
     Push(d3);
}

/* 
 *  --------------------------------------------
 *  Builtin Functions for FE Solution Procedures
 *  --------------------------------------------
 */ 

int Bltin_Mesh()
{
DATUM d1;

   (*(void (*)())*pc)();
   pc = pc + 1;
   Push(d1);
}

int Bltin_Fe_Function()
{
DATUM d1;

   (*(void (*)())*pc)();
   pc = pc + 1;
   Push(d1);
}

int Bltin1_Fe_Function()
{
DATUM d1, d;

   d1 = Pop();
   (*(void (*)())*pc)( d1.m );
   MatrixFree(d1.m);
   pc = pc + 1;
   Push(d);
}

/* 
 *  ------------------------------------------------------
 *  Setup Builtin Element, Section and Material Attributes 
 *  ------------------------------------------------------
 */ 

int Bltin_Element_Attr()
{
DATUM d1, d2, d3, d4, d5; 
ELEMENT_ATTR   *eap;
SYMBOL         *hp;
int  i, length;

     d1 = Pop();     /* element attribute name string , (char *)d1.sym  */
     d2 = Pop();     /* number of iterms in elmt attribute , (int)d2.sym*/  
     length = (int) d2.sym;

     eap = Alloc_Element_Attr_Item();
     eap->map_ldof_to_gdof[0] = (int) 1;
     eap->map_ldof_to_gdof[1] = (int) 2;
     eap->map_ldof_to_gdof[2] = (int) 3;
     eap->map_ldof_to_gdof[3] = (int) 4;
     eap->map_ldof_to_gdof[4] = (int) 5;
     eap->map_ldof_to_gdof[5] = (int) 6;

     for(i = 1; i <= length; i++) {
        d3 = Pop();   /* varible : type/ material /section   , d3.sym        */
        d4 = Pop();   /* names strings: elmt_type/section name/material name,(char *)d4.sym */ 

        if(d3.sym->type ==  MATERIAL ||
           d3.sym->type ==  TYPE     ||
           d3.sym->type ==  SECTION  ||
	   d3.sym->type ==  FIBER    ) {
	   hp = lookup(d3.sym->cpSymName);
	   switch(hp->type) {
               case TYPE:
	            eap->elmt_type = SaveString((char *) d4.sym);
	            break;
               case SECTION:
	            eap->section = SaveString((char *) d4.sym);
	            break;
               case MATERIAL:
	            eap->material = SaveString((char *) d4.sym);
	            break;
               case FIBER:
	            eap->fiber_attr_name = SaveString((char *) d4.sym);
	            break;
               default:
                    break;
	   }
	   free((char *) d4.sym);
        }
        else {
          if(d3.m->iNoRows    != d4.m->iNoRows ||
             d3.m->iNoColumns != d4.m->iNoColumns) {
             printf("local dof matrix and global dof matrix should have same dimensions\n");
             FatalError("In input file: ElementAttr{} ",(char *)NULL);
          }
          else {
             Ldof_to_gdof(eap, d4.m, d3.m);
             PRINT_MAP_DOF = ON;
	     MatrixFree( d3.m );
	     MatrixFree( d4.m );
	  }
        }
     }
     Add_Element_Attr((char *) d1.sym, eap);
     free((char *) d1.sym);
}

int Bltin_Section_Attr()
{
DATUM d1, d2, d3, d4; 
SECTION_ATTR     *sp;
SECTION_ATTR     *sp_temp;
SYMBOL           *hp;
int  i, length;
int  UNITS_SWITCH;

     d1 = Pop(); /* section attribute name, (char *)d1.sym */
     d2 = Pop(); /* no of section property items, (int)d2.sym */
     length = (int) d2.sym;
     UNITS_SWITCH = CheckUnits();

     sp = Alloc_Section_Attr_Item();
     for(i = 1; i <= length; i++) {
        d3 = Pop();   /* variable name, d3.sym */
        d4 = Pop();   /* section property, d4.q */

        /* Check_TYPE */
        if(d3.sym->type != VAR && d3.sym->type != QUAN) 
             FatalError("Syntax error in inputfile: SectionAtr{}",
             " Only Quantities can be put in section attributes",(char *)NULL);

        if(!strcmp(d3.sym->cpSymName , "Ixx")) {
           sp->Ixx.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->Ixx.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
	      UnitsCopy( sp->Ixx.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "Iyy")) {
           sp->Iyy.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->Iyy.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->Iyy.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "Izz")) {
           sp->Izz.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->Izz.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
	      UnitsCopy( sp->Izz.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"Ixz") || !strcmp(d3.sym->cpSymName,"Izx") ){
           sp->Ixz.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->Ixz.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->Ixz.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"Ixy") || !strcmp(d3.sym->cpSymName,"Iyx") ){
           sp->Ixy.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->Ixy.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->Ixy.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"Iyz") || !strcmp(d3.sym->cpSymName,"Izy") ){
           sp->Iyz.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->Iyz.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->Iyz.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "unit_weight")) {
           sp->weight.value  = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->weight.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->weight.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "bf")) {
           sp->bf.value  = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->bf.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->bf.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "tf")) {
           sp->tf.value  = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->tf.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->tf.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "depth")) {
           sp->depth.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->depth.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->depth.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "area")) {
           sp->area.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->area.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->area.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "thickness")) {
           sp->plate_thickness.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->plate_thickness.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->plate_thickness.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"J") ){
           sp->tor_const.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->tor_const.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->tor_const.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"rT") ){
           sp->rT.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->rT.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->rT.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName , "width")) {
           sp->width.value  = d4.q->value;
           sp->bf.value  = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->width.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->width.dimen,  d4.q->dimen );
              sp->bf.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->bf.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"tw") ){
           sp->tw.value = d4.q->value;
           if(UNITS_SWITCH == ON) {
              sp->tw.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy( sp->tw.dimen,  d4.q->dimen );
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
        else if(!strcmp(d3.sym->cpSymName,"shear_factor")) {
           sp->ks = d4.q->value;
	   if(UNITS_SWITCH == ON) {
              free((char *) d4.q->dimen->units_name);
              free((char *) d4.q->dimen);
           }
        }
	free((char *) d4.q);
     }
     Add_Section_Attr((char *) d1.sym, sp);
     free((char *) d1.sym);
}

int Bltin_Material_Attr()
{
DATUM d1, d2, d3, d4; 
MATERIAL_ATTR    *mp;
SECTION_ATTR     *sp_temp;
MATERIAL_ATTR    *mp_temp;
SYMBOL           *hp;
int  i, length;
int  UNITS_SWITCH;

     d1 = Pop();
     d2 = Pop();
     length = (int) d2.sym;
     UNITS_SWITCH = CheckUnits();

     mp = Alloc_Material_Attr_Item();
     for(i = 1; i <= length; i++) {
        d3 = Pop();    /* name ( stored in hash table ) */
        d4 = Pop();   /* scale factor and dimensions */

        switch(d3.sym->type) {
          case TYPE:
            mp->LC_ptr->name = SaveString((char *) d4.sym);
	    free((char *) d4.sym);
            break;
          case QUAN:
          case VAR:
            if(!strcmp(d3.sym->cpSymName,"E")) {
               if(d4.q->value <= 0.0)
                  FatalError("Young's Modulus E cannot be less of equal to 0", (char *)NULL);

               mp->E.value = d4.q->value;
               if(UNITS_SWITCH == ON) {
                  mp->E.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->E.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"Et")) {
               mp->ET.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->ET.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->ET.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"G")) {
               mp->G.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->G.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->G.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"Gt")) {
               mp->Gt.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->Gt.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->Gt.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"poisson")) {
               mp->nu = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"density")) {
               mp->density.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->density.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->density.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"yield")) {
               mp->fy.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->fy.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->fy.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"shear_yield")) {
               mp->fv.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->fv.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->fv.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"ultimate")) {
               mp->fu.value = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  mp->fu.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                  UnitsCopy( mp->fu.dimen,  d4.q->dimen );
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"n")) {
               mp->LC_ptr->n = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"alpha")) {
               mp->LC_ptr->alpha = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"beta")) {
               mp->LC_ptr->beta = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"ialph")) {
               mp->LC_ptr->ialph = (int) d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
            else if(!strcmp(d3.sym->cpSymName,"pen")) {
               mp->LC_ptr->pen = d4.q->value;
	       if(UNITS_SWITCH == ON) {
                  free((char *) d4.q->dimen->units_name);
                  free((char *) d4.q->dimen);
               }
            }
	    free((char *) d4.q);
            break;
          default:
            printf("type = d3.sym->type = %d \n", d3.sym->type);
            FatalError("Incorrect material attribute type ");
            break;
        }
     }
     Add_Material_Attr((char *) d1.sym, mp);
     free((char *) d1.sym);
}

int Bltin_Fiber_Attr()
{
DATUM d1, d2, d3, d4, d5; 
SYMBOL           *hp;
int  i, j, k, length;
DIMENSIONS    *dimen;
int     UNITS_SWITCH;

FIBER_ELMT      *fep;
int         no_fiber;
int        *attr_map;
MATRIX         *attr;
int no_material_attr;

     UNITS_SWITCH = CheckUnits();
     d1 = Pop(); /* number of fibers, d1.q, NoFiber */
     d2 = Pop(); /* fiber attribute name, (char *)d2.sym, "FiberName" */
     d3 = Pop(); /* number of fiber property items, (int)d3.sym */
     length = (int) d3.sym;

     no_fiber = d1.q->value;
     if(UNITS_SWITCH == ON) {
        free((char *) d1.q->dimen->units_name);
        free((char *) d1.q->dimen);
     }
     free((char *) d1.q);

     fep = Alloc_Fiber_Elmt_Attr_Item( no_fiber );

     no_material_attr = 0;
     attr = (MATRIX *)NULL;
     attr_map = (int *)NULL;

     for(i = 1; i <= length; i++) {
        d4 = Pop();   /* variable name, d4.sym */
        d5 = Pop();   /* fiber property matrix, d5.m */

        if(!strcmp(d4.sym->cpSymName , "FiberCoordinate")) {
	   if( d5.m->iNoColumns != no_fiber )
              FatalError("in FiberAttr: FiberCoordinate",
	                 "no. of matrix columns != NoFiber",(char *)NULL);
	   for( j=1 ; j <= no_fiber ; j++ ) {
	      fep->fiber[j-1].y.value = d5.m->uMatrix.daa[0][j-1];
	      if( d5.m->iNoRows > 1 )
	         fep->fiber[j-1].z.value = d5.m->uMatrix.daa[1][j-1];
	      else
	         fep->fiber[j-1].z.value = 0.0;

              if(UNITS_SWITCH == ON) {
                 fep->fiber[j-1].y.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                 fep->fiber[j-1].z.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));

		 dimen = UnitsMult( &(d5.m->spRowUnits[0]), &(d5.m->spColUnits[j-1]) );
                 UnitsCopy( fep->fiber[j-1].y.dimen, dimen );
	         if( d5.m->iNoRows > 1 ) {
		    free((char *) dimen->units_name);
		    free((char *)dimen);
		    dimen = UnitsMult( &(d5.m->spRowUnits[1]), &(d5.m->spColUnits[j-1]) );
                    UnitsCopy( fep->fiber[j-1].z.dimen, dimen );
		    free((char *) dimen->units_name);
		    free((char *)dimen);
	         }
	         else {
                    UnitsCopy( fep->fiber[j-1].z.dimen, dimen );
		    free((char *) dimen->units_name);
		    free((char *) dimen);
	         }
              }
           }
	   MatrixFree(d5.m);
        } /* end of install FiberCoordinate matrix into hash table */

        else if(!strcmp(d4.sym->cpSymName , "FiberArea")) {
	   if( d5.m->iNoColumns != no_fiber )
              FatalError("in FiberAttr: FiberArea",
	                 "no. of matrix columns != NoFiber",(char *)NULL);
	   for( j=1 ; j <= no_fiber ; j++ ) {
	      fep->fiber[j-1].area.value = d5.m->uMatrix.daa[0][j-1];

              if(UNITS_SWITCH == ON) {
                 fep->fiber[j-1].area.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
		 dimen = UnitsMult( &(d5.m->spRowUnits[0]), &(d5.m->spColUnits[j-1]) );
                 UnitsCopy( fep->fiber[j-1].area.dimen, dimen );
		 free((char *) dimen->units_name);
		 free((char *) dimen);
              }
           }
	   MatrixFree(d5.m);
        } /* end of install FiberArea matrix into hash table */

        else if(!strcmp(d4.sym->cpSymName , "FiberMaterialMap")) {
	   if( d5.m->iNoColumns != no_fiber )
              FatalError("in FiberAttr: FiberMaterialMap",
	                 "no. of matrix columns != NoFiber",(char *)NULL);

           /* index array to store fiber material map */
           attr_map = (int *)MyCalloc( no_fiber, sizeof(int) );

	   for( j=1 ; j <= no_fiber ; j++ )
	      attr_map[j-1] = (int) d5.m->uMatrix.daa[0][j-1];

	   MatrixFree(d5.m);
        } /* end of temporyly store FiberMaterialMap matrix into index array attr_map */

        else if(!strcmp(d4.sym->cpSymName , "FiberMaterialAttr")) {
	   if( d5.m->iNoRows != 3 && d5.m->iNoRows != 6 )
              FatalError("in FiberAttr: FiberMaterialAttr",
	                 "no. of matrix rows != 3 or 6",
			 "example: [E1, E2; Et1, Et2; fy1, fy2] for FIBER_*D",
			 "example: [E1,E2; Et1,Et2; fy1,fy2; G1,G2; Gt1,Gt2; fv1,fv2] for FIBER_*DS", (char *)NULL);

	   no_material_attr = d5.m->iNoColumns;
	   attr = MatrixCopy( d5.m );

	   MatrixFree(d5.m);
        } /* end of temporyly store FiberMaterialAttr matrix into attribute matrix attr */
     }

     if( attr == (MATRIX *)NULL )
	FatalError("in FiberAttr: Must give FiberMaterialAttr", (char *)NULL);
     if( attr_map == (int *)NULL )
	FatalError("in FiberAttr: Must give FiberMaterialMap", (char *)NULL);

     /* install the fiber attribution into hash table */
     for( i=1 ; i <= no_fiber ; i++ ) {
	j = attr_map[i-1];
	fep->fiber[i-1].Es.value = attr->uMatrix.daa[0][j-1];
	fep->fiber[i-1].Et.value = attr->uMatrix.daa[1][j-1];
	fep->fiber[i-1].fy.value = attr->uMatrix.daa[2][j-1];
	if( attr->iNoRows == 6 ) {
	   fep->fiber[i-1].Gs.value = attr->uMatrix.daa[3][j-1];
	   fep->fiber[i-1].Gt.value = attr->uMatrix.daa[4][j-1];
	   fep->fiber[i-1].fv.value = attr->uMatrix.daa[5][j-1];
        }
        else {
	   fep->fiber[i-1].Gs.value = 0.0;
	   fep->fiber[i-1].Gt.value = 0.0;
	   fep->fiber[i-1].fv.value = 0.0;
        }

        if(UNITS_SWITCH == ON) {
           fep->fiber[i-1].Es.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           fep->fiber[i-1].Et.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           fep->fiber[i-1].fy.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));

           dimen = UnitsMult( &(attr->spRowUnits[0]), &(attr->spColUnits[j-1]) );
           UnitsCopy( fep->fiber[i-1].Es.dimen, dimen );
           free((char *) dimen->units_name);
           free((char *) dimen);
           dimen = UnitsMult( &(attr->spRowUnits[1]), &(attr->spColUnits[j-1]) );
           UnitsCopy( fep->fiber[i-1].Et.dimen, dimen );
           free((char *) dimen->units_name);
           free((char *) dimen);
           dimen = UnitsMult( &(attr->spRowUnits[2]), &(attr->spColUnits[j-1]) );
           UnitsCopy( fep->fiber[i-1].fy.dimen, dimen );
           free((char *) dimen->units_name);
           free((char *) dimen);

	   if( attr->iNoRows == 6 ) {
              fep->fiber[i-1].Gs.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              fep->fiber[i-1].Gt.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              fep->fiber[i-1].fv.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));

              dimen = UnitsMult( &(attr->spRowUnits[3]), &(attr->spColUnits[j-1]) );
              UnitsCopy( fep->fiber[i-1].Gs.dimen, dimen );
              free((char *) dimen->units_name);
              free((char *) dimen);
              dimen = UnitsMult( &(attr->spRowUnits[4]), &(attr->spColUnits[j-1]) );
              UnitsCopy( fep->fiber[i-1].Gt.dimen, dimen );
              free((char *) dimen->units_name);
              free((char *) dimen);
              dimen = UnitsMult( &(attr->spRowUnits[5]), &(attr->spColUnits[j-1]) );
              UnitsCopy( fep->fiber[i-1].fv.dimen, dimen );
              free((char *) dimen->units_name);
              free((char *) dimen);
           }
	   else {
              fep->fiber[i-1].Gs.dimen  = (DIMENSIONS *) NULL;
              fep->fiber[i-1].Gt.dimen  = (DIMENSIONS *) NULL;
              fep->fiber[i-1].fv.dimen  = (DIMENSIONS *) NULL;
           }
        }
     }
     free( (char *) attr_map );
     MatrixFree( attr );

     Add_Fiber_Elmt_Attr((char *) d2.sym, fep);
     free((char *) d2.sym);
}

int Bltin_Units_Type()
{
DATUM d1;
char  *type;

     d1 = Pop();     /* units type string , (char *)d1.sym  */
     type = (char *)d1.sym;

     if( !strcmp(type,"US") )
         ChangeUnitsType( 1 );
     else  if( !strcmp(type,"SI") )
         ChangeUnitsType( 2 );
     else  if( !strcmp(type,"SI_US") )
         ChangeUnitsType( 3 );
     else
         FatalError("The Units_Type string in input file, SetUnitsType(\"Units_Type\")",
         "In Bltin_Units_Type() : Must give the type of SI, US, or SI_US",(char *)NULL);

     free((char *)d1.sym);
}
