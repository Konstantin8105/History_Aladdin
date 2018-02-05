/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  vector.h : Data structures for vector module                                  
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
 *  Written by: Mark Austin                                             July 1993
 *  ============================================================================= 
 */

#ifndef VECTOR_H
#define VECTOR_H

/* Vector Data Structure */

typedef struct {
        char *cpVectorName;    /*  *name  */
        int        iLength;    /*  length  */
        DATA_TYPE    eType;    /*  type  */
        union {
            int    *ia;        /*  *i  */
            double *da;        /*  *d  */
        } uVector;             /*  array  */
} VECTOR;

#if (__STDC__ == 1)

/* Function Declarations for Standard ANSI C */

VECTOR *VectorAlloc( char * , DATA_TYPE , int );
VECTOR *VectorAdd( VECTOR * , VECTOR * );
VECTOR *VectorSub( VECTOR * , VECTOR * );
void    VectorPrint( VECTOR * );
void    VectorFree( VECTOR * );

void    PrintVectorInteger( VECTOR * );
void    VectorFreeInteger( VECTOR * );
VECTOR *VectorAddInteger( VECTOR * , VECTOR * );
VECTOR *VectorSubInteger( VECTOR * , VECTOR * );
int    *iVectorAlloc( int );

void    PrintVectorDouble( VECTOR * );
void    VectorFreeDouble( VECTOR * );
VECTOR *VectorAddDouble( VECTOR *, VECTOR *);
VECTOR *VectorSubDouble( VECTOR *, VECTOR *);
double *dVectorAlloc( int );

VECTOR *NaiveGaussElimination( MATRIX *, VECTOR *);
VECTOR *GaussElimination( char *, MATRIX *, VECTOR *);
VECTOR *SetupScaleFactors( MATRIX * );
VECTOR *SetupPivotVector( MATRIX * );

MATRIX *LUDecompositionIndirect( MATRIX *, VECTOR *);
MATRIX *LUSubstitutionIndirect( char *, VECTOR *, MATRIX *, MATRIX *);

#else  /* Start case not STDC */

/* Function Declarations for K&R C */

VECTOR *VectorAlloc();
VECTOR *VectorAdd();
VECTOR *VectorSub();
void    VectorPrint();
void    VectorFree();

void    PrintVectorInteger();
void    VectorFreeInteger();
VECTOR *VectorAddInteger();
VECTOR *VectorSubInteger();
int    *iVectorAlloc();

void    PrintVectorDouble();
void    VectorFreeDouble();
VECTOR *VectorAddDouble();
VECTOR *VectorSubDouble();
double *dVectorAlloc();

VECTOR *NaiveGaussElimination();
VECTOR *GaussElimination();
VECTOR *SetupScaleFactors();
VECTOR *SetupPivotVector();

MATRIX *LUDecompositionIndirect();
MATRIX  *LUSubstitutionIndirect();

#endif /* End case not STDC */

#endif /* end case VECTOR_H */
