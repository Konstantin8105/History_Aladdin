/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  matrix.h : Data Structures and Function Declarations for Matrix
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
 *  Written by: Mark Austin and Wane-Jang Lin                       December 1995
 *  ============================================================================= 
 */


#ifndef MATRIX_H
#define MATRIX_H

/* [a] : Data Structures for Matrices of numbers/engineering quantities */

typedef struct {
        double dReal, dImaginary;
} COMPLEX;

typedef enum {
	INTEGER_ARRAY  = 1,
	DOUBLE_ARRAY   = 2,
	COMPLEX_ARRAY  = 3
} DATA_TYPE;

typedef enum {
	SEQUENTIAL = 1,
	INDIRECT   = 2,
	SKYLINE    = 3,
	SPARSE     = 4
} INTERNAL_REP;

typedef struct matrix {
        char      *cpMatrixName;    /*  *name   */
        int             iNoRows;    /*  no_rows   */
        int          iNoColumns;    /*  no_columns   */
	DIMENSIONS  *spRowUnits;    /*  *row_units_buf  */
	DIMENSIONS  *spColUnits;    /*  *col_units_buf  */
        INTERNAL_REP       eRep;
        DATA_TYPE         eType;    /*  type  */
        union {
            int           **iaa;
            double        **daa;    /*  **d  */
            COMPLEX       **caa;
        } uMatrix;
} MATRIX;

/* [b] : Declarations for Root Matrix Functions */

#ifdef __STDC__ 

MATRIX *MatrixAllocate( MATRIX * );
MATRIX *MatrixDiag( MATRIX * );
MATRIX *MatrixZero( MATRIX * );
MATRIX *MatrixOne( MATRIX * );
MATRIX *MatrixAdd( MATRIX * , MATRIX * );
MATRIX *MatrixAddReplace( MATRIX *, MATRIX * );
MATRIX *MatrixSub( MATRIX * , MATRIX * );
MATRIX *MatrixSubReplace( MATRIX *, MATRIX * );
MATRIX *MatrixMult( MATRIX * , MATRIX * );
MATRIX *MatrixPower( MATRIX * , QUANTITY * );
MATRIX *MatrixNegate( MATRIX * );
MATRIX *MatrixNegateReplace( MATRIX *);
MATRIX *MatrixTranspose( MATRIX * );
MATRIX *MatrixCopy( MATRIX * );
MATRIX *MatrixSolve( MATRIX *, MATRIX * );
MATRIX *MatrixLU( MATRIX * );
MATRIX *MatrixFB( MATRIX *, MATRIX * );
MATRIX *MatrixInverse( MATRIX * );
MATRIX *MatrixDimension( MATRIX * );
MATRIX *MatrixScale( MATRIX *, double );
double  MatrixContentScale( MATRIX *, int, int );
MATRIX *MatrixQuanMult( QUANTITY *, MATRIX * );
MATRIX *MatrixQuanDiv( MATRIX * , QUANTITY * );
MATRIX *MatrixZeroUnits( MATRIX *, int, int );
MATRIX *MatrixUnitsLess( MATRIX * );
MATRIX *MatrixUnitsSimplify( MATRIX * );

MATRIX *MatrixPrintVar( MATRIX *, ... );
MATRIX *MatrixPrintCast( MATRIX *, ... );
MATRIX *MatrixColumnUnits( MATRIX *, ... );
MATRIX *MatrixRowUnits( MATRIX *, ... );
MATRIX *MatrixExtract( MATRIX *, ... );
MATRIX *MatrixPut( MATRIX *, ... );

void      MatrixFree( MATRIX * );
QUANTITY *MatrixDet( MATRIX * );
QUANTITY *MatrixMax( MATRIX * );
QUANTITY *MatrixMin( MATRIX * );
QUANTITY *MatrixL2Norm( MATRIX * );
QUANTITY *QuantityCast( MATRIX * );

/* [c] : Matrix Functions with INDIRECT Storage Pattern */

MATRIX   *MatrixAllocIndirect( char *, DATA_TYPE , int , int );
double  **MatrixAllocIndirectDouble( int, int );
int     **MatrixAllocIndirectInteger( int, int );
void      MatrixFreeIndirectDouble( double **, int );
void      MatrixFreeIndirectInteger( int **, int );

void    MatrixPrintIndirectDouble( MATRIX * );
void    MatrixPrintIndirectInteger( MATRIX * );

MATRIX *MatrixCopyIndirectDouble( MATRIX * );
MATRIX *MatrixAddIndirectDouble( MATRIX *, MATRIX *);
MATRIX *MatrixAddReplaceIndirectDouble( MATRIX *, MATRIX *);
MATRIX *MatrixSubIndirectDouble( MATRIX *, MATRIX *);
MATRIX *MatrixSubReplaceIndirectDouble( MATRIX *, MATRIX *);
MATRIX *MatrixMultIndirectDouble( MATRIX *, MATRIX *);
MATRIX *MatrixMultIndirectSkylineDouble( MATRIX *, MATRIX *);
MATRIX *MatrixMultSkylineIndirectDouble( MATRIX *, MATRIX *);
MATRIX *MatrixNegateIndirectDouble( MATRIX * );
MATRIX *MatrixNegateReplaceIndirectDouble( MATRIX * );

MATRIX *MatrixTransposeIndirectDouble( MATRIX * );
MATRIX *MatrixInverseIndirectDouble( MATRIX * );
MATRIX *MatrixScaleIndirectDouble( MATRIX *, double );
double  MatrixContentScaleIndirectDouble( MATRIX *, int, int );

/* [d] : Matrix Functions with SKYLINE Storage Pattern */

MATRIX *MatrixAllocSkyline( char *, DATA_TYPE, int, int, int *);
void    MatrixFreeSkyline( MATRIX * );
void    MatrixPrintSkylineDouble( MATRIX * );

MATRIX *MatrixReallocSkyline( MATRIX * );
MATRIX *MatrixReallocSkylineDouble( MATRIX * );

MATRIX *MatrixAddSkyline( MATRIX *, MATRIX * );
MATRIX *MatrixSubSkyline( MATRIX *, MATRIX * );
MATRIX *MatrixMultSkyline( MATRIX *, MATRIX *);
MATRIX *MatrixNegateSkyline( MATRIX *);
MATRIX *MatrixNegateReplaceSkyline( MATRIX * );
MATRIX *MatrixCopySkyline( MATRIX * );
MATRIX *MatrixTransposeSkyline( MATRIX * );
MATRIX *MatrixInverseSkyline( MATRIX * );
MATRIX *MatrixScaleSkyline( MATRIX *, double );
double  MatrixContentScaleSkyline( MATRIX *, int, int );

MATRIX *LUDecompositionSkyline( MATRIX *);
MATRIX *LUBacksubstitutionSkyline( MATRIX *, MATRIX *);

MATRIX *MatrixIndirectToSkyline( MATRIX * );
MATRIX *MatrixSkylineToIndirect( MATRIX * );

MATRIX *MatrixAssembleSkyline( MATRIX *, MATRIX *, int *, int * );
MATRIX *CholeskyDecompositionIndirect( MATRIX * );

void    MatrixSolveEigen( MATRIX *, MATRIX *, MATRIX *, MATRIX *, int );
MATRIX *Solve_Eigen( MATRIX *, MATRIX *, MATRIX * );
MATRIX *Extract_Eigenvalue( MATRIX * );
MATRIX *Extract_Eigenvector( MATRIX * );
void    Print_Eigen( MATRIX * );

void         dMatrixPrint( char *, double **, int, int );
double     **dMatrixCopy( double **, int, int );
double     **dMatrixCopyRep( double **, double **, int, int );
double     **dVmatrixCrossProduct( double **, double **, int, int, double **, int, int );
double       dVmatrixInnerProduct( double **, int, int, double **, int, int );
double     **dMatrixMult( double **, int, int, double **, int, int );
double     **dMatrixMultRep( double **, double **, int, int, double **, int, int );
double     **dMatrixTranspose( double **, int, int );
double       dVmatrixL2Norm( double **, int, int );
double       dMatrixDet( double **, int, int );

#else  /* start case not STDC */

/* [b] : Declarations for Root Matrix Functions */

MATRIX *MatrixPrintCast();
MATRIX *MatrixPrintVar();
MATRIX *MatrixAllocate();
MATRIX *MatrixDiag();
MATRIX *MatrixZero();
MATRIX *MatrixOne();
MATRIX *MatrixScale();
double  MatrixContentScale();
MATRIX *MatrixCopy();
MATRIX *MatrixTranspose();
MATRIX *MatrixDimension();
MATRIX *MatrixAdd();
MATRIX *MatrixAddReplace();
MATRIX *MatrixSub();
MATRIX *MatrixSubReplace();
MATRIX *MatrixMult();
MATRIX *MatrixPower();
MATRIX *MatrixNegate();
MATRIX *MatrixNegateReplace();
void    MatrixFree();
MATRIX *MatrixSolve();
MATRIX *MatrixLU();
MATRIX *MatrixFB();
void    MatrixSolveEigen();
MATRIX *MatrixInverse();
QUANTITY  *MatrixDet();
QUANTITY  *MatrixL2Norm();
QUANTITY  *MatrixMax();
QUANTITY  *MatrixMin();

/* [b.1] : Operations between MATRIX and QUANTITY */

MATRIX *MatrixQuanMult();
MATRIX *MatrixQuanDiv();

/* [b.2] : Declarations for Matrix Functions about Units */

MATRIX      *MatrixColumnUnits();
MATRIX      *MatrixRowUnits();
MATRIX      *MatrixZeroUnits();
MATRIX      *MatrixUnitsSimplify();
MATRIX      *MatrixUnitsLess();
MATRIX      *MatrixExtract();
MATRIX      *MatrixPut();
QUANTITY    *QuantityCast();


/* [c] : Matrix Functions with INDIRECT Storage Pattern */

MATRIX  *MatrixAllocIndirect();
double **MatrixAllocIndirectDouble();
int    **MatrixAllocIndirectInteger();

void     MatrixFreeIndirectDouble();
void     MatrixFreeIndirectInteger();

MATRIX *MatrixCopyIndirectDouble();
MATRIX *MatrixScaleIndirectDouble();
double  MatrixContentScaleIndirectDouble();

MATRIX *MatrixAddIndirectDouble();
MATRIX *MatrixAddReplaceIndirectDouble();
MATRIX *MatrixSubIndirectDouble();
MATRIX *MatrixSubReplaceIndirectDouble();
MATRIX *MatrixNegateIndirectDouble();
MATRIX *MatrixNegateReplaceIndirectDouble();

MATRIX *MatrixMultIndirectDouble();
MATRIX *MatrixTransposeIndirectDouble();

MATRIX *LUDecompositionIndirect();
MATRIX *LUSubstitutionIndirect();

MATRIX *MatrixInverseIndirectDouble();

/* [d] : Matrix Functions with SKYLINE Storage Pattern */

MATRIX *MatrixAllocSkyline();
void    MatrixPrintSkylineDouble();
void    MatrixFreeSkyline();

MATRIX *ReSkyMatrix();
MATRIX *ReSkyMatrixDouble();

MATRIX *MatrixAddSkyline();
MATRIX *MatrixSubSkyline();
MATRIX *MatrixMultSkyline();
MATRIX *MatrixNegateSkyline();
MATRIX *MatrixNegateReplaceSkyline();

MATRIX *MatrixMultIndirectSkylineDouble();
MATRIX *MatrixMultSkylineIndirectDouble();

MATRIX *MatrixCopySkyline();
MATRIX *MatrixScaleSkyline();
double  MatrixContentScaleSkyline();

MATRIX *LUDecompositionSkyline();
MATRIX *LUBacksubstitutionSkyline();

MATRIX *MatrixTransposeSkyline();
MATRIX *MatrixInverseSkyline();

MATRIX *MatrixIndirectToSkyline();
MATRIX *MatrixSkylineToIndirect();

MATRIX *MatrixAssembleSkyline();
MATRIX *CholeskyDecompositionIndirect();

MATRIX *Solve_Eigen();
MATRIX *Extract_Eigenvalue();
MATRIX *Extract_Eigenvector();
void    Print_Eigen();

/* [e] : Declarations for Double Matrix (without units) Functions */

void         dMatrixPrint();
double     **dMatrixCopy();
double     **dMatrixCopyRep();
double     **dVmatrixCrossProduct();
double       dVmatrixInnerProduct();
double     **dMatrixMult();
double     **dMatrixMultRep();
double     **dMatrixTranspose();
double       dVmatrixL2Norm();
double       dMatrixDet();

#endif /* end case not STDC */

#endif /* end case MATRIX_H */
