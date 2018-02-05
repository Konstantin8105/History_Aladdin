/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  matrix_indirect.c : Functions for Matrices having INDIRECT storage pattern.
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
 *  Written by: Mark Austin                                             1992-1993
 *  ============================================================================= 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"

/* #define DEBUG */
/* #define LU_DEBUG */


/*
 *  =======================================================================
 *  MatrixAllocIndirect() : Allocate memory for Matrix data structure with
 *  INDIRECT storage pattern
 *  
 *  Input  :  char *cpMatrixName -- Pointer to name of matrix.
 *         :  DATA_TYPE eType    -- Data type to be stored in Matrix.
 *         :  int  iNoRows       -- No of Rows in Matrix.
 *         :  int  iNoColumns    -- No of Columns in Matrix.
 *  Output :  MATRIX *Matrix     -- Pointer to matrix data structure.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixAllocIndirect( char *cpMatrixName, DATA_TYPE  eType,
                             int  iNoRows, int   iNoColumns)
#else  /* Start case not STDC */
MATRIX *MatrixAllocIndirect( cpMatrixName, eType, iNoRows, iNoColumns)
char     *cpMatrixName;
DATA_TYPE        eType;
int            iNoRows;
int         iNoColumns;
#endif /* End case not STDC */
{
MATRIX *spMatrix;

   /* [a] : Allocate matrix parent data structure, and cpMatrixName */

      spMatrix = (MATRIX *) MyMalloc( sizeof(MATRIX) );
      if(cpMatrixName != (char *)NULL)
         spMatrix->cpMatrixName = SaveString(cpMatrixName);
      else 
         spMatrix->cpMatrixName = (char *) NULL;

   /* [b] : Set parameters and allocate memory for matrix uMatrix */

      spMatrix->eRep       = INDIRECT;
      spMatrix->eType      = eType;
      spMatrix->iNoRows    = iNoRows;
      spMatrix->iNoColumns = iNoColumns;
      if( CheckUnits() == ON ) {
         spMatrix->spRowUnits 
         = (DIMENSIONS *) MyCalloc(iNoRows, sizeof(DIMENSIONS));
         spMatrix->spColUnits
         = (DIMENSIONS *) MyCalloc(iNoColumns, sizeof(DIMENSIONS));
      }
      else {
         spMatrix->spRowUnits = (DIMENSIONS *)NULL;
         spMatrix->spColUnits = (DIMENSIONS *)NULL;
      }

      switch((int) spMatrix->eType) {
          case DOUBLE_ARRAY:
               spMatrix->uMatrix.daa 
               = MatrixAllocIndirectDouble( iNoRows, iNoColumns);
               break;
          case INTEGER_ARRAY:
               spMatrix->uMatrix.iaa 
               = MatrixAllocIndirectInteger( iNoRows, iNoColumns);
               break;
          case COMPLEX_ARRAY:
               FatalError("In MatrixAllocIndirect() : spMatrix->eType not implemented",
                         (char *) NULL);
               break;
          default:
               FatalError("In MatrixAllocIndirect() : Undefined spMatrix->eType",
                         (char *) NULL);
               break;
      }

      return (spMatrix);
}

/*
 *  =======================================================================
 *  MatrixAllocIndirectDouble() : Allocate INDIRECT storage pattern for
 *  DOUBLE data types
 *  
 *  Input  :  int  iNoRows       -- No of Rows in Matrix.
 *         :  int  iNoColumns    -- No of Columns in Matrix.
 *  Output :  double **Matrix    -- Pointer to matrix data structure.
 *  =======================================================================
 */

#ifdef __STDC__
double **MatrixAllocIndirectDouble( int iNoRows, int iNoColumns)
#else  /* Start case not STDC */
double **MatrixAllocIndirectDouble( iNoRows, iNoColumns)
int    iNoRows;
int iNoColumns;
#endif /* End case not STDC */
{
double **Matrix;
int ii;

      Matrix = (double **) MyCalloc( iNoRows, sizeof(double *));
      for(ii = 1; ii <= iNoRows; ii++) 
          Matrix[ii-1] = (double *) MyCalloc( iNoColumns, sizeof(double));

      return (Matrix);
}

/*
 *  =======================================================================
 *  MatrixAllocIndirectInteger() : Allocate INDIRECT storage pattern for
 *                                 INTEGER data types
 *  
 *  Input  :  int  iNoRows       -- No of Rows in Matrix.
 *         :  int  iNoColumns    -- No of Columns in Matrix.
 *  Output :  double **Matrix    -- Pointer to matrix data structure.
 *  =======================================================================
 */

#ifdef __STDC__
int **MatrixAllocIndirectInteger( int iNoRows, int iNoColumns)
#else  /* Start case not STDC */
int **MatrixAllocIndirectInteger( iNoRows, iNoColumns)
int    iNoRows;
int iNoColumns;
#endif /* End case not STDC */
{
int **Matrix;
int ii;

      Matrix = (int **) MyCalloc( iNoRows, sizeof(int *));
      for(ii = 1; ii <= iNoRows; ii++) 
          Matrix[ii-1] = (int *) MyCalloc( iNoColumns, sizeof(int));

      return (Matrix);
}


/*
 *  =======================================================================
 *  MatrixPrintIndirectDouble() : Print a Matrix [iNoRows][iNoColumns] of
 *                                data type DOUBLE 
 *  
 *  MatrixPrintIndirectInteger() : Print a Matrix [iNoRows][iNoColumns] of
 *                                 data type INTEGER
 *  
 *  Where --  COLUMNS_ACROSS_PAGE  = Number of columns printed across page.
 *            spA->iNoRows    = Number of rows in matrix.
 *            spA->iNoColumns = Number of columns in matrix.
 *            iFirstColumn         = Number of first column in block
 *            iLastColumn          = Number of last  column in block
 *            ib                   = Current No of Matrix Block.
 * 
 *  Input  :  MATRIX *spA         -- Pointer to matrix data structure.
 *  Output :  void
 *  =======================================================================
 */

enum { COLUMNS_ACROSS_PAGE = 6 };                                        /* Item [a] */

#ifdef __STDC__                                                  
void MatrixPrintIndirectDouble( MATRIX *spA )
#else  /* Start case not STDC */
void MatrixPrintIndirectDouble( spA )
MATRIX  *spA;
#endif /* End case not STDC */
{
int ii, ij, ib;    
int iNoBlocks; 
int iFirstColumn, iLastColumn; 
double da;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Compute no of blocks of rows to be printed */               /* Item [c] */

    if( spA->iNoColumns % ((int) COLUMNS_ACROSS_PAGE) == 0 )                          
        iNoBlocks = (spA->iNoColumns/((int) COLUMNS_ACROSS_PAGE));
    else
        iNoBlocks = (spA->iNoColumns/((int) COLUMNS_ACROSS_PAGE)) + 1;

    /* [b] : Loop over blocks of rows */

    for( ib = 1; ib <= iNoBlocks; ib++ ) {                               /* Item [d] */

         iFirstColumn = (ib-1)*((int) COLUMNS_ACROSS_PAGE) + 1;
         iLastColumn  = (int) MIN( ib*((int) COLUMNS_ACROSS_PAGE) , spA->iNoColumns );

         /* [c] : Print title of matrix at top of each block */          /* Item [e] */
       
         if( spA->cpMatrixName != NULL ) 
             printf("\nMATRIX : \"%s\"\n\n", spA->cpMatrixName);
         else 
             printf("\nMATRIX : \"UNTITLED\"\n\n");

         /* [d] : Label row and column nos */

         printf ("row/col         ");
         for( ii = iFirstColumn; ii <= iLastColumn; ii++ )
             printf("       %3d   ", ii);
         printf("\n");

         switch( UNITS_SWITCH) {
           case ON:
               printf ("        units   " );
               for( ii = iFirstColumn; ii <= iLastColumn; ii++ )
                  if(spA->spColUnits[ii-1].units_name != NULL) {
                     printf("%10s   ",spA->spColUnits[ii-1].units_name);
                  }
                  else
                    printf("             ");
               printf("\n");

             /* [e] : Print Contents of Matrix */   /* Item [f] */

               for( ii = 1; ii <= spA->iNoRows; ii++ ) {
                 printf(" %3d ", ii);
                 if(spA->spRowUnits[ii-1].units_name != NULL) {
                    printf("%8s ",spA->spRowUnits[ii-1].units_name);
                 }
                 else
                    printf("         ");
                 for(ij  = iFirstColumn; ij <= iLastColumn; ij++) {
                     da = MatrixContentScaleIndirectDouble(spA, ii, ij);
                     printf(" %12.5e", da );
                 }
                 printf("\n");
               }
               break;
           case OFF:

               /* [e] : Print Contents of Matrix */                            /* Item [f] */

               for( ii = 1; ii <= spA->iNoRows; ii++ ) {
                 printf(" %3d ", ii);
                 printf("         ");
                 for( ij  = iFirstColumn; ij <= iLastColumn; ij++)
                     printf(" %12.5e", spA->uMatrix.daa[ ii-1 ][ ij-1 ]);
                 printf("\n");
               }
               break;
           default:
               break;
         }
    }
}

#ifdef __STDC__                                                  
void MatrixPrintIndirectInteger( MATRIX *spA )
#else  /* Start case not STDC */
void MatrixPrintIndirectInteger( spA )
MATRIX  *spA;
#endif /* End case not STDC */
{
int ii, ij, ib;    
int iNoBlocks; 
int iFirstColumn, iLastColumn; 
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Compute no of blocks of rows to be printed */               /* Item [c] */

    if( spA->iNoColumns % ((int) COLUMNS_ACROSS_PAGE) == 0 )                          
        iNoBlocks = (spA->iNoColumns/((int) COLUMNS_ACROSS_PAGE));
    else
        iNoBlocks = (spA->iNoColumns/((int) COLUMNS_ACROSS_PAGE)) + 1;

    /* [b] : Loop over blocks of rows */

    for( ib = 1; ib <= iNoBlocks; ib++ ) {                               /* Item [d] */

         iFirstColumn = (ib-1)*((int) COLUMNS_ACROSS_PAGE) + 1;
         iLastColumn  = (int) MIN( ib*((int) COLUMNS_ACROSS_PAGE) , spA->iNoColumns );

         /* [c] : Print title of matrix at top of each block */          /* Item [e] */
       
         if( spA->cpMatrixName != NULL ) 
             (void) printf("\nMATRIX : \"%s\"\n\n", spA->cpMatrixName);
         else 
             (void) printf("\nMATRIX : \"UNTITLED\"\n\n");

         /* [d] : Label row and column nos */

         (void) printf ("row/col        ");
         for( ii = iFirstColumn; ii <= iLastColumn; ii++ )
              (void) printf("%3d          ", ii);
         (void) printf("\n");

         switch( UNITS_SWITCH) {
           case ON:
               printf ("      units   " );
               for( ii = iFirstColumn; ii <= iLastColumn; ii++ )
                  if(spA->spColUnits[ii-1].units_name != NULL) {
                     printf("%10s   ",spA->spColUnits[ii-1].units_name);
                  }
               printf("\n");

               /* [e] : Print Contents of Matrix */                            /* Item [f] */

               for( ii = 1; ii <= spA->iNoRows; ii++ ) {
                    (void) printf(" %3d ", ii);
                    if(spA->spRowUnits[ii-1].units_name != NULL) {
                       printf("%8s ",spA->spRowUnits[ii-1].units_name);
                    }
                    else
                       printf("         ");
                    for( ij  = iFirstColumn; ij <= iLastColumn; ij++)
                         (void) printf(" %12d", spA->uMatrix.iaa[ ii-1 ][ ij-1 ]);
                    (void) printf("\n");
               }
               break;
           case OFF:
               /* [e] : Print Contents of Matrix */                            /* Item [f] */

               for( ii = 1; ii <= spA->iNoRows; ii++ ) {
                    printf(" %3d ", ii);
                    printf("         ");
                    for( ij  = iFirstColumn; ij <= iLastColumn; ij++)
                         (void) printf(" %12d", spA->uMatrix.iaa[ ii-1 ][ ij-1 ]);
                    printf("\n");
               }
               break;
           default:
               break;
         }
    }
}


/*
 *  =======================================================================
 *  MatrixFreeIndirect() : Free memory for INDIRECT storage.
 *  
 *  Input  :  MATRIX *Matrix  -- Pointer to matrix data structure.
 *  Output :  void
 *  =======================================================================
 */

#ifdef __STDC__
void MatrixFreeIndirect( MATRIX *spA )
#else  /* Start case not STDC */
void MatrixFreeIndirect( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
int i;

   if( spA == (MATRIX *)NULL )   return;

   /* [a] : Free body of matrix[iNoRows][iNoColumns] */

      switch(spA->eType) {
          case DOUBLE_ARRAY:
               MatrixFreeIndirectDouble( spA->uMatrix.daa, spA->iNoRows);
               break;
          case INTEGER_ARRAY:
               MatrixFreeIndirectInteger( spA->uMatrix.iaa, spA->iNoRows);
               break;
          case COMPLEX_ARRAY:
               FatalError("In MatrixFreeIndirect() : spA->eType not implemented",
                         (char *) NULL);
               break;
          default:
               break;
      }

   /* (b) Free cpMatrixName and parent data structure for matrix */

      if( CheckUnits() == ON ) {
         for( i=1; i<=spA->iNoRows; i++ )
             free ((char *) spA->spRowUnits[i-1].units_name );
         for( i=1; i<=spA->iNoColumns; i++ )
             free ((char *) spA->spColUnits[i-1].units_name );
         free ((char *) spA->spRowUnits);
         free ((char *) spA->spColUnits);
      }
      free ((char *) spA->cpMatrixName);
      free ((char *) spA);
      spA = (MATRIX *)NULL;
}

/*
 *  =======================================================================
 *  MatrixFreeIndirectDouble() : Free memory for INDIRECT storage of data
 *                               type DOUBLE.
 *  
 *  Input  :  double **d  -- Pointer to matrix of data type double.
 *         :  int iNoRows -- No of Rows in Matrix.
 *  Output :  void
 *  =======================================================================
 */

#ifdef __STDC__
void MatrixFreeIndirectDouble( double **d, int iNoRows)
#else  /* Start case not STDC */
void MatrixFreeIndirectDouble( d, iNoRows)
double     **d;
int    iNoRows;
#endif /* End case not STDC */
{
int ii;

   if( d == (double **)NULL )  return;

   for(ii = 1; ii <= iNoRows; ii++)
       free ((char *) d[ii-1]);
   free ((char *) d);
}

/*
 *  =======================================================================
 *  MatrixFreeIndirectInteger() : Free memory for INDIRECT storage of data
 *                                type INTEGER.
 *  
 *  Input  :  int     **i  -- Pointer to matrix of data type double.
 *         :  int iNoRows  -- No of Rows in Matrix.
 *  Output :  void
 *  =======================================================================
 */

#ifdef __STDC__
void MatrixFreeIndirectInteger( int **i, int iNoRows)
#else  /* Start case not STDC */
void MatrixFreeIndirectInteger( i, iNoRows)
int        **i;
int    iNoRows;
#endif /* End case not STDC */
{
int ii;

   if( i == (int **)NULL )  return;

   for(ii = 1; ii <= iNoRows; ii++)
       free ((char *) i[ii-1]);
   free ((char *) i);
}


/*
 *  =========================================================
 *  MatrixCompareSize() : Compare sizes of matrices [A] and [B].
 *                      : If(Size [A] == Size[B])
 *                           return TRUE.
 *                        else
 *                           return FALSE.
 *  =========================================================
 */

#ifdef __STDC__
MatrixCompareSize( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MatrixCompareSize( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{

    if((spA->iNoRows    != spB->iNoRows) ||
       (spA->iNoColumns != spB->iNoColumns)) {
        (void) printf("Problem : spA->iNoRows = %4d spA->iNoColumns = %4d\n",
                                 spA->iNoRows, spB->iNoColumns);
        (void) printf("Problem : spB->iNoRows = %4d spB->iNoColumns = %4d\n",
                                 spB->iNoRows, spB->iNoColumns);

        return (FALSE);
    }
    else
        return (TRUE);

}

/*
 *  =========================================================
 *  MatrixCompareInsideDimensions() : Compare inside
 *  dimensions of matrices [A] and [B].
 *
 *  If(Size [A] == Size[B]) return TRUE.
 *  else                    return FALSE.
 *  =========================================================
 */

#ifdef __STDC__
MatrixCompareInsideDimensions( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MatrixCompareInsideDimensions( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{

    if(spA->iNoColumns != spB->iNoRows) {
        (void) printf("Problem : spA->iNoRows = %4d spA->iNoColumns = %4d\n",
                                 spA->iNoRows, spB->iNoColumns);
        (void) printf("Problem : spB->iNoRows = %4d spB->iNoColumns = %4d\n",
                                 spB->iNoRows, spB->iNoColumns);

        return (FALSE);
    }
    else
        return (TRUE);

}


/*
 *  =======================================================================
 *  MatrixAddIndirectDouble() : Add Matrices [C] = [A] + [B].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to (nxn) matrix A.
 *         :  MATRIX     *spB  -- Pointer to (nxn) matrix B.
 *  Output :  MATRIX     *spC  -- Pointer to (nxn) matrix C.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixAddIndirectDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixAddIndirectDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
int  ii, ij;
int  length, length1, length2;
DIMENSIONS *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Compare Dimensions of Matrices */

       if(MatrixCompareSize( spA, spB ) != TRUE)
          FatalError("Execution halted in MatrixAddIndirectDouble()",
                     "Inconsistent Dimensions",
                     (char *) NULL);

       spC = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, spA->iNoRows, spB->iNoColumns);

       switch( UNITS_SWITCH) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++)
                UnitsCopy( &(spC->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]) );
            for(ij = 1; ij <= spA->iNoColumns; ij++)
                UnitsCopy( &(spC->spColUnits[ij-1]), &(spA->spColUnits[ij-1]) );

         /* [b] : Check Units and Add Matrices : [C] = [A] + [B] */

            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {

                    d1 = UnitsMult( &(spA->spRowUnits[ii-1]), &(spA->spColUnits[ij-1]) );
                    d2 = UnitsMult( &(spB->spRowUnits[ii-1]), &(spB->spColUnits[ij-1]) );

                    if(SameUnits(d1, d2) == TRUE) 
                        spC->uMatrix.daa[ii-1][ij-1] 
                        = spA->uMatrix.daa[ii-1][ij-1] + spB->uMatrix.daa[ii-1][ij-1];
                    else {
                       printf("For row No %d, column No %d \n", ii, ij);
                       FatalError("In MatrixAddIndirectDouble(): Inconsistent Units",
                       (char *)NULL);
                    }
                    free((char *) d1->units_name);
                    free((char *) d1);
                    free((char *) d2->units_name);
                    free((char *) d2);
                }
            }
            break;
         case OFF:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {
                    spC->uMatrix.daa[ii-1][ij-1]
                    = spA->uMatrix.daa[ii-1][ij-1] + spB->uMatrix.daa[ii-1][ij-1];
                }
            }
            break;
         default:
            break;
       }

       return ( spC );
}

/*
 *  =======================================================================
 *  MatrixAddReplaceIndirectDouble() : Matrix Replacement [A] = [A] + [B].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to (nxn) matrix A.
 *         :  MATRIX     *spB  -- Pointer to (nxn) matrix B.
 *  Output :  MATRIX     *spC  -- Pointer to (nxn) matrix A.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixAddReplaceIndirectDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixAddReplaceIndirectDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
int  ii, ij;
int  length1, length2;
DIMENSIONS *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Compare Dimensions of Matrices */

       if(MatrixCompareSize( spA, spB ) != TRUE)
          FatalError("Execution halted in MatrixAddIndirectDouble()",
                     "Inconsistent Dimensions",
                     (char *) NULL);

    /* [b] : Check Units and Add Matrices : [A] = [A] + [B] */

       switch( UNITS_SWITCH) {
          case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {

                    d1 = UnitsMult( &(spA->spRowUnits[ii-1]),&(spA->spColUnits[ij-1]));
                    d2 = UnitsMult( &(spB->spRowUnits[ii-1]),&(spB->spColUnits[ij-1]));

                    if(SameUnits(d1, d2) == TRUE) 
                      spA->uMatrix.daa[ii-1][ij-1] += spB->uMatrix.daa[ii-1][ij-1];
                    else {
                     printf("For row No %d, column No %d \n", ii, ij);
                     FatalError("In MatrixAddIndirectDouble(): Inconsistent Units",(char *)NULL);
                    }
                    free((char *) d1->units_name);
                    free((char *) d1);
                    free((char *) d2->units_name);
                    free((char *) d2);
                }
            }
            break;
         case OFF:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {
                      spA->uMatrix.daa[ii-1][ij-1] += spB->uMatrix.daa[ii-1][ij-1];
                }
            }
            break;
         default:
            break;
       }

       free((char *) spA->cpMatrixName);
       spA->cpMatrixName = (char *)NULL;
       return (spA);
}

/*
 *  =======================================================================
 *  MatrixSubIndirectDouble() : Add Matrices [C] = [A] - [B].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to (nxn) matrix A.
 *         :  MATRIX     *spB  -- Pointer to (nxn) matrix B.
 *  Output :  MATRIX     *spC  -- Pointer to (nxn) matrix C.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixSubIndirectDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixSubIndirectDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
int  ii, ij;
int  length, length1, length2;
DIMENSIONS *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Compare Dimensions of Matrices */

       if(MatrixCompareSize( spA, spB ) != TRUE)
          FatalError("Execution halted in MatrixAddIndirectDouble()",
                     "Inconsistent Dimensions",
                     (char *) NULL);

       spC = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, spA->iNoRows, spB->iNoColumns);

       switch( UNITS_SWITCH) {
          case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++)
                UnitsCopy( &(spC->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]) );
            for(ij = 1; ij <= spA->iNoColumns; ij++)
                UnitsCopy( &(spC->spColUnits[ij-1]), &(spA->spColUnits[ij-1]) );

         /* [b] : Check Units and Sub Matrices : [C] = [A] - [B] */

            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {

                    d1 = UnitsMult(&(spA->spRowUnits[ii-1]),&(spA->spColUnits[ij-1]));
                    d2 = UnitsMult(&(spB->spRowUnits[ii-1]),&(spB->spColUnits[ij-1]));

                    if(SameUnits(d1, d2) == TRUE) 
                      spC->uMatrix.daa[ii-1][ij-1] 
                      = spA->uMatrix.daa[ii-1][ij-1] - spB->uMatrix.daa[ii-1][ij-1];
                    else {
                     printf("For row No %d, column No %d \n", ii, ij);
                     FatalError("In MatrixAddIndirectDouble(): Inconsistent Units",(char *)NULL);
                    }
                    free((char *) d1->units_name);
                    free((char *) d1);
                    free((char *) d2->units_name);
                    free((char *) d2);
                }
            }
            break;
         case OFF:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {
                      spC->uMatrix.daa[ii-1][ij-1] 
                      = spA->uMatrix.daa[ii-1][ij-1] - spB->uMatrix.daa[ii-1][ij-1];
                }
            }
            break;
         default:
            break;
      }

       return ( spC );
}

/*
 *  =======================================================================
 *  MatrixSubReplaceIndirectDouble() : Matrix Replacement [A] = [A] - [B].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to (nxn) matrix A.
 *         :  MATRIX     *spB  -- Pointer to (nxn) matrix B.
 *  Output :  MATRIX     *spC  -- Pointer to (nxn) matrix A.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixSubReplaceIndirectDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixSubReplaceIndirectDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
int ii, ij;
int  length1, length2;
DIMENSIONS *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Compare Dimensions of Matrices */

       if(MatrixCompareSize( spA, spB ) != TRUE)
          FatalError("Execution halted in MatrixAddIndirectDouble()",
                     "Inconsistent Dimensions",
                     (char *) NULL);

    /* [b] : Check Units and Sub Matrices : [A] = [A] - [B] */

       switch( UNITS_SWITCH) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {

                    d1 = UnitsMult(&(spA->spRowUnits[ii-1]),&(spA->spColUnits[ij-1]));
                    d2 = UnitsMult(&(spB->spRowUnits[ii-1]),&(spB->spColUnits[ij-1]));

                    if(SameUnits(d1, d2) == TRUE) 
                      spA->uMatrix.daa[ii-1][ij-1] -= spB->uMatrix.daa[ii-1][ij-1];
                    else {
                     printf("For row No %d, column No %d \n", ii, ij);
                     FatalError("In MatrixAddIndirectDouble(): Inconsistent Units",(char *)NULL);
                    }
                    free((char *) d1->units_name);
                    free((char *) d1);
                    free((char *) d2->units_name);
                    free((char *) d2);
                }
            }
            break;
         case OFF:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                for(ij = 1; ij <= spA->iNoColumns; ij++) {
                      spA->uMatrix.daa[ii-1][ij-1] -= spB->uMatrix.daa[ii-1][ij-1];
                }
            }
            break;
         default:
            break;
      }

       free((char *) spA->cpMatrixName);
       spA->cpMatrixName = (char *)NULL;
       return (spA);
}

/*
 *  =======================================================================
 *  MatrixMultIndirectDouble() : Matrix Multiplication [C] = [A].[B].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to matrix A.
 *         :  MATRIX     *spB  -- Pointer to matrix B.
 *  Output :  MATRIX     *spC  -- Pointer to matrix C.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixMultIndirectDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixMultIndirectDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
MATRIX    *spC;
int ii, ij, ik;
int  length1, length2, length;
DIMENSIONS *d1, *d2, *d3, *d4, *d5;
int UNITS_SWITCH;

   /* [a] : Compare Inside Dimensions of Matrices */

   if(MatrixCompareInsideDimensions( spA, spB ) != TRUE)
      FatalError("Execution halted in MatrixMultIndirectDouble()",
                 "Inconsistent Inside Dimensions",
                 (char *) NULL);

   /* [b] : Allocate memory for result of matrix multiply */

   spC = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, spA->iNoRows, spB->iNoColumns);
       
   /* [c] : Multiply matrices (with units) */

   UNITS_SWITCH = CheckUnits();

   if (UNITS_SWITCH == ON ) {

       for(ii = 1; ii <= spA->iNoRows; ii++) {
       for(ij = 1; ij <= spB->iNoColumns; ij++) {

           d4 = (DIMENSIONS *) MyCalloc( 1, sizeof(DIMENSIONS) );
           ZeroUnits(d4);

           for(ik = 1; ik <= spA->iNoColumns; ik++) {

               d1 = UnitsMult(&(spA->spRowUnits[ii-1]),&(spA->spColUnits[ik-1]));
               d2 = UnitsMult(&(spB->spRowUnits[ik-1]),&(spB->spColUnits[ij-1]));
               d3 = UnitsMult(d1,d2);

               if( ik==1 && d3->units_name!=(char *)NULL )
                   UnitsCopy(d4, d3);

               if( SameUnits(d3, d4) == TRUE) {
                   spC->uMatrix.daa[ii-1][ij-1] 
                        += spA->uMatrix.daa[ii-1][ik-1]*spB->uMatrix.daa[ik-1][ij-1];
               } else
                   FatalError("In MatrixMultIndirectDouble(): Inconsistent Units",
                             (char *)NULL);

               free((char *) d1->units_name);
               free((char *) d1);
               free((char *) d2->units_name);
               free((char *) d2);
               free((char *) d3->units_name);
               free((char *) d3);
          }
          free((char *) d4->units_name);
          free((char *) d4);
      }
      }

      for(ii = 1; ii <=  spC->iNoRows; ii++)
          UnitsMultRep(&(spC->spRowUnits[ii-1]),&(spB->spRowUnits[0]),&(spA->spRowUnits[ii-1]));
      for(ij = 1; ij <= spC->iNoColumns; ij++)
          UnitsMultRep(&(spC->spColUnits[ij-1]),&(spA->spColUnits[0]),&(spB->spColUnits[ij-1]));

   }

   /* [d] : Multiply matrices (without units) */

   if (UNITS_SWITCH == OFF ) {
       for(ii = 1; ii <= spA->iNoRows; ii++)
       for(ij = 1; ij <= spB->iNoColumns; ij++)
           for(ik = 1; ik <= spA->iNoColumns; ik++)
               spC->uMatrix.daa[ii-1][ij-1]
                    += spA->uMatrix.daa[ii-1][ik-1]*spB->uMatrix.daa[ik-1][ij-1];
   }

   return (spC);
}

/*
 *  =======================================================================
 *  MatrixScaleIndirectDouble() : Multiply Matrix by Scalar
 * 
 *  Input  :  MATRIX *spA        -- Pointer to matrix data structure [A].
 *         :  double scale       -- double : scale factor c.
 *  Output :  MATRIX *spB        -- Pointer to scaled matrix [B] = c.[A].
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixScaleIndirectDouble( MATRIX *spA, double dScale)
#else  /* Start case not STDC */
MATRIX *MatrixScaleIndirectDouble( spA, dScale)
MATRIX   *spA;
double dScale;
#endif /* End case not STDC */
{
MATRIX     *spB;
int      ii, ij;

    /* [a] : Check Input matrix [A] */

       if(spA == (MATRIX *) NULL) 
          FatalError("In MatrixScaleIndirectDouble() : [A] == NULL",
                    (char *) NULL);

    /* [b] : Scale Indirect Double Matrix */

       spB = MatrixCopyIndirectDouble( spA );

       for( ii = 1; ii <= spA->iNoRows; ii++) 
           for( ij = 1; ij <= spA->iNoColumns; ij++)
                spB->uMatrix.daa[ii-1][ij-1] = dScale * spA->uMatrix.daa[ii-1][ij-1];

       return ( spB);
}

#ifdef __STDC__
double MatrixContentScaleIndirectDouble(MATRIX *m, int row_no, int col_no )
#else
double MatrixContentScaleIndirectDouble(m, row_no, col_no )
MATRIX  *m;
int row_no;        /* row number    */
int col_no;        /* column number */
#endif
{
int   i, j, UnitsType;
double             da;
DIMENSIONS     *dimen;

#ifdef DEBUG
       printf("*** Enter MatrixContentScaleIndirectDouble() : m->iNoRows    = %4d\n", m->iNoRows);
       printf("                        : m->iNoColumns = %4d\n", m->iNoColumns);
#endif
      if(CheckUnits()==OFF)
         FatalError("You have to set units ON to use this function",
                    "In MatrixContentScaleIndirectDouble",(char *)NULL );
 
      i = row_no;
      j = col_no;

      da = m->uMatrix.daa[i-1][j-1];

#ifdef DEBUG
     printf("\n i = %d, j = %d da = %le \n", i, j, da);
#endif
      if(m->spRowUnits[i-1].scale_factor == 0) {
         printf("==> for Row %d, scale_factor of %s = 0", i, m->spColUnits[i-1].units_name);
         FatalError("Fatal error in MatrixContentScaleIndirectDouble(): ",(char *)NULL);
      }

      if(m->spColUnits[j-1].scale_factor == 0) {
         printf("==> for column %d, scale_factor of %s = 0 \n",j, m->spColUnits[j-1].units_name);
         FatalError("Fatal error in MatrixContentScaleIndirectDouble(): ",(char *)NULL);
      }

      UnitsType = CheckUnitsType();

      if(m->spColUnits[j-1].units_name != NULL) {
         switch(UnitsType) {
           case SI:
             if(!strcmp(m->spColUnits[j-1].units_name, "deg_F") )
                da = ConvertTempUnits(m->spColUnits[j-1].units_name, da, US);
           break;
           case US:
             if(!strcmp(m->spColUnits[j-1].units_name, "deg_C") )
                da = ConvertTempUnits(m->spColUnits[j-1].units_name, da, SI);
           break;
         }
       }

      if(m->spRowUnits[i-1].units_name != NULL) {
         switch(UnitsType) {
           case SI:
             if(!strcmp(m->spRowUnits[i-1].units_name, "deg_F") )
                da = ConvertTempUnits(m->spRowUnits[i-1].units_name, da, US);
           break;
           case US:
             if(!strcmp(m->spRowUnits[i-1].units_name, "deg_C") )
                da = ConvertTempUnits(m->spRowUnits[i-1].units_name, da, SI);
           break;
         }
       }

     /* [a] Scale for Column of matrix m */

          da = da / m->spColUnits[j-1].scale_factor;

     /* [b] Scale for Row of matrix m */

          da = da / m->spRowUnits[i-1].scale_factor;

#ifdef DEBUG
       printf("*** Leave MatrixContentScaleIndirectDouble()\n");
#endif
   return (da);
}

/*
 *  =======================================================================
 *  MatrixNegateIndirectDouble() : Matrix Negation [B] = -[A].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to matrix A.
 *  Output :  MATRIX     *spB  -- Pointer to matrix B.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixNegateIndirectDouble( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixNegateIndirectDouble( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX *spB;
int  ii, ij;
int  length;

     spB = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, spA->iNoRows, spA->iNoColumns);

     for(ii = 1; ii <= spA->iNoRows; ii = ii + 1)
     for(ij = 1; ij <= spA->iNoColumns; ij = ij + 1)
         spB->uMatrix.daa[ii-1][ij-1] = -spA->uMatrix.daa[ii-1][ij-1];

     if( CheckUnits() == ON ) {
        for(ii = 1; ii <= spA->iNoRows; ii++)
            UnitsCopy( &(spB->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]) );
        for(ij = 1; ij <= spA->iNoColumns; ij++)
            UnitsCopy( &(spB->spColUnits[ij-1]), &(spA->spColUnits[ij-1]) );
     }
         
     return (spB);
}

/*
 *  =======================================================================
 *  MatrixNegateReplaceIndirectDouble() : Matrix Negation [A] = -[A].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to matrix A.
 *  Output :  MATRIX     *spB  -- Pointer to matrix A.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixNegateReplaceIndirectDouble( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixNegateReplaceIndirectDouble( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
int  ii, ij;

     for(ii = 1; ii <= spA->iNoRows; ii = ii + 1)
     for(ij = 1; ij <= spA->iNoColumns; ij = ij + 1)
         spA->uMatrix.daa[ii-1][ij-1] *= -1.0;

     return (spA);
}


/*
 *  =======================================================================
 *  MatrixTransposeIndirectDouble() : Matrix Transpose [B] = [A]^T.
 *  
 *  Input  :  MATRIX *spA  -- Pointer to matrix A.
 *  Output :  MATRIX *spB  -- Pointer to matrix B.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixTransposeIndirectDouble( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixTransposeIndirectDouble( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX  *spB;
int   ii, ij;
int   length;

    /* [a] : Transpose Matrix */

       spB = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, spA->iNoColumns, spA->iNoRows);

       for(ii = 1; ii <= spA->iNoRows; ii++)
       for(ij = 1; ij <= spA->iNoColumns; ij++)
           spB->uMatrix.daa[ ij-1 ][ ii-1 ] = spA->uMatrix.daa[ ii-1 ][ ij-1 ];

       if( CheckUnits() == ON ) {
          for(ii = 1; ii <= spA->iNoRows; ii++)
              UnitsCopy(&(spB->spColUnits[ii-1]), &(spA->spRowUnits[ii-1]));
          for(ij = 1; ij <= spA->iNoColumns; ij++)
              UnitsCopy(&(spB->spRowUnits[ij-1]), &(spA->spColUnits[ij-1]));
       }
           
       return (spB);
}

/*
 *  ===================================================
 *  MatrixCopyIndirectDouble() : Matrix Copy [B] = [A].
 *  
 *  Input  :  MATRIX *spA  -- Pointer to matrix A.
 *  Output :  MATRIX *spB  -- Pointer to matrix B.
 *  ===================================================
 */

#ifdef __STDC__
MATRIX *MatrixCopyIndirectDouble( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixCopyIndirectDouble( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX *spB;
int  ii, ij;
int  length;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Make Copy of Matrix */

       spB = MatrixAllocIndirect( spA->cpMatrixName, DOUBLE_ARRAY, spA->iNoRows, spA->iNoColumns);

       switch( UNITS_SWITCH) {
         case ON:
             for(ii = 1; ii <= spA->iNoRows; ii++) {
                 UnitsCopy(&(spB->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]));
                 for(ij = 1; ij <= spA->iNoColumns; ij++) {
                    spB->uMatrix.daa[ii-1][ij-1] = spA->uMatrix.daa[ii-1][ij-1];
                 }
             }
             for(ij = 1; ij <= spA->iNoColumns; ij++)
                 UnitsCopy(&(spB->spColUnits[ij-1]), &(spA->spColUnits[ij-1]));
             break;
         case OFF:
             for(ii = 1; ii <= spA->iNoRows; ii++)
                 for(ij = 1; ij <= spA->iNoColumns; ij++)
                    spB->uMatrix.daa[ii-1][ij-1] = spA->uMatrix.daa[ii-1][ij-1];
             break;
         default:
             break;
       }

       return ( spB );
}


#ifdef __STDC__
MATRIX *MatrixInverseIndirectDouble ( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixInverseIndirectDouble ( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX     *spAinv;
MATRIX *lu, *F, *X;
double          *d;
VECTOR       *indx;
int           i, j;
int         length;
int   UNITS_SWITCH;
     
     UNITS_SWITCH = CheckUnits();
     /* [0] : Check whether matrix spA is a square matrix    */
        
        if(spA->iNoRows != spA->iNoColumns) {
           printf("******matrix %s is not a square matrix :\n", spA->cpMatrixName);
           printf("******       %s->iNoRows    = %d\n", spA->cpMatrixName, spA->iNoRows); 
           printf("******       %s->iNoColumns = %d\n", spA->cpMatrixName, spA->iNoColumns); 
           printf("******check input file \n");
           FatalError("****** Only the square matrix has inverse. In MatrixInverse()",(char *)NULL);
        }

     /* [a] : Set up the column vector of ones */

     F    = MatrixAllocIndirect("F", DOUBLE_ARRAY, spA->iNoRows, 1);
     for(j = 1; j<= spA->iNoRows; j++)
            F->uMatrix.daa[j-1][0] = 0.0;

     if( UNITS_SWITCH == ON ) {
        ZeroUnits(&(F->spColUnits[0]));
        for(i = 1; i <= spA->iNoRows; i++)
            ZeroUnits(&(F->spRowUnits[i-1]));
     }

     /* [b] : Compute [L][U] decomposition of spA (copy) */

     indx = SetupPivotVector(spA);
     lu   = LUDecompositionIndirect(spA,indx);

     /* [c] : Solve series of equations for inverse */

     spAinv = MatrixAllocIndirect("Inverse(spA)", DOUBLE_ARRAY, spA->iNoRows, spA->iNoColumns);
     for(i = 1; i <= spA->iNoColumns; i++) {
         F->uMatrix.daa[i-1][0] = 1.0;
         X = LUSubstitutionIndirect((char *)NULL,indx,lu,F);
         for(j = 1; j<= spA->iNoRows; j++)
             spAinv->uMatrix.daa[j-1][i-1] = X->uMatrix.daa[j-1][0];
         MatrixFreeIndirect(X);
         F->uMatrix.daa[i-1][0] = 0.0;
     }

     /* [d] : Assign units to spAinv */
     if( UNITS_SWITCH == ON ) {
        for(i = 1; i <= spA->iNoRows; i++){
           UnitsPowerRep( &(spAinv->spColUnits[i-1]), &(spA->spRowUnits[i-1]), -1.0, YES );
           UnitsPowerRep( &(spAinv->spRowUnits[i-1]), &(spA->spColUnits[i-1]), -1.0, YES );
        }
     }

     /* [e] : Free memory, return result */

     VectorFree(indx);
     MatrixFree( lu );
     MatrixFree( F );

     return (spAinv);
}

/* 
 *  =====================================================================
 *  LUDecompositionDouble() : Use method of Crout Reduction with pivoting
 *  and scaled equtions to decompose a (nxn) matrix [A] into [L][U].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to (nxn) matrix A.
 *         :  VECTOR *spPivot  -- Pointer to (nx1) pivot vector.
 *  Output :  MATRIX     *spA  -- Pointer to (nxn) [L][U] product.
 *  =====================================================================
 */

#ifdef __STDC__
MATRIX *LUDecompositionIndirect( MATRIX *spA, VECTOR *spPivot)
#else  /* Start case not STDC */
MATRIX *LUDecompositionIndirect( spA, spPivot )
MATRIX      *spA;
VECTOR *spPivot;
#endif /* End case not STDC */
{
MATRIX         *spLU;
VECTOR      *spScale;
double          dSum;
double    dNumerator;
double  dDenominator;
int iOrder1, iOrder2;
int       ii, ij, ik;

    /* [a] : Make sure that Matrix spA is square    */

       if((spA->iNoRows != spA->iNoColumns)) {
           (void) printf("Problem : spA->iNoRows = %4d spA->iNoColumns = %4d\n",
                                    spA->iNoRows,      spA->iNoColumns);
           FatalError("Execution halted in LUDecomposition()",
                      "Matrix spA must be square",
                      (char *) NULL);
       }

    /* [b] : Setup Matrix of Scale Factors */

       spLU    = MatrixCopyIndirectDouble( spA );
       spScale = SetupScaleFactors( spLU );
    
    /* [c] : Loop over Columns and Compute Crout Reduction */

       ij = 1;
       Pivot( spLU, spScale, spPivot, ij );
       for (ij = 2; ij <= spLU->iNoColumns; ij = ij + 1) {
            iOrder1 = spPivot->uVector.ia[ 0 ];

            spLU->uMatrix.daa[iOrder1-1][ij - 1] =
                 spLU->uMatrix.daa[iOrder1-1][ij - 1]/spLU->uMatrix.daa[iOrder1-1][ 0 ];
       }

       for (ij = 2; ij <= spLU->iNoColumns - 1; ij = ij + 1) {
            for (ii = ij; ii <= spLU->iNoColumns; ii = ii + 1) {

                 iOrder1 = spPivot->uVector.ia[ ii - 1 ];

                 dSum = 0.0;
                 for (ik = 1; ik <= ij - 1; ik = ik + 1) {
                      iOrder2 = spPivot->uVector.ia[ ik - 1 ];
                      dSum += spLU->uMatrix.daa[ iOrder1-1 ][ ik - 1 ] *
                              spLU->uMatrix.daa[ iOrder2-1 ][ ij - 1 ];
                 }

                 spLU->uMatrix.daa[ iOrder1-1 ][ ij - 1 ] -= dSum;
            }

            Pivot( spLU, spScale, spPivot, ij );

            iOrder1 = spPivot->uVector.ia[ ij - 1 ];

            for (ik = ij + 1; ik <= spLU->iNoColumns; ik = ik + 1) {

                 dSum = 0.0;
                 for (ii = 1; ii <= ij - 1; ii = ii + 1) {
                      iOrder2 = spPivot->uVector.ia[ ii - 1 ];
                      dSum += spLU->uMatrix.daa[ iOrder1-1 ][ ii - 1 ] *
                              spLU->uMatrix.daa[ iOrder2-1 ][ ik - 1 ];
                 }

                 dNumerator   = spLU->uMatrix.daa[ iOrder1-1 ][ ik-1 ] - dSum;
                 dDenominator = spLU->uMatrix.daa[ iOrder1-1 ][ ij-1 ];

                 spLU->uMatrix.daa[ iOrder1-1 ][ ik - 1 ] = dNumerator/dDenominator;
            }
       }

       iOrder1 = spPivot->uVector.ia[ spLU->iNoColumns - 1 ];

       dSum = 0.0;
       for (ik = 1; ik <= spLU->iNoColumns - 1; ik = ik + 1) {
            iOrder2 = spPivot->uVector.ia[ ik - 1 ];
            dSum += spLU->uMatrix.daa[ iOrder1-1 ][ ik - 1 ] *
                    spLU->uMatrix.daa[ iOrder2-1 ][ spLU->iNoColumns - 1 ];
       }

       spLU->uMatrix.daa[ iOrder1 - 1 ][ spLU->iNoColumns - 1 ] -= dSum;

       VectorFree( spScale );
       return( spLU );
}


/* 
 *  =====================================================================
 *  LUSubstitutionDouble() : Compute Forward and Backward Substitution.
 *  
 *  Input  :  cpMatrixNameOfSolution -- Name of Solution Matrix
 *         :  VECTOR *spPivot        -- Pointer to (nx1) pivot vector.
 *         :  MATRIX    *spLU        -- Pointer to (nxn) matrix [L][U].
 *         :  MATRIX     *spB        -- Pointer to (nx1) matrix [B].
 *  Output :  MATRIX      spX        -- Pointer to (nx1) solution matrix.
 *  =====================================================================
 */

#ifdef __STDC__
MATRIX *LUSubstitutionIndirect( char *cpMatrixNameOfSolution, VECTOR *spPivot, MATRIX *spLU, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *LUSubstitutionIndirect( cpMatrixNameOfSolution, spPivot, spLU, spB)
char *cpMatrixNameOfSolution;
VECTOR              *spPivot;
MATRIX                 *spLU;
MATRIX                  *spB;
#endif /* End case not STDC */
{
MATRIX            *spDisp;
double       dDenominator;
double         dNumerator;
double               dSum;
int       iOrder1, ii, ij;
int                length;
DIMENSIONS   *d, *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
   /* [a] : Allocate Memory for Displacement Vector */

      spDisp = MatrixAllocIndirect( cpMatrixNameOfSolution, DOUBLE_ARRAY, spLU->iNoRows, 1 );

   /* [0] : Calculate the units of displacements */

   if( UNITS_SWITCH == ON ) {
      for (ii = 1; ii <= spLU->iNoRows; ii++) {
         d1 = UnitsMult( &(spB->spRowUnits[ii-1]), &(spB->spColUnits[0]) );
         d2 = UnitsMult( &(spLU->spRowUnits[ii-1]), &(spLU->spColUnits[ii-1]) );
         UnitsDivRep( &(spDisp->spRowUnits[ii-1]), d1, d2, YES); 

         free((char *) d1->units_name);
         free((char *) d1);
         free((char *) d2->units_name);
         free((char *) d2);
      }
      ZeroUnits(&(spDisp->spColUnits[0]));
   }

   /* [b] : Forward Substitution with Unscrambling of Permutations */

      iOrder1 = spPivot->uVector.ia[ 0 ];
      spDisp->uMatrix.daa[0][0] = spB->uMatrix.daa[ iOrder1 - 1 ][ 0 ] /
                                  spLU->uMatrix.daa[ iOrder1 - 1 ][ 0 ];

      for( ii = 2; ii <= spLU->iNoRows; ii = ii + 1) {

           iOrder1 = spPivot->uVector.ia[ ii - 1 ];

           dSum = 0.0;
           for( ij = 1; ij <= ii - 1; ij = ij + 1) {
                dSum += spLU->uMatrix.daa[ iOrder1 - 1 ][ ij - 1] *
                        spDisp->uMatrix.daa[ ij - 1][0];
           }               

           dNumerator   = spB->uMatrix.daa[ iOrder1-1 ][0] - dSum;
           dDenominator = spLU->uMatrix.daa[ iOrder1-1 ][ ii-1 ];

           spDisp->uMatrix.daa[ ii - 1 ][0] = dNumerator/dDenominator;
      }               

   /* [c] : Backsubstitution with Unscrambling of Permutations */

      for( ii = spLU->iNoRows - 1; ii >= 1; ii = ii - 1) {

           iOrder1 = spPivot->uVector.ia[ ii - 1 ];

           dSum = 0.0;
           for( ij = ii + 1; ij <= spLU->iNoRows; ij = ij + 1) {
                dSum += spLU->uMatrix.daa[ iOrder1 - 1 ][ ij - 1] *
                        spDisp->uMatrix.daa[ ij - 1][0];
           }               

           spDisp->uMatrix.daa[ ii - 1 ][0] -= dSum;
      }               

      return(spDisp);
}


/*
 *  ================================================================
 *  CholeskyDecompositionIndirect() : [L][L]^T Cholesky Factorization 
 *                                    of Symmetric Indirect Matrix.
 * 
 *  Input  :  MATRIX *spA   -- Pointer to symmetric matrix [A].
 *  Output :  MATRIX *spA   -- Pointer to [L][L]^T.
 *  ================================================================
 */

#ifdef __STDC__
MATRIX *CholeskyDecompositionIndirect( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *CholeskyDecompositionIndirect( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
double     dSum;
int  ii, ij, ik;

    /* [a] : Check Input matrix [A] */

       if( spA == NULL )
           FatalError("In CholeskyDecompositionIndirect() : Pointer to matrix is NULL",
                     (char *) NULL);

       if( spA->eType != DOUBLE_ARRAY ) 
           FatalError("In CholeskyDecompositionIndirect() : spA->eType != DOUBLE_ARRAY",
                     (char *) NULL);

       if( spA->eRep != INDIRECT ) 
           FatalError("In CholeskyDecompositionIndirect() : spA->eRep != INDIRECT",
                     (char *) NULL);

    /* [c] : Compute [A] = [L][L]^T */

       for (ij = 1; ij <= spA->iNoColumns; ij = ij + 1) {
            for (ii = ij; ii <= spA->iNoRows; ii = ii + 1) {

            if( ii == ij) {
                dSum  = 0.0;
                for (ik = 1; ik < ii; ik = ik + 1)
                     dSum += spA->uMatrix.daa[ ii-1 ][ ik-1 ] *
                             spA->uMatrix.daa[ ii-1 ][ ik-1 ];

                if(spA->uMatrix.daa[ ii-1 ][ ii-1 ] > dSum)
                   spA->uMatrix.daa[ii-1][ii-1] = sqrt(spA->uMatrix.daa[ ii-1 ][ ii-1 ] - dSum);
                else
                   FatalError("In CholeskyDecompositionIndirect()",
                              "Matrix is not positive definite",
                              (char *) NULL);
            }

            if( ii > ij) {
                dSum  = 0.0;
                for (ik = 1; ik < ij; ik = ik + 1)
                     dSum += spA->uMatrix.daa[ ii-1 ][ ik-1 ] *
                             spA->uMatrix.daa[ ij-1 ][ ik-1 ];

                spA->uMatrix.daa[ii-1][ij-1] = (spA->uMatrix.daa[ ii-1 ][ ij-1 ] - dSum) /
                                                spA->uMatrix.daa[ij-1][ij-1];

                spA->uMatrix.daa[ij-1][ii-1] = spA->uMatrix.daa[ii-1][ij-1];
            }

            }
       }

       return ( spA );
}


/* 
 *  =================================================================
 *  SetupScaleFactors() : Initialize vector for scale factors. If [A]
 *  is a (nxn) matrix, then spScale will be a (nx1) vector containing
 *  the maximum absolute value in each row of [A].
 *  
 *  Input :   MATRIX spA       -- Pointer to (nxn) matrix [A].
 *  Output :  VECTOR spScale   -- Pointer to (nx1) vector spScale.
 *  =================================================================
 */ 

#ifdef __STDC__
VECTOR *SetupScaleFactors( MATRIX *spA )
#else  /* Start case not STDC */
VECTOR *SetupScaleFactors( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
VECTOR *spScale;
double    dTemp;
int      ii, ij;

    spScale = VectorAlloc("Scale Factor",  DOUBLE_ARRAY, spA->iNoRows );

    for (ii = 1; ii <= spA->iNoRows; ii = ii + 1) {
         spScale->uVector.da[ii-1] = ABS(spA->uMatrix.daa[ii-1][0]);
         for (ij = 2; ij <= spA->iNoColumns; ij = ij + 1) {
              dTemp = ABS(spA->uMatrix.daa[ii-1][ij-1]);
              spScale->uVector.da[ii-1] = (double) MAX( spScale->uVector.da[ii-1], dTemp );
         }

         if( spScale->uVector.da[ii-1] == 0)
             FatalError("Execution halted in SetupScaleFactors()",
                        "Matrix A is Singular !!",
                        (char *) NULL);
    }

    return ( spScale );
}

/* 
 *  ==================================================================
 *  SetupPivotVector() : Initialize (nx1) vector for row permutations.
 *  
 *  Input :   MATRIX spA       -- Pointer to (nxn) matrix [A].
 *  Output :  VECTOR spPivot   -- Pointer to (nx1) vector spPivot.
 *  ==================================================================
 */ 

#ifdef __STDC__
VECTOR *SetupPivotVector( MATRIX *spA )
#else  /* Start case not STDC */
VECTOR *SetupPivotVector( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
VECTOR *spPivot;
int          ii;

    spPivot = VectorAlloc(  "Row Pivots", INTEGER_ARRAY, spA->iNoRows );

    for (ii = 1; ii <= spA->iNoRows; ii = ii + 1)
         spPivot->uVector.ia[ii-1] = ii;

    return ( spPivot );
}


/* 
 *  ==================================================================
 *  Pivot() : Find largest element below diagonal and pivot rows
 *  
 *  Input  :   MATRIX spA       -- Pointer to (nxn) matrix [A].
 *         :   VECTOR spScale   -- Pointer to (nx1) vector spScale.
 *         :   VECTOR spPivot   -- Pointer to (nx1) vector spPivot.
 *         :   int iColumnNo    -- Column No for Pivoting
 *  Output :   MATRIX spA       -- [A] and [spPivot] are rearranged by
 *                                 the pivoting procedure.
 *  ==================================================================
 */ 

#ifdef __STDC__
Pivot( MATRIX *spA, VECTOR *spScale, VECTOR *spPivot, int iColumnNo )
#else  /* Start case not STDC */
Pivot( spA, spScale, spPivot, iColumnNo )
VECTOR *spScale, *spPivot;
MATRIX     *spA;
int   iColumnNo;
#endif /* End case not STDC */
{
double dLargest,     dTemp;
int  ii , iTemp,    iOrder;
int      iPivot;

    /* [a] : Find Row No of Largest Scaled Equation  */

       iPivot  = iColumnNo;
       iOrder  = spPivot->uVector.ia[ iColumnNo - 1 ];

       dLargest = spA->uMatrix.daa[ iOrder-1 ][ iColumnNo - 1 ] /
                  spScale->uVector.da[ iOrder-1 ];
       dLargest = ABS( dLargest );

       for (ii = iColumnNo + 1; ii <= spA->iNoColumns; ii = ii + 1) {
            iOrder = spPivot->uVector.ia[ ii - 1 ];

            dTemp  = spA->uMatrix.daa[ iOrder-1 ][ iColumnNo - 1 ] /
                     spScale->uVector.da[ iOrder-1 ];
            dTemp  = ABS( dTemp );

            if( dTemp > dLargest ) {
                dLargest = dTemp;
                iPivot   = ii;
            }
       }

       if( dLargest == 0 )
           FatalError("Execution halted in Pivot()",
                      "Matrix A is Singular !!",
                      (char *) NULL);

    /* [b] : Swap Row Nos in Order Vector */

       iTemp  = spPivot->uVector.ia[ iPivot - 1 ];
       spPivot->uVector.ia[ iPivot - 1 ]    = spPivot->uVector.ia[ iColumnNo - 1 ];
       spPivot->uVector.ia[ iColumnNo - 1 ] = iTemp;
}


/*
 *  =================================================================
 *  GeneralizedEigen() : Solve Eigen Problem [A][x] = [B][x][Lambda].
 *                       where [A] and [B] are small fully populated  
 *                       symmetric matrices.  
 *  
 *  We convert [A][x] = [B][x][Lambda] into [A*][y] = [y][Lambda] in
 *  a four-step procedure:  
 *  
 *        (a) : Compute [B] = [L].[L]^T  
 *        (b) : Solve   [L][C] = [A].
 *        (c) : Compute [C]^T.
 *        (d) : Solve   [L][A*]^T = [C]^T; matrix [A*] is symmetric.  
 *  
 *  Use Householder method to solve [A*][y] = [y][Lamda], followed by 
 *  backward substitutution for [L]^T.[x] = [y].
 *  
 *  Written By : M. Austin                              November 1993 
 *  =================================================================
 */

GeneralizedEigen( spA, spB, spEigenvalue, spEigenvector )
MATRIX   *spA,    *spB;
MATRIX   *spEigenvalue;
MATRIX  *spEigenvector;
{
MATRIX *spEigenvectorY;
MATRIX       *spCwork2;
MATRIX        *spAwork;
MATRIX        *spCwork;
MATRIX          *spLLT;
int        iSize;
int   ii, ij, ik;
double       sum;

#ifdef MYDEBUG
       printf("\n*** Enter GeneralizedEigen()\n");
       MatrixPrint ( spA );
       MatrixPrint ( spB );
#endif

    /* [a] : Check Input and Allocate Working Arrays */

       iSize = spA->iNoRows;
       spCwork = MatrixAllocIndirect("Working [C]", DOUBLE_ARRAY, iSize, iSize );
       spAwork = MatrixAllocIndirect("Working [A]", DOUBLE_ARRAY, iSize, iSize );

    /* [b] : Compute Cholesky Decomposition [B] = [L].[L]^T */

       spLLT = MatrixCopy( spB );
       spLLT = CholeskyDecompositionIndirect( spLLT );

    /* [c] : Solve [L][C] = [A] via Forward Substitution */

       for (ij = 1; ij <= spA->iNoRows; ij = ij + 1) {
            for (ii = 1; ii <= spA->iNoColumns; ii = ii + 1) {

            sum  = 0.0;
            for (ik = 1; ik < ii; ik = ik + 1)
                 sum += spLLT->uMatrix.daa[ ii-1 ][ ik-1 ] *
                               spCwork->uMatrix.daa[ ik-1 ][ij-1];

            spCwork->uMatrix.daa[ii-1][ij-1] = (spA->uMatrix.daa[ ii-1 ][ ij-1 ] - sum) /
                                                spLLT->uMatrix.daa[ ii-1 ][ ii-1 ];

            }
       }

    /* [d] : Replace [C] by [C]^T. */

       spCwork2 = MatrixTranspose ( spCwork );

    /* [e] : Solve [L][A*] = [C]^T via Forward Substitution */

       for (ij = 1; ij <= spA->iNoRows; ij = ij + 1) {
            for (ii = 1; ii <= spA->iNoColumns; ii = ii + 1) {

            sum  = 0.0;
            for (ik = 1; ik < ii; ik = ik + 1)
                 sum += spLLT->uMatrix.daa[ ii-1 ][ ik-1 ] *
                               spAwork->uMatrix.daa[ ik-1 ][ij-1];

            spAwork->uMatrix.daa[ii-1][ij-1] =
                    (spCwork2->uMatrix.daa[ ii-1 ][ ij-1 ] - sum) /
                     spLLT->uMatrix.daa[ ii-1 ][ ii-1 ];
            }
       }

    /* [f] : Use Householder Method to Solve Standard Eigenvalue Problem */

       iSize = spAwork->iNoRows;
       spEigenvectorY = MatrixAllocIndirect("Eigenvector [Y]", DOUBLE_ARRAY, iSize, iSize);

    /* [e] : Use Householder Transformation to convert [spStiff] to Triangular form */

       Householder( spAwork, spEigenvalue, spEigenvectorY );

    /* [g] : Compute [L]^T.[X] = [Y] via Back Substitution */

       for (ij = spAwork->iNoRows; ij >= 1 ; ij = ij - 1) {
            for (ii = 1; ii <= spAwork->iNoColumns; ii = ii + 1) {

            sum  = 0.0;
            for (ik = ij + 1; ik <= spAwork->iNoRows; ik = ik + 1)
                 sum += spLLT->uMatrix.daa[ ij-1 ][ ik-1 ] *
                               spEigenvector->uMatrix.daa[ ik-1 ][ii-1];

            spEigenvector->uMatrix.daa[ij-1][ii-1] =
                          (spEigenvectorY->uMatrix.daa[ ij-1 ][ ii-1 ] - sum) /
                           spLLT->uMatrix.daa[ ij-1 ][ ij-1 ];

            }
       }

    /* [h] : Cleanup before leaving */

       MatrixFree( spAwork );
       MatrixFree( spCwork );
       MatrixFree( spCwork2 );
       MatrixFree( spEigenvectorY );
       MatrixFree( spLLT );

#ifdef MYDEBUG
       printf("\n*** Leave GeneralizedEigen()\n");
#endif

}


/*
 *  ========================================================================
 *  Householder() : Compute Eigenvalues/Eigenvectors of Symmetric Matrix Via
 *  Householder's Transformation and QL-Algorithm.
 *  
 *  Input  :  MATRIX *spK           -- Pointer to large matrix [K].
 *         :  MATRIX *spEigenvalue  -- Pointer to (nx1) eigenvalue matrix.
 *         :  MATRIX *spEigenvector -- Pointer to (nxn) eigenvector matrix.
 *  Output :  MATRIX *spEigenvalue  -- Pointer to (nx1) eigenvalue matrix.
 *         :  MATRIX *spEigenvector -- Pointer to (nxn) eigenvector matrix.
 *  
 *  Written By: Mark Austin                                    November 1993
 *  ========================================================================
 */

Householder( spK, spEigenvalue, spEigenvector )
MATRIX            *spK;
MATRIX   *spEigenvalue;
MATRIX  *spEigenvector;
{
double  *t1, *t2;
int         stat;
int       ii, ij;
int iNoEquations;

 /* [a] : Setup Working Matrices for Eigenvalues and Eigenvectors */

    iNoEquations = spK->iNoRows;
    t1 = (double *) MyCalloc( iNoEquations, sizeof(double));
    t2 = (double *) MyCalloc( iNoEquations, sizeof(double));

    for( ii = 1; ii <= iNoEquations; ii = ii + 1)
    for( ij = 1; ij <= iNoEquations; ij = ij + 1)
         spEigenvector->uMatrix.daa[ ii-1 ][ ij-1 ] = spK->uMatrix.daa[ ii-1 ][ ij-1 ];

 /* [b] : Convert [K] to Tridiagonal Form with Householder Transformation */

    HouseHolderTransformation( spEigenvector, t1, t2, iNoEquations);

 /* [c] : Use ql Algorithm to compute Eigenvalues and Eigenvectors */

    stat = QlAlgorithm( spEigenvector, t1, t2, iNoEquations);

 /* [d] : Transfer Working Matrix to Eigenvalues */

    for( ii = 1; ii <= iNoEquations; ii = ii + 1)
         spEigenvalue->uMatrix.daa[ ii-1 ][0] = t1[ ii-1 ];

    free((char *) t1);
    free((char *) t2);

    return stat;
}


/*
 *  ==============================================================
 *  HouseHolderTrans() : Use HouseHolder Transformation to convert
 *                       Symmetric Matrix to tridiagonal form.
 *  ==============================================================
 */

static double TOLERANCE = 1.3e-16;

HouseHolderTransformation( spK , d, e, iNoEquations)
MATRIX      *spK;
double    *d, *e;
int iNoEquations;
{
int    ii, ij, ik, im;
double dTempf, dTempg;
double dTemph, dTempk;

    ii = iNoEquations;
    while (--ii) {

        im = ii - 2;
        dTempf = spK->uMatrix.daa[ii][ii-1];
        dTempg = 0.0;
        for (ij = 0; ij <= im; ij++)
             dTempg += spK->uMatrix.daa[ii][ij] * spK->uMatrix.daa[ii][ij];

        dTemph = dTempg + dTempf * dTempf;
        if (dTempg < TOLERANCE) {
                 e[ii]  = dTempf;
                 dTemph = 0.0;
        } else {
             im++;
             if (dTempf >= 0)
                 dTempg = -sqrt( dTemph );
             else
                 dTempg =  sqrt( dTemph );

             e[ii] = dTempg;
             dTemph -= dTempf * dTempg;
             spK->uMatrix.daa[ii][ii-1] = dTempf - dTempg;

             dTempf = 0;
             for (ij = 0; ij <= im; ij++) {
                  spK->uMatrix.daa[ij][ii] = spK->uMatrix.daa[ii][ij] / dTemph;
                  dTempg = 0;
                  for (ik = 0; ik <= ij; ik++)
                       dTempg += spK->uMatrix.daa[ij][ik] * spK->uMatrix.daa[ii][ik];
                  for (ik = ij+1; ik <= im; ik++)
                       dTempg += spK->uMatrix.daa[ik][ij] * spK->uMatrix.daa[ii][ik];
                  e[ij] = dTempg / dTemph;
                  dTempf += dTempg * spK->uMatrix.daa[ij][ii];
             }

             dTempk = dTempf / (dTemph + dTemph);
             for (ij = 0; ij <= im; ij++) {
                  dTempf = spK->uMatrix.daa[ii][ij];
                  dTempg = e[ij] - dTempk * dTempf;
                  e[ij]  = dTempg;
                  for (ik=0; ik <= ij; ik++)
                       spK->uMatrix.daa[ij][ik] = spK->uMatrix.daa[ij][ik] - dTempf * e[ik] -
                                                  dTempg * spK->uMatrix.daa[ii][ik];
             }
        }

        d[ii] = dTemph;

    }

    d[0] = 0.0;
    e[0] = 0.0;
    for (ii=0; ii < iNoEquations; ii++) {
         im = ii - 1;
         if (d[ii] != 0.0) 
         for (ij=0; ij <= im; ij++) {
              dTempg = 0.0;
              for (ik=0; ik <= im; ik++)
                   dTempg += spK->uMatrix.daa[ii][ik] * spK->uMatrix.daa[ik][ij];
              for (ik=0; ik <= im; ik++)
                   spK->uMatrix.daa[ik][ij] -= dTempg * spK->uMatrix.daa[ik][ii];
         }

         d[ii] = spK->uMatrix.daa[ii][ii];
         spK->uMatrix.daa[ii][ii] = 1.0;
         for (ij=0; ij <= im; ij++)
              spK->uMatrix.daa[ii][ij] = spK->uMatrix.daa[ij][ii] = 0.0;
    }
}


/* 
 *  ================================================
 *  QL() : QL Algorithm : Adapted from EISPACK......
 *  ================================================
 */

QlAlgorithm( spK, d, e, iNoEquations )
MATRIX      *spK;
double    *d, *e;
int iNoEquations;
{
int ii, ij, ik, im, in;
int imm,         iTemp;
int               flag;
double b,c,f,g,h,p,r,s;

   for (ii = 2; ii <= iNoEquations; ii = ii + 1)
        e[ ii-2 ] = e[ ii-1 ];
   e[ iNoEquations - 1 ] = 0.0;
   b = 0.0;
   f = 0.0;

   for (im = 1; im <= iNoEquations; im = im + 1) {

        ij = 0;
        h = EPS * (ABS(d[im-1]) + ABS(e[im-1]));
        if (b < h)
            b = h;

        for(in = im; in <= iNoEquations; in++)
            if(ABS(e[ in-1 ]) <= b)
               break;

        if (in != im) {

            flag = 1;
            while (flag) {
               if (ij == 30) return 1;

               ij  = ij + 1;
               imm = im + 1;

               g = d[ im-1 ];
               p = (d[ imm - 1 ] - g)/(e[ im-1 ] + e[ im-1 ]);

               r = sqrt((double) (p*p + 1.0));
               if (p < 0.0)
                   d[ im-1 ] = e[ im-1 ] / (p - r);
               else
                   d[ im-1 ] = e[ im-1 ] / (p + r);

               h = g - d[ im-1 ];

               for(ii = imm; ii <= iNoEquations; ii = ii + 1) 
                   d[ ii-1 ] = d[ ii-1 ] - h;

               f = f + h;
               p = d[ in-1 ];
               c = 1.0;
               s = 0.0;

               for (ii = 1; ii <= (in - im); ii = ii + 1) {

                    iTemp = in - ii;
                    g     = c * e[( iTemp - 1 )];
                    h     = c * p;

                    if (ABS(p) >= ABS(e[ iTemp-1 ])) {
                        c = e[ iTemp - 1 ] / p;
                        r = sqrt((double) (c * c + 1.0));
                        e[ iTemp ] = s * p * r;
                        s = c / r;
                        c = 1.0 / r;
                    } else {
                        c = p / e[ iTemp - 1 ];
                        r = sqrt((double )(c * c + 1.0));
                        e[ iTemp ] = s * e[ iTemp - 1 ] * r;
                        s = 1.0 / r;
                        c = c * s;
                    }

                    p = c * d[ iTemp - 1 ] - s * g;
                    d[ iTemp ] = h + s * (c * g + s * d[ iTemp - 1 ]);

                    for (ik = 1; ik <= iNoEquations; ik = ik + 1) {

                         h = spK->uMatrix.daa[ ik-1 ][ iTemp ];
                         spK->uMatrix.daa[ ik-1 ][ iTemp ] =
                              s * spK->uMatrix.daa[ ik-1 ][ iTemp-1 ] + c * h;
                         spK->uMatrix.daa[ ik-1 ][ iTemp - 1] =
                              c * spK->uMatrix.daa[ ik-1 ][ iTemp-1 ] - s * h;
                    }
               }

               e[ im-1 ] = s * p;
               d[ im-1 ] = c * p;
               if(ABS(e[ im-1 ]) <= b) flag = 0;
            }
        }

        d[ im-1 ] = d[ im-1 ] + f;
    }

    for (ii = 2; ii <= iNoEquations; ii = ii + 1) {
         iTemp = ii - 1;
         ik    = iTemp;
         p     = d[ iTemp-1 ];

         for (ij = ii; ij <= iNoEquations; ij = ij + 1) {
              if (ABS(d[ ij-1 ]) > ABS(p)) {
                  ik = ij;
                  p  = d[ ij - 1 ];
              }
         }

         if (ik != iTemp) {
             d[ ik-1 ] = d[ iTemp-1 ];
             d[ iTemp-1 ] = p;
             for (ij = 1; ij <= iNoEquations; ij = ij + 1) {
                  p = spK->uMatrix.daa[ ij-1 ][ iTemp-1 ];
                  spK->uMatrix.daa[ ij-1 ][ iTemp-1 ] = spK->uMatrix.daa[ ij-1 ][ ik-1 ];
                  spK->uMatrix.daa[ ij-1 ][ ik-1 ]    = p;
             }
         }
    }

    return 0;
}

