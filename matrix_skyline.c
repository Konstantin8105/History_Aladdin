/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  matrix_skyline.c : Functions for Operations on Skyline Matrices.
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
 *  Written by: Mark Austin and Lanheng Jin                              May 1994
 *  ============================================================================= 
 */

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"

/* #define MYDEBUG */

/*
 *  =======================================================================
 *  MatrixAllocSkyline() : Allocate (and free) memory for SKYLINE matrix
 *                         data structure. 
 *  
 *  Input  :  char *cpMatrixName -- Pointer to name of matrix.
 *         :  DATA_TYPE eType    -- Data type to be stored in Matrix.
 *         :  int  iNoRows       -- No of Rows in Matrix.
 *         :  int  iNoColumns    -- No of Columns in Matrix.
 *         :  int  *ld           -- Length of each column in skyline profile
 *  Output :  MATRIX *spA        -- Pointer to matrix data structure.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixAllocSkyline( char *cpMatrixName,  DATA_TYPE eType,
                            int  iNoRows,  int  iNoColumns,  int  *ld)
#else  /* start case not STDC */
MATRIX *MatrixAllocSkyline( cpMatrixName,  eType,  iNoRows, iNoColumns, ld)
char     *cpMatrixName;
DATA_TYPE        eType;
int            iNoRows;
int         iNoColumns;
int                *ld;
#endif /* end case not STDC */
{
MATRIX *spA;
int      ii;

   /* [a] : Check input parameters/arrays */

      if(ld == (int *)NULL) 
         FatalError("In MatrixAllocSkyline() : ld[] = NULL",
                   (char *) NULL);

   /* [b] : Allocate matrix parent data structure, and name */

      spA = (MATRIX *) MyMalloc( sizeof(MATRIX) );
      if(cpMatrixName != (char *)NULL)
         spA->cpMatrixName = SaveString(cpMatrixName);
      else 
         spA->cpMatrixName = (char *) NULL;

   /* [c] : Set parameters and allocate memory for matrix array */

      spA->eRep       = SKYLINE;
      spA->eType      = eType;
      spA->iNoRows    = iNoRows;
      spA->iNoColumns = iNoColumns;

      if( CheckUnits() == ON ) {
         spA->spRowUnits = (DIMENSIONS *) MyCalloc( iNoRows, sizeof(DIMENSIONS) );
         spA->spColUnits = (DIMENSIONS *) MyCalloc( iNoColumns, sizeof(DIMENSIONS) );
      }
      else {
         spA->spRowUnits = (DIMENSIONS *)NULL;
         spA->spColUnits = (DIMENSIONS *)NULL;
      }

      switch((int) spA->eType) {
          case DOUBLE_ARRAY:
               spA->uMatrix.daa = (double **)
                                  MyCalloc( iNoColumns, sizeof(double*));
               for(ii = 0; ii < iNoColumns; ii = ii + 1) {
                   spA->uMatrix.daa[ ii ] = (double *)
                                            MyCalloc((ld[ii]+1), sizeof(double));
                   spA->uMatrix.daa[ ii ][0] = ld[ ii ];
               }
               break;
          default:
               FatalError("In MatrixAllocSkyline() : Undefined spA->eType",
                         (char *) NULL);
               break;
      }

      return ( spA );
}

/*
 *  =======================================================================
 *  MatrixFreeSkyline() : Free memory in Skyline Matrix
 *  
 *  Input  :  MATRIX *Matrix     -- Pointer to matrix data structure.
 *  Output :  void
 *  =======================================================================
 */

#ifdef __STDC__
void MatrixFreeSkyline( MATRIX *spA )
#else  /* Start case not STDC */
void MatrixFreeSkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
int ii;

   if( spA == (MATRIX *)NULL )  return;

   /* [a] : Free body of skyline matrix[][] */

      switch( spA->eType ) {
          case DOUBLE_ARRAY:
               for(ii = 1; ii <= spA->iNoColumns; ii++)
                   free ((char *) (spA->uMatrix.daa[ ii-1 ]));
               free ((char *) (spA->uMatrix.daa));
               break;
          default:
               FatalError("In MatrixFreeSkyline() : Undefined spA->eType",
                         (char *) NULL);
               break;
      }

   /* [b] : Free name and parent data structure for matrix */

      if( CheckUnits() == ON ) {
         for( ii=1 ; ii<=spA->iNoRows ; ii++ )
             free ((char *) spA->spRowUnits[ii-1].units_name );
         for( ii=1 ; ii<=spA->iNoColumns ; ii++ )
             free ((char *) spA->spColUnits[ii-1].units_name );
         free ((char *) spA->spRowUnits);
         free ((char *) spA->spColUnits);
      }
      free ((char *) spA->cpMatrixName);
      free ((char *) spA);
      spA = (MATRIX *)NULL;
}


/*
 *  =======================================================================
 *  Print a Matrix [iNoRows][iNoColumns] of data type DOUBLE 
 *  Dump Profile of Skyline Matrix.
 *  
 *  Where --  spA->iNoRows         = Number of rows in matrix.
 *            spA->iNoColumns      = Number of columns in matrix.
 *            COLUMNS_ACROSS_PAGE  = Number of columns printed across page.
 * 
 *            ib                   = Current No of Matrix Block.
 *            iFirstColumn         = Number of first column in block
 *            iLastColumn          = Number of last  column in block
 * 
 *  Input  :  MATRIX *spA        -- Pointer to matrix data structure.
 *  Output :  void
 *  =======================================================================
 */

enum { COLUMNS_ACROSS_PAGE = 6 };                                        /* Item [a] */

#ifdef __STDC__                                                  
void MatrixPrintSkylineDouble( MATRIX *spA )
#else  /* Start case not STDC */
void MatrixPrintSkylineDouble( spA )
MATRIX  *spA;
#endif /* End case not STDC */
{
int ii, ij, ik, im, ib;    
int iNoBlocks; 
int iFirstColumn, iLastColumn; 
double  da;
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
         iLastColumn  = (int) MIN( ib*((int) COLUMNS_ACROSS_PAGE) ,
                             spA->iNoColumns );

         /* [c] : Print title of matrix at top of each block */          /* Item [e] */
       
         if( spA->cpMatrixName != NULL ) 
             printf("\nSKYLINE MATRIX : \"%s\"\n\n", spA->cpMatrixName);
         else 
             printf("\nSKYLINE MATRIX : \"UNTITLED\"\n\n");

         /* [d] : Label row and column nos */

         printf ("row/col        ");
         for( ii = iFirstColumn; ii <= iLastColumn; ii++ )
             printf("        %3d   ", ii);
         printf("\n");

         switch( UNITS_SWITCH ) {
           case ON:
              printf ("        units  " );
              for( ii = iFirstColumn; ii <= iLastColumn; ii++ )
                 if(spA->spColUnits[ii-1].units_name != NULL) {
                    printf(" %10s   ",spA->spColUnits[ii-1].units_name);
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
                   for( ij  = iFirstColumn; ij <= iLastColumn; ij++) {
                       da = MatrixContentScaleSkyline(spA, ii, ij);
                       printf(" %12.5e ", da );
                   }
                   printf("\n");
              }
              break;
           case OFF:
              /* [e] : Print Contents of Matrix */                            /* Item [f] */

              for( ii = 1; ii <= spA->iNoRows; ii++ ) {
                   printf("%3d ", ii);
                   printf("         ");
                   for( ij  = iFirstColumn; ij <= iLastColumn; ij++) {

                        ik = (int) MIN( ii, ij );
                        im = (int) MAX( ii, ij );

                        if((im-ik+1) <= spA->uMatrix.daa[im-1][0])
                             printf(" %12.5e ", spA->uMatrix.daa[ im-1 ][ im-ik+1 ]);
                        else
                             printf(" %12.5e ", 0.0);

                   }
                   printf("\n");
              }
              break;
           default:
              break;
         }
    }
}


/*
 *  ================================================================
 *  Matrix Operations : Addition, Subtraction, Multiplication, Copy.
 *  ================================================================
 */

/*
 *  ========================================================================
 *  MatrixReallocSkyline() : Reallocate Memory for Skyline Matrix. Eliminate
 *                           matrix elements that are smaller than MINREAL
 * 
 *  Input  :  MATRIX *spA  -- Pointer to matrix data structure.
 *  Output :  MATRIX *spB  -- Pointer to new skyline matrix.
 *  ========================================================================
 */

#ifdef __STDC__
MATRIX *MatrixReallocSkyline ( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixReallocSkyline ( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX     *spB;
int ii, ij, *ld;
int      length;

    if( spA == NULL )
        FatalError("In MatrixReallocSkyline() : spA = NULL",
                  (char *) NULL);

    ld = iVectorAlloc( spA->iNoColumns ); 
    for(ii=0; ii < spA->iNoColumns; ii = ii + 1) {
        for (ij = spA->uMatrix.daa[ii][0]; (ij > 1) &&
             (ABS(spA->uMatrix.daa[ii][ij]) < MINREAL); ij = ij - 1);
        ld[ii] = ij;
    }

    spB = MatrixAllocSkyline( spA->cpMatrixName, DOUBLE_ARRAY,
                              spA->iNoRows, spA->iNoColumns, ld);

    for(ii = 0; ii < spA->iNoColumns; ii++)
        for(ij = 1; ij <= spB->uMatrix.daa[ii][0]; ij++)
            spB->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij];

    if( CheckUnits() == ON ) {
       for(ii = 1; ii <= spA->iNoRows; ii++)
           UnitsCopy( &(spB->spRowUnits[ii-1]),  &(spA->spRowUnits[ii-1]) );
       for(ij = 1; ij <= spA->iNoColumns; ij++)
           UnitsCopy( &(spB->spColUnits[ij-1]),  &(spA->spColUnits[ij-1]) );
    }

    MatrixFreeSkyline(spA);
    free((char *) ld);
    return (spB); 
}

/*
 *  =======================================================================
 *  MatrixAddSkyline() : Add Skyline Matrices
 * 
 *  Input  :  MATRIX *spA        -- Pointer to matrix data structure.
 *         :  MATRIX *spB        -- Pointer to matrix data structure.
 *  Output :  MATRIX *spC        -- Pointer to matrix sum
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixAddSkyline( MATRIX *spA , MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixAddSkyline( spA, spB )
MATRIX *spA, *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
int *ld, min_ldi;
int  ii, ij;
int  length, length1, length2;
DIMENSIONS *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Check for compatible matrix dimensions and data types */

       if((spA == NULL) || (spB == NULL))
          FatalError("In MatrixAddSkyline() : Either spA = NULL or spB = NULL",
                    (char *) NULL);

       if(spA->eRep != spA->eRep) 
          FatalError("In MatrixAddSkyline() : spA->eRep != spB->eRep",
                    (char *) NULL );

       if(spA->eType != spB->eType) 
          FatalError("In MatrixAddSkyline() : spA->eType != spB->eType",
                    (char *) NULL);

    /* [b] : Compute Skyline Profile for [A] + [B] */

       ld = iVectorAlloc( spA->iNoColumns ); 
       for(ii = 0; ii < spA->iNoColumns; ii++) 
           ld[ii] = (int) MAX (spA->uMatrix.daa[ii][0], spB->uMatrix.daa[ii][0]);

    /* [c] : Allocate Skyline Matrix for [C] : [C] = [A] + [B] */

       spC = MatrixAllocSkyline((char *)NULL,DOUBLE_ARRAY,spA->iNoRows, spA->iNoColumns, ld);

    /* [d] : Compute [C] = [A] + [B] */

       switch( UNITS_SWITCH) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++)
                UnitsCopy( &(spC->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]) );
            for(ij = 1; ij <= spA->iNoColumns; ij++)
                UnitsCopy( &(spC->spColUnits[ij-1]), &(spA->spColUnits[ij-1]) );

         /* [b] : Check Units and Add Matrices : [C] = [A] + [B] */

            for(ii = 0; ii < spA->iNoColumns; ii++) {
                min_ldi = (int) MIN ( spA->uMatrix.daa[ii][0], spB->uMatrix.daa[ii][0]);
                for(ij = 1; ij <= spA->iNoColumns; ij++) {

                    d1 = UnitsMult( &(spA->spRowUnits[ii]), &(spA->spColUnits[ij-1]) );
                    d2 = UnitsMult( &(spB->spRowUnits[ii]), &(spB->spColUnits[ij-1]) );

                    if(SameUnits(d1, d2) == TRUE) {
                       if( ij <= min_ldi )
                          spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij] 
                                                   + spB->uMatrix.daa[ii][ij];
                       else if( (spA->uMatrix.daa[ii][0] > min_ldi) && (ij <= ld[ii]) ) 
                          spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij];
                       else if( (spB->uMatrix.daa[ii][0] > min_ldi) && (ij <= ld[ii]) ) 
                          spC->uMatrix.daa[ii][ij] = spB->uMatrix.daa[ii][ij];
                    }
                    else {
                       printf("For row No %d, column No %d \n", ii, ij);
                       FatalError(" In MatrixAddIndirectDouble(): Inconsistent Units",
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
            for(ii = 0; ii < spA->iNoColumns; ii++) {
                min_ldi = (int) MIN ( spA->uMatrix.daa[ii][0], spB->uMatrix.daa[ii][0]);
                for(ij = 1; ij <= min_ldi; ij++)
                    spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij] +
                                               spB->uMatrix.daa[ii][ij];

                if( spA->uMatrix.daa[ii][0] > min_ldi) 
                    for(ij = min_ldi+1; ij <= ld[ii]; ij++) 
                        spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij];
                 
                if( spB->uMatrix.daa[ii][0] > min_ldi) 
                    for(ij = min_ldi+1; ij <= ld[ii]; ij++) 
                        spC->uMatrix.daa[ii][ij] = spB->uMatrix.daa[ii][ij];
            }
            break;
         default:
            break;
       }
       spC = MatrixReallocSkyline ( spC );

       free((char *) ld);
       return ( spC );
}

/*
 *  =======================================================================
 *  MatrixSubSkyline() : Subtract two Skyline Matrices
 * 
 *  Input  :  MATRIX *spA        -- Pointer to matrix data structure.
 *         :  MATRIX *spB        -- Pointer to matrix data structure.
 *  Output :  MATRIX *spC        -- Pointer to matrix differnce.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixSubSkyline( MATRIX *spA , MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixSubSkyline( spA, spB )
MATRIX *spA, *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
int *ld, min_ldi;
int  ii, ij;
int  length, length1, length2;
DIMENSIONS *d1, *d2;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();
    /* [a] : Check for compatible matrix dimensions and data types */

       if((spA == NULL) || (spB == NULL))
          FatalError("In MatrixSubSkyline() : Either spA = NULL or spB = NULL",
                    (char *) NULL);

       if(spA->eRep != spA->eRep) 
          FatalError("In MatrixSubSkyline() : spA->eRep != spB->eRep",
                    (char *) NULL );

       if(spA->eType != spB->eType) 
          FatalError("In MatrixSubSkyline() : spA->eType != spB->eType",
                    (char *) NULL);

    /* [b] : Compute Skyline Profile for [A] - [B] */

       ld = iVectorAlloc( spA->iNoRows ); 
       for(ii = 0; ii < spA->iNoRows; ii++) 
           ld[ii] = (int) MAX (spA->uMatrix.daa[ii][0], spB->uMatrix.daa[ii][0]);

    /* [c] : Allocate Skyline Matrix for [C] : [C] = [A] - [B] */

       spC = MatrixAllocSkyline((char *)NULL, DOUBLE_ARRAY,
                        spA->iNoRows, spA->iNoColumns, ld);

    /* [d] : Compute [C] = [A] - [B] */

       switch( UNITS_SWITCH) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++)
                UnitsCopy( &(spC->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]) );
            for(ij = 1; ij <= spA->iNoColumns; ij++)
                UnitsCopy( &(spC->spColUnits[ij-1]), &(spA->spColUnits[ij-1]) );

         /* [b] : Check Units and Sub Matrices : [C] = [A] - [B] */

            for(ii = 0; ii < spA->iNoRows; ii++) {
                min_ldi = (int) MIN ( spA->uMatrix.daa[ii][0], spB->uMatrix.daa[ii][0]);
                for(ij = 1; ij <= spA->iNoColumns; ij++) {

                    d1 = UnitsMult( &(spA->spRowUnits[ii]), &(spA->spColUnits[ij-1]) );
                    d2 = UnitsMult( &(spB->spRowUnits[ii]), &(spB->spColUnits[ij-1]) );

                    if(SameUnits(d1, d2) == TRUE) {
                       if( ij <= min_ldi )
                          spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij]
                                                   - spB->uMatrix.daa[ii][ij];
                       else if( (spA->uMatrix.daa[ii][0] > min_ldi) && (ij <= ld[ii]) ) 
                          spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij];
                       else if( (spB->uMatrix.daa[ii][0] > min_ldi) && (ij <= ld[ii]) ) 
                          spC->uMatrix.daa[ii][ij] = -spB->uMatrix.daa[ii][ij];
                    }
                    else {
                       printf("For row No %d, column No %d \n", ii, ij);
                       FatalError("In MatrixSubIndirectDouble(): Inconsistent Units",
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
            for(ii = 0; ii < spA->iNoRows; ii++) {
                min_ldi = (int) MIN ( spA->uMatrix.daa[ii][0], spB->uMatrix.daa[ii][0]);
                for(ij = 1; ij <= min_ldi; ij++)
                    spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij] -
                                               spB->uMatrix.daa[ii][ij];

                if( spA->uMatrix.daa[ii][0] > min_ldi) 
                    for(ij = min_ldi+1; ij <= ld[ii]; ij++) 
                        spC->uMatrix.daa[ii][ij] = spA->uMatrix.daa[ii][ij];
                 
                if( spB->uMatrix.daa[ii][0] > min_ldi) 
                    for(ij = min_ldi+1; ij <= ld[ii]; ij++) 
                        spC->uMatrix.daa[ii][ij] = -spB->uMatrix.daa[ii][ij];
            }
            break;
         default:
            break;
       }
       spC = MatrixReallocSkyline ( spC );

       free((char *) ld);
       return ( spC );
}

/*
 *  =======================================================================
 *  MatrixMultSkyline() : Multiply two Skyline Matrices.
 * 
 *  Input  :  MATRIX *spA -- Pointer to matrix data structure [A].
 *         :  MATRIX *spB -- Pointer to matrix data structure [B].
 *  Output :  MATRIX *spC -- Pointer to matrix product [C] = [A]*[B]
 * 
 *  Note : If [A] and [B] are symmetric, then [A]*[B] is symmetric only iff
 *         [A]*[B] = [B]*[A]. In general, this condition will not hold -- so
 *         we store matrix product [ma]*[mb] in INDIRECT form.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixMultSkyline( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixMultSkyline( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
int ii, ij, ik, init;
int im, in, kka, kkb;
int  length1, length2, length;
DIMENSIONS *d1, *d2, *d3, *d4, *d5;
int UNITS_SWITCH;

#define DEBUG

#ifdef DEBUG
       printf("*** Enter MatrixMultSkyline()\n");
#endif
     UNITS_SWITCH = CheckUnits();

    /* [a] : Check for compatible matrix dimensions and data types etc... */

       if((spA == NULL) || (spB == NULL))
          FatalError("In MatrixMultSkyline() : Either spA = NULL or spB = NULL",
                    (char *) NULL);

       if(spA->eRep != spB->eRep) 
          FatalError("In MatrixMultSkyline() : spA->eRep != spB->eRep",
                    (char *) NULL);

       if(spA->eType != spB->eType) 
          FatalError("In MatrixMultSkyline() : spA->eType != spB->eType",
                    (char *) NULL);

       if(spA->iNoColumns != spB->iNoRows) 
          FatalError("In MatrixMultSkyline() : spA->iNoColumns != spB->iNoRows",
                    (char *) NULL);

    /* [b] : Allocate Memory for Matrix Product */

       spC = MatrixAllocIndirect( (char *) NULL, DOUBLE_ARRAY,
                               spA->iNoRows, spB->iNoColumns); 

    /* [c] : Compute elements of matrix product */

       switch( UNITS_SWITCH) {
         case ON: 
            for(ii = 1; ii <= spA->iNoRows; ii++) {
               for(ij = 1; ij <= spB->iNoColumns; ij++) {

                  d4 = (DIMENSIONS *) MyCalloc( 1, sizeof(DIMENSIONS) );
                  ZeroUnits(d4);

                  init = (int) MAX((ii - spA->uMatrix.daa[ii-1][0]+1),
                                   (ij - spB->uMatrix.daa[ij-1][0]+1));

                  for(ik = 1; ik <= spA->iNoColumns; ik++) {

                     d1 = UnitsMult(&(spA->spRowUnits[ii-1]),&(spA->spColUnits[ik-1]));
                     d2 = UnitsMult(&(spB->spRowUnits[ik-1]),&(spB->spColUnits[ij-1]));
                     d3 = UnitsMult(d1,d2);

                     if( ik == 1  &&  d3->units_name != (char *)NULL )
                         UnitsCopy(d4, d3);

                     if( ik >= init ) {
                        if( SameUnits(d3,d4) == TRUE ) {

                           im  = (int) MAX(ik,ii);
                           kka = (int) MIN(ik,ii);

                           in  = (int) MAX(ik,ij);
                           kkb = (int) MIN(ik,ij);

                           if((im-kka+1) <= spA->uMatrix.daa[ im-1 ][0] &&
                              (in-kkb+1) <= spB->uMatrix.daa[ in-1 ][0] )
                               spC->uMatrix.daa[ ii-1 ][ ij-1 ] +=
                                  spA->uMatrix.daa[im-1][im-kka+1]
                                 *spB->uMatrix.daa[in-1][in-kkb+1];
                        }
                        else
                          FatalError("In MatrixMultSkyline(): Inconsistent Units",(char *)NULL);
                     }
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
               UnitsMultRep( &(spC->spRowUnits[ii-1]), &(spB->spRowUnits[0]),&(spA->spRowUnits[ii-1]));

           for(ij = 1; ij <= spC->iNoColumns; ij++)
               UnitsMultRep( &(spC->spColUnits[ij-1]), &(spA->spColUnits[0]),&(spB->spColUnits[ij-1]));
           break;
        case OFF:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
               for(ij = 1; ij <= spB->iNoColumns; ij++) {
                  init = (int) MAX((ii - spA->uMatrix.daa[ii-1][0]+1),
                                   (ij - spB->uMatrix.daa[ij-1][0]+1));
                  for(ik = init; ik <= spA->iNoColumns; ik++) {
                     im  = (int) MAX(ik,ii);
                     kka = (int) MIN(ik,ii);

                     in  = (int) MAX(ik,ij);
                     kkb = (int) MIN(ik,ij);

                     if((im-kka+1) <= spA->uMatrix.daa[im-1][0] &&
                        (in-kkb+1) <= spB->uMatrix.daa[in-1][0] )
                         spC->uMatrix.daa[ ii-1 ][ ij-1 ] +=
                         spA->uMatrix.daa[im-1][im-kka+1]*spB->uMatrix.daa[in-1][in-kkb+1];
                  }
               }
            }
           break;
        default:
           break;
      }

#ifdef DEBUG
       printf("*** Leaving MatrixMultSkyline()\n");
#endif

      return (spC);
}
#undef DEBUG

/*
 *  =======================================================================
 *  MatrixNegateSkyline() : Matrix Negation [B] = -[A].
 *  
 *  Input  :  MATRIX     *spA  -- Pointer to matrix A.
 *  Output :  MATRIX     *spB  -- Pointer to matrix B.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixNegateSkyline( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixNegateSkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX *spB;
int  ii, ij;
int     *ld;
int  length;

    /* [a] : Check Input matrix [A] */

       if(spA == (MATRIX *) NULL) 
          FatalError("In MatrixNegateSkyline() : [A] == NULL", (char *) NULL);

    /* [b] : Negate Skyline Matrix */
    
       ld = iVectorAlloc ( spA->iNoRows); 
       for(ii = 1; ii <= spA->iNoRows; ii++)
           ld[ii-1] = spA->uMatrix.daa[ii-1][0];
    
       spB = MatrixAllocSkyline( (char *)NULL, DOUBLE_ARRAY,
                                    spA->iNoRows, spA->iNoColumns, ld);

       for(ii = 1; ii <= spA->iNoRows; ii++)
           for(ij = 1; ij <= ld[ii-1]; ij++)
               spB->uMatrix.daa[ii-1][ij] = -spA->uMatrix.daa[ii-1][ij];

       if( CheckUnits() == ON ) {
          for(ii = 1; ii <= spA->iNoRows; ii++)
              UnitsCopy(&(spB->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]));
          for(ij = 1; ij <= spA->iNoColumns; ij++)
              UnitsCopy(&(spB->spColUnits[ij-1]), &(spA->spColUnits[ij-1]));
       }
           
     free((char *) ld);
     return (spB);
}

#ifdef __STDC__
MATRIX *MatrixNegateReplaceSkyline( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixNegateReplaceSkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
int  ii, ij;

       for(ii = 1; ii <= spA->iNoRows; ii++)
           for(ij = 1; ij <= spA->uMatrix.daa[ii-1][0]; ij++)
               spA->uMatrix.daa[ii-1][ij] = -spA->uMatrix.daa[ii-1][ij];

     return (spA);
}

#ifdef __STDC__
MATRIX *MatrixTransposeSkyline( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixTransposeSkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX  *spB;
int   ii, ij;
int      *ld;
int   length;

    /* [a] : Transpose Matrix */

       ld = iVectorAlloc ( spA->iNoColumns); 
       for(ii = 1; ii <= spA->iNoColumns; ii++)
           ld[ii-1] = spA->uMatrix.daa[ii-1][0];
    
       spB = MatrixAllocSkyline( (char *)NULL, DOUBLE_ARRAY,
                                    spA->iNoRows, spA->iNoColumns, ld);

       for(ii = 1; ii <= spA->iNoRows; ii++)
           for(ij = 1; ij <= ld[ii-1]; ij++)
              spB->uMatrix.daa[ii-1][ij] = spA->uMatrix.daa[ii-1][ij];

       if( CheckUnits() == ON ) {
          for(ii = 1; ii <= spA->iNoRows; ii++)
              UnitsCopy(&(spB->spColUnits[ii-1]), &(spA->spRowUnits[ii-1]));
          for(ij = 1; ij <= spA->iNoColumns; ij++)
              UnitsCopy(&(spB->spRowUnits[ij-1]), &(spA->spColUnits[ij-1]));
       }
      
       free((char *) ld);
       return (spB);
}


/*
 *  ===============================================================================
 *  MatrixMultSkylineIndirectDouble() : Compute Product of Skyline-Indirect Matrix.
 * 
 *  Input  :  MATRIX *spA -- Pointer to matrix data structure [A].
 *         :  MATRIX *spB -- Pointer to matrix data structure [B].
 *  Output :  MATRIX *spC -- Pointer to matrix product [C] = [A]*[B]
 * 
 *  Note : If [A] and [B] are symmetric, then [A]*[B] is symmetric only iff
 *         [A]*[B] = [B]*[A]. In general, this condition will not hold -- so
 *         we store matrix product [ma]*[mb] in INDIRECT form.
 *  ===============================================================================
 */

#ifdef __STDC__
MATRIX *MatrixMultSkylineIndirectDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixMultSkylineIndirectDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
double       dAvalue;
int ii, ij, ik, init;
int   im,   in;
int  length1, length2, length;
DIMENSIONS *d1, *d2, *d3, *d4, *d5;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();

    /* [a] : Check for compatible matrix dimensions and data types etc... */

       if((spA == NULL) || (spB == NULL))
          FatalError("In MatrixMultSkyline() : Either spA = NULL or spB = NULL",
                    (char *) NULL);

       if(spA->iNoColumns != spB->iNoRows) 
          FatalError("In MatrixMultSkyline() : spA->iNoColumns != spB->iNoRows",
                    (char *) NULL);

    /* [b] : Allocate Memory for Matrix Product */

       spC = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, spA->iNoRows, spB->iNoColumns); 

    /* [c] : Compute elements of matrix product */

       switch( UNITS_SWITCH == ON ) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
               for(ij = 1; ij <= spB->iNoColumns; ij++) {
                  d4 = (DIMENSIONS *) MyCalloc( 1, sizeof(DIMENSIONS) );
                  ZeroUnits(d4);

                  init = ii - spA->uMatrix.daa[ii-1][0] + 1;

                  for(ik = 1; ik <= spA->iNoColumns; ik++) {

                     d1 = UnitsMult( &(spA->spRowUnits[ii-1]), &(spA->spColUnits[ik-1]) );
                     d2 = UnitsMult( &(spB->spRowUnits[ik-1]), &(spB->spColUnits[ij-1]) );
                     d3 = UnitsMult(d1,d2);

                     if( ik==1 && d3->units_name!=(char *)NULL )
                        UnitsCopy(d4, d3);

                     if( ik >= init ) {
                         if( SameUnits(d3, d4) == TRUE) {

                           im = (int) MIN( ii, ik );
                           in = (int) MAX( ii, ik );

                           if((in-im+1) <= spA->uMatrix.daa[in-1][0])
                               dAvalue = spA->uMatrix.daa[ in-1 ][ in-im+1 ];
                           else
                               dAvalue = 0.0;

                           spC->uMatrix.daa[ ii-1 ][ ij-1 ] +=
                                dAvalue * spB->uMatrix.daa[ ik - 1 ][ ij - 1 ];
                         } 
                         else
                           FatalError(" In MatrixMultSkylineIndirectDouble(): Inconsistent Units",(char *)NULL);
                     }
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
               UnitsMultRep( &(spC->spRowUnits[ii-1]), &(spB->spRowUnits[0]), &(spA->spRowUnits[ii-1]) );
           for(ij = 1; ij <= spC->iNoColumns; ij++)
               UnitsMultRep( &(spC->spColUnits[ij-1]), &(spA->spColUnits[0]), &(spB->spColUnits[ij-1]) );

           break;
        case OFF:
           for( ii=1; ii <= spA->iNoRows; ii++) {
           for( ij=1; ij <= spB->iNoColumns; ij++) {

              init = ii - spA->uMatrix.daa[ii-1][0] + 1;
              for(ik = init; ik <= spA->iNoColumns; ik++) {
        
                  im = MIN( ii, ik );
                  in = MAX( ii, ik );

                  if((in-im+1) <= spA->uMatrix.daa[in-1][0])
                      dAvalue = spA->uMatrix.daa[ in-1 ][ in-im+1 ];
                  else
                      dAvalue = 0.0;

                  spC->uMatrix.daa[ ii-1 ][ ij-1 ] +=
                       dAvalue * spB->uMatrix.daa[ ik - 1 ][ ij - 1 ];
            
             }
          }
          }
           break;
        default:
           break;
      }

#ifdef DEBUG
       printf("*** Leaving MatrixMultSkylineIndirectDouble()\n");
#endif

      return (spC);
}

#ifdef __STDC__
MATRIX *MatrixMultIndirectSkylineDouble( MATRIX *spA, MATRIX *spB )
#else  /* Start case not STDC */
MATRIX *MatrixMultIndirectSkylineDouble( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
MATRIX *spC;
double       dBvalue;
int ii, ij, ik, init;
int   im,   in;
int  length1, length2, length;
DIMENSIONS *d1, *d2, *d3, *d4, *d5;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();

    /* [a] : Check for compatible matrix dimensions and data types etc... */

       if((spA == NULL) || (spB == NULL))
          FatalError("In MatrixMultSkylineIndirect() : Either spA = NULL or spB = NULL",
                    (char *) NULL);

       if(spA->iNoColumns != spB->iNoRows) 
          FatalError("In MatrixMultSkylineIndirect() : spA->iNoColumns != spB->iNoRows",
                    (char *) NULL);

    /* [b] : Allocate Memory for Matrix Product */

       spC = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, spA->iNoRows, spB->iNoColumns); 

    /* [c] : Compute elements of matrix product */

       switch( UNITS_SWITCH == ON ) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
               for(ij = 1; ij <= spB->iNoColumns; ij++) {
                  d4 = (DIMENSIONS *) MyCalloc( 1, sizeof(DIMENSIONS) );
                  d4 = ZeroUnits(d4);

                  init = ii - spB->uMatrix.daa[ii-1][0] + 1;

                  for(ik = 1; ik <= spA->iNoColumns; ik++) {

                     d1 = UnitsMult(&(spA->spRowUnits[ii-1]),&(spA->spColUnits[ik-1]));
                     d2 = UnitsMult(&(spB->spRowUnits[ik-1]),&(spB->spColUnits[ij-1]));
                     d3 = UnitsMult(d1,d2);

                     if( ik == 1  &&  d3->units_name != (char *)NULL )
                         UnitsCopy(d4, d3);

                     if( ik >= init ) {
                         if( SameUnits(d3, d4) == TRUE) {

                           im = (int) MIN( ik, ij );
                           in = (int) MAX( ik, ij );

                           if((in-im+1) <= spB->uMatrix.daa[in-1][0])
                               dBvalue = spB->uMatrix.daa[ in-1 ][ in-im+1 ];
                           else
                               dBvalue = 0.0;

                           spC->uMatrix.daa[ ii-1 ][ ij-1 ] +=
                                spA->uMatrix.daa[ ii - 1 ][ ik - 1 ]* dBvalue;
                         } 
                         else
                           FatalError(" In MatrixMultIndirectSkylineDouble(): Inconsistent Units",(char *)NULL);
                     }
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

           break;
        case OFF:
            for( ii=1; ii <= spA->iNoRows; ii++) {
            for( ij=1; ij <= spB->iNoColumns; ij++) {

               init = ii - spB->uMatrix.daa[ii-1][0] + 1;
               for(ik = init; ik <= spB->iNoRows; ik++) {
        
                   im = MIN( ik, ij );
                   in = MAX( ik, ij );

                   if((in-im+1) <= spB->uMatrix.daa[in-1][0])
                       dBvalue = spB->uMatrix.daa[ in-1 ][ in-im+1 ];
                   else
                       dBvalue = 0.0;

                   spC->uMatrix.daa[ ii-1 ][ ij-1 ] +=
                        spA->uMatrix.daa[ ii - 1 ][ ik - 1 ]* dBvalue;
            
              }
           }
           }
           break;
        default:
           break;
      }

#ifdef DEBUG
       printf("*** Leaving MatrixMultIndirectSkylineDouble()\n");
#endif

      return (spC);
} 

/*
 *  =======================================================================
 *  MatrixCopySkyline() : Copy Skyline Matrices
 * 
 *  Input  :  MATRIX *spA        -- Pointer to matrix data structure.
 *  Output :  MATRIX *spCopy     -- Pointer to matrix copy
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixCopySkyline( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixCopySkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX   *spCopy;
int  *ld, ii, ij;
int       length;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();

    /* [a] : Check input for [ma] */

       if(spA == NULL)
          FatalError("In MatrixCopySkyline() : [A] = NULL",
                    (char *) NULL);

    /* [b] : Copy Matrix [spCopy] = [spA] */

       ld = iVectorAlloc ( spA->iNoColumns); 
       for(ii = 1; ii <= spA->iNoColumns; ii++)
           ld[ii-1] = spA->uMatrix.daa[ii-1][0];

       spCopy = MatrixAllocSkyline( spA->cpMatrixName, DOUBLE_ARRAY,
                                    spA->iNoRows, spA->iNoColumns, ld);

       switch( UNITS_SWITCH == ON ) {
         case ON:
            for(ii = 1; ii <= spA->iNoRows; ii++) {
                UnitsCopy(&(spCopy->spRowUnits[ii-1]), &(spA->spRowUnits[ii-1]));
                for(ij = 1; ij <= spA->iNoColumns; ij++) {
                    if( ij <= ld[ii-1] )
                        spCopy->uMatrix.daa[ii-1][ij] = spA->uMatrix.daa[ii-1][ij];
                }
            }
            for(ij = 1; ij <= spA->iNoColumns; ij++)
                UnitsCopy(&(spCopy->spColUnits[ij-1]), &(spA->spColUnits[ij-1]));
            break;
        case OFF:
            for(ii = 1; ii <= spA->iNoRows; ii++)
                for(ij = 1; ij <= ld[ii-1]; ij++)
                   spCopy->uMatrix.daa[ii-1][ij] = spA->uMatrix.daa[ii-1][ij];
           break;
        default:
           break;
      }

       free((char *) ld);
       return ( spCopy ); 
}


/* 
 *  ===========================================================
 *  Functions to Cast Matrix Storage Patterns 
 * 
 *     [a] : Indirect matrix to skyline matrix.
 *     [b] : Skyline matrix to indirect matrix.
 * 
 *  Input  :  MATRIX *spA  -- Pointer to matrix data structure.
 *  Output :  MATRIX *spB  -- Pointer to recast matrix.
 *  ===========================================================
 */

#ifdef __STDC__
MATRIX *MatrixIndirectToSkyline( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixIndirectToSkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
double **dpB;
int      *ld;
int   ii, ij;
static double REALTOL = 1E-7;

    /* [a] : Check input matrix [A] is defined and symmetric */

       if(spA == NULL)
          FatalError("In MatrixIndirectToSkyline() : spA == NULL",
                    (char *) NULL);

       for(ii = 1; ii <= spA->iNoColumns-1; ii++) {
          for(ij = 0; ij <= ii-1; ij++){
            if(ABS( spA->uMatrix.daa[ij][ii] - spA->uMatrix.daa[ii][ij]) > REALTOL) {
                printf("In MatrixIndirectToSkyline() : %s[%d][%d]-%s[%d][%d] = %le\n",
                       spA->cpMatrixName, ij, ii, spA->cpMatrixName,ii, ij, spA->uMatrix.daa[ij][ii] - spA->uMatrix.daa[ii][ij]);
                printf("In MatrixIndirectToSkyline() : [%s] is not symmetric\n", spA->cpMatrixName);
                FatalError("In MatrixIndirectToSkyline() :", (char *) NULL);
            }
          }
       }

    /* [c] : Compute required Ld[] matrix */

       ld = iVectorAlloc( spA->iNoColumns );
       for(ii = 1; ii <= spA->iNoColumns; ii++) {
       for(ij = 1; (ij < ii)&&(ABS(spA->uMatrix.daa[ij-1][ii-1]) < MINREAL); ij++);
           ld[ii-1] = ii-ij+1;
       }

       switch((int) spA->eType) {
          case DOUBLE_ARRAY:
               dpB = (double **)MyCalloc(spA->iNoRows,sizeof(double*));
               for(ii = 1; ii <= spA->iNoRows; ii++ ) {
                   dpB[ii-1] = (double *) MyCalloc((ld[ii-1]+1), sizeof(double));
                   dpB[ii-1][0] = ld[ii-1];
                   for(ij = 1; ij <= ld[ii-1]; ij++)
                       dpB[ii-1][ij] = spA->uMatrix.daa[ii-ij][ii-1];
               }
               for(ii = 1; ii <= spA->iNoRows; ii++)
                   free ((char *) (spA->uMatrix.daa[ ii-1 ]));
               free ((char *) (spA->uMatrix.daa));
               break;
          default:
               FatalError("In MatrixIndirectToSkyline() : Undefined spA->eType", (char *) NULL);
               break;
      }
      spA->uMatrix.daa = dpB;
      spA->eRep = SKYLINE;
      free((char *) ld);

       return ( spA );
}

/*
 *  ===================================================================
 *  MatrixSkylineToIndirect() : Cast Skyline matrix into Indirect Form.
 *  ===================================================================
 */

#ifdef __STDC__
MATRIX *MatrixSkylineToIndirect( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixSkylineToIndirect( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
int  ii, ij, ik;
double    **dpB;

    dpB = MatrixAllocIndirectDouble( spA->iNoRows, spA->iNoColumns );
    for( ii = 1; ii <= spA->iNoColumns; ii++) 
        for( ij = 1; ij <= spA->uMatrix.daa[ii][0]; ij++)
             ik = ABS(ii-ij)+1;
             dpB[ii-1][ik-1] = dpB[ik-1][ii-1] = spA->uMatrix.daa[ii-1][ij];

    for( ii = 1; ii <= spA->iNoColumns; ii++) 
        free((char *) spA->uMatrix.daa[ii-1]);
    free(spA->uMatrix.daa);

    spA->uMatrix.daa = dpB;
    spA->eRep = INDIRECT;

    return (spA);
}

/*
 *  =======================================================================
 *  MatrixScaleSkyline() : Multiply Matrix by Scalar
 * 
 *  Input  :  MATRIX *spA        -- Pointer to matrix data structure [A].
 *         :  double scale       -- double : scale factor c.
 *  Output :  MATRIX *spB        -- Pointer to scaled matrix [B] = c.[A].
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *MatrixScaleSkyline( MATRIX *spA, double dScale)
#else  /* Start case not STDC */
MATRIX *MatrixScaleSkyline( spA, dScale)
MATRIX   *spA;
double dScale;
#endif /* End case not STDC */
{
MATRIX     *spB;
int *ld, ii, ij;

    /* [a] : Check Input matrix [A] */

       if(spA == (MATRIX *) NULL) 
          FatalError("In MatrixScaleToSkyline() : [A] == NULL",
                    (char *) NULL);

    /* [b] : Scale Skyline Matrix */

       spB = MatrixCopySkyline( spA );

       for( ii = 1; ii <= spA->iNoRows; ii++) 
           for( ij = 1; ij <= spA->uMatrix.daa[ii-1][0]; ij++)
                spB->uMatrix.daa[ii-1][ij] = dScale * spA->uMatrix.daa[ii-1][ij];

       spB = MatrixReallocSkyline ( spB);

       return ( spB);
}

#ifdef __STDC__
double MatrixContentScaleSkyline(MATRIX *spA, int row_no, int col_no )
#else
double MatrixContentScaleSkyline(spA, row_no, col_no )
MATRIX  *spA;
int row_no;        /* row number    */
int col_no;        /* column number */
#endif
{
int iMin, iMax;
int     ii, ij;
double      da;

#ifdef DEBUG
       printf("*** Enter MatrixContentScale() : spA->iNoRows    = %4d\n", spA->iNoRows);
       printf("                        : spA->iNoColumns = %4d\n", spA->iNoColumns);
#endif
      if(CheckUnits()==OFF)
         FatalError("You have to set units ON to use this function","In MatrixScale",(char *)NULL );
 
      iMin = MIN( row_no, col_no );
      iMax = MAX( row_no, col_no );
      ii = iMax - 1;
      ij = iMax - iMin + 1;

     /* [a] Scale for Column of matrix spA */

      if( ij <= spA->uMatrix.daa[ii][0] )
            da = spA->uMatrix.daa[ii][ij];
      else
            da = 0.0;

      if(spA->spColUnits[col_no-1].scale_factor != 0) 
         da = da / spA->spColUnits[col_no-1].scale_factor;
      else {
         printf("==> for column %d, scale_factor of %s = 0 \n",
               col_no, spA->spColUnits[col_no-1].units_name);
         FatalError("Fatal error in MatrixContentScaleSkyline(): ",(char *)NULL);
      }

     /* [b] Scale for Row of matrix spA */

      if(spA->spRowUnits[row_no-1].scale_factor != 0) 
         da = da / spA->spRowUnits[row_no-1].scale_factor;
      else {
         printf("==> for Row %d, scale_factor of %s = 0",
                row_no, spA->spColUnits[row_no-1].units_name);
         FatalError("Fatal error in MatrixScaleSkyline(): ",(char *)NULL);
      }
      return (da);

#ifdef DEBUG
       printf("*** Leave MatrixContentScaleSkyline()\n");
#endif
}
                      
/*
 *  =======================================================================
 *  LUDecompositionSkyline() : [L][D][L]^T Factorization of Skyline Matrix
 * 
 *  Input  :  MATRIX *spA   -- Pointer to skyline matrix [A].
 *  Output :  MATRIX *spLU  -- Pointer to decomposed matrix [L][D][L]^T.
 *  =======================================================================
 */

#ifdef __STDC__
MATRIX *LUDecompositionSkyline( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *LUDecompositionSkyline( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX   *spLU;
int ii, ij, ik;
int    min_col;

    /* [a] : Check Input matrix [A] */

       if(spA == NULL)
          FatalError("In LUDecompositionSkyline() : Pointer to matrix is NULL",
                    (char *) NULL);

    /* [b] : Compute [L][D][L]^T decomposition */

       spLU = MatrixCopySkyline( spA ); 
       switch((int) spLU->eType ) {
           case DOUBLE_ARRAY:
                for(ii = 0; ii < spLU->iNoColumns; ii++) {
                    for(ij = spLU->uMatrix.daa[ii][0]; ij > 1; ij = ij - 1) {

                        min_col = (int) MIN((spLU->uMatrix.daa[ii][0]),
                                      (spLU->uMatrix.daa[ii-ij+1][0]+ij-1));        

                        for(ik=ij+1; ik <= min_col; ik++) 
                            spLU->uMatrix.daa[ii][ij] -=
                                spLU->uMatrix.daa[ii][ik] *
                                spLU->uMatrix.daa[ii-ij+1][ik-ij+1] *
                                spLU->uMatrix.daa[ii-ik+1][1];

                        spLU->uMatrix.daa[ii][ij] /= spLU->uMatrix.daa[ii-ij+1][1];
                    }

                    for(ij = spLU->uMatrix.daa[ii][0]; ij > 1; ij--) 
                        spLU->uMatrix.daa[ii][1] -=
                            spLU->uMatrix.daa[ii][ij] *
                            spLU->uMatrix.daa[ii][ij] *
                            spLU->uMatrix.daa[ii-ij+1][1];

                    if(ABS(spLU->uMatrix.daa[ii][1]) <= MINREAL) {
                       printf(" In LUDecompositionSkyline() : %s[%d][1] = %le \n", spLU->cpMatrixName,ii+1,
                                spLU->uMatrix.daa[ii][1]);
                       FatalError("In LUDecompositionSkyline() : Singular Matrix ",
                                 (char *) NULL);
                    }
                }
                break;
           default:
                FatalError("In LUDecompositionSkyline() : Undefined Matrix Data Type",
                          (char *) NULL);
        }
    
    return ( spLU );
}

/*
 *  ==================================================================================
 *  LUBacksubstitutionSkyline() : Forward/Back Substitution of Factored Skyline Matrix
 * 
 *  Input  :  MATRIX *spA -- Pointer to matrix data structure.
 *  Output :  MATRIX *spX -- Pointer to matrix solution
 *  ==================================================================================
 */

#ifdef __STDC__
MATRIX *LUBacksubstitutionSkyline( MATRIX *spA, MATRIX *spB)
#else  /* Start case not STDC */
MATRIX *LUBacksubstitutionSkyline( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* End case not STDC */
{
int            ii, ij, ik;
DIMENSIONS        *d1,*d2;
MATRIX            *spDisp;


#ifdef MYDEBUG
       printf("*** Enter LUBacksubstitutionSkyline()\n");
#endif

    if(spA == NULL) {
       FatalError("In LUBacksubstitutionSkyline() : spA == NULL",
                 (char *) NULL);
    }

    if(spA->iNoColumns != spB->iNoRows) {
       FatalError("In LUBacksubstitutionSkyline() : Dimensions of",
                  "A[][] and b[] are incompatible                ",
                 (char *) NULL);
    }

    spDisp = MatrixCopyIndirectDouble(spB);

    switch((int) spA->eType ) {
        case DOUBLE_ARRAY:

           /* [0] : Calculate the units of displacements */

         if( CheckUnits() == ON ) {
           for (ii = 1; ii <= spDisp->iNoRows; ii++) {
              d1 = UnitsMult( &(spB->spRowUnits[ii-1]), &(spB->spColUnits[0]) );
              d2 = UnitsMult( &(spA->spRowUnits[ii-1]), &(spA->spColUnits[ii-1]) );
              UnitsDivRep(&(spDisp->spRowUnits[ii-1]), d1, d2, YES); 

              free((char *) d1->units_name);
              free((char *) d1);
              free((char *) d2->units_name);
              free((char *) d2);
           }
              ZeroUnits(&(spDisp->spColUnits[0]));

         }
           /* [a] : Forward Substitution */

           for(ij=0; ij < spDisp->iNoColumns; ij++)
               for(ii = 0; ii < spDisp->iNoRows; ii++)
                   for(ik = spA->uMatrix.daa[ii][0]; ik > 1; ik--) 
                       spDisp->uMatrix.daa[ii][ij] -= spDisp->uMatrix.daa[ii-ik+1][ij] *
                                                   spA->uMatrix.daa[ii][ik];

           for(ij = 0; ij < spDisp->iNoColumns; ij++)
               for(ii=0; ii < spDisp->iNoRows; ii++)
                   spDisp->uMatrix.daa[ii][ij] /= spA->uMatrix.daa[ii][1];
   
           /* [b] : Back Substitution */

           for(ij = 0; ij < spDisp->iNoColumns; ij++)
               for(ii = spDisp->iNoRows-1; ii >= 0; ii--)
                   for(ik = ii+1; ik <= spDisp->iNoRows-1; ik++)
                       if((ik+1-ii) <= spA->uMatrix.daa[ik][0] )
                           spDisp->uMatrix.daa[ii][ij] -= spDisp->uMatrix.daa[ik][ij] *
                                                       spA->uMatrix.daa[ik][ik+1-ii];

           break;
        default:
             FatalError("In LUBacksubstitutionSkyline() : Undefined Matrix Data Type",
                       (char *) NULL);
    }

#ifdef MYDEBUG
       printf("*** Leave LUBacksubstitutionSkyline()\n");
#endif

    return ( spDisp );
}

/*
 *  ===================================================================
 *  MatrixInverseSkyline() : Compute inverse of Skyline Double Matrix
 * 
 *  Input  :  MATRIX *spA -- Pointer to matrix data structure [A].
 *  Output :  MATRIX *spB -- Pointer to matrix inverse [A]^{-1}.
 *  ===================================================================
 */

#ifdef __STDC__
MATRIX *MatrixInverseSkyline ( MATRIX *spA )
#else  /* Start case not STDC */
MATRIX *MatrixInverseSkyline ( spA )
MATRIX *spA;
#endif /* End case not STDC */
{
MATRIX     *spLU;
MATRIX        *F;
MATRIX   *spAinv;
int       ii, ij;
int UNITS_SWITCH;

     UNITS_SWITCH = CheckUnits();

    /* [a] : Check Input */

       if( spA == NULL )
           FatalError("In MatrixInverseSkyline () : [ma] == NULL \n",
                     (char *) NULL);

    /* [b] : Compute LU Decomposition and column vectors of inverse */

       F = MatrixAllocIndirect( spA->cpMatrixName, DOUBLE_ARRAY,
                                spA->iNoRows, spA->iNoColumns);

       for(ii = 0; ii < spA->iNoRows; ii++) {
           for(ij = 0; ij < spA->iNoColumns; ij++) {
               if(ii == ij)
                  F->uMatrix.daa[ii][ij] = 1.0;
               else 
                  F->uMatrix.daa[ii][ij] = 0.0;
           }
       }

       if( UNITS_SWITCH == ON ) {
          for(ii = 1; ii <= spA->iNoRows; ii++)
              ZeroUnits(&(F->spRowUnits[ii-1]));
          for(ii = 1; ii <= spA->iNoColumns; ii++)
              ZeroUnits(&(F->spColUnits[ii-1]));
       }

       spLU   = LUDecompositionSkyline(spA);
       spAinv = LUBacksubstitutionSkyline( spLU, F );

       spAinv = MatrixIndirectToSkyline(spAinv);
       spAinv = MatrixReallocSkyline(spAinv);

       /* Assign units to spAinv */

       if( UNITS_SWITCH == ON ) {
          for(ii = 1; ii <= spA->iNoRows; ii++){
              UnitsPowerRep( &(spAinv->spColUnits[ii-1]), &(spA->spRowUnits[ii-1]), -1.0, YES );
              UnitsPowerRep( &(spAinv->spRowUnits[ii-1]), &(spA->spColUnits[ii-1]), -1.0, YES );
          }
       }

       MatrixFreeSkyline(spLU);
       MatrixFreeIndirect(F);

       return (spAinv);
} 


/*
 *  =================================================================
 *  Subspace() : Solve [A][x] = [Lambda][B][x] by Subspace Iteration.
 *             : Assume [A] and [B] are stored in Skyline Form, and 
 *             : that we will compute "p = iNoVectors" eigenvectors.
 * 
 *             : [A] and [B] are large (nxn) matrices
 *             : Size of [X_{k}] = [X_{k+1}]   = (nxm) matrix
 *             : Size of [Y_{k}] = [Y_{k+1}]   = (nxm) matrix
 *             : Size of [A_{k+1}] = [B_{k+1}] = (mxm) matrix
 *             : Size of [Eigenvalue_{k+1}] = (mxm) diagonal matrix
 *             : Size of [Q_{k+1}]          = (mxm) matrix
 * 
 *  Note : [A] and [B] are stored in skyline form. All other working  
 *         matrices have INDIRECT storage patterns.
 * 
 *  Input :  spEigenvalue  = Allocated Vector for Eigenvalues
 *           spEigenvector = Allocated Vector for Eigenvectors
 *           iNoVectors    = No Eigenvectors to be computed. 
 *  Output : spEigenvalue  = Allocated Vector for Eigenvalues
 *           spEigenvector = Allocated Vector for Eigenvectors
 * 
 *  Written By : M. Austin                             November 1993. 
 *  =================================================================
 */

enum { MAX_SUBSPACE_ITERATIONS = 50 };

#define TOLERANCE  0.01

Subspace( spA, spB, spEigenvalue, spEigenvector , iNoEigen)
MATRIX   *spA,    *spB;
MATRIX   *spEigenvalue;
MATRIX  *spEigenvector;
int           iNoEigen;
{
MATRIX    *spEigenvectorQ;
MATRIX *spAwork, *spBwork;
MATRIX           *spTemp1;
MATRIX            *spA_LU;
MATRIX               *spY;
MATRIX            *spYnew;
MATRIX   *spEigenvalueOld;
MATRIX               *spV;
MATRIX              *spVk;
int     iSize, ii, ij, ik;
int   iMin, iMax,   iLoop;
int  iSubspaceConvergence;
double dEigenvalueOld;
double dConvergence;
double dEigenvalue;
double dMaxValue;

       printf("\n*** Enter Subspace() : Compute %3d Eigenvalues and Eigenvectors\n",
              iNoEigen);

    /* [a] : Check Input and Allocate Working Matrices */

       if( spA == NULL || spB == NULL )
           FatalError("In Subpace() : Pointer spA == NULL or spB == NULL",
                     (char *) NULL);

       if( spA->eType != DOUBLE_ARRAY || spB->eType != DOUBLE_ARRAY ) 
           FatalError("In Subspace() : spA->eType != DOUBLE_ARRAY or spB->eType != DOUBLE_ARRAY",
                     (char *) NULL);

       if( spA->eRep != SKYLINE || spB->eRep != SKYLINE ) 
           FatalError("In Subspace() : spA->eRep != SKYLINE or spB->eRep != SKYLINE",
                     (char *) NULL);

    /* [b] : Initialize Starting Eigenvectors [X_{k}] */

       iSize    = spA->iNoRows;
       for(ii = 1; ii <= iNoEigen; ii = ii + 1)
           spEigenvector->uMatrix.daa[ii-1][ii-1] = 1.0;

    /* [c] : Loops for Subspace Iteration */

       spAwork        = MatrixAllocIndirect("[A_{k+1}]", DOUBLE_ARRAY,
                                             iNoEigen, iNoEigen);
       spBwork        = MatrixAllocIndirect("[B_{k+1}]", DOUBLE_ARRAY,
                                             iNoEigen, iNoEigen);

       spEigenvectorQ = MatrixAllocIndirect("Eigenvector [Q]", DOUBLE_ARRAY,
                                             iNoEigen, iNoEigen);

       spVk   = MatrixAllocIndirect("[V_{k}]",   DOUBLE_ARRAY, iSize, iNoEigen);
       spY    = MatrixAllocIndirect("[Y_{k}]",   DOUBLE_ARRAY, iSize, iNoEigen);
       spYnew = MatrixAllocIndirect("[Y_{k+1}]", DOUBLE_ARRAY, iSize, iNoEigen);

       spEigenvalueOld = MatrixCopy( spEigenvalue );

       iLoop = 0; iSubspaceConvergence = FALSE;
       while( iLoop <= (int) MAX_SUBSPACE_ITERATIONS &&
              iSubspaceConvergence == FALSE ) {
           iLoop = iLoop + 1;

           printf("*** In Subspace() : Start Iteration %3d \n", iLoop);

           /* [c.1] : Compute [Y_{k}] = [B]*[X_{k}] */

           for(ii = 1; ii <= iSize; ii = ii + 1) 
           for(ij = 1; ij <= iNoEigen; ij = ij + 1) {
               spY->uMatrix.daa[ ii-1 ][ ij-1 ] = 0.0;
               for(ik = 1; ik <= iSize; ik = ik + 1)  {

                   iMin = (int) MIN( ii, ik );
                   iMax = (int) MAX( ii, ik );

                   if((iMax-iMin+1) <= spB->uMatrix.daa[iMax-1][0]) {
                       spY->uMatrix.daa[ ii-1 ][ ij-1 ] += 
                            spB->uMatrix.daa[ iMax-1 ][ iMax-iMin+1 ] *
                            spEigenvector->uMatrix.daa[ ik-1 ][ ij-1 ] ; 
                   }
               }
           }

           /* [c.2] : Solve [A][V_{k+1}] = [Y_{k}] */

           if(iLoop == 1) {
              spA_LU  = LUDecompositionSkyline(spA);
           }

           for(ii = 1; ii <= iSize; ii = ii + 1) 
           for(ij = 1; ij <= iNoEigen; ij = ij + 1)
               spVk->uMatrix.daa[ ii-1 ][ ij-1 ] = spY->uMatrix.daa[ ii-1 ][ ij-1 ];

           if(iLoop != 1) {
              MatrixFree (spTemp1);
              MatrixFree (spV);
           }
           spV = LUBacksubstitutionSkyline( spA_LU , spVk );

           /* [c.3] : Compute [Y_{new}] = [A] . [V_{k+1}] */

           for(ii = 1; ii <= iSize; ii = ii + 1) 
           for(ij = 1; ij <= iNoEigen; ij = ij + 1) {
               spYnew->uMatrix.daa[ ii-1 ][ ij-1 ] = 0; 
               for(ik = 1; ik <= iSize; ik = ik + 1)  {

                   iMin = (int) MIN( ii, ik );
                   iMax = (int) MAX( ii, ik );

                   if((iMax-iMin+1) <= spA->uMatrix.daa[iMax-1][0]) {
                       spYnew->uMatrix.daa[ ii-1 ][ ij-1 ] += 
                               spA->uMatrix.daa[ iMax-1 ][ iMax-iMin+1 ] *
                               spV->uMatrix.daa[ ik-1 ][ ij-1 ] ; 
                   }
               }
           }

           /* [c.4] : Compute [A_{k+1}] = [V_{k+1}]^T . [Y_{new}] */

           for(ii = 1; ii <= iNoEigen; ii = ii + 1) 
           for(ij = 1; ij <= iNoEigen; ij = ij + 1) {
               spAwork->uMatrix.daa[ ii-1 ][ ij-1 ] = 0.0; 
               for(ik = 1; ik <= iSize; ik = ik + 1) 
                   spAwork->uMatrix.daa[ ii-1 ][ ij-1 ] += 
                            spV->uMatrix.daa[ ik-1 ][ ii-1 ] * 
                            spYnew->uMatrix.daa[ ik-1 ][ ij-1 ];
           }

           /* [c.5] : Compute [Y_{new}] = [B] . [V_{k+1}] */

           for(ii = 1; ii <= iSize; ii = ii + 1) 
           for(ij = 1; ij <= iNoEigen; ij = ij + 1) {
               spYnew->uMatrix.daa[ ii-1 ][ ij-1 ] = 0; 
               for(ik = 1; ik <= iSize; ik = ik + 1)  {

                   iMin = (int) MIN( ii, ik );
                   iMax = (int) MAX( ii, ik );

                   if((iMax-iMin+1) <= spB->uMatrix.daa[iMax-1][0]) {
                       spYnew->uMatrix.daa[ ii-1 ][ ij-1 ] += 
                               spB->uMatrix.daa[ iMax-1 ][ iMax-iMin+1 ] *
                               spV->uMatrix.daa[ ik-1 ][ ij-1 ] ; 
                   }
               }
           }

           /* [c.6] : Compute [B_{k+1}] = [V_{k+1}]^T . [Y_{k+1}] */

           for(ii = 1; ii <= iNoEigen; ii = ii + 1) 
           for(ij = 1; ij <= iNoEigen; ij = ij + 1) {
               spBwork->uMatrix.daa[ ii-1 ][ ij-1 ] = 0.0; 
               for(ik = 1; ik <= iSize; ik = ik + 1) 
                   spBwork->uMatrix.daa[ ii-1 ][ ij-1 ] += 
                            spV->uMatrix.daa[ ik-1 ][ ii-1 ] * 
                            spYnew->uMatrix.daa[ ik-1 ][ ij-1 ];
           }

           /* [c.7] : Solve [A_{k+1}].[Q_{k+1}] = [B_{k+1}].[Q_{k+1}].[Lambda_{k+1}] */

           GeneralizedEigen( spAwork, spBwork, spEigenvalue, spEigenvectorQ );

           /* [c.8] : Check Convergence Criteria */

           iSubspaceConvergence = TRUE;
           for( ii = 1; ii <= spEigenvalue->iNoRows; ii = ii + 1) {
                dEigenvalue    = spEigenvalue->uMatrix.daa[ii-1][0];
                dEigenvalueOld = spEigenvalueOld->uMatrix.daa[ii-1][0];

                dConvergence = ABS(dEigenvalueOld - dEigenvalue)/dEigenvalue;
                if(dConvergence > TOLERANCE) {
                   iSubspaceConvergence = FALSE;
                }
           }

           MatrixFree ( spEigenvalueOld );
           spEigenvalueOld = MatrixCopy( spEigenvalue );

           /* [c.9] : Update [X_{k+1}] = [V_{k+1}].[Q_{k+1}] */

              spTemp1 = MatrixMult( spV, spEigenvectorQ );

           /* [c.10] : Scale [X_{k+1}] so that Max Value is 1 */

              for(ij = 1; ij <= iNoEigen; ij = ij + 1)  {
                  dMaxValue = 0.0;
                  for(ii = 1; ii <= iSize; ii = ii + 1) 
                      if(ABS(spTemp1->uMatrix.daa[ ii-1 ][ ij-1 ]) >= ABS(dMaxValue))
                         dMaxValue = spTemp1->uMatrix.daa[ ii-1 ][ ij-1 ];

                  for(ii = 1; ii <= iSize; ii = ii + 1) 
                      spEigenvector->uMatrix.daa[ ii-1 ][ ij-1 ] =
                                     spTemp1->uMatrix.daa[ ii-1 ][ ij-1 ]/dMaxValue;
              }

       }

       /* [d] : Summarize Algorithm Performance */

       printf("\n*** SUBSPACE ITERATION CONVERGED IN %2d ITERATIONS \n", iLoop);

       /* [e] : Clean Up before Leaving */

       MatrixFree ( spTemp1 );
       MatrixFree ( spV );
       MatrixFree ( spA_LU );
       MatrixFree ( spAwork );
       MatrixFree ( spBwork );
       MatrixFree ( spEigenvectorQ );
       MatrixFree ( spEigenvalueOld );

       MatrixFree ( spYnew );
       MatrixFree ( spY );
       MatrixFree ( spVk );
}

