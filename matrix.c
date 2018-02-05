/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  matrix.c : High-Level Functions for Matrix Operations and Printing
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

#ifdef  __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"

/*
#define DEBUG 
*/


/*
 *  ===============================================================================
 *  MatrixPrintCast(matrix, ...) : Print a matrix casting appropriate row or column
 *                                 units to those listed in [ units_m ].
 *                                                                 
 *  Usage :  In Aladdin Input File                                      
 *           MatrixPrintCast(m)                                     
 *        or MatrixPrintCast(m1, m2, units_m)                      
 *        or MatrixPrintCast(m1, m2,...,m4, units_m)               
 *                                                                 
 *        Column units are scaled with [ units_m ] = [ a , b , c ]
 *        Row units are scaled with    [ units_m ] = [ a ; b ; c ]
 *                                                                 
 *  Usage :  In c-code                                       
 *           MatrixPrint(m, (MATRIX *) NULL)                  
 *        or MatrixPrint(m1, m2, units_m, (MATRIX *) NULL)    
 *  ===============================================================================
 */

/* #define MATRIX_DEBUG */

#ifdef  __STDC__
MATRIX *MatrixPrintCast( MATRIX *first_matrix, ... ) {
#else
MATRIX *MatrixPrintCast(va_alist)
va_dcl
{
MATRIX *first_matrix;
#endif

va_list arg_ptr;
MATRIX      *p, *q;
DIMENSIONS *d1,*d2;
int col_counter, *COL_COUNT; 
int row_counter, *ROW_COUNT; 
int NO_COL_WITH_UNITS, NO_ROW_WITH_UNITS; 
int NO_COL_WITH_NO_UNITS, NO_ROW_WITH_NO_UNITS; 
int UNITS_NO_q, UNITS_NO_p, UNITS_SWITCH;
int length, length1, length2; 
int i, j, k, ii, jj, COL_ROW; 
int iFLAG  = FALSE;

#ifdef MATRIX_DEBUG
       printf("Enter MatrixPrintCast() function\n");
#endif

    /* [a] : Get units_matrix */

#ifdef  __STDC__
    va_start(arg_ptr, first_matrix);
#else
    va_start(arg_ptr);
    first_matrix = va_arg(arg_ptr, MATRIX *);
#endif

    UNITS_SWITCH = CheckUnits();
    i = 0;
    q = (MATRIX *) 1;  /* a NON NULL pointer */

    /* [b] : Compute number of matrices in argument list */

    for(p = first_matrix; q != NULL; q = va_arg(arg_ptr, MATRIX *)) {
        if( i != 0) p = q;
        i = i + 1;
    }
    va_end(arg_ptr);

    /* [c] : One matrix in argument list (call standard matrixprint function) */

    if( i == 1) {
        MatrixPrintVar( (MATRIX *) first_matrix, (MATRIX *) NULL );
        return;
    }

    /* [d] : Check that units are on for rest of function */

    if( UNITS_SWITCH == OFF) {
        FatalError("In MatrixPrintCast() : Units must be on to use this function",
                  (char *)NULL);
    }

    /* [e] : Walk along Matrix Argument List and print contents */

#ifdef __STDC__
    va_start(arg_ptr, first_matrix);
#else
    va_start(arg_ptr);
    first_matrix = va_arg(arg_ptr, MATRIX *);
#endif

    for(q = first_matrix; q != p;  q = va_arg(arg_ptr, MATRIX *)) {

        /* [e.1] : Allocate memory for row and columns ... */ 

        ROW_COUNT = (int *) MyCalloc(q->iNoRows,    sizeof(int));
        COL_COUNT = (int *) MyCalloc(q->iNoColumns, sizeof(int));

        /* [e.2] : Count number of columns in matrix q with units          */ 
        /*                                                                 */
        /*         Array COL_COUNT[i] stores the no of the i-th column in  */
        /*         matrix q having units.                                  */
        /*                                                                 */

        col_counter = 0;
        for(k = 1; k <= q->iNoColumns; k++) {
            /* Column units */
            d2 = &(q->spColUnits[k-1]);
            if(d2->length_expnt != 0 || d2->time_expnt   != 0 || d2->mass_expnt != 0 ||
               d2->temp_expnt   != 0 || d2->radian_expnt != 0 ) {
               col_counter = col_counter + 1;
               COL_COUNT[col_counter-1] = k;  /* record column with units */
            } 
        }

        /* [e.3] : Count number of rows in matrix q having units        */ 
        /*                                                              */
        /*         Array ROW_COUNT[i] stores the no of the i-th row in  */
        /*         matrix q having units.                               */
        /*                                                              */

        row_counter = 0;
        for(j = 1; j<= q->iNoRows; j++) {
            d2 = &(q->spRowUnits[j-1]);
            if(d2->length_expnt != 0 || d2->time_expnt   != 0 || d2->mass_expnt != 0 ||
               d2->temp_expnt   != 0 || d2->radian_expnt != 0 ) {
               row_counter = row_counter + 1;
               ROW_COUNT[row_counter-1] = j;  /* record row  with units */
            }
        }

        /* Save cols and rows with (and without) units.            */

        NO_ROW_WITH_UNITS = row_counter;
        NO_ROW_WITH_NO_UNITS = q->iNoRows - NO_ROW_WITH_UNITS;
        NO_COL_WITH_UNITS = col_counter;
        NO_COL_WITH_NO_UNITS = q->iNoColumns - NO_COL_WITH_UNITS;

        /* If matrix "q" is dimensionless, print error message and */
        /* terminate program execution                             */

        if(NO_COL_WITH_NO_UNITS == q->iNoColumns &&
           NO_ROW_WITH_NO_UNITS == q->iNoRows) { 
           FatalError("In MatrixPrintCase()"
                      "Trying to print out a non-dimensional matrix with units",
                       (char *) NULL);
        } 

        /* Units matrix must be larger than (1x1) matrix */

        if(p->iNoRows == 1 && p->iNoColumns == 1) {
           FatalError("In MatrixPrintCase()"
                      "Units matrix must have more than one column or row",
                       (char *) NULL);
        }

        /* p->iNoRows == 1 implies scaling of column units */

        if(p->iNoRows == 1) {
           UNITS_NO_p = p->iNoColumns;
           UNITS_NO_q = NO_COL_WITH_UNITS;
           COL_ROW    = COLUMN;
        }

        /* p->iNoColumns == 1 implies scaling of column units */

        if(p->iNoColumns == 1) {
           UNITS_NO_p = p->iNoRows;
           UNITS_NO_q = NO_ROW_WITH_UNITS;
           COL_ROW    = ROW;
        }

        /* If matrix "p" contains more rows/columns of units than matrix "q" */
        /* print error message and terminate program execution               */

        if(UNITS_NO_q < UNITS_NO_p) {
           FatalError("In MatrixPrintCase()"
                      "The units matrix contains too many rows/columns",
                      (char *) NULL); 
        }

        /* [e.5] : No of columns with units in matrix q are larger   */
        /*         than or equal to no of columns in units matrix p  */

        /* Case 1 : Scale Column Units */

        if( COL_ROW == COLUMN ) {
           for(ii = 1; ii <= UNITS_NO_q; ii++) {
               j  = COL_COUNT[ii-1];
               d1 = &(q->spColUnits[j-1]);
               for(k = 1; k <= p->iNoColumns; k++) {
                   d2 = &(p->spColUnits[k-1]);

                   if((d1->length_expnt  == d2->length_expnt) &&
                      (d1->mass_expnt    == d2->mass_expnt)   &&
                      (d1->time_expnt    == d2->time_expnt)   &&
                      (d1->temp_expnt    == d2->temp_expnt)   &&
                      (d1->radian_expnt  == d2->radian_expnt))  {
                      UnitsCopy( &(q->spColUnits[j-1]), &(p->spColUnits[k-1]) );
                      iFLAG = TRUE;
                   }
               }
            }
            if(iFLAG == FALSE)
               FatalError("In MatrixPrintCast(): Inconsistent Units",
                         (char *) NULL);
        }

        /* Case 2 : Scale Row Units */

        if( COL_ROW == ROW ) {
            for(j = 1; j <= UNITS_NO_q; j++) {
               k  = ROW_COUNT[j-1];
               d1 = &(q->spRowUnits[k-1]);
               for(jj = 1; jj <= p->iNoRows; jj++) {
                   d2 = &(p->spRowUnits[jj-1]);
                   if((d1->length_expnt  == d2->length_expnt) &&
                      (d1->mass_expnt    == d2->mass_expnt)   &&
                      (d1->time_expnt    == d2->time_expnt)   &&
                      (d1->temp_expnt    == d2->temp_expnt)   &&
                      (d1->radian_expnt  == d2->radian_expnt))  {
                         UnitsCopy( &(q->spRowUnits[k-1]), &(p->spRowUnits[jj-1]) );
                         iFLAG = TRUE;
                   }
               }
           }
           if(iFLAG == FALSE)
              FatalError("In MatrixPrintCast(): Inconsistent Units",
                        (char *) NULL );
        }

        /* Now print matrix (without simplification of units) */

        switch((int)p->eRep) {
          case  SEQUENTIAL:
          case  INDIRECT:
                switch((int) p->eType) {
                    case DOUBLE_ARRAY:
                         MatrixPrintIndirectDouble(q);
                         break;
                    case INTEGER_ARRAY:
                         MatrixPrintIndirectInteger(q);
                         break;
                    default:
                         FatalError("In MatrixPrintVar() : Undefined m->eType",
                                   (char *) NULL);
                         break;
                }
                break;
          case  SKYLINE:
                MatrixPrintSkylineDouble(q);
                break;
          default:
                FatalError("In MatrixPrintVar() : Undefined m->eRep",
                          (char *) NULL);
                break;
       } 

       free( (char *) ROW_COUNT );
       free( (char *) COL_COUNT );
   }

   va_end(arg_ptr);

#ifdef MATRIX_DEBUG
       printf("Leave MatrixPrintCast() function\n");
#endif

}


/*
 *  =====================================================================
 *  MatrixPrintVar() : Print a variable number of matrices
 *
 *  In Aladdin input file this function is called with PrintMatrix( .. );
 *  =====================================================================
 */

#ifdef  __STDC__
MATRIX *MatrixPrintVar(MATRIX *first_matrix, ...) {
#else
MATRIX *MatrixPrintVar(va_alist)
va_dcl
{
MATRIX  *first_matrix;
#endif

va_list      arg_ptr;
MATRIX            *p;
int     UNITS_SWITCH;

    /* [a] : Get first function argument */

#ifdef __STDC__
    va_start(arg_ptr, first_matrix);
#else
    va_start(arg_ptr);
    first_matrix = va_arg(arg_ptr, MATRIX *);
#endif

    /* [b] : Print matrices : units switch == ON */

    UNITS_SWITCH = CheckUnits();

    if( UNITS_SWITCH == ON ) {
    for(p = first_matrix; p != NULL;  p = va_arg(arg_ptr, MATRIX *)) {
         
        switch((int)p->eRep) {
            case  SEQUENTIAL:
            case  INDIRECT:
                  switch((int) p->eType) {
                      case DOUBLE_ARRAY:
                           MatrixUnitsSimplify(p);
                           MatrixPrintIndirectDouble(p);
                           break;
                      case INTEGER_ARRAY:
                           MatrixUnitsSimplify(p);
                           MatrixPrintIndirectInteger(p);
                           break;
                      default:
                           FatalError("In MatrixPrintVar() : Undefined m->eType",
                                     (char *) NULL);
                           break;
                  }
                  break;
            case  SKYLINE:
                  MatrixUnitsSimplify(p);
                  MatrixPrintSkylineDouble(p);
                  break;
            default:
                  FatalError("In MatrixPrintVar() : Undefined m->eRep",
                            (char *) NULL);
                  break;
        } 
    }
    }

    /* [c] : Print matrices : units switch == OFF */

    if( UNITS_SWITCH == OFF ) {
    for(p = first_matrix; p != NULL;  p = va_arg(arg_ptr, MATRIX *)) {
         
        switch((int)p->eRep) {
            case  SEQUENTIAL:
            case  INDIRECT:
                  switch((int) p->eType) {
                      case DOUBLE_ARRAY:
                           MatrixPrintIndirectDouble(p);
                           break;
                      default:
                           FatalError("In MatrixPrintVar() : Undefined matrix->eType",
                                     (char *) NULL);
                           break;
                  }
                  break;
            case  SKYLINE:
                  MatrixPrintSkylineDouble(p);
                  break;
            default:
                  FatalError("In MatrixPrintVar() : Undefined matrix->eRep",
                            (char *) NULL);
                  break;
       }
    }
    }

    va_end(arg_ptr);
}



/*
 *  ==============================================================
 *  Allocate Matrix data structure, and array[iNoRows][iNoColumns]
 *  ==============================================================
 */

#ifdef __STDC__
MATRIX *MatrixAllocate(MATRIX *m)
#else
MATRIX *MatrixAllocate(m)
MATRIX *m;
#endif
{
MATRIX         *m1;
int        iNoRows;
int     iNoColumns;
int              i;

   /* [a] : Check size and values of argument matrix m */

   if((m->iNoRows < 1) || (m->iNoColumns < 2))
       FatalError("In MatrixAllocate() : Argument Matrix Too Small",(char *)NULL);

   /* [b] : Allocate and return Matrix */

   iNoRows    = (int) m->uMatrix.daa[0][0];
   iNoColumns = (int) m->uMatrix.daa[0][1];
 
   m1 = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, iNoRows, iNoColumns);

   /* [c] : Zero out units in row/column units vectors */

   if( CheckUnits() == ON ) {
      for( i=1 ; i<= iNoRows ; i++ )
        ZeroUnits(&(m1->spRowUnits[i-1]));
      for( i=1 ; i<= iNoColumns ; i++ )
        ZeroUnits(&(m1->spColUnits[i-1]));
   }

   return(m1);
}

#ifdef __STDC__
MATRIX *MatrixDiag(MATRIX *m)
#else
MATRIX *MatrixDiag(m)
MATRIX *m;
#endif
{
MATRIX              *m1;
int    matrix_dimension;
int           length, i;

   /* [a] : Check size and values of argument matrix m */

   if((m->iNoRows < 1) || (m->iNoColumns < 2))
       FatalError("In MatrixDiag() : Argument Matrix Too Small",(char *)NULL);

   /* [b] : Check dimensions of matrix */

   matrix_dimension = (int) m->uMatrix.daa[0][0];
   if(matrix_dimension < 1)
      FatalError("In MatrixDiag() : Matrix dimension < 1",(char *)NULL);

   /* [c] : Instantiate and return Matrix                    */
   /*                                                        */
   /* Note : UNITS OF MATRIX ELEMENTS ARE STORED             */
   /*        IN COLUMN_UNITS_BUF AS AN DEFAULT               */

   m1 = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, matrix_dimension, matrix_dimension);

   for(i = 1; i <= matrix_dimension; i++) {
       m1->uMatrix.daa[i-1][i-1]= (double) m->uMatrix.daa[0][1];

       if( CheckUnits() == ON ) {
           UnitsCopy(&(m1->spColUnits[i-1]), &(m->spColUnits[1]));
           ZeroUnits(&(m1->spRowUnits[i-1]));
       }
   }

   return(m1);
}

/*==================================*/
/* Usage :                          */
/* MatrixZero([n]):                 */
/*   generate a nxn matrix of zeros */
/* MatrixZero([n,m]):               */
/*   generate a nxm matrix of zeros */
/*                                  */
/*==================================*/

#ifdef __STDC__
MATRIX *MatrixZero(MATRIX *m)
#else
MATRIX *MatrixZero(m)
MATRIX *m;
#endif
{
MATRIX         *m1;
int        iNoRows;
int     iNoColumns;
int              i;
int              j;


#ifdef DEBUG
      printf("Enter MatrixZero() : \n");
#endif 
   /* Check size and values of argument matrix m */

      if((m->iNoRows < 1) || (m->iNoColumns < 1))
         FatalError("In MatrixZero() : Argument Matrix Too Small",(char *)NULL);

      iNoRows    = (int) m->uMatrix.daa[0][0];

      switch((int) m->iNoColumns) {
        case 1:
          iNoColumns = (int) m->uMatrix.daa[0][0];
        break;
        case 2:
          iNoColumns = (int) m->uMatrix.daa[0][1];
        break;
        default:
          printf("*** too many columns in the matrix. \n Check Zero() in input file \n");
          FatalError("*** Syntax error: two columns at most in Zero([]), (char *) NULL");
        break;
      }

      m1 = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, iNoRows, iNoColumns);
      if(CheckUnits()==ON) {
         for(i = 1; i <= iNoRows; i++)
             ZeroUnits(&(m1->spRowUnits[i-1]));
         for(i = 1; i <= iNoColumns; i++)
             ZeroUnits(&(m1->spColUnits[i-1]));
      }

#ifdef DEBUG
      printf("In MatrixZero():  iNoRows    = %d \n", iNoRows);
      printf("In MatrixZero():  iNoColumns = %d \n", iNoColumns);
#endif

      for(i = 1; i <= iNoRows; i++)
         for(j = 1; j <= iNoColumns; j++)
             m1->uMatrix.daa[i-1][j-1] = 0.0;

#ifdef DEBUG
      printf("Leaving MatrixZero() :  iNoRows = %d \n", iNoRows);
#endif 
      return(m1);
}

#ifdef __STDC__
MATRIX *MatrixOne(MATRIX *m)
#else
MATRIX *MatrixOne(m)
MATRIX *m;
#endif
{
MATRIX         *m1;
int        iNoRows;
int     iNoColumns;
int            i,j;

   /* Check size and values of argument matrix m */

      if((m->iNoRows < 1) || (m->iNoColumns < 1))
         FatalError("In MatrixOne() : Argument Matrix Too Small",(char *)NULL);

      iNoRows    = (int) m->uMatrix.daa[0][0];

      switch((int) m->iNoColumns) {
        case 1:
          iNoColumns = (int) m->uMatrix.daa[0][0];
        break;
        case 2:
          iNoColumns = (int) m->uMatrix.daa[0][1];
        break;
        default:
          printf("*** too many columns in the matrix. \n Check Zero() in input file \n");
          FatalError("*** Syntax error: two columns at most in Zero([]), (char *) NULL");
        break;
      }

      m1 = MatrixAllocIndirect((char *) NULL, DOUBLE_ARRAY, iNoRows, iNoColumns);
      if(CheckUnits()==ON) {
         for(i = 1; i <= iNoRows; i++)
             ZeroUnits(&(m1->spRowUnits[i-1]));
         for(i = 1; i <= iNoColumns; i++)
             ZeroUnits(&(m1->spColUnits[i-1]));
      }

      for(i = 1; i <= iNoRows; i++)
         for(j = 1; j <= iNoColumns; j++)
             m1->uMatrix.daa[i-1][j-1] = 1.0;

      return(m1);
}


/*
 *  ===========================================
 *  MatrixFree() : Free Matrix Storage
 *  
 *  Input  : Matrix spA -- pointer to matrix A.  
 *  Output : void.
 *  ===========================================
 */

#ifdef __STDC__
void MatrixFree( MATRIX *spA )
#else  /* start case not STDC */
void MatrixFree( spA )
MATRIX *spA;
#endif /* end case not STDC */
{
      if ( spA==(MATRIX *)NULL )   return;

      switch((int) spA->eRep) {
          case INDIRECT:
               MatrixFreeIndirect( spA );
               break;
          case SKYLINE:
               MatrixFreeSkyline( spA );
               break;
          default:
               FatalError("In MatrixFree() : Undefined spA->eRep",
                         (char *) NULL);
               break;
      }
}


/*
 *  ===================================================
 *  MatrixAdd() : Matrix Add Operation [c] = [a] + [b].
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *             MATRIX spB       -- Pointer to Matrix B
 *  Output :   MATRIX spC       -- Pointer to Matrix C 
 *  ===================================================
 */

#ifdef __STDC__
MATRIX *MatrixAdd( MATRIX *spA, MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixAdd( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* end case not STDC */
{
MATRIX *spC;

    /* [a] : Check compatibility of matrix storage schemes, data types, and sizes */

    if(spA->eRep != spB->eRep) {
       printf("WARNING : Incompatible Matrix Representations in MatrixAdd()\n");
       printf("WARNING : *** spA->eRep = %4d spB->eRep = %4d\n",
               spA->eRep, spB->eRep);
    }

    if(spA->eType != spB->eType) {
       printf("WARNING : Incompatible Matrix Types in MatrixAdd()\n");
       printf("WARNING : *** spA->eType = %4d spB->eType = %4d\n",
               spA->eType, spB->eType);
    }

    /* [b] : Compute matrix operation */

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spC = MatrixAddIndirectDouble( spA, spB );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixAdd() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spC = MatrixAddSkyline( spA, spB );
             break;
        case SPARSE:
             break;
        default:
             FatalError("In MatrixAdd() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spC );
}

/*
 *  ==============================================================
 *  MatrixAddReplace() : Compute matrix assignment [A] = [A] + [B]
 *  
 *  Input :    MATRIX spA -- Pointer to Matrix A
 *             MATRIX spB -- Pointer to Matrix B
 *  Output :   MATRIX spA -- Pointer to Replacement Matrix A 
 *  ==============================================================
 */

#ifdef __STDC__
MATRIX *MatrixAddReplace( MATRIX *spA, MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixAddReplace( spA , spB )
MATRIX *spA;
MATRIX *spB;
#endif /* end case not STDC */
{
MATRIX *spC;

    /* [a] : Check that matrix types are compatible */

    if(spA->eRep != spB->eRep) {
       printf("WARNING : Incompatible Matrix Representations in MatrixAddReplace()\n");
       printf("WARNING : *** m1->eRep = %4d m2->eRep = %4d\n", spA->eRep, spB->eRep);
    }

    if(spA->eType != spB->eType) {
       printf("WARNING : Incompatible Matrix Types in MatrixAdd()\n");
       printf("WARNING : *** m1->eType = %4d m2->eType = %4d\n", spA->eType, spB->eType);
    }

    /* [b] : Compute matrix operation */

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spA = MatrixAddReplaceIndirectDouble( spA, spB );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixAddReplace() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spC = MatrixAddSkyline( spA, spB );
             MatrixFreeSkyline( spA );
             return( spC );
             break;
        default:
             FatalError("In MatrixAddReplace() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spA );
}

/*
 *  ==========================================================
 *  MatrixSub() : Matrix Subtraction Operation [c] = [a] - [b]
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *             MATRIX spB       -- Pointer to Matrix B
 *  Output :   MATRIX spC       -- Pointer to Matrix C 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixSub( MATRIX *spA, MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixSub( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* end case not STDC */
{
MATRIX *spC;

    /* [a] : Check that matrix types are compatible */

    if( spA->eRep != spB->eRep ) {
        printf("WARNING : Incompatible Matrix Representations in MatrixSub()\n");
        printf("WARNING : *** m1->eRep = %4d m2->eRep = %4d\n",
                spA->eRep, spB->eRep);
    }

    if( spA->eType != spB->eType ) {
        printf("WARNING : Incompatible Matrix Types in MatrixSub()\n");
        printf("WARNING : *** m1->eType = %4d m2->eType = %4d\n",
                spA->eType, spB->eType);
    }

    /* [b] : Compute Matrix Operation */

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spC = MatrixSubIndirectDouble( spA, spB );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixSub() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spC = MatrixSubSkyline( spA, spB );
             break;
        case SPARSE:
             break;
        default:
             FatalError("In MatrixSub() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spC );
}

/*
 *  ==============================================================
 *  MatrixSubReplace() : Compute matrix assignment [A] = [A] - [B]
 *  
 *  Input :    MATRIX spA -- Pointer to Matrix A
 *             MATRIX spB -- Pointer to Matrix B
 *  Output :   MATRIX spA -- Pointer to Replacement Matrix A 
 *  ==============================================================
 */

#ifdef __STDC__
MATRIX *MatrixSubReplace( MATRIX *spA, MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixSubReplace( spA , spB )
MATRIX *spA;
MATRIX *spB;
#endif /* end case not STDC */
{
MATRIX *spC;

    /* [a] : Check that matrix types are compatible */

    if(spA->eRep != spB->eRep) {
       printf("WARNING : Incompatible Matrix Representations in MatrixAddReplace()\n");
       printf("WARNING : *** spA->eRep = %4d spB->eRep = %4d\n",
               spA->eRep, spB->eRep);
    }

    if(spA->eType != spB->eType) {
       printf("WARNING : Incompatible Matrix Types in MatrixAdd()\n");
       printf("WARNING : *** spA->eType = %4d spB->eType = %4d\n",
               spA->eType, spB->eType);
    }

    /* [b] : Compute matrix operation */

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spA = MatrixSubReplaceIndirectDouble( spA, spB );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixAddReplace() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spC = MatrixSubSkyline( spA, spB );
             MatrixFreeSkyline( spA );
             return( spC );
             break;
        default:
             FatalError("In MatrixSubReplace() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spA );
}

/*
 *  ==========================================================
 *  MatrixMult() : Matrix Multiplication [c] = [a].[b]
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *             MATRIX spB       -- Pointer to Matrix B
 *  Output :   MATRIX spC       -- Pointer to Matrix C 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixMult( MATRIX *spA , MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixMult( spA, spB )
MATRIX *spA;
MATRIX *spB;
#endif /* end case not STDC */
{
MATRIX *spC;

    /* [a] : Check that matrix data types are compatible */

    if(spA->eType != spB->eType) {
       printf("WARNING : Incompatible Matrix Types in MatrixMult()\n");
       printf("WARNING : *** spA->eType = %4d spB->eType = %4d\n",
               spA->eType, spB->eType);
    }

    /* [b] : Compute Matrix Multiplication for Indirect Storage */

    if((int) spA->eRep == INDIRECT && (int) spB->eRep == INDIRECT) {
       switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                spC = MatrixMultIndirectDouble( spA, spB );
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixMult() : Undefined spA->eType",
                          (char *) NULL);
                break;
       }
    }

    /* [c] : Compute Matrix Multiplication for Skyline-Indirect Storage */

    if((int) spA->eRep == SKYLINE && (int) spB->eRep == INDIRECT) {
       switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                spC = MatrixMultSkylineIndirectDouble( spA, spB );
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixMult() : Undefined spA->eType",
                          (char *) NULL);
                break;
       }
    }

    /* [d] : Compute Matrix Multiplication for Skyline-Indirect Storage */

    if((int) spA->eRep == INDIRECT && (int) spB->eRep == SKYLINE) {
       switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                spC = MatrixMultIndirectSkylineDouble( spA, spB );
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixMult() : Undefined spA->eType",
                          (char *) NULL);
                break;
       }
    }

    /* [e] : Compute Matrix Multiplication for Skyline-Skyline Storage */

    if((int) spA->eRep == SKYLINE && (int) spB->eRep == SKYLINE) {
       switch((int) spA->eType) {
           case DOUBLE_ARRAY:
       printf(" case : SKYLINE : DOUBLE \n");
                spC = MatrixMultSkyline( spA, spB );
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixMult() : Undefined spA->eType",
                          (char *) NULL);
                break;
       }
    }

    return ( spC );
}

/*
 *  ==========================================================
 *  MatrixPower() : Matrix Power function [c] = [a]^n
 *  
 *  Input :    MATRIX   spA     -- Pointer to Matrix   A
 *             QUANTITY spQ     -- Pointer to Qauntity Q
 *  Output :   MATRIX   spC     -- Pointer to Matrix   C 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixPower( MATRIX *spA , QUANTITY *spQ )
#else  /* start case not STDC */
MATRIX *MatrixPower( spA, spQ )
MATRIX   *spA;
QUANTITY *spQ;
#endif /* end case not STDC */
{
MATRIX       *spC;
double       dexp; 
int       i, j, n; 
static double EPS = 1.0E-10;

    /* [a] : Check that matrix are square matrix */

    if((int) spA->iNoRows != (int) spA->iNoColumns) {
       printf(" *** In MatrixPower(): \n");
       printf(" *** No of rows of %s    = %d \n", spA->iNoRows, spA->cpMatrixName);
       printf(" *** No of columns of %s = %d \n", spA->iNoColumns, spA->cpMatrixName);
       FatalError(" Matrix must be a quare matrix to use MatrixPower() function",
                 (char *) NULL);
    }

    /* [b] : Check that exponent is an integer */

    dexp = spQ->value;
    n    = (int) dexp;
    if(abs(dexp/((double) n) - 1.0) > EPS) {
       printf(" *** In MatrixPower(): \n");
       printf(" *** ratio  = %lf\n", abs(dexp/(double) n));
       printf(" *** Exponent of = %lf\n", dexp);
       FatalError(" Exponent must be an integer to use MatrixPower() function",
                  (char *) NULL);
    }

    /* [c] : Check that exponent larger than -1 */

    if(dexp < 0.0 && n != (int) -1) {
       printf(" *** In MatrixPower(): \n");
       printf(" *** Exponent of %lf = \n", dexp);
       FatalError(" Exponent must be larger -1 to use MatrixPower() function",
                  (char *) NULL);
    }
     
    /* [d] : Exponent is zero */

    if(n == (int) 0 ) {
       spC = MatrixAllocIndirect((char *) NULL,DOUBLE_ARRAY,spA->iNoRows,spA->iNoColumns);
      if(CheckUnits()==ON) {
         for(i = 1; i <= spA->iNoRows; i++)
             ZeroUnits( &(spC->spRowUnits[i-1]) );
         for(i = 1; i <= spA->iNoColumns; i++)
             ZeroUnits( &(spC->spColUnits[i-1]) );
      }
      for(i = 1; i <= spA->iNoRows; i++)
         for(j = 1; j <= spA->iNoColumns; j++)
             spC->uMatrix.daa[i-1][j-1] = 1.0;
      return(spC);
     }

    /* [e] : Calculate the inverse for Indirect Storge */

     if((int) spA->eRep == INDIRECT && n == (int) -1) {
         switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                spC = MatrixInverseIndirectDouble(spA);
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixPower() : Undefined spA->eType",
                          (char *) NULL);
                break;
         }
     }

    /* [f] : Calculate the inverse for Skyline Storge */

     if((int) spA->eRep == SKYLINE && n == (int) -1) {
         switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                spC = MatrixInverseSkyline(spA);
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixPower() : Undefined spA->eType",
                          (char *) NULL);
                break;
         }
     }


    /* [g] : Compute Matrix Power for Indirect Storage */

    if((int) spA->eRep == INDIRECT && n > 0) {
       switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                if(n == (int) 1)
                   spC = MatrixCopyIndirectDouble(spA);
                else{
                   spC = MatrixCopyIndirectDouble(spA);
                   for(i = 1; i <= n-1; i++) 
                       spC = MatrixMultIndirectDouble( spC, spA );
                }
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixPower() : Undefined spA->eType",
                          (char *) NULL);
                break;
       }
    }

    /* [c] : Compute Matrix Multiplication for Skyline Storage */

    if((int) spA->eRep == SKYLINE && n > 0) {
       switch((int) spA->eType) {
           case DOUBLE_ARRAY:
                if(n == (int) 1)
                   spC = MatrixCopySkyline(spA);
                else{
                   spC = MatrixCopySkyline(spA);
                   for(i = 1; i <= n-1; i++) 
                       spC = MatrixMultSkyline(spC, spA );
                }
                break;
           case INTEGER_ARRAY:
           case COMPLEX_ARRAY:
           default:
                FatalError("In MatrixPower() : Undefined spA->eType",
                          (char *) NULL);
                break;
       }
    }

    return ( spC );
}

/*
 *  ==========================================================
 *  MatrixNegate() : Compute Negative Matrix [B] = -[A];
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *  Output :   MATRIX spB       -- Pointer to Matrix B 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixNegate( MATRIX *spA )
#else  /* start case not STDC */
MATRIX *MatrixNegate( spA )
MATRIX *spA;
#endif /* end case not STDC */
{
MATRIX *spB;

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spB = MatrixNegateIndirectDouble( spA );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixNegate() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spB = MatrixNegateSkyline( spA );
             break;
        default:
             FatalError("In MatrixNegate() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spB );
}

/*
 *  ==========================================================
 *  MatrixNegateReplace() : Matrix Replacement [A] = -[A];
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *  Output :   MATRIX spA       -- Pointer to Matrix -A 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixNegateReplace( MATRIX *spA )
#else  /* start case not STDC */
MATRIX *MatrixNegateReplace( spA )
MATRIX *spA;
#endif /* end case not STDC */
{

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spA = MatrixNegateReplaceIndirectDouble(spA);
                      break;
                 default:
                      FatalError("In MatrixNegateReplace() : Undefined m->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spA = MatrixNegateReplaceSkyline( spA );
             break;
        default:
             FatalError("In MatrixNegateReplace() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return (spA);
}

/*
 *  ==========================================================
 *  MatrixTranspose() : Compute Matrix Transpose [B] = [A]^T.
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *  Output :   MATRIX spB       -- Pointer to Matrix B 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixTranspose( MATRIX *spA )
#else  /* start case not STDC */
MATRIX *MatrixTranspose( spA )
MATRIX *spA;
#endif /* end case not STDC */
{
MATRIX *spB;

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spB = MatrixTransposeIndirectDouble( spA );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixTranspose() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spB = MatrixTransposeSkyline( spA );
             break;
        case SPARSE:
             break;
        default:
             FatalError("In MatrixTranspose() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spB );
}

/*
 *  ==========================================================
 *  MatrixCopy() : Make Copy of Matrix [B] = [A].
 *  
 *  Input :    MATRIX spA       -- Pointer to Matrix A
 *  Output :   MATRIX spB       -- Pointer to Matrix B 
 *  ==========================================================
 */

#ifdef __STDC__
MATRIX *MatrixCopy( MATRIX *spA )
#else  /* start case not STDC */
MATRIX *MatrixCopy( spA )
MATRIX *spA;
#endif /* end case not STDC */
{
MATRIX *spB;

    switch((int) spA->eRep) {
        case INDIRECT:
             switch((int) spA->eType) {
                 case DOUBLE_ARRAY:
                      spB = MatrixCopyIndirectDouble( spA );
                      break;
                 case INTEGER_ARRAY:
                 case COMPLEX_ARRAY:
                 default:
                      FatalError("In MatrixCopy() : Undefined spA->eType",
                                (char *) NULL);
                      break;
             }
             break;
        case SKYLINE:
             spB = MatrixCopySkyline( spA );
             break;
        default:
             FatalError("In MatrixCopy() : Undefined spA->eRep",
                       (char *) NULL);
             break;
    }

    return ( spB );
}


/*
 *  ==========================================================================
 *  MatrixSolve() : Solve Systems of Linear Equations [A][X] = [B].
 * 
 *  Input :    MATRIX spA -- Pointer to Matrix A
 *             MATRIX spB -- Pointer to r-h-s Matrix B
 *  Output :   MATRIX spX -- Pointer to solution matrix X
 * 
 *  Also calls functions : MatrixLU() - LU decomposition.
 *                         MatrixFB() - Forward/Back substitution.
 *
 *  WARNING : This high-level equation solver works for only one decomposition
 *            followed by multiple forward/backward substitutions.
 *  ==========================================================================
 */

static VECTOR *spPivot = (VECTOR *)NULL;

#ifdef __STDC__
MATRIX *MatrixSolve( MATRIX *spA , MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixSolve( spA, spB )
MATRIX      *spA, *spB;
#endif /* end case not STDC */
{
MATRIX    *spLU;
MATRIX     *spX;

       /* [a] LUP Decomposition followed by Forward/Backward Substitution */

       switch((int) spA->eRep ) {
           case INDIRECT:
                switch((int) spA->eType) {
                    case DOUBLE_ARRAY:
                         if(spPivot != (VECTOR *)NULL) {
                            VectorFreeDouble(spPivot);
                            spPivot = (VECTOR *)NULL;
                         }
                         spPivot = SetupPivotVector( spA );
                         spLU    = LUDecompositionIndirect( spA , spPivot);
                         spX     = LUSubstitutionIndirect( (char *)NULL, spPivot, spLU, spB);
                         break;
                    default:
                         FatalError("In MatrixSolve() : Undefined A->eType", (char *) NULL);
                         break;
                }
                break;
           case SKYLINE:
                spLU = LUDecompositionSkyline( spA );
                spX  = LUBacksubstitutionSkyline( spLU, spB );
                break;
           case SPARSE:
                break;
          default:
                FatalError("In MatrixSolve() : Undefined spA->eRep",
                          (char *) NULL);
                break;
       }

       MatrixFree( spLU );
       return ( spX );
}

/*
 *  ========================================================================
 *  MatrixLU() : Compute LU decompostion with Pivoting of Rows and U_ii = 1.
 * 
 *  Input :    MATRIX spA  -- Pointer to Matrix A
 *  Output :   MATRIX spLU -- Pointer to decomposed matrix LU
 *  ========================================================================
 */


#ifdef __STDC__
MATRIX *MatrixLU( MATRIX *spA )
#else  /* start case not STDC */
MATRIX *MatrixLU( spA )
MATRIX *spA;
#endif /* end case not STDC */
{
MATRIX    *spLU;

       /* [a] : LUP Decomposition */

       switch((int) spA->eRep) {
           case INDIRECT:
                switch((int) spA->eType) {
                    case DOUBLE_ARRAY:
                         if(spPivot != (VECTOR *)NULL) {
                            VectorFreeDouble(spPivot);
                            spPivot = (VECTOR *)NULL;
                         }
                         spPivot = SetupPivotVector( spA );
                         spLU    = LUDecompositionIndirect( spA , spPivot);
                         break;
                    default:
                         FatalError("In MatrixLU() : Undefined A->eType", (char *) NULL);
                         break;
                }
                break;
           case SKYLINE:
                spLU = LUDecompositionSkyline( spA );
                break;
           case SPARSE:
                break;
          default:
                FatalError("In MatrixLU() : Undefined spA->eRep",
                          (char *) NULL);
                break;
       }

       return ( spLU );
}

/*
 *  =======================================================
 *  MatrixFB() : Compute forward and backward substitution.
 * 
 *  Input :    MATRIX spLU -- Pointer to Matrix A
 *             MATRIX spB  -- Pointer to r.h.s. matrix B
 *  Output :   MATRIX spX  -- Pointer to solution matrix X
 *  =======================================================
 */

#ifdef __STDC__
MATRIX *MatrixFB( MATRIX *spLU, MATRIX *spB )
#else  /* start case not STDC */
MATRIX *MatrixFB( spLU, spB )
MATRIX *spLU;
MATRIX  *spB;
#endif /* end case not STDC */
{
MATRIX *spX;

       /* [a] : Forward/Backward Substitution */

       switch((int) spLU->eRep) {
           case INDIRECT:
                switch((int) spLU->eType) {
                    case DOUBLE_ARRAY:
                         spX = LUSubstitutionIndirect((char *) NULL, spPivot, spLU, spB);
                         break;
                    default:
                         FatalError("In MatrixFB() : Undefined lu->eType",
                                   (char *) NULL);
                         break;
                }
                break;
           case SKYLINE:
                spX = LUBacksubstitutionSkyline( spLU, spB );
                break;
          default:
                FatalError("In MatrixFB() : Undefined spLU->eRep",
                          (char *) NULL);
                break;
       }

       return ( spX );
}

/*
 *  =============================================================
 *  MatrixInverse() : Compute matrix inverse
 * 
 *  Input :    MATRIX spA    -- Pointer to Matrix A
 *  Output :   MATRIX spAinv -- Pointer to inversed matrix spAinv
 *  =============================================================
 */

#ifdef __STDC__
MATRIX *MatrixInverse( MATRIX *spA )
#else  /* start case not STDC */
MATRIX *MatrixInverse( spA )
MATRIX  *spA;
#endif /* end case not STDC */
{
MATRIX  *spAinv;

       switch((int) spA->eRep) {
           case INDIRECT:
                switch((int) spA->eType) {
                    case DOUBLE_ARRAY:
                         spAinv = MatrixInverseIndirectDouble( spA );
                         break;
                    default:
                         FatalError("In MatrixInverse() : Undefined spA->eType",
                                   (char *) NULL);
                         break;
                }
                break;
           case SKYLINE:
                spAinv = MatrixInverseSkyline( spA );
                break;
          default:
                FatalError("In MatrixInverse() : Undefined spA->eRep",
                          (char *) NULL);
                break;
       }

       return ( spAinv );
}

/*
 *  =======================================================
 *  MatrixDet() : Compute determinant of Matrix
 * 
 *  Input :    MATRIX     m  -- Pointer to Matrix m
 *  Output :   QUANTITY   q  -- Matrix Determinant.
 *  =======================================================
 */

#ifdef __STDC__
QUANTITY  *MatrixDet(MATRIX *m)
#else
QUANTITY  *MatrixDet(m)
MATRIX    *m;
#endif
{
MATRIX   *m1;
QUANTITY  *q;

#ifdef DEBUG
       printf("*** Enter MatrixDet() : m->iNoRows    = %4d\n", m->iNoRows);
       printf("                      : m->iNoColumns = %4d\n", m->iNoColumns);
#endif
      q         = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
      if(CheckUnits() == ON) {
         q->dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(q->dimen);
      }
      else
         q->dimen = (DIMENSIONS *)NULL;

      m1 = MatrixCopy(m);
      if( m1->eRep==SKYLINE )   m1 = MatrixSkylineToIndirect(m1);
          
      q->value = dMatrixDet( m1->uMatrix.daa, m1->iNoRows, m1->iNoColumns );

      MatrixFree(m1);

#ifdef DEBUG
       printf("*** Leave MatrixDet()\n");
#endif

   return (q);
}

/*
 *  ===========================================================
 *  MatrixMax() : Get the maximum number in the matrix elements
 *  MatrixMin() : Get the minimum number in the matrix elements
 * 
 *  Input :    MATRIX     m  -- Pointer to Matrix m
 *  Output :   QUANTITY   q  -- Max or Min Number
 *  ===========================================================
 */

#ifdef __STDC__
QUANTITY  *MatrixMax(MATRIX *m)
#else
QUANTITY  *MatrixMax(m)
MATRIX    *m;
#endif
{
int     i, j;
QUANTITY  *q;

      q  = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
      if(CheckUnits() == ON) {
         q->dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(q->dimen);
      }
      else
         q->dimen = (DIMENSIONS *)NULL;

      q->value = m->uMatrix.daa[0][0];
      for( i=1 ; i <= m->iNoRows ; i++ ) {
          for( j=1 ; j <= m->iNoColumns ; j++ ) {
              q->value = MAX( q->value, m->uMatrix.daa[i-1][j-1] );
          }
      }

   return (q);
}

#ifdef __STDC__
QUANTITY  *MatrixMin(MATRIX *m)
#else
QUANTITY  *MatrixMin(m)
MATRIX    *m;
#endif
{
int     i, j;
QUANTITY  *q;

      q = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
      if(CheckUnits() == ON) {
         q->dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(q->dimen);
      }
      else
         q->dimen = (DIMENSIONS *)NULL;

      q->value = m->uMatrix.daa[0][0];
      for( i=1 ; i <= m->iNoRows ; i++ ) {
          for( j=1 ; j <= m->iNoColumns ; j++ ) {
              q->value = MIN( q->value, m->uMatrix.daa[i-1][j-1] );
          }
      }

   return (q);
}

QUANTITY *MatrixL2Norm(m)
MATRIX   *m;
{
QUANTITY *q;
double   sum = 0.0;
int i;

#ifdef DEBUG
       printf("*** Enter MatrixL2Norm() : m->iNoRows    = %4d\n", m->iNoRows);
       printf("                         : m->iNoColumns = %4d\n", m->iNoColumns);
#endif

   /* Check that matrix is either a column (or row) vector */

      if((m->iNoRows != 1) && (m->iNoColumns != 1)) {
          FatalError("In MatrixL2Norm() : Matrix is not a row (or column) vector",(char *)NULL);
      }

   /* Compute L2 Norm for column/row vector */

   sum       =  dVmatrixL2Norm(m->uMatrix.daa, m->iNoRows, m->iNoColumns);
   q         = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
   q->value  = sum;
   if(CheckUnits() == ON) {
      q->dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
      ZeroUnits(q->dimen);
   } else
      q->dimen  = (DIMENSIONS *)NULL;

   return (q);
}

/* 
 *  =================
 *  Matrix Dimensions 
 *  ================= 
 */ 

#ifdef __STDC__
MATRIX *MatrixDimension(MATRIX *m1)
#else
MATRIX *MatrixDimension(m1)
MATRIX *m1;
#endif
{
MATRIX *m2;

      m2 = MatrixAllocIndirect("Dimensions", DOUBLE_ARRAY, 1, 2);
      m2->uMatrix.daa[0][0]  = m1->iNoRows;
      m2->uMatrix.daa[0][1]  = m1->iNoColumns;
      if(CheckUnits()==ON) {
         ZeroUnits(&(m2->spRowUnits[0]));
         ZeroUnits(&(m2->spColUnits[0]));
         ZeroUnits(&(m2->spColUnits[1]));
      }

      return (m2);
}

/* 
 *  ================================================= 
 *  MatrixScale() : Scale a matrix by a double factor
 *  ================================================= 
 */ 

#ifdef __STDC__
MATRIX *MatrixScale(MATRIX *spA, double factor )
#else
MATRIX *MatrixScale(spA, factor )
MATRIX  *spA;
double  factor;
#endif
{
MATRIX *spB;

#ifdef DEBUG
       printf("*** Enter MatrixScale() : spA->iNoRows    = %4d\n", spA->iNoRows);
       printf("                        : spA->iNoColumns = %4d\n", spA->iNoColumns);
#endif
 
       switch((int) spA->eRep) {
           case INDIRECT:
                switch((int) spA->eType) {
                    case DOUBLE_ARRAY:
                         spB = MatrixScaleIndirectDouble( spA, factor );
                         break;
                    default:
                         FatalError("In MatrixScale() : Undefined spA->eType",
                                   (char *) NULL);
                         break;
                }
                break;
           case SKYLINE:
                spB = MatrixScaleSkyline( spA, factor );
                break;
          default:
                FatalError("In MatrixScale() : Undefined spA->eRep",
                          (char *) NULL);
                break;
       }

#ifdef DEBUG
       printf("*** Leave MatrixScale()\n");
#endif
       return ( spB );
}

#ifdef __STDC__
double MatrixContentScale(MATRIX *spA, int row_no, int col_no )
#else
double MatrixContentScale(spA, row_no, col_no )
MATRIX  *spA;
int row_no;        /* row number    */
int col_no;        /* column number */
#endif
{
double  da;

#ifdef DEBUG
       printf("*** Enter MatrixContentScale() : spA->iNoRows    = %4d\n", spA->iNoRows);
       printf("                        : spA->iNoColumns = %4d\n", spA->iNoColumns);
#endif
      if(CheckUnits()==OFF) 
         FatalError("You have to set units ON to use this function","In MatrixScale",(char *)NULL );
 
     /* [a] Scale for Column of matrix spA */

       switch((int) spA->eRep) {
           case INDIRECT:
                switch((int) spA->eType) {
                    case DOUBLE_ARRAY:
                         da = MatrixContentScaleIndirectDouble( spA, row_no, col_no );
                         break;
                    default:
                         FatalError("In MatrixContentScale() : Undefined spA->eType",
                                   (char *) NULL);
                         break;
                }
                break;
           case SKYLINE:
                da = MatrixContentScaleSkyline( spA, row_no, col_no );
                break;
          default:
                FatalError("In MatrixContentScale() : Undefined spA->eRep",
                          (char *) NULL);
                break;
       }

#ifdef DEBUG
       printf("*** Leave MatrixContentScale()\n");
#endif

       return ( da );
}

/* ======================================
   Operations between QUANTITY and MATRIX
   ====================================== */

#ifdef __STDC__
MATRIX *MatrixQuanMult(QUANTITY *q1, MATRIX *m2)
#else
MATRIX *MatrixQuanMult(q1, m2)
QUANTITY *q1;
MATRIX   *m2;
#endif
{
MATRIX     *m3;
int      i,j,k;
int     length;
DIMENSIONS *d1;

   /* [a] : Multiply a quantity to a matrix */

   m3 = MatrixScale( m2, q1->value );
       
   if( CheckUnits() == ON ) {
       for(i = 1; i <= m2->iNoRows; i++) {
           UnitsCopy( &(m3->spRowUnits[i-1]), &(m2->spRowUnits[i-1]) );
           for(j = 1; j <= m2->iNoColumns; j++)
               UnitsMultRep( &(m3->spColUnits[j-1]), q1->dimen, &(m2->spColUnits[j-1]) );
       }
   }

   return (m3);
}

#ifdef __STDC__
MATRIX *MatrixQuanDiv(MATRIX *m2, QUANTITY *q1)
#else
MATRIX *MatrixQuanDiv(m2, q1)
QUANTITY *q1;
MATRIX   *m2;
#endif
{
MATRIX         *m3;
DIMENSIONS     *d1;
int          i,j,k;
int         length;

    /* [a] : matrix is divdied by a quantity */

    m3 = MatrixScale( m2, 1.0/q1->value );
       
    if( CheckUnits()==ON ) {
        for(i = 1; i <= m2->iNoRows; i++) {
            UnitsCopy( &(m3->spRowUnits[i-1]), &(m2->spRowUnits[i-1]) );
            for(j = 1; j <= m2->iNoColumns; j++)
                UnitsDivRep( &(m3->spColUnits[j-1]), &(m2->spColUnits[j-1]), q1->dimen, YES);
        }
    }

    return (m3);
}

/* 
 *  =================================================
 *  QuantityCast() : convert 1x1 matrix into quantity 
 *  =================================================
 */ 

#ifdef __STDC__
QUANTITY *QuantityCast( MATRIX *m )
#else
QUANTITY *QuantityCast(m)
MATRIX   *m;
#endif
{
QUANTITY         *q;
int          length;

   /* [a] : Check that matrix is either a column (or row) vector */

   if((m->iNoRows != 1) || (m->iNoColumns != 1)) {
       FatalError("In Quantity_Cast() : Matrix is not a 1x1 Matrix",(char *)NULL);
   }

   q            = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
   q->value     = m->uMatrix.daa[0][0];

   switch(CheckUnits()) {
      case ON:
           q->dimen = UnitsMult( &(m->spRowUnits[0]), &(m->spColUnits[0]) );
           break;
      case OFF:
           q->dimen = (DIMENSIONS *)NULL;
           break;
   }
     
   return (q);
}


/* 
 *  ==========================================================
 *  MatrixColumnUnits() : Specify Column Units in a Matrix
 *  
 *  Usage:  1st Argument is pointer to matrix concerned.
 *          2nd Argument is matrix of units to be applied.
 *          3rd Argument column number for units.
 *  
 *  Examples : matrix = ColumnUnits( matrix, [ kg ] );
 *             matrix = ColumnUnits( matrix, [ kg ], [2] );
 *  ==========================================================
 */ 

#ifdef  __STDC__
MATRIX *MatrixColumnUnits(MATRIX *first_matrix, ...) {
#else
MATRIX *MatrixColumnUnits(va_alist)
va_dcl
{
MATRIX     *first_matrix;
#endif

va_list          arg_ptr;
MATRIX                *p;
MATRIX *m, *m1, *units_m;
int     i, j, k, counter;
int    NO_COL_WITH_UNITS;
int NO_COL_WITH_NO_UNITS;
int               *COUNT;
int           ID, COL_NO;
int               length;
int       iFLAG  = FALSE;
int       iFLAG1 = FALSE;

   /* [a] : Check units flags and get appropriate units */

   if( CheckUnits() == OFF) {
       FatalError("In MatrixColumnUnits () : To call this function, units must be set to ON ",
                 (char *)NULL );
   }

   /* [b] : Get appropriate units */

#ifdef  __STDC__
    va_start(arg_ptr, first_matrix);
#else
    va_start(arg_ptr);
    first_matrix = va_arg(arg_ptr, MATRIX *);
#endif

    m1 = first_matrix;
    units_m = va_arg(arg_ptr, MATRIX *);
    p       = va_arg(arg_ptr, MATRIX *);
    m       = MatrixCopyIndirectDouble(m1);
    va_end(arg_ptr);

   /* [c] : Retrieve and check appropriate units */

   if( p != NULL) {
       ID = 1; 
       COL_NO  = p->uMatrix.daa[0][0]; 
   } else
       ID = 2;

   if(units_m->iNoRows != 1) 
      FatalError("Fatal Error in MatrixColumnUnits(): Units_matrix should be a 1xn vector",
                (char *) NULL); 

   /* [d] : Case 1 : Update "column units" in one column only.              */
   /*                                                                       */
   /*       Two cases exist -- (1.1) Column i of matrix m has no units      */
   /*                          (1.2) Column i of matrix m has units.        */
   /*                                Check compatibility with units_m.      */
   /*                                Copy unmits_m to m.                    */
   /*                                                                       */
   /*       Examples : stiff = Matrix([2,2]);                               */
   /*                  stiff = ColumnUnits ( stiff, [  sec], [1] );         */
   /*                  stiff = ColumnUnits ( stiff, [sec^2], [2] );         */
   /*                                                                       */

   if (ID == 1) {

      if(units_m->iNoColumns != 1) 
         FatalError("In MatrixColumnUnits(): too many units for one column ",(char *)NULL);

      i = COL_NO; /* column no i is to be assigned an unit  */

      if(m->spColUnits[i-1].length_expnt == 0 &&
         m->spColUnits[i-1].time_expnt == 0   &&
         m->spColUnits[i-1].mass_expnt == 0   &&
         m->spColUnits[i-1].temp_expnt == 0) {

         /* Case (1.1) : Assign units from units_m to m */
          
         UnitsCopy( &(m->spColUnits[i-1]), &(units_m->spColUnits[0]) );
         for(j = 1; j <= m->iNoRows; j++)
             m->uMatrix.daa[j-1][i-1] = 
                units_m->spColUnits[0].scale_factor*m->uMatrix.daa[j-1][i-1];
      } else {

         /* Case (1.2) : Check with units_m, and copy units_m to units of m */
           
         if( SameUnits( &(m->spColUnits[i-1]), &(units_m->spColUnits[0])) == TRUE ){
             UnitsCopy( &(m->spColUnits[i-1]), &(units_m->spColUnits[0]) );
         } else {
             FatalError("In MatrixColumnUnits(): Inconsistent units",
                       (char *) NULL );
         }
      }
   }

   /* [e] : Case 2 : Update "units" in multiple matrix columns.                      */
   /*                                                                                */
   /*       Cases -- (2.1) Matrix m and units_m have the same number of columns.     */
   /*                (2.1.1) : Matrix m has no units. Assign units from units_m to m */
   /*                (2.1.2) : Matrix m has units. Copy units_m to units of m.       */
   /*                (2.2) Matrices m and units_m have different number of columns.  */
   /*                (2.2.1) : matrix m rows have no units. Assign all columns of m  */
   /*                          with units of units_m                                 */
   /*                (2.2.2) : Assign unit_m to columns of m where they match.       */
   /*                                                                                */
   /*       Examples : stiff = Matrix([3,3]);                                        */
   /*                  stiff = ColumnUnits ( stiff, [cm/sec^2], [1] );               */
   /*                  stiff = ColumnUnits ( stiff, [sec] );                         */
   /*                  stiff = ColumnUnits ( stiff, [in/sec^2] );                    */
   /*                                                                                */

   if (ID == 2) {

      /* Case (2.1) : Matrices m and units_m have same number of columns */

      if(units_m->iNoColumns == m->iNoColumns) {
         for (i = 1; i <= m->iNoColumns; i++) { 

            if(m->spColUnits[i-1].length_expnt == 0 && 
               m->spColUnits[i-1].time_expnt == 0   &&
               m->spColUnits[i-1].mass_expnt == 0   &&
               m->spColUnits[i-1].temp_expnt == 0) {

               /* Case (2.1.1) : Matrix m has no units -- Assign units from units_m to m */

               UnitsCopy( &(m->spColUnits[i-1]), &(units_m->spColUnits[i-1]) );
               for(j = 1; j <= m->iNoRows; j++)
                   m->uMatrix.daa[j-1][i-1] = 
                      units_m->spColUnits[i-1].scale_factor*m->uMatrix.daa[j-1][i-1];

            } else {

               /* Case (2.1.2) : Matrix m has units. Check, and copy units_m to units of m */
           
               if(SameUnits(&(m->spColUnits[i-1]), &(units_m->spColUnits[i-1])) == TRUE ){
                  UnitsCopy( &(m->spColUnits[i-1]), &(units_m->spColUnits[i-1]) );
               } else {
                  FatalError("Fatal Error in MatrixColumnUnits(): inconsistent units ",
                            (char *)NULL);
               }
            }
         }
      }

      /* Case (2.2) : Matrices m and units_m have different number of columns */

      if(units_m->iNoColumns != m->iNoColumns) {

         k = 0;
         COUNT = (int *) MyCalloc(m->iNoColumns, sizeof(int));
         counter = 0;

         /* Find and record column numbers with matching units */

         for (i = 1; i <= m->iNoColumns; i++) {
            if(m->spColUnits[i-1].length_expnt != 0 ||
               m->spColUnits[i-1].time_expnt   != 0 ||
               m->spColUnits[i-1].mass_expnt   != 0 ||
               m->spColUnits[i-1].temp_expnt   != 0) {
                k = k + 1;
                COUNT[k-1] = i; 
            }
         }

         NO_COL_WITH_UNITS = k;
         NO_COL_WITH_NO_UNITS = m->iNoColumns - NO_COL_WITH_UNITS;

         if( NO_COL_WITH_NO_UNITS == m->iNoColumns) {

            /* Case (2.2.1) : matrix m rows have no units. Assign all columns of m */
            /*                with units of units_m                                */

            if(units_m->iNoColumns != 1)
               FatalError(" in MatrixColumnUnits(): Too many/few units ",
                         (char *)NULL);

            if(units_m->iNoColumns == 1) { 
               for(i = 1; i <= m->iNoColumns; i++) {
                   UnitsCopy( &(m->spColUnits[i-1]), &(units_m->spColUnits[0]) );
                   for(j = 1; j <= m->iNoRows; j++)
                       m->uMatrix.daa[j-1][i-1] *= units_m->spColUnits[0].scale_factor;
               }
            }
         }

         if( NO_COL_WITH_NO_UNITS != m->iNoColumns) {

            /* Case (2.2.2) : Assign unit_m to columns of m where they match */

            for( k = 1; k <=  NO_COL_WITH_UNITS; k++) {
               i = COUNT[k-1];
               for(j = 1; j <= units_m->iNoColumns; j++) {
                   if(SameUnits(&(m->spColUnits[i-1]), &(units_m->spColUnits[j-1])) == TRUE) {
                      UnitsCopy( &(m->spColUnits[i-1]), &(units_m->spColUnits[j-1]) );
                      iFLAG1 = TRUE;
                   } else 
                      iFLAG1 = MAX(iFLAG1, FALSE);

                   if(iFLAG1 == TRUE)
                      break;
               }
               if(iFLAG1 == TRUE)
                  iFLAG = TRUE;
               iFLAG1 = FALSE;
            }

            if(iFLAG == FALSE)
               FatalError("In MatrixColumnUnits() : Inconsistent units", (char *)NULL);
         }

         free((char *) COUNT);
      }
   }

   /* [f] : Initialize non_units row buffer */

   for( i = 1; i <= m->iNoRows; i++)
        if(m->spRowUnits[i-1].units_name == (char *)NULL)
           ZeroUnits(&(m->spRowUnits[i-1]));

   return (m);
}


/* 
 *  ==========================================================
 *  MatrixRowUnits() : Specify Row Units in a Matrix
 *  
 *  Usage:  1st Argument is pointer to matrix concerned.
 *          2nd Argument is matrix of units to be applied.
 *          3rd Argument row number for units.
 *  
 *  Examples : matrix = RowUnits( matrix, [ kg ] );
 *             matrix = RowUnits( matrix, [ kg ], [2] );
 *  ==========================================================
 */ 

#ifdef  __STDC__
MATRIX *MatrixRowUnits(MATRIX *first_matrix, ...) {
#else
MATRIX *MatrixRowUnits(va_alist)
va_dcl
{
MATRIX     *first_matrix;
#endif

va_list          arg_ptr;
MATRIX                *p;
MATRIX *m, *m1, *units_m;
int     i, j, k, counter;
int    NO_ROW_WITH_UNITS;
int NO_ROW_WITH_NO_UNITS;
int               *COUNT;
int           ID, ROW_NO;
int               length;
int       iFLAG1 = FALSE;
int       iFLAG  = FALSE;

   /* [a] : Check units flags and get appropriate units */

   if( CheckUnits() == OFF) {
       FatalError("In MatrixRowUnits () : To call this function, units must be set to ON ",
                 (char *)NULL );
   }

   /* [b] : Get appropriate units */

#ifdef  __STDC__
    va_start(arg_ptr, first_matrix);
#else
    va_start(arg_ptr);
    first_matrix = va_arg(arg_ptr, MATRIX *);
#endif

   /* [c] : Retrieve and check appropriate units */

   m1 = first_matrix;
   units_m = va_arg(arg_ptr, MATRIX *);
   p       = va_arg(arg_ptr, MATRIX *);
   va_end(arg_ptr);

   m = MatrixCopyIndirectDouble(m1);

   if( p != NULL) {
       ID = 1;
       ROW_NO  = p->uMatrix.daa[0][0];
   } else
       ID = 2;

   if(units_m->iNoRows < 1) 
      FatalError("In MatrixRowUnits() : Units_matrix should be at least 1xn matrix",
                (char *) NULL); 

   /* [d] : Case 1 : Update "row units" in one row only.                    */
   /*                                                                       */
   /*       Two cases exist -- (1.1) Row i of matrix m has no units         */
   /*                          (1.2) Row i of matrix m has units.           */
   /*                                Check compatibility with units_m.      */
   /*                                Copy unmits_m to m.                    */
   /*                                                                       */
   /*       Examples : stiff = Matrix([2,2]);                               */
   /*                  stiff = RowUnits ( stiff, [  sec], [1] );            */
   /*                  stiff = RowUnits ( stiff, [sec^2], [2] );            */
   /*                                                                       */

   if (ID == 1) {

      if(units_m->iNoColumns > 1) 
         FatalError("In MatrixRowUnits() : Too many units for one column ",(char *)NULL);

      i = ROW_NO;  /* Row no i is to be assigned an unit  */

      /* Case [1] : row_no i of matrix m has no units  */
      /*            assign units from units_m to m     */

      if( m->spRowUnits[i-1].length_expnt == 0 && m->spRowUnits[i-1].time_expnt == 0 &&
          m->spRowUnits[i-1].mass_expnt == 0   && m->spRowUnits[i-1].temp_expnt == 0) {
          UnitsCopy( &(m->spRowUnits[i-1]), &(units_m->spColUnits[0]) );
          for(j = 1;j <= m->iNoColumns; j++) 
                   m->uMatrix.daa[i-1][j-1]  *=  units_m->spColUnits[0].scale_factor; 
      } else {

          /* Case [2] : row no i of  matrix m has units. check with  */
          /*            units_m and copy units_m to units of m       */ 

          if( SameUnits(&(m->spRowUnits[i-1]), &(units_m->spColUnits[0])) == TRUE) {
              UnitsCopy( &(m->spRowUnits[i-1]), &(units_m->spColUnits[0]) );
          } else {
              FatalError("In MatrixRowUnits() : Inconsistent units ",
                        (char *) NULL);
          }
      }
   }

   /* [e] : Case 2 : Update "units" in multiple matrix rows.                         */
   /*                                                                                */
   /*       Cases -- (2.1) Matrix m and units_m have the same number of rows.        */
   /*                (2.1.1) : Matrix m has no units. Assign units from units_m to m */
   /*                (2.1.2) : Matrix m has units. Copy units_m to units of m.       */
   /*                (2.2) Matrices m and units_m have different number of rows.     */
   /*                (2.2.1) : matrix m rows have no units. Assign all rows of m     */
   /*                          with units of units_m                                 */
   /*                (2.2.2) : Assign unit_m to rows of m, where they match.         */
   /*                                                                                */
   /*       Examples : stiff = Matrix([3,3]);                                        */
   /*                  stiff = RowUnits ( stiff, [cm/sec^2], [1] );                  */
   /*                  stiff = RowUnits ( stiff, [sec] );                            */
   /*                  stiff = RowUnits ( stiff, [in/sec^2] );                       */
   /*                                                                                */

   if (ID == 2) {

       /* CASE[a]: Matrix m's row number is the same as Units Matrix */
       /*          units_m's column number                           */

       if(units_m->iNoColumns == m->iNoRows) {

          /* Case [a.1] : matrix m has no units          */
          /*              assign units from units_m to m */

          for(i = 1; i <= m->iNoRows; i++) { 
              if(m->spRowUnits[i-1].length_expnt == 0 &&
                   m->spRowUnits[i-1].time_expnt == 0   &&
                   m->spRowUnits[i-1].mass_expnt == 0   &&
                   m->spRowUnits[i-1].temp_expnt == 0) {
                   UnitsCopy( &(m->spRowUnits[i-1]), &(units_m->spColUnits[i-1]) );
                   for(j = 1; j <= m->iNoColumns; j++)
                       m->uMatrix.daa[i-1][j-1] *= units_m->spColUnits[i-1].scale_factor;
              } else {

                  /* Case [a.2] : matrix m has units. check with units_m     */
                  /*              and copy units_m to units of m             */ 

                  if(SameUnits(&(m->spRowUnits[i-1]), &(units_m->spColUnits[i-1])) == TRUE) {
                     UnitsCopy( &(m->spRowUnits[i-1]), &(units_m->spColUnits[i-1]) );
                  } else {
                     FatalError("In MatrixRowUnits() : Inconsistent units ",
                               (char *) NULL);
                  }
              }
          }
       }

       /* CASE [b] : No of rows in matrix m and No of */
       /*            columns in units_m are different */

       if( units_m->iNoColumns != m->iNoRows) {
           k = 0;
           COUNT = (int *) MyCalloc(m->iNoRows, sizeof(int));
           counter = 0;

           for(i = 1; i <= m->iNoRows; i++) {
               if(m->spRowUnits[i-1].length_expnt != 0 ||
                  m->spRowUnits[i-1].time_expnt   != 0 ||
                  m->spRowUnits[i-1].mass_expnt   != 0 ||
                  m->spRowUnits[i-1].temp_expnt   != 0) {
                  k = k + 1;
                  COUNT[k-1] = i;     /* record rows with units */
                }
           }

           NO_ROW_WITH_UNITS = k;
           NO_ROW_WITH_NO_UNITS = m->iNoRows - NO_ROW_WITH_UNITS;

           /* all rows of m have no units */

           if(NO_ROW_WITH_NO_UNITS == m->iNoRows) {

               /* Case [b.1] : matrix m rows have no units */
               if(units_m->iNoColumns == 1){
                  for(i = 1; i <= m->iNoRows; i++) {
                      UnitsCopy( &(m->spRowUnits[i-1]), &(units_m->spColUnits[0]) );
                      for(j = 1; j <= m->iNoColumns; j++)
                          m->uMatrix.daa[i-1][j-1] *= units_m->spColUnits[0].scale_factor;
                  }
               }
               else  
                   FatalError(" in MatrixRowUnits(): Too many units ",(char *)NULL);
           }

           /* Some row of matrix m have units */

           if(NO_ROW_WITH_NO_UNITS != m->iNoRows) {

              /* Case [b.2] : Matrix m has units. Search & check with units_m   */
              /*              Selectively copy units of units_m to units of m   */

              for(k = 1; k <=  NO_ROW_WITH_UNITS; k++){
                   i = COUNT[k-1];
                   for(j = 1; j <= units_m->iNoColumns; j++) {
                       if(SameUnits(&(m->spRowUnits[i-1]), &(units_m->spColUnits[j-1])) == TRUE){
                          UnitsCopy( &(m->spRowUnits[i-1]), &(units_m->spColUnits[j-1]) );
                          iFLAG1 = TRUE;
                       }
                       else
                          iFLAG1 = MAX(iFLAG1, FALSE);
                       if(iFLAG1 == TRUE) break;
                   }
                   if(iFLAG1 == TRUE) iFLAG = TRUE;
                   iFLAG1 = FALSE;
              }
              if(iFLAG == FALSE)
                 FatalError(" Fatal Error in MatrixRowUnits(): inconsistent units ",
                           (char *) NULL);
           }
           free((char *) COUNT);
        }
   }

   /* [f] : Initialize non_units column buffer */

   for( i = 1; i <= m->iNoColumns; i++) { 
        if(m->spColUnits[i-1].units_name == (char *) NULL)
           ZeroUnits(&(m->spColUnits[i-1]));
   }

   return (m);
}

/*
 *  ==================================================
 *  MatrixZeroUnits() : Initialize matrix units buffer
 *  
 *  Input :  
 *  Ouput :  
 *  ==================================================
 */

#ifdef __STDC__
MATRIX *MatrixZeroUnits(MATRIX *m, int iNoRows, int no_cols)
#else
MATRIX *MatrixZeroUnits(m, iNoRows, no_cols)
MATRIX *m;
int iNoRows, no_cols;
#endif
{

      m = (MATRIX *) MyMalloc(sizeof(MATRIX));

      m->cpMatrixName = "Zero_units_Matrix";
      m->iNoRows      = iNoRows;
      m->iNoColumns   = no_cols;
      m->spRowUnits   = BufferInit(m->iNoRows);
      m->spColUnits   = BufferInit(m->iNoColumns);
      m->eType        = DOUBLE_ARRAY;
      m->uMatrix.daa  = MatrixAllocIndirectDouble(m->iNoRows, m->iNoColumns);
       
      return (m);
}

/*
 *  =======================================================
 *  Change a matrix with units into a matrix without units.
 * 
 *  Input : 
 *  Output: 
 *  =======================================================
 */

#ifdef __STDC__
MATRIX *MatrixUnitsLess(MATRIX *m)
#else
MATRIX *MatrixUnitsLess(m)
MATRIX *m;
#endif
{
int  i,j,k;
MATRIX *m1;

   for(i = 1; i <= m->iNoRows; i++)
       ZeroUnits(&(m->spRowUnits[i-1]));

   for(j = 1; j <= m->iNoColumns; j++) 
       ZeroUnits(&(m->spColUnits[j-1]));

   m1 = MatrixCopy(m);

   return (m1); 
}

/*
 *  ====================================================
 *  MatrixUnitsSimplify() : Simplify matrix units buffer
 * 
 *  Input  : 
 *  Output : 
 *  ====================================================
 */

#ifdef __STDC__
MATRIX *MatrixUnitsSimplify(MATRIX *m)
#else
MATRIX *MatrixUnitsSimplify(m)
MATRIX *m;
#endif
{
DIMENSIONS *d;
int   i, j, k;
int    length;

   /* check whether the m is a vector */

   if(m->iNoColumns == 1 || m->iNoRows == 1) {
      if(m->iNoRows == 1) { /* row vector */
         for(i = 1; i <= m->iNoColumns; i++) {
             UnitsMultRep( &(m->spColUnits[i-1]), &(m->spRowUnits[0]), &(m->spColUnits[i-1]) );
             UnitsSimplify(&(m->spColUnits[i-1]));
         }
         ZeroUnits( &(m->spRowUnits[0]) );
      }

      if(m->iNoColumns == 1 && m->iNoRows != 1) { /* column vector */
         for(i = 1; i <= m->iNoRows; i++) {
             UnitsMultRep( &(m->spRowUnits[i-1]), &(m->spRowUnits[i-1]), &(m->spColUnits[0]) );
             UnitsSimplify(&(m->spRowUnits[i-1]));
         }
         ZeroUnits( &(m->spColUnits[0]) );
      }
   } else {

      /* m is not a vector */

      for(i = 1; i <= m->iNoRows; i++) {
          UnitsSimplify( &(m->spRowUnits[i-1]) );
      }
      for(j = 1; j <= m->iNoColumns; j++) {
          UnitsSimplify( &(m->spColUnits[j-1]) );
      }
   }
       
   return (m);
}


/* 
 *  ======================================================
 *  MatrixExtract():  Extract part of matrix and assign  
 *                    to another matrix                  
 *  Usage:                                               
 *      Copy sub-matrix of m2, bounded by top_left corner A
 *      m1 = MatrixExtract(m1, m2, A)                       
 *      and A->uMatrix.daa = (i,j);               
 *  ======================================================
 */ 

#ifdef __STDC__
MATRIX *MatrixExtract(MATRIX *m1, ...)
{
#else
MATRIX *MatrixExtract(va_alist)
va_dcl
{
MATRIX *m1;
#endif

va_list         arg_ptr;
MATRIX          *m2, *A;
MATRIX               *m;
int          i, j, k, n;
int  ii, jj, iMin, iMax;

#ifdef DEBUG
     printf("\n Enter MatrixExtract()\n");
#endif 

#ifdef  __STDC__
    va_start(arg_ptr, m1);
#else
    va_start(arg_ptr);
    m1 = va_arg(arg_ptr, MATRIX *);
#endif

    m2 = va_arg(arg_ptr, MATRIX *);
    A  = va_arg(arg_ptr, MATRIX *);
    va_end(arg_ptr);

    m = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, m1->iNoRows, m1->iNoColumns );

    /* [a] : matrix bounds check */

    k = m->iNoRows    + (int) A->uMatrix.daa[0][0] - 1;
    n = m->iNoColumns + (int) A->uMatrix.daa[0][1] - 1;

    if((int) A->uMatrix.daa[0][0] < 1 || (int) A->uMatrix.daa[0][1] < 1 ||
        k > m2->iNoRows || n > m2->iNoColumns ) {
        printf("\n **** The corners are out of the bounds of matrix %s \n\n",
                m2->cpMatrixName);
        FatalError("Check input file at Extract()",
                "Matrix rows/columns are outbounds in MatrixExtract()",(char *)NULL);
    }

    switch( (int)m2->eRep ) {
       case INDIRECT :
          i = (int) A->uMatrix.daa[0][0]; 
          for(ii = 1; ii <= m->iNoRows; ii++) {
             j = (int) A->uMatrix.daa[0][1]; 
             for(jj = 1; jj <= m->iNoColumns; jj++) {
                m->uMatrix.daa[ii-1][jj-1] = m2->uMatrix.daa[i-1][j-1];
                j++;
             }
             i++;
          }
       break;

       case SKYLINE :
          i = (int) A->uMatrix.daa[0][0]; 
          for(ii = 1; ii <= m->iNoRows; ii++) {
             j = (int) A->uMatrix.daa[0][1]; 
             for(jj = 1; jj <= m->iNoColumns; jj++) {
                iMin = MIN( i, j );
                iMax = MAX( i, j );
                if( (iMax-iMin+1) <= m2->uMatrix.daa[iMax-1][0] )
                   m->uMatrix.daa[ii-1][jj-1] = m2->uMatrix.daa[iMax-1][iMax-iMin+1];
                else
                   m->uMatrix.daa[ii-1][jj-1] = 0.0;
                j++;
             }
             i++;
          }
       break;
       default :
       break;
    }

    if( CheckUnits() == ON ) {
       for(ii=1, i=(int)A->uMatrix.daa[0][0]; ii <= m->iNoRows; ii++, i++)
          UnitsCopy(&(m->spRowUnits[ii-1]), &(m2->spRowUnits[i-1]));
       for(jj=1, j=(int)A->uMatrix.daa[0][1]; jj <= m->iNoColumns; jj++, j++)
          UnitsCopy(&(m->spColUnits[jj-1]), &(m2->spColUnits[j-1]));
    }

    return (m);
}


/*
 *  ======================================================
 *  MatrixPut():  Put a matrix into a bigger matrix     
 *                so that the former matrix becomes the 
 *                sub-matrix of the latter matrix       
 *                the content of sub-matrix of the      
 *                bigger matrix is replaced by the new  
 *                submatrix                             
 *  Usage:                                              
 *     Copy matrix m2 into matrix m1 at location bounded 
 *     by top_left corner A;
 *     m1 = MatrixPut(m1, m2, A)                        
 *          and A->uMatrix.daa = (i,j);                  
 *  ======================================================
 */

#ifdef __STDC__
MATRIX *MatrixPut(MATRIX *m1, ...)
{
#else
MATRIX *MatrixPut(va_alist)
va_dcl
{
MATRIX *m1;
#endif

va_list         arg_ptr;
MATRIX          *m2, *A;
MATRIX               *m;
int          i, j, k, n;
int  ii, jj, iMin, iMax;
int        UNITS_SWITCH;

#ifdef DEBUG
     printf("\n Enter MatrixPut()\n");
#endif 

#ifdef  __STDC__
    va_start(arg_ptr, m1);
#else
    va_start(arg_ptr);
    m1 = va_arg(arg_ptr, MATRIX *);
#endif

    m2 = va_arg(arg_ptr, MATRIX *);
    A  = va_arg(arg_ptr, MATRIX *);
    va_end(arg_ptr);
    UNITS_SWITCH = CheckUnits();
    
    m = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, m1->iNoRows, m1->iNoColumns );

    switch( (int)m1->eRep ) {
       case INDIRECT:
          for(ii = 1; ii <= m->iNoRows; ii++)
             for(jj = 1; jj <= m->iNoColumns; jj++)
                m->uMatrix.daa[ii-1][jj-1] = m1->uMatrix.daa[ii-1][jj-1];
       break;
       case SKYLINE:
          for(ii = 1; ii <= m->iNoRows; ii++) {
             for(jj = 1; jj <= m->iNoColumns; jj++) {
                iMin = MIN( ii, jj );
                iMax = MAX( ii, jj );
                if( (iMax-iMin+1) <= m1->uMatrix.daa[iMax-1][0] )
                   m->uMatrix.daa[ii-1][jj-1] = m1->uMatrix.daa[iMax-1][iMax-iMin+1];
                else
                   m->uMatrix.daa[ii-1][jj-1] = 0.0;
             }
          }
       break;
       default:
       break;
    }
    if( UNITS_SWITCH == ON ) {
       MatrixUnitsSimplify( m1 );
       for(ii = 1; ii <= m->iNoRows; ii++)
          UnitsCopy(&(m->spRowUnits[ii-1]), &(m1->spRowUnits[ii-1]));
       for(jj = 1; jj <= m->iNoColumns; jj++)
          UnitsCopy(&(m->spColUnits[jj-1]), &(m1->spColUnits[jj-1]));
    }

    /* [a] matrix bounds check */

    k = m2->iNoRows    + (int) A->uMatrix.daa[0][0] - 1;
    n = m2->iNoColumns + (int) A->uMatrix.daa[0][1] - 1;

    if((int) A->uMatrix.daa[0][0] < 1 || (int) A->uMatrix.daa[0][1] < 1 ||
        k > m1->iNoRows ||  n > m1->iNoColumns ) {
        printf("\n **** The corners are out of the bounds of matrix %s \n\n",m1->cpMatrixName);
        FatalError("Check input file at Put()",
                   "Matrix rows/columns are outbounds in MatrixPut()",(char *)NULL);
    }

    switch( (int)m2->eRep ) {
       case INDIRECT :
          i = (int) A->uMatrix.daa[0][0]; 
          for(ii = 1; ii <= m2->iNoRows; ii++) {
             j = (int) A->uMatrix.daa[0][1]; 
             for(jj = 1; jj <= m2->iNoColumns; jj++) {
                m->uMatrix.daa[i-1][j-1] = m2->uMatrix.daa[ii-1][jj-1];
                j++;
             }
             i++;
          }
       break;

       case SKYLINE :
          i = (int) A->uMatrix.daa[0][0]; 
          for(ii = 1; ii <= m2->iNoRows; ii++) {
             j = (int) A->uMatrix.daa[0][1]; 
             for(jj = 1; jj <= m2->iNoColumns; jj++) {
                iMin = MIN( ii, jj );
                iMax = MAX( ii, jj );
                if( (iMax-iMin+1) <= m2->uMatrix.daa[iMax-1][0] )
                   m->uMatrix.daa[i-1][j-1] = m2->uMatrix.daa[iMax-1][iMax-iMin+1];
                else
                   m->uMatrix.daa[i-1][j-1] = 0.0;
                j++;
             }
             i++;
          }
       break;
       default :
       break;
    }

    if( UNITS_SWITCH == ON ) {
       for(ii=1, i=(int)A->uMatrix.daa[0][0]; ii <= m2->iNoRows; ii++, i++)
          UnitsCopy(&(m->spRowUnits[i-1]), &(m2->spRowUnits[ii-1]));
       for(jj=1, j=(int)A->uMatrix.daa[0][1]; jj <= m2->iNoColumns; jj++, j++)
          UnitsCopy(&(m->spColUnits[j-1]), &(m2->spColUnits[jj-1]));
    }

    return (m);
}


/*
 *  ===================================================================
 *  void MatrixSolveEigen() : Use method of Subspace iteration to
 *                            solve [A][x] = [Lambda][B][x].
 *
 *  Algorithm : [A] and [B] are large (nxn) matrices
 *            : Size of [X_{k}] = [X_{k+1}]   = (nxm) matrix
 *            : Size of [Y_{k}] = [Y_{k+1}]   = (nxm) matrix
 *            : Size of [A_{k+1}] = [B_{k+1}] = (mxm) matrix
 *            : Size of [Eigenvalue_{k+1}] = (mxm) diagonal matrix
 *            : Size of [Q_{k+1}]          = (mxm) matrix
 * 
 *  Input :  spA, spB
 *  Output : spEigenvector
 * 
 *  Written By : M. Austin                             November 1993. 
 *  ===================================================================
 */

enum { MAX_SUBSPACE_ITERATIONS = 50 };

#define TOLERANCE  0.00001

#ifdef __STDC__
void MatrixSolveEigen( MATRIX *spA, MATRIX *spB , MATRIX *spEigenvalue, MATRIX *spEigenvector, int iNoEigen ) 
#else
void MatrixSolveEigen( spA, spB , spEigenvalue, spEigenvector, iNoEigen ) 
MATRIX *spA, *spB;
MATRIX *spEigenvalue;
MATRIX *spEigenvector;
int     iNoEigen;
#endif
{
MATRIX   *spEigenvalueOld;
MATRIX    *spEigenvectorQ;
MATRIX *spAwork, *spBwork;
MATRIX            *spA_LU;
MATRIX               *spY;
MATRIX            *spYnew;
MATRIX    *spV, *spVtrans;
MATRIX     *spX, *spTemp1;
int     iSize, ii, ij, ik;
int                 iLoop;
int  iSubspaceConvergence;
double     dEigenvalueOld;
double       dConvergence;
double        dEigenvalue;
double          dMaxValue;
int          UNITS_SWITCH;

#ifdef DEBUG
       printf("*** In MatrixSolveEigen()\n");
#endif 

    /* [c] : Loops for Subspace Iteration */

       spEigenvalueOld = MatrixCopyIndirectDouble( spEigenvalue );

       iSize = spA->iNoRows;
       iLoop = 0; iSubspaceConvergence = FALSE;
       while( iLoop <= (int) MAX_SUBSPACE_ITERATIONS && iSubspaceConvergence == FALSE ) {

           iLoop = iLoop + 1;

#ifdef DEBUG
       printf("*** In MatrixSolveEigen() : Start Iteration %3d \n", iLoop);
#endif

           /* [c.1] : Compute [Y_{k}] = [B]*[X_{k}] */

           spY = MatrixMult( spB, spEigenvector );

           /* [c.2] : Solve [A][V_{k+1}] = [Y_{k}] */

           if(iLoop == 1) {
              spA_LU = MatrixLU( spA );
              spV = MatrixAllocIndirect((char *)NULL,DOUBLE_ARRAY, iSize, iNoEigen);
           }

           spX = MatrixAllocIndirect((char *)NULL,DOUBLE_ARRAY, iSize, 1);
           for( ii=1 ; ii<=iNoEigen ; ii++ ) {
               for( ij=1 ; ij<=iSize ; ij++ )
                   spX->uMatrix.daa[ij-1][0] = spY->uMatrix.daa[ij-1][ii-1];

               spTemp1 = MatrixFB( spA_LU, spX );
               for( ij=1; ij<=iSize; ij++ )
                   spV->uMatrix.daa[ij-1][ii-1] = spTemp1->uMatrix.daa[ij-1][0];

               MatrixFree( spTemp1 );
           }
           MatrixFree( spX );

           spVtrans = MatrixTranspose( spV );

           /* [c.3] : Compute [Y_{new}] = [A] . [V_{k+1}] */

           spYnew = MatrixMult( spA, spV );

           /* [c.4] : Compute [A_{k+1}] = [V_{k+1}]^T . [Y_{new}] */

           spAwork  = MatrixMult( spVtrans, spYnew );

           /* [c.5] : Compute [Y_{new}] = [B] . [V_{k+1}] */

           MatrixFree( spYnew );
           spYnew = MatrixMult( spB, spV );

           /* [c.6] : Compute [B_{k+1}] = [V_{k+1}]^T . [Y_{new}] */

           spBwork = MatrixMult( spVtrans, spYnew );

           /* [c.7] : Solve [A_{k+1}].[Q_{k+1}] = [B_{k+1}].[Q_{k+1}].[Lambda_{k+1}] */

           spEigenvectorQ  = MatrixAllocIndirect("Eigenvector [Q]", DOUBLE_ARRAY, iNoEigen, iNoEigen);
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
           spEigenvalueOld = MatrixCopyIndirectDouble( spEigenvalue );

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

              MatrixFree ( spTemp1 );
              MatrixFree ( spY );
              MatrixFree ( spYnew );
              MatrixFree ( spAwork );
              MatrixFree ( spBwork );
              MatrixFree ( spEigenvectorQ );
       }

       /* [d] : Summarize Algorithm Performance */

       printf("\n*** SUBSPACE ITERATION CONVERGED IN %2d ITERATIONS \n", iLoop);

       /* [e] : Clean Up before Leaving */

       MatrixFree ( spEigenvalueOld );
       MatrixFree ( spV );
       MatrixFree ( spVtrans );

#ifdef DEBUG
    printf("*** Leaving MatrixSolveEigen()\n");
#endif 

}

#ifdef __STDC__
MATRIX *Solve_Eigen( MATRIX *spA, MATRIX *spB, MATRIX *spN )
#else
MATRIX *Solve_Eigen( spA, spB, spN )
MATRIX *spA, *spB, *spN;
#endif
{
int              iNoEigen;
int     ii, ij, ik, iSize;
MATRIX           *spEigen;
MATRIX      *spEigenvalue;
MATRIX     *spEigenvector;
int                length;
DIMENSIONS  *d,*d1,*dimen;
int          UNITS_SWITCH;
int             UnitsType;

    /* [a] : Check Input and Allocate Working Matrices */

    if( spA == NULL || spB == NULL )
        FatalError("In Solve_Eigen() : Pointer spA == NULL or spB == NULL",
                  (char *) NULL);

    if( spA->eType != DOUBLE_ARRAY || spB->eType != DOUBLE_ARRAY ) 
        FatalError("In Solve_Eigen() : spA->eType != DOUBLE_ARRAY or spB->eType != DOUBLE_ARRAY",
                  (char *) NULL);

    /* [a.1]  Set up size of Eigenvectors */

    iSize    = spA->iNoRows;
    iNoEigen = (int) spN->uMatrix.daa[0][0];

    /* [a.2]  Allocate space for eigenvectors and eigenvalues */

    spEigen       = MatrixAllocIndirect("eigen",  DOUBLE_ARRAY, iSize+1, iNoEigen);
    spEigenvalue  = MatrixAllocIndirect("Eigenvalues",  DOUBLE_ARRAY, iNoEigen, 1);
    spEigenvector = MatrixAllocIndirect("Eigenvectors", DOUBLE_ARRAY, iSize, iNoEigen);

    /* [b] : Initialize Starting Eigenvectors [X_{k}] */

    for(ii = 1; ii <= iNoEigen; ii = ii + 1)
        spEigenvector->uMatrix.daa[ii-1][ii-1] = 1.0;

    if( CheckUnits()==ON ) {
        UNITS_SWITCH = ON;
        UnitsType = CheckUnitsType();
        SetUnitsOff();
    }

    MatrixSolveEigen( spA, spB, spEigenvalue, spEigenvector, iNoEigen );

    for( ij=1 ; ij<=iNoEigen ; ij++ ) {
         spEigen->uMatrix.daa[0][ij-1] = spEigenvalue->uMatrix.daa[ij-1][0];
         for( ii=1 ; ii<=iSize ; ii++ )
              spEigen->uMatrix.daa[ii][ij-1] = spEigenvector->uMatrix.daa[ii-1][ij-1];
    }

    if( UNITS_SWITCH == ON ) {
        d     = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
        d1    = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
        dimen = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
 
        /* Compute units for eigenvalues */

       for( ii=1; ii<=iNoEigen; ii++ ) {
            UnitsMultRep( dimen, &(spA->spRowUnits[ii-1]), &(spA->spColUnits[ii-1]) );
            UnitsMultRep(    d1, &(spB->spRowUnits[ii-1]), &(spB->spColUnits[ii-1]) );
            UnitsDivRep( &(spEigen->spColUnits[ii-1]), dimen, d1, NO );
            spEigen->spColUnits[ii-1].radian_expnt = 0;
            spEigen->spColUnits[ii-1].units_type   = UnitsType;
       }

       /* Compute units for Eigenvectors */

       ZeroUnits( &(spEigen->spRowUnits[0]) );
       ZeroUnits( d1 );

       for( ii=1; ii<=iSize; ii++ ) {
           UnitsMultRep( dimen, &(spA->spRowUnits[ii-1]), &(spA->spColUnits[ii-1]) );
           UnitsMultRep(    d1, &(spB->spRowUnits[ii-1]), &(spB->spColUnits[ii-1]) );
           UnitsMultRep(     d, &(spEigen->spColUnits[0]), d1 );
           UnitsDivRep( &(spEigen->spRowUnits[ii]), d, dimen, NO );
           spEigen->spRowUnits[ii].units_type = UnitsType;
       }

       SetUnitsOn();

       free((char *)d->units_name);
       free((char *)d);
       free((char *)d1->units_name);
       free((char *)d1);
       free((char *)dimen->units_name);
       free((char *)dimen);
    }

#ifdef DEBUG
       printf("*** In Solve_Eigen() : End Units Calculation \n");
#endif 

    MatrixFree( spEigenvalue );
    MatrixFree( spEigenvector );

    return( spEigen );
}

/*
 *  ====================================================
 *  Extract_Eigenvalue() : Extract Eigenvalues from spA.
 *
 *  ====================================================
 */

#ifdef __STDC__
MATRIX *Extract_Eigenvalue( MATRIX *spA )
#else
MATRIX *Extract_Eigenvalue( spA )
MATRIX *spA;
#endif
{
MATRIX *spB;
int      ii;
int  length;

    if( spA == NULL )
        FatalError("In Extract_Eigenvalue() : run Eigen() before using this function",(char *)NULL );

    spB = MatrixAllocIndirect("EigenValue",DOUBLE_ARRAY,spA->iNoColumns,1);
    for( ii=1 ; ii<=spA->iNoColumns ; ii++ )
         spB->uMatrix.daa[ii-1][0] = spA->uMatrix.daa[0][spA->iNoColumns-ii];

    if(CheckUnits()==ON) {
       for( ii=1 ; ii<=spB->iNoRows ; ii++ ) {
            UnitsCopy(&(spB->spRowUnits[ii-1]), &(spA->spColUnits[spA->iNoColumns-ii]));
       }
       ZeroUnits( &(spB->spColUnits[0]) );
    }

    return( spB );
}

/*
 *  ======================================================
 *  Extract_Eigenvector() : Extract Eigenvectors from spA.
 *
 *  ======================================================
 */

#ifdef __STDC__
MATRIX *Extract_Eigenvector( MATRIX *spA )
#else
MATRIX *Extract_Eigenvector( spA )
MATRIX *spA;
#endif
{
MATRIX *spB;
int   ii,jj;
int  length;

    if( spA == NULL )
        FatalError("In Extract_Eigenvector() : run Eigen() before using this function",(char *)NULL );

    spB = MatrixAllocIndirect("EigenVector",DOUBLE_ARRAY,(spA->iNoRows-1),spA->iNoColumns);
    for( ii=1 ; ii<=spB->iNoRows ; ii++ )
    for( jj=1 ; jj<=spB->iNoColumns ; jj++ )
         spB->uMatrix.daa[ii-1][jj-1] = spA->uMatrix.daa[ii][spB->iNoColumns-jj];

    if(CheckUnits()==ON) {
       for( ii=1 ; ii<=spB->iNoRows ; ii++ )
	    UnitsCopy(&(spB->spRowUnits[ii-1]), &(spA->spRowUnits[ii]));

       for( ii=1 ; ii<=spB->iNoColumns ; ii++ )
            ZeroUnits( &(spB->spColUnits[ii-1]) );
    }

    return( spB );
}
