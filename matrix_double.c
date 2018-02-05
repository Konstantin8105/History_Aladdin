/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  matrix_double.c : Functions for Matrices of data type double
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
 *  Written by: Mark Austin                              July 1992 - October 1993
 *  ============================================================================= 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"

/*
#define DEBUG 
*/

/*
 ********************************************************
      Function: dMatrixPrint
      Function to print a matrix (nxm) of type double
      Usage: dMatrixPrint(array_name, array, n, m);
      where  array_name: name of the array.
             **array: pointer to a pointer to an array.
             n: number of rows.
             m: number of columns.
             NCOL: number of columns printed across page.
 ********************  **********************************
 */

#ifdef __STDC__
void dMatrixPrint(char *array_name, double **array, int n, int m)
#else
void dMatrixPrint(array_name, array, n, m)
char *array_name;
double **array;
int n,m;
#endif
{
int i,j,k=0;
int index, max_cr;
int mcr,l;

#define NCOL 4

  max_cr = (m%NCOL==0) ? (m/NCOL)-1 : m/NCOL;

  if (m<=NCOL) {
      mcr = m;
      max_cr = 0;
  }
  else mcr = NCOL;

  for (index=0;index<=max_cr;index++) {
     printf ("\nMatrix \"%s\"\n\n",array_name);
     printf ("row/col        ");

     for(i=(index*NCOL),l=0;i<=(index*mcr)+NCOL-1;i++,l++) {
        if (((index*mcr)+l)>=m) break;
        printf("%3d          ",i+1);
     }
     printf("\n");

     k = 0;
     for (i=1;i<=n;i++) {
        printf(" %3d ",++k);
        for (j=(index*NCOL)+1,l=0;j<=((index*NCOL)+NCOL);j++,l++) {
           if (((index*NCOL)+l)>=m) break;
           printf(" %12.5e",array[i-1][j-1]);
        }
        printf("\n");
     }
  }
}

/*
 *  =================
 *  Matrix Operations
 *  =================
 */

#ifdef __STDC__
double **dMatrixMult(double **m1, int m1_row, int m1_colum, double **m2 ,int m2_row, int m2_colum )
#else
double **dMatrixMult(m1, m1_row, m1_colum, m2 ,m2_row, m2_colum )
double **m1;
double **m2;
int m1_colum, m1_row;
int m2_colum, m2_row;
#endif
{
double **m3;
int i,j,k;
DIMENSIONS d1, d2;

#ifdef DEBUG
       printf("*** In dMatrixMult()\n");
#endif

    /* [a] : Check for Compatible Matrices */

       if(m1_colum != m2_row) {
           printf("FATAL ERROR >> Execution halted in dMatrixMult()\n");
           printf("FATAL ERROR >> Problem : m1_colum = %4d\n", m1_colum);
           printf("                       : m2->row  = %4d \n", m2_row);
           FatalError("In dMatrixMult() : Inconsistent Dimensions",(char *)NULL);
       }

    /* [b] : Multiply Matrices */
     
       m3  = MatrixAllocIndirectDouble(m1_row, m2_colum);

       for(i = 1; i <= m1_row; i++) {
          for(j = 1; j <= m2_colum; j++) {
             for(k = 1; k <= m1_colum; k++) {
                m3[i-1][j-1] += m1[i-1][k-1]*m2[k-1][j-1];
             }
          }
       }

#ifdef DEBUG
       printf("*** leaving dMatrixMult()\n");
#endif

   return (m3);
}

#ifdef __STDC__
double **dMatrixMultRep(double **m3, double **m1, int m1_row, int m1_colum, double **m2, int m2_row, int m2_colum )
#else
double **dMatrixMultRep(m3, m1, m1_row, m1_colum, m2, m2_row, m2_colum )
double **m3;
double **m1;
double **m2;
int m1_colum, m1_row;
int m2_colum, m2_row;
#endif
{
int         i,j,k;
DIMENSIONS d1, d2;

#ifdef DEBUG
       printf("*** In dMatrixMultRep()\n");
#endif

    /* [0] : Initialize m3                */
   
       for(i = 1; i <= m1_row; i++) {
          for(j = 1; j <= m2_colum; j++) {
                m3[i-1][j-1] = 0.0;
          }
       }

    /* [a] : Check for Compatible Matrices */

       if(m1_colum != m2_row) {
           printf("FATAL ERROR >> Execution halted in dMatrixMultRep()\n");
           printf("FATAL ERROR >> Problem : m1_colum = %4d\n", m1_colum);
           printf("                       : m2->row  = %4d \n", m2_row);
           FatalError("In dMatrixMult() : Inconsistent Dimensions",(char *)NULL);
       }

    /* [b] : Multiply Matrices */
     
       for(i = 1; i <= m1_row; i++) {
          for(j = 1; j <= m2_colum; j++) {
             for(k = 1; k <= m1_colum; k++) {
                m3[i-1][j-1] += m1[i-1][k-1]*m2[k-1][j-1];
             }
          }
       }

#ifdef DEBUG
       printf("*** leaving dMatrixMultRep()\n");
#endif

   return (m3);
}

#ifdef __STDC__
double **dMatrixTranspose(double **m1, int iNoRows, int iNoColumns)
#else
double **dMatrixTranspose(m1, iNoRows, iNoColumns)
double **m1;
int iNoRows, iNoColumns;
#endif
{
double **m2;
int i,j,k;

#ifdef DEBUG
       printf("*** In dMatrixTranspose()\n");
#endif

    /* [a] : Transpose Matrix */

       m2 = MatrixAllocIndirectDouble(iNoColumns, iNoRows);

       for(i = 1; i <= iNoRows; i++){
         for(j = 1; j <= iNoColumns; j++) {
           m2[j-1][i-1] = m1[i-1][j-1];
         }
       }
#ifdef DEBUG
       printf("*** Leave dMatrixTranspose()\n");
#endif
       return (m2);
}


#ifdef __STDC__
double **dMatrixCopy(double **m1, int iNoRows, int iNoColumns)
#else
double **dMatrixCopy(m1, iNoRows, iNoColumns)
double **m1;
int iNoRows, iNoColumns;
#endif
{
double **m2;
int i,j,k;

#ifdef DEBUG
       printf("*** In dMatrixCopy()\n");
#endif
       m2  = MatrixAllocIndirectDouble(iNoRows, iNoColumns);

       for(i = 1; i <= iNoRows; i++) {
           for(j = 1; j <= iNoColumns; j++) {
              m2[i-1][j-1] = m1[i-1][j-1];
           }
       }

#ifdef DEBUG
       printf("*** Leave dMatrixCopy()\n");
#endif
       return (m2);
}

#ifdef __STDC__
double **dMatrixCopyRep(double **m2, double **m1, int iNoRows, int iNoColumns)
#else
double **dMatrixCopyRep(m2, m1, iNoRows, iNoColumns)
double **m2;
double **m1;
int iNoRows, iNoColumns;
#endif
{
int i,j,k;

#ifdef DEBUG
       printf("*** In dMatrixCopyRep()\n");
#endif
       
       for(i = 1; i <= iNoRows; i++) {
           for(j = 1; j <= iNoColumns; j++) {
              m2[i-1][j-1] = m1[i-1][j-1];
           }
       }

#ifdef DEBUG
       printf("*** Leave dMatrixCopyRep()\n");
#endif
       return (m2);
}



/* ================================================================ */
/* Vector cross product a = b x c;                                  */
/* output :             a = a[n1][1];                               */
/* input  :             b = b[n2][1] or b[1][n2]                    */
/*        :             c = c[n3][1] or c[1][n3]                    */
/*        :             n1,n2,n3 = 1, 2, 3                          */
/* ================================================================ */

#ifdef __STDC__
double **dVmatrixCrossProduct(double **a, double **b, int b_rows, int b_cols, double **c, int c_rows, int c_cols)
#else
double **dVmatrixCrossProduct(a, b, b_rows, b_cols, c, c_rows, c_cols)
int b_rows, b_cols;
int c_rows, c_cols;
double **a, **b, **c;
#endif
{
double b1, b2,b3;
double c1, c2,c3;

#ifdef DEBUG
     printf(" enter dVmatrixCrossProduct(): \n");
#endif

    /* Check that matrix is either a 3x1 (or 1x3)  column (or row) vector */

      if((b_rows != 1) && (b_cols != 1) || 
         (b_rows != 3) && (b_cols != 3) ) {
          FatalError("In dVmatrixCrossProduct() : Matrix b is not a 1x3 row (or 3x1 column) vector",(char *)NULL);
      }

      if((c_rows != 1) && (c_cols != 1) ||
         (c_rows != 3) && (c_cols != 3)) {
          FatalError("In dVmatrixCrossProduct() : Matrix c is not a 1x3 row (or 3x1 column) vector",(char *)NULL);
      }

    /* calculate the vector cross product */

     if(b_rows == 1) {
        b1 = b[0][0];
        b2 = b[0][1];
        b3 = b[0][2];
     }
     else if(b_cols == 1) {
        b1 = b[0][0];
        b2 = b[1][0];
        b3 = b[2][0];
     }

     if(c_rows == 1) {
        c1 = c[0][0];
        c2 = c[0][1];
        c3 = c[0][2];
     }
     else if(c_cols == 1) {
        c1 = c[0][0];
        c2 = c[1][0];
        c3 = c[2][0];
     }

     a[0][0] = b2*c3 - c2*b3;
     a[1][0] = b3*c1 - c3*b1;
     a[2][0] = b1*c2 - c1*b2;

#ifdef DEBUG
     printf(" leaving dVmatrixCrossProduct(): \n");
#endif

     return(a);
}

/* ================================================================ */
/* Vector inner (dot) product a = b . c;                            */
/* output :             a = a number ;                              */
/* input  :             b = b[n2][1] or b[1][n2]                    */
/*        :             c = c[n3][1] or c[1][n3]                    */
/* ================================================================ */

#ifdef __STDC__
double dVmatrixInnerProduct(double **b, int b_rows, int b_cols, double **c, int c_rows, int c_cols)
#else
double dVmatrixInnerProduct(b, b_rows, b_cols, c, c_rows, c_cols)
int b_rows, b_cols;
int c_rows, c_cols;
double    **b, **c;
#endif
{
double sum = 0.0;
double b1, b2,b3;
double c1, c2,c3;
int      i, j, k;

#ifdef DEBUG
     printf(" enter dVmatrixInnerProduct(): \n");
#endif

    /* Check that matrix is either a column (or row) vector */

      if((b_rows != 1) && (b_cols != 1) ) { 
          FatalError("In dVmatrixInnerProduct() : Matrix b is not a row (or column) vector",(char *)NULL);
      }

      if((c_rows != 1) && (c_cols != 1) ) {
          FatalError("In dVmatrixInnerProduct() : Matrix c is not a row (or column) vector",(char *)NULL);
      }

    /* Check vectors  b and c have same dimensions */

      if((c_rows == 1) && 
         ((c_cols != b_rows) &&(c_cols != b_cols))   ) {
          FatalError("In dVmatrixInnerProduct() : Vectors b and c do not have same dimensions",(char *)NULL);
      }

      if((c_cols == 1) && 
         ((c_rows != b_rows) &&(c_rows != b_cols))   ) {
          FatalError("In dVmatrixInnerProduct() : Vectors b and c do not have same dimensions",(char *)NULL);
      }

    /* calculate the vector inner (dot) product */

     if(b_rows == 1) {
        if(c_rows == 1) {
           for (i = 1; i <= b_cols; i++)   
              sum += b[0][i-1]*c[0][i-1];
        }
        else if(c_cols == 1) {
           for (i = 1; i <= b_cols; i++)   
              sum += b[0][i-1]*c[i-1][0];
     
        }
     }
     else if(b_cols == 1) {
        if(c_rows == 1) {
           for (i = 1; i <= b_rows; i++)   
              sum += b[i-1][0]*c[0][i-1];
        }
        else if(c_cols == 1) {
           for (i = 1; i <= b_rows; i++)  { 
              sum += b[i-1][0]*c[i-1][0];
           }
        }
     }

#ifdef DEBUG
     printf(" leaving dVmatrixInnerProduct(): \n");
#endif

     return(sum);
}

#ifdef __STDC__
double dVmatrixL2Norm(double **m, int iNoRows, int iNoCols)
#else
double dVmatrixL2Norm(m, iNoRows, iNoCols)
double **m;
int iNoRows, iNoCols;
#endif
{
double sum = 0.0;
int    i;
#ifdef DEBUG
       printf("*** Enter dVmatrixL2Norm() \n");
#endif

   /* Compute L2 Norm for column/row vector */

      if(iNoRows == 1) {
         for(i = 1; i <= iNoCols; i++)
             sum = sum + (m[0][i-1]*m[0][i-1]);
      }
      else if (iNoCols == 1) {
         for(i = 1; i <= iNoRows; i++)
             sum = sum + (m[i-1][0]*m[i-1][0]);
      }

      sum = sqrt(sum);

#ifdef DEBUG
       printf("*** Leaving dVmatrixL2Norm() \n");
#endif

     return(sum);
}

#ifdef __STDC__
double dMatrixDet(double **m, int iNoRows, int iNoCols)
#else
double dMatrixDet(m, iNoRows, iNoCols)
double **m;
int iNoRows, iNoCols;
#endif
{
double   **temp_m;
double      value;
int           i,j;
int   ii,ij,im,in;

   if( iNoRows!=iNoCols )
       FatalError("In dMatrixDet() : matrix is not a square matrix",(char *)NULL);

   if( iNoRows==1 ) {
       value = m[0][0];
       return( value );
   }
   if( iNoRows==2 ) {
       value = m[0][0]*m[1][1]-m[0][1]*m[1][0];
       return( value );
   }
   if( iNoRows==3 ) {
       value = m[0][0]*m[1][1]*m[2][2]+m[0][1]*m[1][2]*m[2][0]+m[0][2]*m[2][1]*m[1][0]
              -m[0][2]*m[1][1]*m[2][0]-m[0][0]*m[2][1]*m[1][2]-m[0][1]*m[1][0]*m[2][2];
       return( value );
   }

   value = 0.0;
   for( i=0 ; i<iNoRows ; i++ ) {
       temp_m = MatrixAllocIndirectDouble( iNoRows-1, iNoCols-1 );
       for( ii=0, im=1 ; ii<(iNoRows-1) ; ii++, im++ ) {
           for( ij=0, in=0 ; ij<(iNoRows-1) ; ij++, in++ ) {
               if( in==i )  in++;
               temp_m[ii][ij] = m[im][in];
           }
       }

       value += pow(-1.0,i)*m[0][i]*dMatrixDet( temp_m, iNoRows-1, iNoRows-1 );
       MatrixFreeIndirectDouble( temp_m, iNoRows-1 );
   }

   return(value);
}


