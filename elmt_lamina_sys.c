/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_lamina_sys.c : SHELL_4N Element (Linear-Elastic/Elastic-Plastic)
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
 *  Written by: Xiaoguang Chen                                       January 1993
 *  ============================================================================= 
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "fe_functions.h"
#include "elmt.h"

/*
#define DEBUG 
*/


/* ============================================================== */
/*   Element SHELL_BELYTSCHKO                                     */
/*   3D   Shell Element:                                          */
/*        material properties array                               */
/*        Input Properties:                                       */
/* ============================================================== */
/*    p->work_material[0] = E;
      p->work_material[1] = G;
      p->work_material[2] = fy;
      p->work_material[3] = ET;
      p->work_material[4] = nu;
      p->work_material[5] = density;

      p->work_section[0] = Ixx;
      p->work_section[1] = Iyy;
      p->work_section[2] = Izz;
      p->work_section[3] = Ixy;
      p->work_section[4] = Ixz;
      p->work_section[5] = Iyz;
      p->work_section[6] = weight;
      p->work_section[7] = bf;
      p->work_section[8] = tf;
      p->work_section[9] = depth;                                  
      p->work_section[10] = area;                                 
      p->work_section[11] = thickness;                            */
/* ============================================================== */

/* =======================================================*/
/* Calculation of Direction Matrix; Transform coordinate  */
/* & velocity into lamina coordiante system               */
/* =======================================================*/

#ifdef __STDC__
void Lamina_Sys( ARRAY *p, double **Direction_Matrix, double **co_coord, double **co_velocity)
#else
void Lamina_Sys(p, Direction_Matrix, co_coord, co_velocity)
ARRAY  *p;
double **Direction_Matrix, **co_coord, **co_velocity;
#endif
{
double cs, sn;
double dof, size, temp, temp11, temp12, temp21, temp22;

double x21, x31, x42;
double y21, y31, y42;
double z21, z31, z42;

double co_x31, co_y31, co_z31;
double co_x42, co_y42, co_z42;

double co_vel_x31, co_vel_y31, co_vel_z31;
double co_vel_x42, co_vel_y42, co_vel_z42;

double co_rot_x31, co_rot_y31, co_rot_z31;
double co_rot_x42, co_rot_y42, co_rot_z42;

double Phi_dot, delta_phi;    /* the rate of the angle phi, Phi is the angle between x-axis and x^ axis */
double Omega_z;    /* spin rate about z-axis         */
double Omega_z21;  /* angular velocity of side 1-2   */
double Jac;          /* Determinant of Jacobian matrix */
double s_norm;

double **r21_ptr, **r31_ptr, **r42_ptr;
double **s_ptr, **e1_ptr, **e2_ptr, **e3_ptr;
double **A_matrix, **coord;

double **displ;
double delta_t = 0.0;

int i, j, k, ii, jj, kk, n, n1, n2, k1, k2;
int counter;
    
#ifdef DEBUG
    printf(" Enter Lamina_Sys() \n");
#endif

    dof  = p->dof_per_node;
    size = p->size_of_stiff;

    e1_ptr  = MatrixAllocIndirectDouble(3, 1);
    e2_ptr  = MatrixAllocIndirectDouble(3, 1);
    e3_ptr  = MatrixAllocIndirectDouble(3, 1);


 /* --------------------------------------------------- */
 /* [1] Orientatuions of local base vectors             */ 
 /* --------------------------------------------------- */

 /* [a] Compute e3_ptr                      */
 /*     pointer normal to the shell surface */
     
#ifdef DEBUG
    printf(" In elmt_shell_Belytschko_explicit(): \n");
    printf("    : computing e3_ptr \n");
#endif

    x21 =  p->coord[0][1].value - p->coord[0][0].value;
    y21 =  p->coord[1][1].value - p->coord[1][0].value;
    z21 =  p->coord[2][1].value - p->coord[2][0].value;

    x31 =  p->coord[0][2].value - p->coord[0][0].value;
    y31 =  p->coord[1][2].value - p->coord[1][0].value;
    z31 =  p->coord[2][2].value - p->coord[2][0].value;

    x42 =  p->coord[0][3].value - p->coord[0][1].value;
    y42 =  p->coord[1][3].value - p->coord[1][1].value;
    z42 =  p->coord[2][3].value - p->coord[2][1].value;
  
    r42_ptr = MatrixAllocIndirectDouble(3, 1);
    r42_ptr[0][0] = x42;
    r42_ptr[1][0] = y42;
    r42_ptr[2][0] = z42;

    r31_ptr = MatrixAllocIndirectDouble(3, 1);
    r31_ptr[0][0] = x31;
    r31_ptr[1][0] = y31;
    r31_ptr[2][0] = z31;

    r21_ptr = MatrixAllocIndirectDouble(3, 1);
    r21_ptr[0][0] = x21;
    r21_ptr[1][0] = y21;
    r21_ptr[2][0] = z21;

#ifdef DEBUG
    dMatrixPrint("r31_ptr", r31_ptr, 3,1);
    dMatrixPrint("r42_ptr", r42_ptr, 3,1);
#endif

    s_ptr   = MatrixAllocIndirectDouble(3, 1);
    s_ptr   = dVmatrixCrossProduct(s_ptr, r31_ptr, 3, 1, r42_ptr, 3,1);

#ifdef DEBUG
    dMatrixPrint("s_ptr", s_ptr, 3,1);
#endif

    s_norm  = (double) dVmatrixL2Norm(s_ptr, 3, 1);

    e3_ptr[0][0] = s_ptr[0][0]/s_norm;
    e3_ptr[1][0] = s_ptr[1][0]/s_norm;
    e3_ptr[2][0] = s_ptr[2][0]/s_norm;

#ifdef DEBUG
    dMatrixPrint("e3_ptr", e3_ptr, 3,1);
#endif

 /* [b] Procedure 1: Compute e1_ptr  for x^ axis is embedded in the */
 /*                  element between nodes 1 and 2, side 1-2.       */
 /*     Note :       procedure is quite accurate if the shear       */
 /*                  strains are less than 10%                      */
 /*     e1_ptr, e2_ptr: pointers tangent to shell surface           */

#ifdef DEBUG
    printf(" In elmt_shell_Belytschko_explicit(): \n");
    printf("                   : computing e1_ptr \n");
#endif

    temp = dVmatrixInnerProduct(r21_ptr, 3, 1, e3_ptr, 3, 1);

    s_ptr[0][0]  = x21 - temp*e3_ptr[0][0]; 
    s_ptr[1][0]  = y21 - temp*e3_ptr[1][0]; 
    s_ptr[2][0]  = z21 - temp*e3_ptr[2][0]; 

    s_norm       = dVmatrixL2Norm(s_ptr, 3, 1);

    e1_ptr[0][0] = s_ptr[0][0]/s_norm;
    e1_ptr[1][0] = s_ptr[1][0]/s_norm;
    e1_ptr[2][0] = s_ptr[2][0]/s_norm;
    
#ifdef DEBUG
       dMatrixPrint("e1_ptr", e1_ptr, 3,1);
#endif

 /* [c] Procedure 1: Compute e2_ptr  for x^ axis is embedded in the element */
 /*                  between nodes 1 and 2, side 1-2.                       */
 /*                  Lamina coordinate                                      */

#ifdef DEBUG
    printf(" In elmt_shell_Belytschko_explicit(): \n");
    printf("    : computing e2_ptr \n");
#endif

    e2_ptr = dVmatrixCrossProduct(e2_ptr, e3_ptr, 3, 1, e1_ptr, 3, 1);
    
#ifdef DEBUG
    dMatrixPrint("e2_ptr", e2_ptr, 3,1);
    printf(" In elmt_shell_Belytschko_explicit(): \n");
    printf("    : computing direction matrix  \n");
#endif

    MatrixFreeIndirectDouble(r21_ptr, 3);
    MatrixFreeIndirectDouble(r31_ptr, 3);
    MatrixFreeIndirectDouble(r42_ptr, 3);

    for ( i = 1; i <= 6; i++) {
        if(i <= 3) {
          Direction_Matrix[0][i-1] = e1_ptr[i-1][0];
          Direction_Matrix[1][i-1] = e2_ptr[i-1][0];
          Direction_Matrix[2][i-1] = e3_ptr[i-1][0];
          Direction_Matrix[3][i-1] = 0.0;
          Direction_Matrix[4][i-1] = 0.0;
          Direction_Matrix[5][i-1] = 0.0;
        }
        else {
          j = i - 3;
          Direction_Matrix[0][i-1] = 0.0;
          Direction_Matrix[1][i-1] = 0.0;
          Direction_Matrix[2][i-1] = 0.0;
          Direction_Matrix[3][i-1] = e1_ptr[j-1][0];
          Direction_Matrix[4][i-1] = e2_ptr[j-1][0];
          Direction_Matrix[5][i-1] = e3_ptr[j-1][0];
        }
    }

    MatrixFreeIndirectDouble(e1_ptr, 3);
    MatrixFreeIndirectDouble(e2_ptr, 3);
    MatrixFreeIndirectDouble(e3_ptr, 3);
    MatrixFreeIndirectDouble(s_ptr, 3);

    A_matrix = MatrixAllocIndirectDouble(3, 3);
    for ( i = 1; i <= 3; i++) 
       for (j = 1; j <= 3; j++)
            A_matrix[i-1][j-1] = Direction_Matrix[i-1][j-1]; 

    counter = 1;

    coord = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);


START:

 /* Transform velocity and coordinates x,y z */
 /* from global to local co-rational system  */ 

    for( i = 1; i <= 3; i++ ) {
       for ( j = 1; j <= p->nodes_per_elmt; j++) {
          coord[i-1][j-1] = co_coord[i-1][j-1]; 
       }
    }

#ifdef DEBUG
    dMatrixPrint("coord", coord, 3,4);
    dMatrixPrint("A_matrix", A_matrix, 3,3);
#endif

    co_coord = dMatrixMultRep(co_coord, A_matrix, 3, 3, coord, 3, p->nodes_per_elmt); 

#ifdef DEBUG
       dMatrixPrint("co_coord", co_coord, 3,4);
#endif
    
    for(k = 1; k <= p->nodes_per_elmt; k++) { 
    for(i = 1; i <= p->dof_per_node; i++) {
        co_velocity[i-1][k-1] = 0.0;
    }
    } 

#ifdef DEBUG
       dMatrixPrint(" co_velocity",  co_velocity, p->dof_per_node, p->nodes_per_elmt);
#endif

 /* [d] Procedure 2: Compute e1_ptr, e2_ptr for x^ axis is not embedded in */
 /*                  the element along the side 1-2.                       */

    /* angular velocity of side 1-2 */

       temp      = sqrt(x21*x21 + y21*y21 + z21*z21);

       Omega_z21 = (co_velocity[1][1] - co_velocity[1][0])/temp;
       
       co_x31 =  co_coord[0][2] - co_coord[0][0];
       co_y31 =  co_coord[1][2] - co_coord[1][0];
       co_z31 =  co_coord[2][2] - co_coord[2][0];

       co_x42 =  co_coord[0][3] - co_coord[0][1];
       co_y42 =  co_coord[1][3] - co_coord[1][1];
       co_z42 =  co_coord[2][3] - co_coord[2][1];
    
       co_vel_x31 =  co_velocity[0][2] - co_velocity[0][0];
       co_vel_y31 =  co_velocity[1][2] - co_velocity[1][0];
       co_vel_z31 =  co_velocity[2][2] - co_velocity[2][0];

       co_vel_x42 =  co_velocity[0][3] - co_velocity[0][1];
       co_vel_y42 =  co_velocity[1][3] - co_velocity[1][1];
       co_vel_z42 =  co_velocity[2][3] - co_velocity[2][1];

       Jac = (x31*y42-x42*y31);

       Omega_z   = (co_y42*co_vel_y31+co_y31*co_vel_y42 - co_x42*co_vel_x31 - co_x31*co_vel_x42)/Jac;
       Phi_dot   = Omega_z - Omega_z21;

       delta_phi = Phi_dot*delta_t;

       cs        = p->direc_cos[0];
       sn        = p->direc_cos[1];

       p->direc_cos[0] = cs - delta_phi * sn - 0.5*delta_phi*delta_phi*cs; 
       p->direc_cos[1] = sn - delta_phi * cs - 0.5*delta_phi*delta_phi*sn;
       
      /* modify the direction cosine vectors e1, e2, e3 of loacal coord */

       cs = p->direc_cos[0];
       sn = p->direc_cos[1];

       temp11 =  cs*A_matrix[0][0] + sn*A_matrix[1][0];
       temp12 =  cs*A_matrix[0][1] + sn*A_matrix[1][1];

       temp21 = -sn*A_matrix[0][0] + cs*A_matrix[1][0];
       temp22 = -sn*A_matrix[0][1] + cs*A_matrix[1][1];
         
       A_matrix[0][0] = temp11;
       A_matrix[0][1] = temp12;
       A_matrix[0][2] = 0.0;

       A_matrix[1][0] = temp21;
       A_matrix[2][1] = temp22;
       A_matrix[2][2] = 0.0;

       A_matrix[2][0] =  0.0;
       A_matrix[2][1] =  0.0;
       A_matrix[2][2] =  1.0;
       
        /* need a two-pass procedure to obtain the */
        /* correct direction cosine and local      */
        /* coordindate and loocal velocity         */

        if(counter == 1) {
           counter++; 
           goto START;
        }

END:
        MatrixFreeIndirectDouble(coord, 3);
        MatrixFreeIndirectDouble(A_matrix, 3);

#ifdef DEBUG
        dMatrixPrint("co_coord", co_coord, 3,4);
        dMatrixPrint("co_velocity", co_velocity, 3,4);
        printf(" Leaving Lamina_Sys()\n");
#endif

}
