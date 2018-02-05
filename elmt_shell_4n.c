/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_shell_4n.c : STATIC/DYNAMIC SHELL_4N Element
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
 *  ------------------------------------------------------------------- 
 *
 *  The Algorithms are rederived in order to suit for Implicit Algorithms.
 *  For implicit integration: linear strain and displacement relation is 
 *  employed
 *                                                                    
 *  ------------------------------------------------------------------- 
 *                                                                    
 *  Written by: Xiaoguang Chen                                       January 1995
 *  ============================================================================= 
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"
#include "fe_database.h"
#include "fe_functions.h"
#include "elmt.h"

/* function declarations */

ARRAY *sld04( ARRAY *, int );


/* ============================================================== */
/*   Element SHELL_FOUR_NODES                                     */
/*        Implicit code                                           */
/*        Shell Element:                                          */
/*        material properties array                               */
/*        Input Properties:                                       */
/* ============================================================== */
/*    p->work_material[0] = E;
      p->work_material[1] = G;
      p->work_material[2] = fy;
      p->work_material[3] = ET;
      p->work_material[4] = nu;
      p->work_material[5] = density;
      p->work_material[6] = fu;

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

ARRAY *elmt_shell_4n(p, isw)
ARRAY    *p;
int     isw;
{

#ifdef DEBUG
       printf(" enter elmt_shell_4n() \n");
#endif

   p = elmt_shell_4nodes_implicit(p, isw);

#ifdef DEBUG
       printf(" Leaving elmt_shell_4n() \n");
#endif

   return(p);
}

ARRAY *elmt_shell_4nodes_implicit(p, isw)
ARRAY    *p;
int     isw;
{
static double                                   nu; 
static QUANTITY      fy, G, E, ET, density, lambda; 
static double         Ixx, Iyy, Izz, Ixy, Ixz, Iyz;
static double        bf, tf, A, depth, h, EA, EIzz;

static DIMENSIONS *dp_length, *dp_force, *dp_moment;
static DIMENSIONS *dp_stress, *dp_degree, *dp_temperature;

double                               x21, x31, x42;
double                               y21, y31, y42;
double                               z21, z31, z42;

double                      co_x31, co_y31, co_z31;
double                      co_x42, co_y42, co_z42;

double             EE, Jac,jacobian, k_bar = 0.833;     /* k_bar is shear factor */
double                   elmt_length, aspect_ratio;

double                                  **co_coord;     /*  co_coord is a coordinate rotated with material deformation */
double                   **T_matrix, **T_Transpose;
double                          **Direction_Matrix;

static double                *integ_coord, *weight;

double                           **stress, **displ;
double             **nodal_load, **nodal_load_temp;
double                             diff, sum, temp;
double                      **Cep, **stiff, **mass; 
int         i, j, k, ii, jj, kk, n, n1, n2, k1, k2;
int      surface_pts, surf_integ_pts, no_integ_pts;
int            size, dof, length, length1, length2;
int                        UNITS_SWITCH, UnitsType;


#ifdef DEBUG
       printf("*** Enter elmt_shell_4nodes_implicit() : isw = %4d\n", isw);
#endif


/* =======================================================*/
/* Calculation of Direction Matrix; Transform coordinate  */
/* into lamina coordiante system                          */
/* =======================================================*/
     
       UNITS_SWITCH = CheckUnits();
       UnitsType    = CheckUnitsType();

       co_coord         = MatrixAllocIndirectDouble(p->no_dimen, p->nodes_per_elmt);
       Direction_Matrix = MatrixAllocIndirectDouble(6,6);

       for( i = 1; i <= 3; i++ ) {
           for ( j = 1; j <= p->nodes_per_elmt; j++) {
               co_coord[i-1][j-1] = p->coord[i-1][j-1].value; /* set initial co_coord */
           }
       }

       Lamina_Sys_Implicit(p, Direction_Matrix, co_coord);

#ifdef DEBUG
        dMatrixPrint("co_coord", co_coord, 3,4);
        dMatrixPrint("Direction_Matrix", Direction_Matrix, 6,6);
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

      co_x31 =  co_coord[0][2] - co_coord[0][0];
      co_y31 =  co_coord[1][2] - co_coord[1][0];
      co_z31 =  co_coord[2][2] - co_coord[2][0];

      co_x42 =  co_coord[0][3] - co_coord[0][1];
      co_y42 =  co_coord[1][3] - co_coord[1][1];
      co_z42 =  co_coord[2][3] - co_coord[2][1];

      Jac = (x31*y42-x42*y31);
      no_integ_pts = p->integ_ptr->thickness_pts;
      dof          = p->dof_per_node;
      size         = p->size_of_stiff;

/*****************************************************************/
#ifdef DEBUG
    printf(" ***** In elmt_shell_4nodes_implicit(): \n");
    printf("                                      : enter main switch \n");
#endif

    switch(isw) {
        case PROPTY: 

#ifdef DEBUG
    printf(" In elmt_shell_4nodes_implicit(): \n");
    printf("    : enter case of PROPTY\n");
#endif

             /* material properties:  elastic */

             E.value       =  p->work_material[0].value;
             fy.value      =  p->work_material[2].value;
             ET.value      =  p->work_material[3].value;
             nu            =  p->work_material[4].value;
             density.value =  p->work_material[5].value;
             lambda.value  =  nu*E.value/(1.0+nu)/(1.0-2.0*nu);

           /* (1)   check  possion_ratio value */

            if( nu == 0.0 || nu > 0.5 ) {
                printf("WARNING >> ... In shell_Belytschko element() -  nu value = %9.4f,reset to 0.3 !\n",
                       nu);
                nu = 0.3;    /* default poi_ratio value */
            }

            /* (2)   calculate  G value */
            
            G.value = p->work_material[1].value = E.value/(1.0 - 2.0*nu) ;

            if(E.value/((1.0 - 2.0*nu)) != p->work_material[1].value) {
                    printf(" elmt_shell(): WARNING: G is not equal to E/(1-2nu), check G for homogeneous material \n");
                    printf("             : ignore this message for non-homogeneous materials \n");
            }


            if( UNITS_SWITCH == ON ) {
               E.dimen       =  p->work_material[0].dimen;
               fy.dimen      =  p->work_material[2].dimen;
               ET.dimen      =  p->work_material[3].dimen;
               density.dimen =  p->work_material[5].dimen;
               lambda.dimen  =  E.dimen;
               G.dimen = E.dimen;
            }

             Ixx    = p->work_section[0].value;
             Iyy    = p->work_section[1].value;
             Izz    = p->work_section[2].value;
             Ixy    = p->work_section[3].value;
             Ixz    = p->work_section[4].value;
             Iyz    = p->work_section[5].value;
             bf     = p->work_section[7].value;
             tf     = p->work_section[8].value;
             depth  = p->work_section[9].value;
             A      = p->work_section[10].value;
             h      = p->work_section[11].value;    /* thickness of the shell */

             EA     = E.value*A;
             EIzz   = E.value*Izz;
         
#ifdef DEBUG
    printf(" In elmt_shell_4nodes_implicit(): \n");
    printf("    : leaving case of PROPTY\n");
#endif
             break;
             case STIFF: /* form element stiffness */

#ifdef DEBUG
       printf("*** In elmt_shell() : start case STIFF\n");
       printf("                    : Density         = %14.4e\n", density.value);
       printf("                    : shell_thickness = %8.2f\n", h);
       printf("                    : Jac             = %8.2f\n", Jac);
       printf("                    : E               = %lf\n", E.value);
       printf("                    : nu              = %lf\n", nu);
       printf("                    : no_integ_pts in z_direction = %d\n", no_integ_pts);
#endif

                 /* Integration loop over thickness */
      
                 stiff        = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
                 T_matrix     = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
                 T_Transpose  = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
                 integ_coord  = dVectorAlloc(no_integ_pts + 1);
                 weight       = dVectorAlloc(no_integ_pts + 1);

                 size         = p->size_of_stiff;
                 gauss(integ_coord,weight,no_integ_pts);

                 /* Integration loop over surface */

                 for (ii = 1; ii <= no_integ_pts; ii++) {

                      Shell_Stiff_Plane_4node(stiff, p, co_coord, integ_coord[ii], ii, E.value, nu); 	

                    for ( i = 1; i <= size; i++) {
                        for ( j = 1; j <= size; j++) {
                            p->stiff->uMatrix.daa[i-1][j-1] += stiff[i-1][j-1]*h*0.5*weight[ii];     
                        }
                    }
                  
                 }  /* end of gaussian integration through thickness */

#ifdef DEBUG
     printf(" End of integration through the thickness \n");
#endif

#ifdef DEBUG
           dMatrixPrint(" stiff before rotation", p->stiff->uMatrix.daa, p->size_of_stiff, p->size_of_stiff); 
#endif
                 free((char *) integ_coord);
                 free((char *) weight);

                 /***************************************************/
                 /* Transform stiffness matrix from local lamina    */
                 /* coordinate system to global coordinate system   */
                 /***************************************************/

                /* store the direction matrix : T_matrix */

                for (i = 1; i <= p->size_of_stiff; i++)  {
                    for(j = 1; j <= p->nodes_per_elmt; j++) {
                        for( k = 1; k <= p->dof_per_node; k++) {
                             n1  = p->dof_per_node*(j-1);
                             n2  = p->dof_per_node*j;
                             n   = n1 + k;
                             if(i > n1 && i <= n2) {
                                ii = i - n1;
                                T_matrix[i-1][n-1] = Direction_Matrix[ii-1][k-1];
                             }
                             else 
                                T_matrix[i-1][n-1] = 0.0;
                        }
                    }
                }

                stiff =  dMatrixMultRep(stiff, p->stiff->uMatrix.daa, size, size, T_matrix, size, size);

                for(i = 1; i <= size; i++){
                    for(j = 1; j <= size; j++) {
                        T_Transpose[i-1][j-1] = T_matrix[j-1][i-1];
                    }
                }
                p->stiff->uMatrix.daa = dMatrixMultRep(p->stiff->uMatrix.daa,
                                 T_Transpose, size, size, stiff, size, size);
    
#ifdef DEBUG
               /* check the symmetry of the stiffness matrix */

                for(i = 1; i <= size; i++){
                    for(j = 1; j <= size; j++) {
                        diff = p->stiff->uMatrix.daa[i-1][j-1] -p->stiff->uMatrix.daa[j-1][i-1];

                        if(diff > 1E-7 && 
                           (ABS(p->stiff->uMatrix.daa[i-1][j-1])) > 1E-7 &&
                           (ABS(diff/p->stiff->uMatrix.daa[i-1][j-1])) > 1E-1) {

                           printf("K[%d][%d] = %lf \n",
                                   i, j, p->stiff->uMatrix.daa[i-1][j-1]);
                           printf("K[%d][%d] = %lf \n",
                                   j, i, p->stiff->uMatrix.daa[j-1][i-1]);
                           printf("diff = %le \n", 
                                  p->stiff->uMatrix.daa[i-1][j-1] -
                                  p->stiff->uMatrix.daa[j-1][i-1]);
                           printf("elmtNo = %d Stiffness matrix IS NOT SYMMETRIC \n",
                                  p->elmt_no);
                           break;
                       }
                       else {
                           if( i == size && j == size)
                              printf("elmtNo = %d Stiffness matrix IS SYMMETRIC \n",
                                      p->elmt_no);
                       }
                    }
                }
#endif
                MatrixFreeIndirectDouble(T_matrix, size);
                MatrixFreeIndirectDouble(T_Transpose, size);
                MatrixFreeIndirectDouble(stiff, size);

                /**************************************************/
                /* Assign Units to Stiffness Matrix               */
                /**************************************************/

       /* Initiation of Stiffness Units Buffer                      */

       switch( UNITS_SWITCH ) {
         case ON:
           if(UnitsType == SI || UnitsType == SI_US ) {
              dp_stress = DefaultUnits("Pa");
              dp_length = DefaultUnits("m");
           }
           else {
              dp_stress = DefaultUnits("psi");
              dp_length = DefaultUnits("in");
           }

          /* node 1 */
           UnitsMultRep( &(p->stiff->spColUnits[0]), dp_stress, dp_length );
           UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );
           UnitsCopy( &(p->stiff->spColUnits[2]), &(p->stiff->spColUnits[0]) );
           UnitsMultRep( &(p->stiff->spColUnits[3]), &(p->stiff->spColUnits[0]) , dp_length );
           UnitsCopy( &(p->stiff->spColUnits[4]), &(p->stiff->spColUnits[3]) );

           ZeroUnits( &(p->stiff->spRowUnits[0]) );
           ZeroUnits( &(p->stiff->spRowUnits[1]) );
           ZeroUnits( &(p->stiff->spRowUnits[2]) );
           UnitsCopy( &(p->stiff->spRowUnits[3]), dp_length );
           UnitsCopy( &(p->stiff->spRowUnits[4]), dp_length );

          /* node i  i > 1*/
           for ( i = 2; i <= p->nodes_per_elmt; i++) {
                kk = p->dof_per_node*(i-1) + 3; 
                for( j = 1; j <= p->dof_per_node; j++) {
                     k  = p->dof_per_node*(i-1) + j;
                     if( k <= kk) {
                        UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                        UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[0]) );
                     }
                     if(k > kk) {
                        UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[3]) );
                        UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[3]) );
                     }
                 }
           }
           free((char *) dp_stress->units_name);
           free((char *) dp_stress);
           free((char *) dp_length->units_name);
           free((char *) dp_length);

          break;
          case OFF:
          break;
          default:
          break;
        }
             break;
	case PRESSLD:
	case EQUIV_NODAL_LOAD:

#ifdef DEBUG
       printf("*** In elmt_shell() : enter case EQUIV_NODAL_LOAD\n");
#endif
         /* Form external nodal load vector */
         /* due to distributed loading      */
            
       if(p->elmt_load_ptr != (ELEMENT_LOADS *) NULL) {
           p = sld04(p, PRESSLD); /* Equivalent Load in local_coordinate */
            
       /***************************************************/
       /* Transform nodal_load vector from local lamina   */
       /* coordinate system to global coordinate system   */
       /***************************************************/

       T_matrix    = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
       T_Transpose = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
       nodal_load  = MatrixAllocIndirectDouble(p->size_of_stiff, 1);

      /* store the direction matrix in T_matrix */
      /* and transfer vector p->nodal_load into */
      /* nodal_load as a matrix form            */

         for (i = 1; i <= p->size_of_stiff; i++)  {
             for(j = 1; j <= p->nodes_per_elmt; j++) {
                for( k = 1; k <= p->dof_per_node; k++) {
                   n1  = p->dof_per_node*(j-1);
                   n2  = p->dof_per_node*j;
                   n   = n1 + k;
                   if(i > n1 && i <= n2) {
                      ii = i - n1;
                      T_matrix[i-1][n-1] = Direction_Matrix[ii-1][k-1];
                   }
                   else
                      T_matrix[i-1][n-1] = 0.0;
                }
            }
         }

      /* =====================================*/
      /* Calculate the inverse of T_matrix    */
      /*     T_matrix inverse [T]^-1 =        */
      /*        = T_matrix trnspose : [T]^t   */
      /* =====================================*/

         for( i = 1; i <= p->size_of_stiff; i++) {
             for( j = 1; j <= p->size_of_stiff; j++) 
                T_Transpose[i-1][j-1] = T_matrix[j-1][i-1];
         }
	
         nodal_load = dMatrixMultRep(nodal_load,
                      T_Transpose, p->size_of_stiff,p->size_of_stiff,
                      p->equiv_nodal_load->uMatrix.daa, p->size_of_stiff,1);

         for(i = 1; i <= p->size_of_stiff; i++)
            p->equiv_nodal_load->uMatrix.daa[i-1][0]
            = nodal_load[i-1][0];

         MatrixFreeIndirectDouble(T_matrix, p->size_of_stiff);
         MatrixFreeIndirectDouble(T_Transpose, p->size_of_stiff);
         MatrixFreeIndirectDouble(nodal_load, p->size_of_stiff);
       }

/* ------------------------ UNITS -----------------------------------*/

      UNITS_SWITCH = CheckUnits();
      switch( UNITS_SWITCH ) {
        case ON:
         if(UnitsType == SI) {
            dp_length = DefaultUnits("m");
            dp_force  = DefaultUnits("N");
         }
         if(UnitsType == US) {
            dp_length = DefaultUnits("in");
            dp_force  = DefaultUnits("lbf");
         }
         dp_moment = UnitsMult( dp_force, dp_length );
             
         for(i= 1; i<= p->dof_per_node; i++) {
            for(j = 1; j <= p->nodes_per_elmt; j++) {
                k  = p->dof_per_node*(j-1)+i;
                k1 =  p->dof_per_node*j-3;
                if ( k <= k1) {      /* force units */
                   UnitsCopy( &(p->equiv_nodal_load->spRowUnits[k-1]), dp_force );
                }
                else {  /* k > k1 moment units */
                   UnitsCopy( &(p->equiv_nodal_load->spRowUnits[k-1]), dp_moment );
                }
             }
          }

          free((char *) dp_force->units_name);
          free((char *) dp_force);
          free((char *) dp_length->units_name);
          free((char *) dp_length);
          free((char *) dp_moment->units_name);
          free((char *) dp_moment);
      break;
      case OFF:
      break;
      default:
      break;
   }

#ifdef DEBUG
       printf("*** In elmt_shell() : print equiv_nodal_load\n\n");
       MatrixPrintIndirectDouble(p->equiv_nodal_load);
       printf("*** In elmt_shell() : end case EQUIV_NODAL_LOAD\n");
#endif

             break;

        case STRESS_UPDATE:

        /* update stress for given displacement */
        /* and displacement incremental         */
             break;
        case STRESS:
        case STRESS_LOAD:
        case LOAD_MATRIX:
        /* ====================================== */
        /* Form internal nodal load vector due to */
        /* stress at previous load step           */
        /* ====================================== */

#ifdef DEBUG
       printf("****** In elmt_shell_4nodes_implicit() : \n");
       printf("       enter cases: STRESS_UPDATE, LOAD_MATRIX, STRESS_LOAD \n");     
       printf(" elemt_state = %d \n", p->elmt_state);
#endif
          nodal_load      = MatrixAllocIndirectDouble(p->size_of_stiff, 1); 
          nodal_load_temp = MatrixAllocIndirectDouble(p->size_of_stiff, 1); 
          T_matrix        = MatrixAllocIndirectDouble(p->size_of_stiff, 
                                                      p->size_of_stiff);
          T_Transpose     = MatrixAllocIndirectDouble(p->size_of_stiff,
                                                      p->size_of_stiff);

          h  = p->work_section[11].value;    /* thickness of the shell */
          for(i = 1; i <= p->nodes_per_elmt; i++) {
              elmt_length = ABS(p->coord[0][i-1].value - p->coord[0][i].value);

              if(elmt_length != 0)
                 break;
           }

           aspect_ratio = h/elmt_length;
           k_bar = 2.0*(1+nu)*aspect_ratio*aspect_ratio; /* reduced shear factor */

/*         EE  = (p->mater_matrix->uMatrix.daa[0][0]
                 +p->mater_matrix->uMatrix.daa[1][1])*0.5; */

           EE  = E.value;

        /* =================================================*/
        /* nodal load due to the inital stress or stress at */
        /* previous step                                    */
        /* =================================================*/

      /* Integration loop over thickness */

          integ_coord = dVectorAlloc(no_integ_pts+1);
          weight      = dVectorAlloc(no_integ_pts+1);

          gauss(integ_coord,weight,no_integ_pts);

          for(ii = 1; ii <= no_integ_pts; ii++) {
              nodal_load_temp 
              = Shell_Nodal_Load_Plane(nodal_load_temp, p, co_coord,
                             integ_coord[ii], ii, E.value, nu, isw);
              for(i = 1; i<= p->size_of_stiff; i++) {
                  nodal_load[i-1][0] += nodal_load_temp[i-1][0]*weight[ii]*0.5*h;  
              }

          } /* gaussian integration ends */

          free((char *) integ_coord);
          free((char *) weight);

          for(i= 1; i<= p->size_of_stiff; i++) {
              p->nodal_loads[i-1].value = nodal_load[i-1][0];  
          }

       /***************************************************/
       /* Transform nodal_load vector from local lamina   */
       /* coordinate system to global coordinate system   */
       /***************************************************/

      /* =====================================*/
      /* Calculate the inverse of T_matrix    */
      /*     T_matrix inverse [T]^-1 =        */
      /*        = T_matrix trnspose : [T]^t   */
      /* =====================================*/

      /* store the direction matrix : T_matrix */

   if( isw == LOAD_MATRIX ) {
      for (i = 1; i <= p->size_of_stiff; i++)  {
           for(j = 1; j <= p->nodes_per_elmt; j++) {
               for( k = 1; k <= p->dof_per_node; k++) {
                    n1  = p->dof_per_node*(j-1);
                    n2  = p->dof_per_node*j;
                    n   = n1 + k;
                    if(i > n1 && i <= n2) {
                       ii = i - n1;
                       T_matrix[i-1][n-1] = Direction_Matrix[ii-1][k-1];
                    }
                    else 
                       T_matrix[i-1][n-1] = 0.0;
                }
           }
       }

      for( i = 1; i <= p->size_of_stiff; i++) {
         for( j = 1; j <= p->size_of_stiff; j++) 
           T_Transpose[i-1][j-1] = T_matrix[j-1][i-1];
      }
     
      size = p->size_of_stiff;	
      nodal_load_temp = dMatrixMultRep(nodal_load_temp,
                        T_Transpose, size, size, nodal_load, size,1);

      for(i= 1; i<= p->size_of_stiff; i++) {
          p->nodal_loads[i-1].value = nodal_load_temp[i-1][0];  
      }
    }

      MatrixFreeIndirectDouble(nodal_load,p->size_of_stiff);
      MatrixFreeIndirectDouble(nodal_load_temp,p->size_of_stiff);
      MatrixFreeIndirectDouble(T_Transpose,p->size_of_stiff);
      MatrixFreeIndirectDouble(T_matrix,p->size_of_stiff);


   /* ------------NODAL LOAD UNITS ------------------------*/
   /* The units type is determined by the SetUnitsType()   */
   /* ---------------------------------------------------- */

     switch( UNITS_SWITCH ) {
       case ON:
          if(UnitsType == SI || UnitsType == SI_US ) {
              dp_force  = DefaultUnits("N");
              dp_length = DefaultUnits("m");
          }
          else {
              dp_force  = DefaultUnits("lbf");
              dp_length = DefaultUnits("in");
          }

          /* node no 1 */
          UnitsCopy( p->nodal_loads[0].dimen, dp_force );
          UnitsCopy( p->nodal_loads[1].dimen, dp_force );
          UnitsCopy( p->nodal_loads[2].dimen, dp_force );
          UnitsMultRep( p->nodal_loads[3].dimen, dp_force, dp_length );
          UnitsCopy( p->nodal_loads[4].dimen, p->nodal_loads[3].dimen );

          /* node no > 1 */
          for(i = 2; i <= p->nodes_per_elmt; i++) {    
              for(j = 1; j <= p->dof_per_node; j++) {
                  k = p->dof_per_node*(i-1)+j;
                  if(j <= 3) 
                     UnitsCopy( p->nodal_loads[k-1].dimen, p->nodal_loads[0].dimen );
                  else 
                     UnitsCopy( p->nodal_loads[k-1].dimen, p->nodal_loads[3].dimen );
              }
          }
          free((char *) dp_length->units_name);
          free((char *) dp_length);
          free((char *) dp_force->units_name);
          free((char *) dp_force);
          break;
       case OFF:
          break;
       default:
          break;
     }

#ifdef DEBUG
printf("*** In elmt_shell_4nodes_implicit() : end case LOAD_MATRIX, STRESS_UPDATE, STRESS_LOAD \n");
#endif
	     break;

        case MASS_MATRIX:  /* form mass matrix */

#ifdef DEBUG
       printf("*** In elmt_shell_4nodes_implicit() : start case MASS\n");
       printf("                : Density = %14.4e\n", density.value);
       printf("                : shell_thickness = %8.2f\n", h);
       printf("                : Jac             = %8.2f\n", Jac);
#endif
   /*====================================================*/
   /*  CALCULATING MASS MATRIX IN CO_COORDINATE SYSTEM   */ 
   /*====================================================*/

      Shell_4Node_Mass(p, p->stiff, co_coord, density.value, h);

#ifdef DEBUG
       printf("*** In elmt_shell() : end case MASS\n");
#endif
             break;
        default:
             break;
    }

    MatrixFreeIndirectDouble(co_coord, p->no_dimen);
    MatrixFreeIndirectDouble(Direction_Matrix,6);

#ifdef DEBUG
       printf("*** leaving elmt_shell() \n");
#endif

    return(p);
}

/* Print SHELL_4N Element Properties */
#ifdef __STDC__
void print_property_shell_4n(EFRAME *frp, int i)
#else
void print_property_shell_4n(frp, i)
EFRAME    *frp;
int          i;                 /* elmt_attr_no */
#endif
{
int     UNITS_SWITCH;
ELEMENT_ATTR    *eap;

#ifdef DEBUG
       printf("*** Enter print_property_shell_4n()\n");
#endif

     UNITS_SWITCH = CheckUnits();
     eap = &frp->eattr[i-1];

     if( PRINT_MAP_DOF == ON ) {
        if(frp->no_dof == 3 || frp->no_dof == 2) { 
           printf("             ");
           printf("         : gdof [0] = %4d : gdof[1] = %4d : gdof[2] = %4d\n",
                           eap->map_ldof_to_gdof[0],
                           eap->map_ldof_to_gdof[1],
                           eap->map_ldof_to_gdof[2]);
        }

        if(frp->no_dof == 6) { /* 3d analysis */
           printf("             ");
           printf("         : dof-mapping : gdof[0] = %4d : gdof[1] = %4d : gdof[2] = %4d\n",
                           eap->map_ldof_to_gdof[0],
                           eap->map_ldof_to_gdof[1],
                           eap->map_ldof_to_gdof[2]);
           printf("             ");
           printf("                         gdof[3] = %4d : gdof[4] = %4d : gdof[5] = %4d\n",
                           eap->map_ldof_to_gdof[3],
                           eap->map_ldof_to_gdof[4],
                           eap->map_ldof_to_gdof[5]);
        } 
     }

     switch(UNITS_SWITCH) {
       case ON:
        UnitsSimplify( eap->work_material[0].dimen );
        UnitsSimplify( eap->work_material[2].dimen );
        UnitsSimplify( eap->work_material[5].dimen );
        UnitsSimplify( eap->work_section[2].dimen );
        UnitsSimplify( eap->work_section[10].dimen );
        if( eap->work_material[0].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Young's Modulus =  E = %16.3e %s\n",
                           eap->work_material[0].value/eap->work_material[0].dimen->scale_factor,
                           eap->work_material[0].dimen->units_name);
        }
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's ratio = nu = %16.3e   \n", eap->work_material[4].value);
        }
        if( eap->work_material[2].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Yielding Stress = fy = %16.3e %s\n",
                           eap->work_material[2].value/eap->work_material[2].dimen->scale_factor,
                           eap->work_material[2].dimen->units_name);
        }
	if( eap->work_material[5].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Density         = %16.3e %s\n",
                           eap->work_material[5].value/eap->work_material[5].dimen->scale_factor,
                           eap->work_material[5].dimen->units_name);
	}
	if( eap->work_section[2].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Inertia Izz     = %16.3e %s\n",
                           eap->work_section[2].value/eap->work_section[2].dimen->scale_factor,
                           eap->work_section[2].dimen->units_name);
	}
	if( eap->work_section[10].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Area            = %16.3e %s\n",
                           eap->work_section[10].value/eap->work_section[10].dimen->scale_factor,
                           eap->work_section[10].dimen->units_name);
	}
       break;
       case OFF:
         printf("             ");
         printf("         : Young's Modulus =  E = %16.3e\n",
                            eap->work_material[0].value);
         printf("             ");
         printf("         : Yielding Stress = fy = %16.3e\n",
                            eap->work_material[2].value);
         printf("             ");
         printf("         : Poisson's ratio = nu = %16.3e   \n", eap->work_material[4].value);
         printf("             ");
         printf("         : Density         = %16.3e\n",
                            eap->work_material[5].value);
         printf("             ");
         printf("         : Inertia Izz     = %16.3e\n",
                            eap->work_section[2].value);
         printf("             ");
         printf("         : Area            = %16.3e\n",
                            eap->work_section[10].value);
        break;
        default:
        break;
     }
#ifdef DEBUG
       printf("*** Leave print_property_shell_4n()\n");
#endif
}

/* ============================================*/
/* Equivalent Load due to distributed pressure */
/* for shell Belytschko element                */
/* in local coordinate                         */
/* ============================================*/

ARRAY *sld04(p,task)
ARRAY *p;
int task;
{
ELEMENT_LOADS               *elsptr;
ELOAD_LIB                      *elp;

double                  *nodal_load;
double      Jac, x31, y31, x42, y42; 
double px,py,pz, mx,my,mz, bx,by,bz;

int  i,j,k, ii, jj, kk, n, n1,n2,n3;

#ifdef DEBUG
      printf("**** enter sld04(): \n");
#endif

    /* Initialize total load */

    nodal_load = dVectorAlloc(p->size_of_stiff);

    for(i =1; i <= p->size_of_stiff ;i++)
        nodal_load[i-1] = 0.0;

#ifdef DEBUG
      printf(" **** In sld04(): begin switch()  \n");
#endif

    switch(task){
       case PRESSLD:
#ifdef DEBUG
      printf(" **** In sld04(): case PRESSLD\n");
#endif

       /* ==========================================*/
       /* For distrbuted surface loading:           */
       /* PRESSURE OR TRACTION                      */
       /* nodal equivalent force is calculated by : */
       /* force = integral N^T*traction/pressure dA */
       /* ==========================================*/
       
       /*================================================*/
       /* use one-point integration rule in surface area */
       /* Ni = 1/4.0 [N] = 1/4 [I, I, I, I]              */
       /* [I] = 5x5 unit_matrix                          */
       /* Let traction = [p]5x1                          */
       /* Integ [N]^t*traction dA  = (1/4)*[p]*Jac*4.0   */
       /*                                  [p]           */
       /*                                  [p]           */
       /*                                  [p]           */
       /*                                  [p]           */
       /*================================================*/
         elsptr  =  p->elmt_load_ptr;
         for (j=1; j<= elsptr->no_loads_faces;j++) {
              elp = &elsptr->elib_ptr[j-1];
       
              x31 =  p->coord[0][2].value - p->coord[0][0].value;
              y31 =  p->coord[1][2].value - p->coord[1][0].value;

              x42 =  p->coord[0][3].value - p->coord[0][1].value;
              y42 =  p->coord[1][3].value - p->coord[1][1].value;

              Jac = 0.5*(x31*y42-x42*y31);
 
              for(i = 1; i <= p->size_of_stiff; i++) {
                  n = i/p->dof_per_node;
                  k = i - n*p->dof_per_node;
                  nodal_load[i-1] = 0.25*elp->traction[k-1].value*4.0*Jac; 
              }
         }
       break;
    default:
       break;
    }

    for(i = 1; i <= p->size_of_stiff; i++) {
       p->equiv_nodal_load->uMatrix.daa[i-1][0] = nodal_load[i-1];
    }
    
    free(nodal_load);

#ifdef DEBUG
      printf("**** leaving sld04(): \n");
#endif
    return(p);
}


/* ============================*/
/* Material matrix             */
/* ============================*/

double **MATER_MAT_SHELL(m1, E, nu)
double  **m1;
double E, nu;
{
double           temp;
double       kk = 1.2;       
    
#ifdef DEBUG
     printf(" enter MATER_MAT_SHELL() \n");
     printf(" E = %lf , nu = %lf \n", E, nu);
#endif
    /* Stress : stress_x, stress_y,  stress_xy,  stress_yz,  stress_zx  */
    /* Strain : strain_x, strain_y, 2strain_xy, 2strain_yz, 2strain_zx  */

    
    temp =  E/(1-nu*nu);

    m1[0][0]= m1[1][1] = temp;
    m1[0][1]= m1[1][0] = nu*temp;
    m1[2][2]= E/2.0/(1.0+nu);
    m1[3][3] = m1[4][4] = E/2.0/(1.0+nu)/kk;

    m1[0][2] = m1[2][0] = m1[1][2] = m1[2][1] = 0.0;
    m1[3][0] = m1[3][1] = m1[3][2] = m1[3][4] = 0.0;
    m1[4][0] = m1[4][1] = m1[4][2] = m1[4][3] = 0.0;


#ifdef DEBUG
     dMatrixPrint("m1", m1, 5, 5);
     printf("\n leaving MATER_MAT_SHELL() \n");
#endif
  return (m1);
}

void MATER_SHELL_UPDATE(p, co_coord, E, nu, integ_pt)
ARRAY                              *p;
double                     **co_coord;
double                          E, nu;
int                          integ_pt; /* integration point */
{
double             **m1, **mater_temp;
double **stress_dev, **Norm_Transpose;
double         A_2, mean_stress, R, A;
double                           temp;
int            i, j, k, ii, jj, kk, n;
int                      dof, size, G;

#ifdef DEBUG
    printf(" enter MATER_SHELL_UPDATE() \n");
#endif

  /* Given i -state displ, stress, [C] and strain */
  /* calculate i+1 state displ, stress, these i+1 */
  /* state stress can be used for i+1 state [Cep] */
  /* matrix                                       */

#ifdef DEBUG
    printf(" In MATER_SHELL_UPDATE(): calculating MATER_MAT_SHELL() \n");
#endif
      dof             = p->dof_per_node;
      size            = p->size_of_stiff;
      ii              = integ_pt;
      G               = E/2.0/(1.0+nu);


#ifdef DEBUG
    printf(" In MATER_SHELL_UPDATE(): update p->material = %s \n", p->material_name);
    printf(" In MATER_SHELL_UPDATE(): integ_pts          = %d \n", integ_pt);
    printf(" In MATER_SHELL_UPDATE(): dof                = %d \n", dof);
    printf(" In MATER_SHELL_UPDATE(): deformation state  = %d \n", p->elmt_state);
#endif

  /* Elastic material */

  if(p->material_name != NULL && !strcmp(p->material_name, "ELASTIC")) {

     p->mater_matrix->uMatrix.daa 
     = MATER_MAT_SHELL(p->mater_matrix->uMatrix.daa,E, nu);
   }

  /* Elastic Perfectly Plastic Materials */


  if(p->material_name != NULL &&
     !strcmp(p->material_name,"ELASTIC_PERFECTLY_PLASTIC")) {

     m1              = MatrixAllocIndirectDouble(p->dof_per_node, p->dof_per_node);
     stress_dev      = MatrixAllocIndirectDouble(dof, 1);
     Norm_Transpose  = MatrixAllocIndirectDouble(1,dof);
     mater_temp      = MatrixAllocIndirectDouble(dof, dof);

     /* Elastic Deformation */

     switch(p->elmt_state) {

       /* Elastic deformation */
       case 0: 
         p->mater_matrix->uMatrix.daa 
         = MATER_MAT_SHELL(p->mater_matrix->uMatrix.daa,E, nu);
       break;
       case 1:
       /* plastic deformation */

       /* Step 1: Calculate the spherical stress and */
       /*         deviatric stress of trial stress   */
       /*         also the radiaus square A_2        */

            mean_stress = (p->stress->uMatrix.daa[0][ii-1] - p->LC_ptr->back_stress[0][ii-1]
                          +p->stress->uMatrix.daa[1][ii-1] - p->LC_ptr->back_stress[1][ii-1])/3.0;
            A_2 = 0.0;
            for(i = 1; i <= dof; i++) {
               if(i <= 2)           
                 stress_dev[i-1][0] = p->stress->uMatrix.daa[i-1][ii-1] - mean_stress - p->LC_ptr->back_stress[i-1][ii-1];
               else
                 stress_dev[i-1][0] = p->stress->uMatrix.daa[i-1][ii-1] - p->LC_ptr->back_stress[i-1][ii-1];

               A_2 += stress_dev[i-1][0]*stress_dev[i-1][0]; 
            }

       /* Step 2: comparision                        */

            R = p->LC_ptr->R[ii-1];

            if(A_2 <= R*R) { /* ELASTIC DEFORMATION */
               p->mater_matrix->uMatrix.daa 
               = MATER_MAT_SHELL(p->mater_matrix->uMatrix.daa,E, nu);
            }
            else {           /* PLASTIC DEFORMATION */
        
              m1  = MATER_MAT_SHELL(m1, E, nu);

              /* Step 3: normal direction on yield surface */ 
              /*         stored temporarily on stress_dev  */
              /* Norm_Direc */

              for(i = 1; i <= dof; i++) {
                 stress_dev[i-1][0] /= R; 
                 Norm_Transpose[0][i-1] = stress_dev[i-1][0];
              }
                 
              /* Step 4: Calculate [Cep]                   */ 
              /*         [Cep] = [C] -[C]*[n][n]^t         */
              
               p->mater_matrix->uMatrix.daa 
               = dMatrixMultRep(p->mater_matrix->uMatrix.daa,
                 stress_dev, dof,1, Norm_Transpose, 1, dof);
               mater_temp
               = dMatrixMultRep(mater_temp, m1, dof, dof,
                 p->mater_matrix->uMatrix.daa, dof, dof);

               for(i = 1; i <= dof; i++) {
                  for(j = 1; j <= dof; j++) {
                     p->mater_matrix->uMatrix.daa[i-1][j-1]
                     = m1[i-1][j-1] - mater_temp[i-1][j-1];
                  }
               }
         
            } /* end of else : plastic deformation */

       break;
       default:
         printf(" In MATER_SHELL_UPDATE(): elmt_no = %d\n", p->elmt_no);
         printf(" elmt_state = 0 : Elastic_deformation \n");
         printf(" elmt_state = 1 : Elastic-plastic_deformation \n");
         printf(" elmt_state = %d: p->elmt_state \n");
         FatalError(" Unknown elmt state ",(char *)NULL);
       break;
     }

     MatrixFreeIndirectDouble(stress_dev, p->dof_per_node);
     MatrixFreeIndirectDouble(Norm_Transpose, 1);
     MatrixFreeIndirectDouble(m1, p->dof_per_node);
     MatrixFreeIndirectDouble(mater_temp, p->dof_per_node);
  }

  /* Elastic Plastic Materials      */
  /* with combined strain hardening */

  if(p->material_name != NULL &&
     !strcmp(p->material_name,"ELASTIC_PLASTIC")) {

     m1              = MatrixAllocIndirectDouble(p->dof_per_node, p->dof_per_node);
     stress_dev      = MatrixAllocIndirectDouble(dof, 1);
     Norm_Transpose  = MatrixAllocIndirectDouble(1,dof);
     mater_temp      = MatrixAllocIndirectDouble(dof, dof);
    
     switch(p->elmt_state) {
       /* Elastic deformation */
       case 0: 
         p->mater_matrix->uMatrix.daa 
         = MATER_MAT_SHELL(p->mater_matrix->uMatrix.daa, E, nu);
       break;
       case 1:
       /* plastic deformation */

       /* Step 1: Calculate the spherical stress and */
       /*         deviatric stress of trial stress   */
       /*         also the radiaus square A_2        */

         mean_stress = (p->stress->uMatrix.daa[0][ii-1]- p->LC_ptr->back_stress[0][ii-1]
                       +p->stress->uMatrix.daa[1][ii-1]- p->LC_ptr->back_stress[1][ii-1])/3.0;
         
         A_2 = 0.0;
         for(i = 1; i <= dof; i++) {
             if(i <= 2)           
               stress_dev[i-1][0] = p->stress->uMatrix.daa[i-1][ii-1]
                                  - mean_stress - p->LC_ptr->back_stress[i-1][ii-1];
             A_2 += stress_dev[i-1][0]*stress_dev[i-1][0]; 
         }
         A = sqrt(A_2);

       /* Step 2: comparision   */

         R = p->LC_ptr->R[ii-1];

         if(A_2 <= R*R) { /* ELASTIC DEFORMATION */
            p->mater_matrix->uMatrix.daa 
            = dMatrixCopyRep(p->mater_matrix->uMatrix.daa,
              m1, p->dof_per_node, p->dof_per_node);
         }
         else {  /* PLASTIC DEFORMATION */
        
           m1  = MATER_MAT_SHELL(m1, E, nu);

           /* Step 3: normal direction on yield surface */ 
           /*         stored temporarily on stress_dev  */
           for(i = 1; i <= dof; i++) {
               stress_dev[i-1][0] /= R;  /* Norm_Direc */
               Norm_Transpose[0][i-1] = stress_dev[i-1][0];
           }

       /* Step 4: Calculate [Cep]                         */ 
       /*      [Cep] = [C] -[C]*[n][n]^t*[C]/(2G+H*2/3) */
       /*      [Cep] = [C] -[C]*[n][n]^t*2G/(2G+H*2/3)  */
              
             p->mater_matrix->uMatrix.daa = dMatrixMultRep(p->mater_matrix->uMatrix.daa,
                                            stress_dev, dof,1, Norm_Transpose, 1, dof);

             mater_temp = dMatrixMultRep(mater_temp, m1, dof, dof,
                          p->mater_matrix->uMatrix.daa, dof, dof);

/*
             p->mater_matrix->uMatrix.daa = dMatrixMultRep(p->mater_matrix->uMatrix.daa,
                                            mater_temp, dof, dof, m1, dof, dof);
*/
             for(i = 1; i <= dof; i++) {
                for(j = 1; j <= dof; j++) {
                    p->mater_matrix->uMatrix.daa[i-1][j-1] 
                    = m1[i-1][j-1]- p->mater_matrix->uMatrix.daa[i-1][j-1]
                      *2.0*G/(2.0*G+p->LC_ptr->H[ii-1]*2.0/3.0);
                }
             }
         
         } /* end of else : plastic deformation */

       break;
       default:
         printf(" In MATER_SHELL_UPDATE():  elmt_no  = %d\n", p->elmt_no);
         printf(" elmt_state = 0 : Elastic_deformation \n");
         printf(" elmt_state = 1 : Elastic-plastic_deformation \n");
         printf(" elmt_state = %d: p->elmt_state \n");
         FatalError(" Unknown elmt state ",(char *)NULL);
       break;
     }
     MatrixFreeIndirectDouble(stress_dev, p->dof_per_node);
     MatrixFreeIndirectDouble(Norm_Transpose, 1);
     MatrixFreeIndirectDouble(m1, p->dof_per_node);
     MatrixFreeIndirectDouble(mater_temp, p->dof_per_node);
  }

#ifdef DEBUG
    printf(" Leaving MATER_SHELL_UPDATE() \n");
#endif
}


void DISPL_UPDATE(p)
ARRAY    *p;
{
int i,j, k;

   for(i = 1; i <= p->dof_per_node; i++) {
      for(j = 1; j <= p->nodes_per_elmt; j++) {
          p->displ->uMatrix.daa[i-1][j-1] += p->displ_incr->uMatrix.daa[i-1][j-1];
      }
   } 

}

void Stress_Update(p, co_coord, B_matrix, integ_pt)
ARRAY                   *p;
double          **co_coord;
double          **B_matrix;
int               integ_pt;
{
double                      H, R, fy, E, Et, nu;
double          **displ_incr,**m1,**strain_incr;
double          **back_stress_incr, **strain_pl;
double    **stress, **stress_dev, **stress_incr;
double     temp, G, A_2, A, Lambda, mean_stress;
double         eff_pl_strainTemp, eff_pl_strain;
double       Beta1 = 0.05, temp1, effect_stress;
int                    iNo_iter_step, dof, size;
int                      i, j, k, ii, jj, kk, n;
DIMENSIONS                               *dimen;

  /* Given i -state displ, stress mater_matrix, [Cep], */
  /* and strain,  calculate i+1 state displ_incr       */
  /* these i+1 state stress can be used for            */
  /* i+1 state [Cep] matrix                             */


#ifdef DEBUG
     printf("\n enter Stress_Update(): \n");
#endif

  kk   = integ_pt;
  dof  = p->dof_per_node;
  size = p->size_of_stiff;
  E    = p->work_material[0].value;
  nu   = p->work_material[4].value;
  fy   = p->work_material[2].value;
  G    = E/(1.0+nu)/2.0;

  m1         = MatrixAllocIndirectDouble(p->dof_per_node, p->dof_per_node);
  m1         = MATER_MAT_SHELL(m1, E, nu);

  /* Step 1: calculate strain_incremental  */
  /*         from displacement incremental */
  /*         at given integration point    */

   displ_incr       = MatrixAllocIndirectDouble(p->size_of_stiff, 1);
   strain_incr      = MatrixAllocIndirectDouble(dof,1);
   strain_pl        = MatrixAllocIndirectDouble(dof,1);
   stress           = MatrixAllocIndirectDouble(dof,1); /* stress or stress incremental   */
   stress_dev       = MatrixAllocIndirectDouble(dof,1); /* deviatric component of stress  */
   stress_incr      = MatrixAllocIndirectDouble(dof,1); /* stress incremental             */
   back_stress_incr = MatrixAllocIndirectDouble(dof,1); /* back_stress stress incremental */

#ifdef DEBUG
      printf(" In Stress_Update(): materials_name = %s \n",p->material_name);
#endif

  /* ================= */
  /* Elastic material  */
  /* ================= */

  if(p->material_name != NULL && !strcmp(p->material_name,"ELASTIC")) {

     /* Transfer p->displ into vector form */
     /* -----------------------------------*/
     /* There is no needed to calculate    */
     /* incrementally for the case of      */
     /* linear elasticity                  */
     /* -----------------------------------*/

     for(i = 1; i <= p->dof_per_node; i++) {
         for(j = 1; j <= p->nodes_per_elmt; j++) {
             k = p->dof_per_node*(j-1)+i;
             displ_incr[k-1][0] = p->displ->uMatrix.daa[i-1][j-1];
         }
     }

     strain_incr = dMatrixMultRep(strain_incr, B_matrix, dof,
                                  size, displ_incr, size, 1); 
     stress      = dMatrixMultRep(stress, m1, dof, dof,
                                  strain_incr, dof, 1);

     for (i = 1; i <= dof; i++){
         p->stress->uMatrix.daa[i-1][kk-1] += stress[i-1][0];
     }
     SaveRespondBuffer(p, kk);
  }

  /* ====================================*/
  /* Elastic Plastic & Elastic Perfectly */
  /* Plastic Materials                   */
  /* ====================================*/

  if(p->material_name != NULL &&
      (!strcmp(p->material_name, "ELASTIC_PLASTIC") ||
       !strcmp(p->material_name, "ELASTIC_PERFECTLY_PLASTIC")) ) {

     /* Transfer p->displ_incr into vector form */

     switch(p->elmt_state) {
       case 0 :       /* elastic state */
         for(i = 1; i <= p->dof_per_node; i++) {
             for(j = 1; j <= p->nodes_per_elmt; j++) {
                 k = p->dof_per_node*(j-1)+i;
                 displ_incr[k-1][0] = p->displ->uMatrix.daa[i-1][j-1];
             }
         }
       break;
       case 1 :       /* inelastic state */
         for(i = 1; i <= p->dof_per_node; i++) {
             for(j = 1; j <= p->nodes_per_elmt; j++) {
                 k = p->dof_per_node*(j-1)+i;
                 displ_incr[k-1][0] = p->displ_incr->uMatrix.daa[i-1][j-1];
             }
         }
       break;
       default:
         printf("**** In Stress_Update(): elmt_no = %d\n", p->elmt_no);
         printf(" elmt_state = %d: p->elmt_state \n");
         FatalError("*****Undefine element state ()", (char *) NULL);
       break;
     }

     /* Calculate stress increment */

     strain_incr = dMatrixMultRep(strain_incr, B_matrix, dof,
                                  size, displ_incr, size, 1); 

     stress      = dMatrixMultRep(stress, m1, dof, dof,
                                  strain_incr, dof, 1);

#ifdef DEBUG
     dMatrixPrint("displ_incr", displ_incr, size, 1);
     dMatrixPrint("stress", stress, dof, 1);
     dMatrixPrint("strain", strain_incr, dof, 1);
     printf(" \n stress before incremented \n");
     for(i = 1; i <= dof; i++)
         printf(" p->stress->uMatrix.daa[%d][%d]= %lf \n",
               i, kk, p->stress->uMatrix.daa[i-1][kk-1]);
#endif

     switch(p->elmt_state) {

       /* Elastic deformation */
       case 0: 
         for(i = 1; i <= dof; i++)
             p->stress->uMatrix.daa[i-1][kk-1] = stress[i-1][0];

         SaveRespondBuffer(p, kk);

       break;
       case 1:
         
       /* plastic deformation */
       /* Step 1: Calculate the spherical stress and */
       /*         deviatric stress of trial stress   */
       /*         also the radiaus square A_2        */

#ifdef DEBUG
     for (i = 1; i <= dof; i++)
         printf("incrmental stress[%d] = %lf \n", i, stress[i-1][0]);
#endif
         for (i = 1; i <= dof; i++) {
           stress_incr[i-1][0] = stress[i-1][0];
           stress[i-1][0] += p->stress->uMatrix.daa[i-1][kk-1];
         }

#ifdef DEBUG
      for (i = 1; i <= dof; i++)
         printf(" total stress[%d] = %lf \n", i, stress[i-1][0]);
#endif
         /* Note : sigma_z == 0 for shell elmt */
         
         mean_stress    = (stress[0][0] - p->LC_ptr->back_stress[0][kk-1]
                          +stress[1][0] - p->LC_ptr->back_stress[1][kk-1])/3.0;
         A_2  = 0.0;
         for(i = 1; i <= dof; i++) {
             if(i <= 2)
               stress_dev[i-1][0] = stress[i-1][0] - mean_stress
                                    - p->LC_ptr->back_stress[i-1][kk-1];
             else
               stress_dev[i-1][0] = stress[i-1][0] -
                                    p->LC_ptr->back_stress[i-1][kk-1];

             A_2  += stress_dev[i-1][0]*stress_dev[i-1][0];
         }
         A = sqrt(A_2);
         R = p->LC_ptr->R[kk-1];
         eff_pl_strain = p->effect_pl_strain[kk-1];

       /* Step 2: comparision                        */

        if(A <= R) { /* ELASTIC DEFORMATION */
          for(i = 1; i <= dof; i++) 
             p->stress->uMatrix.daa[i-1][kk-1] = stress[i-1][0];

           SaveRespondBuffer(p, kk);

#ifdef DEBUG1
    printf(" +++++++ ELASTIC DEFORMATION A = %lf R = %lf \n", A, R);
    printf(" at elmt_no = %d, integ_pt = %d\n", p->elmt_no, kk);
#endif
        }
        else {       /* PLASTIC DEFORMATION */
        

       /* Step 3 Estimate number of sub-incrementations needed     */

          if( ABS(p->LC_ptr->beta) < 1E-10){ /* Only for beta = 0 kinematic hardening case */

       /* Step 3.1 Estimate the effective plastic strain increment */

              temp = sqrt(3.0/2.0)*(A-R)/(p->LC_ptr->H[kk-1]+3.0*G);

       /* Step 3.2 Estimate the back stress increment              */
           
              /* Estimate H' */
              if(!strcmp(p->material_name, "ELASTIC_PERFECTLY_PLASTIC")) {
                  H = 0.0;
              }else {
                  if(!strcmp(p->LC_ptr->name, "Ramberg-Osgood")) {
                      effect_stress = A*sqrt(3.0/2.0);
                      Load_Curve(p, &H, effect_stress, E,fy);
                  }
                  if(!strcmp(p->LC_ptr->name, "Bi-Linear")) 
                      H = p->LC_ptr->H[kk-1];
             }
          /* calculate the effective plastic strain incremental */

             temp1 = sqrt(2.0/3.0)*(H + 3.0*G);
             temp  = A/temp1 - R/temp1;            /* eff_pl_strain_incr */

       /* Step 3.3: Compute Lambda  and pl_strain_incr   */
       /*         Lambda = sqrt(3/2)*eff_pl_strain_incr  */
       /* plastic strain incr is now stored in strain_incr */

             Lambda = sqrt(3.0/2.0)*temp;
             for(i = 1; i <= dof; i++) 
                 strain_incr[i-1][0] = Lambda*stress_dev[i-1][0]/A;

          /* Step 3.4:  calculate back stress increment  */
          /* beta = 0 for kinematic hardening */
          /* beta = 1 for isotropic hardening */

              for(i = 1; i <= dof; i++) {
                  back_stress_incr[i-1][0] = H*strain_incr[i-1][0]*2.0/3.0;
              }
          }

          mean_stress = (stress_incr[0][0] - back_stress_incr[0][0]
                        +stress_incr[1][0] - back_stress_incr[1][0])/3.0;
          temp   = 0.0;
          for(i = 1; i <= dof; i++) {
              if(i <= 2)
                 stress_dev[i-1][0] = stress_incr[i-1][0] - mean_stress - back_stress_incr[i-1][0];
              else
                 stress_dev[i-1][0] = stress_incr[i-1][0] - back_stress_incr[i-1][0];

              temp += stress_dev[i-1][0]*stress_dev[i-1][0];
          }
          temp = sqrt(temp);
          iNo_iter_step = (int) (2*temp/R/Beta1) + 1 ;

#ifdef DEBUG1
    printf(" ******Plastic DEFORMATION A  = %lf, R = %lf \n", A, R);
    printf(" ******Plastic DEFORMATION dA = %lf, R = %lf \n", temp, R);
    printf(" at elmt_no = %d, integ_pt = %d\n", p->elmt_no, kk);
    printf(" No of sub-Incremental steps = %d \n",iNo_iter_step);
#endif
#undef DEBUG1
          /* Step 4 Start sub-incrementation  */

         /* copy the plastic strain before sub-incrementation */

         eff_pl_strainTemp = eff_pl_strain;
         for(i = 1; i <= dof; i++) 
             strain_pl[i-1][0] = p->strain_pl->uMatrix.daa[i-1][kk-1];

         switch(iNo_iter_step) {
             case 1: 
               ii = 1;
               Plastic_Deform(p, &H, &R, &eff_pl_strain, stress, stress_dev,
                              strain_incr, E, fy, G, A, iNo_iter_step, dof, ii, kk);
               SaveRespondBuffer(p, kk);
#ifdef DEBUG
   printf(" A= %lf R = %lf H = %lf eff_pl_strain = %lf\n", A, R, H, eff_pl_strain);
#endif
             break;
             default: 
               for(i = 1; i <= p->dof_per_node; i++) {
                   for(j = 1; j <= p->nodes_per_elmt; j++) {
                       k = p->dof_per_node*(j-1)+i;
                       displ_incr[k-1][0] 
                       = p->displ_incr->uMatrix.daa[i-1][j-1]/((double) iNo_iter_step);
                   }
               }
               /* Calculate stress increment */

               strain_incr = dMatrixMultRep(strain_incr,B_matrix, dof, 
                                            size, displ_incr, size, 1);
               stress_incr = dMatrixMultRep(stress_incr, m1, dof, dof,
                                            strain_incr, dof, 1);
               for(ii = 1; ii <= iNo_iter_step; ii++) {
                   /* Trial stress */
                   if(ii == 1) {
                      for(i = 1; i <= dof; i++) 
                          stress[i-1][0] = stress_incr[i-1][0] + 
                                           p->stress->uMatrix.daa[i-1][kk-1];
                   }
                   else
                      for(i = 1; i <= dof; i++) 
                          stress[i-1][0] += stress_incr[i-1][0];

                   mean_stress = (stress[0][0] - p->LC_ptr->back_stress[0][kk-1]
                                 +stress[1][0] - p->LC_ptr->back_stress[1][kk-1])/3.0;
                   A_2 = 0.0;
                   for(i = 1; i <= dof; i++) {
                       if(i <= 2)
                          stress_dev[i-1][0] = stress[i-1][0] - mean_stress
                                             - p->LC_ptr->back_stress[i-1][kk-1];
                       else
                          stress_dev[i-1][0] = stress[i-1][0]
                                             - p->LC_ptr->back_stress[i-1][kk-1];
                       A_2 += stress_dev[i-1][0]*stress_dev[i-1][0];
                   }
                   A = sqrt(A_2);
                   if(A <= R) { /* ELASTIC DEFORMATION */
                      if(i == iNo_iter_step){
                         for(i = 1; i <= dof; i++)
                             p->stress->uMatrix.daa[i-1][kk-1] = stress[i-1][0];
                      }
                   /* go to next sub incremental iteration */

                   }else {   /* PLASTIC DEFORMATION */

                      Plastic_Deform(p, &H, &R, &eff_pl_strain, stress,stress_dev,
                                     strain_incr, E, fy, G, A, iNo_iter_step, dof, ii, kk);
                   }
               } /* end of sub-incremental iteration */
             break;
          } /* end of switch for sub incrementation */
          
          p->LC_ptr->R[kk-1] = R;
          p->LC_ptr->H[kk-1] = H;
          p->effect_pl_strain[kk-1]   = eff_pl_strain;
          p->eff_pl_strain_incr[kk-1] = eff_pl_strain - eff_pl_strainTemp;
          for(i = 1; i <= dof; i++) {
              p->strain_pl_incr->uMatrix.daa[i-1][kk-1]
              = p->strain_pl->uMatrix.daa[i-1][kk-1] - strain_pl[i-1][0];
          }
          SaveRespondBuffer(p, kk);
        }
       break;
       default:
         printf(" In Stress_Update(): elmt_no \n", p->elmt_no);
         printf(" elmt_state = 0 : Elastic_deformation \n");
         printf(" elmt_state = 1 : plastic_deformation \n");
         printf(" elmt_state = %d: p->elmt_state \n");
         FatalError(" Unknown elmt state ",(char *)NULL);
       break;
     }
  }

  /* ASSIGN UNITS TO p ARRAY */
  
   if(CheckUnits() == ON) {
         switch(CheckUnitsType()) {
           case SI:
             dimen = DefaultUnits("Pa");
           break;
           case US:
             dimen = DefaultUnits("psi");
           break;
       }
       for(i = 1; i <= dof; i++)
           p->stress->spRowUnits[i-1] = *DefaultUnits("psi");

       free((char *) dimen->units_name);
       free((char *) dimen);
   }

   MatrixFreeIndirectDouble(m1, dof);
   MatrixFreeIndirectDouble(strain_incr, dof);
   MatrixFreeIndirectDouble(stress, dof); 
   MatrixFreeIndirectDouble(stress_dev,dof);
   MatrixFreeIndirectDouble(stress_incr,dof);
   MatrixFreeIndirectDouble(back_stress_incr,dof);

#ifdef DEBUG
      dMatrixPrint("p->stress in Stress_Update() leaving ", p->stress->uMatrix.daa, p->dof_per_node, 12);
     printf(" Leaving Stress_Update() \n");
#endif
}

void
Plastic_Deform(p, H, R, eff_pl_strain, stress, stress_dev,
               strain_incr, E, fy, G, A, iNo_iter_step, dof, ii, kk)
ARRAY                              *p;
double            *H, *R, fy, E, G, A;
double         **stress, **stress_dev;
double                  **strain_incr;
double                 *eff_pl_strain;
int        iNo_iter_step, dof, ii, kk;
{
double                           temp;
double                         Lambda;
double             eff_pl_strain_incr;
double           temp1, effect_stress;
int                           i, j, k;

#ifdef DEBUG
    printf(" enter Plastic_Deform() \n");
#endif

  /* Step 5 Compute effective incremental */
  /*        plastic strain within each    */
  /*        sub-incrementation            */
              
  /* Estimate H' */
  if(!strcmp(p->material_name, "ELASTIC_PERFECTLY_PLASTIC")) {
     *H = 0.0;
  }else {
     if(!strcmp(p->LC_ptr->name, "Ramberg-Osgood")) {
         effect_stress = A*sqrt(3.0/2.0);

#ifdef DEBUG
    printf(" before H = %le \n", *H);
    printf(" fy = %lf \n", fy);
    printf(" E  = %lf \n", E);
    printf(" effect_stress= %le \n", effect_stress);
#endif
         Load_Curve(p, H, effect_stress, E,fy);
     }
     if(!strcmp(p->LC_ptr->name, "Bi-Linear")) 
         *H = p->LC_ptr->H[kk-1];
  }
  /* calculate the effective plastic strain incremental */

  temp1 = sqrt(2.0/3.0)*(*H + 3.0*G);
  eff_pl_strain_incr  = A/temp1 - (*R)/temp1;
 *eff_pl_strain      += eff_pl_strain_incr;

#ifdef DEBUG3
     printf("========== In Plastic_Deform() : H = %le\n", *H);
     printf("========== In Plastic_Deform() : eff_pl_strain_incr = %le\n", eff_pl_strain_incr);
#endif

  /* Step 6: Compute Lambda  and pl_strain_incr       */
  /*         Lambda = sqrt(3/2)*eff_pl_strain_incr/A  */
  /* plastic strain incr is now stored in strain_incr */

  Lambda = sqrt(3.0/2.0)*eff_pl_strain_incr;
  for(i = 1; i <= dof; i++) {
      strain_incr[i-1][0] = Lambda*stress_dev[i-1][0]/A;
      p->strain_pl->uMatrix.daa[i-1][kk-1] += strain_incr[i-1][0];
  }
  /* Step 7:  Update back stress and Stress */

  /* beta = 0 for kinematic hardening */
  /* beta = 1 for isotropic hardening */

  if( ABS(p->LC_ptr->beta) < 1E-10){
      for(i = 1; i <= dof; i++) {
          temp = (*H)*strain_incr[i-1][0]*2.0/3.0;
          stress[i-1][0] = stress[i-1][0]-2.0*G*strain_incr[i-1][0] + temp;
          p->LC_ptr->back_stress[i-1][kk-1] += temp;
      }
  }else {

  /* Step 8:  Update stress */
  for(i = 1; i <= dof; i++) 
      stress[i-1][0] = stress[i-1][0]-2.0*G*strain_incr[i-1][0];
 }

  if(ii == iNo_iter_step) 
     for(i = 1; i <= dof; i++) 
      p->stress->uMatrix.daa[i-1][kk-1] = stress[i-1][0];

  /* Step 9:  Update R */
#ifdef DEBUG
    printf(" before R = %lf \n", *R);
    printf(" H       = %le \n",  *H);
    printf(" beta    = %lf \n", p->LC_ptr->beta);
    printf(" eff_pl_strain_incr = %le \n", eff_pl_strain_incr);
#endif
    if( ABS(p->LC_ptr->beta -1.0) < 1E-10)
       (*R) += sqrt(2.0/3.0)*(*H)*eff_pl_strain_incr;

#ifdef DEBUG
    printf(" after R = %lf \n", *R);
#endif

#ifdef DEBUG
    printf(" Leaving Plastic_Deform() \n");
#endif

}

double **B_MATRIX_4Node(B_matrix, p, shp, z_coord)
double       **B_matrix;
ARRAY                *p;
double            **shp;
double          z_coord;
{
int       i,j, k, ii, n;
double                h;


/* ======================================================= */
/* B_matrix : strain = Transpose(B_matrix)* nodal_displ.   */
/* ======================================================= */

#ifdef DEBUG 
     printf(" enter B_MATRIX_4Node() \n");
#endif

      /* B_matrix = [B', B"] */

      for(j = 1; j <= p->nodes_per_elmt; j++) {
          k = p->dof_per_node*(j-1)-1;

      /* ------------------------------------------ */
      /* Bi' mattrix is independent of z-coordinate */
      /* Bi' mattrix is estimated first here        */
      /* Bi'  = []5x3                               */
      /* ------------------------------------------ */

          B_matrix[0][k+1] = shp[0][j-1]; 
          B_matrix[0][k+2] = B_matrix[0][k+3] = 0.0;

          B_matrix[1][k+2] = shp[1][j-1]; 
          B_matrix[1][k+1] = B_matrix[1][k+3] = 0.0;

          B_matrix[2][k+1] = shp[1][j-1]; 
          B_matrix[2][k+2] = shp[0][j-1]; 
          B_matrix[2][k+3] = 0.0;

          B_matrix[3][k+3] = shp[1][j-1]; 
          B_matrix[3][k+1] = B_matrix[3][k+2] = 0.0;

          B_matrix[4][k+3] = shp[0][j-1]; 
          B_matrix[4][k+1] = B_matrix[4][k+2] = 0.0;

          h  = p->work_section[11].value;    /* thickness of the shell */

      /* ----------------------------------*/
      /*  Calculate Bi" matrix             */
      /*  Bi" = [ ] 5x2                    */
      /* ----------------------------------*/

          B_matrix[0][k+4] =  0.0;
          B_matrix[0][k+5] =  shp[0][j-1]*h*0.5*z_coord; 

          B_matrix[1][k+4] = -shp[1][j-1]*h*0.5*z_coord; 
          B_matrix[1][k+5] =  0.0;

          B_matrix[2][k+4] = -shp[0][j-1]*h*0.5*z_coord; 
          B_matrix[2][k+5] =  shp[1][j-1]*h*0.5*z_coord; 

          B_matrix[3][k+4] = -shp[2][j-1];
          B_matrix[3][k+5] =  0.0;

          B_matrix[4][k+4] =  0.0;
          B_matrix[4][k+5] =  shp[2][j-1];
      }

#ifdef DEBUG 
     printf(" Leaving B_MATRIX_4Node() \n");
#endif

      return(B_matrix);
}


void Load_Curve(p, H, stress, E, fy)
ARRAY         *p;
double    stress;
double        *H;
double     E, fy;
{
double temp1, temp2;
    
#ifdef DEBUG
    printf(" enter Loac_Curve() \n");
#endif
    temp1 = stress/fy;
    temp2 = p->LC_ptr->n -1.0;

    *H = E/(p->LC_ptr->n*p->LC_ptr->alpha*pow(temp1, temp2));

#ifdef DEBUG
    printf(" Inside Load_Curve(): H = %le\n", *H);
    printf(" Leaving Load_Curve() \n");
#endif

 /* =================================================*/
 /* for Ramberg-Osgood stress-strain relations       */
 /* : strain/strain0 =                               */
 /*    stress/stress0+alpha*(stress/stress0)^n       */
 /*  pl_strain = strain0a*alpha*(stress/stress0)^n   */
 /* =================================================*/
 /* H = d(stress)/d(pl_strain) =                     */
 /*    = E/(alpha*n)*(stress/stress0)^(1-n)          */
 /* =================================================*/
}


void BB_Vector(co_coord, BB1, BB2) 
double **co_coord;
double *BB1, *BB2;
{
double co_x31, co_y31, co_z31;
double co_x42, co_y42, co_z42;
int i, j, k;

#ifdef DEBUG
      printf(" enter BB_Vector() \n");
      dMatrixPrint("co_coord", co_coord, 3, 4);
#endif

      /* BB_Vector */

      co_x31 =  co_coord[0][2] - co_coord[0][0];
      co_y31 =  co_coord[1][2] - co_coord[1][0];
      co_z31 =  co_coord[2][2] - co_coord[2][0];

      co_x42 =  co_coord[0][3] - co_coord[0][1];
      co_y42 =  co_coord[1][3] - co_coord[1][1];

      BB1[0] = -co_y42/4.0;
      BB1[1] =  co_y31/4.0;
      BB1[2] =  co_y42/4.0;
      BB1[3] = -co_y31/4.0;

      BB2[0] =  co_x42/4.0;
      BB2[1] = -co_x31/4.0;
      BB2[2] = -co_x42/4.0;
      BB2[3] =  co_x31/4.0;

 /* in Belyschenko's paper                    */
 /* B1I = (1/2) [...]  instead of (1/4) [...] */
 /* B2I = (1/2) [...]  instead of (1/4) [...] */
    

#ifdef DEBUG
  for(i = 1; i <= 4; i++)
      printf(" BB1[%d] = %lf, BB2[%d] = %lf \n", i, BB1[i-1], i, BB2[i-1]);
      printf(" leaving BB_Vector() \n");
#endif
  
}


double **Shell_Nodal_Load_Plane(nodal_load, p, co_coord, z_coord, z_integ_pt, EE, nu, task)
double                 **nodal_load;
ARRAY                            *p;
double                   **co_coord;
double                      z_coord;
int                      z_integ_pt;      /* integration point in z-direction */
double                       EE, nu;
int                            task;
{
int           i, j, k,n, ii, jj, kk;
int                       dof, size;
int                    UNITS_SWITCH;
int       surface_pts, no_integ_pts;
double               *x_integ_coord;
double               *y_integ_coord;
double         weight[16], jacobian;
double                   sum, **shp;
double                   **B_matrix;
double                **B_Transpose;
double                     **stress; 
double            **nodal_load_temp;


#ifdef DEBUG 
     printf(" Enter Shell_Nodal_Load_Plane() \n");
#endif

    UNITS_SWITCH = CheckUnits();

    shp            = MatrixAllocIndirectDouble(3,4);
    x_integ_coord  = dVectorAlloc(16);
    y_integ_coord  = dVectorAlloc(16);

    B_matrix         = MatrixAllocIndirectDouble(p->dof_per_node, p->size_of_stiff);
    B_Transpose      = MatrixAllocIndirectDouble(p->size_of_stiff, p->dof_per_node);
    nodal_load_temp  = MatrixAllocIndirectDouble(p->size_of_stiff, 1);
    stress           = MatrixAllocIndirectDouble(p->dof_per_node, 1);

    size           = p->size_of_stiff;
    dof            = p->dof_per_node;
    surface_pts    = p->integ_ptr->surface_pts;
    no_integ_pts   = 0;
    pgauss(surface_pts, &no_integ_pts, x_integ_coord, y_integ_coord, weight);
   
    for(i = 1; i <= size; i++) {
      nodal_load[i-1][0]      = 0.0;
      nodal_load_temp[i-1][0] = 0.0;
    }

       if(task == STRESS && z_integ_pt == 1){ /* print elmement stress */
          if(p->elmt_no == 1)
             printf(" Element : %s \n Material : %s \n\n", p->elmt_type, p->material_name);

          printf("\n STRESS in  Element No  %d \n",p->elmt_no);
          printf(" =============================================================================================================== \n");
          printf(" Gaussion    xi       eta        gamma    Stre-xx         Stre-yy         Stre-xy         Stre-yz        Stre-zx \n");
          if(UNITS_SWITCH == OFF)
             printf("  Points \n");
       }


    for( ii = 1; ii <= no_integ_pts; ii++) {

       elmt_shell_shape_4node(p, co_coord, x_integ_coord[ii-1],y_integ_coord[ii-1],shp,&jacobian,STIFF);

       B_matrix = B_MATRIX_4Node(B_matrix, p, shp, z_coord);

       for(i = 1; i <= p->dof_per_node; i++) {
           for(j = 1; j <= p->size_of_stiff; j++) {
               B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];
           }
       }
       /* update the stress in array pointer */

       kk = no_integ_pts*(z_integ_pt-1)+ii;

#ifdef DEBUG
       printf("jj = %d, kk = %d z_integ_pt = %d, no_integ_pts = %d \n",
               ii, kk, z_integ_pt, no_integ_pts);
       printf(" In Shell_Nodal_Load_Plane(): begin Stress_Update(): \n");
#endif
       Stress_Update(p, co_coord, B_matrix, kk);

#ifdef DEBUG
       printf(" end of Stress_Update() \n");
#endif

       /* IF THE TASK IS TO DO STRESS UPDATE or PRINT STRESS, */
       /* STOP HERE AND CONTINUE FOR THE NEXT LOOP            */
      
       if(task == STRESS_UPDATE)
          goto STRESS_UPDATE_END; 

       if(task == STRESS){ /* print elmement stress */
          if(UNITS_SWITCH == ON) {
             if(kk == 1) {
                printf("  Points                                    %s             %s             %s             %s            %s \n",
                     p->stress->spRowUnits[0].units_name,
                     p->stress->spRowUnits[1].units_name,
                     p->stress->spRowUnits[2].units_name,
                     p->stress->spRowUnits[3].units_name,
                     p->stress->spRowUnits[4].units_name);
                printf(" ---------------------------------------------------------------------------------------------------------------\n \n");
             }
             printf("   %d  %10.4f %10.4f %10.4f", kk, x_integ_coord[ii-1], y_integ_coord[ii-1], z_coord);
             printf("\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",
                    p->stress->uMatrix.daa[0][kk-1]/p->stress->spRowUnits[0].scale_factor,
                    p->stress->uMatrix.daa[1][kk-1]/p->stress->spRowUnits[1].scale_factor,
                    p->stress->uMatrix.daa[2][kk-1]/p->stress->spRowUnits[2].scale_factor,
                    p->stress->uMatrix.daa[3][kk-1]/p->stress->spRowUnits[3].scale_factor,
                    p->stress->uMatrix.daa[4][kk-1]/p->stress->spRowUnits[4].scale_factor);
          }
          if(UNITS_SWITCH == OFF) {
             if(kk == 1) 
                printf(" ---------------------------------------------------------------------------------------------------------------\n \n");
             printf("   %d  %10.4f %10.4f %10.4f", kk, x_integ_coord[ii-1], y_integ_coord[ii-1], z_coord);
             printf("\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",
                    p->stress->uMatrix.daa[0][kk-1],
                    p->stress->uMatrix.daa[1][kk-1],
                    p->stress->uMatrix.daa[2][kk-1],
                    p->stress->uMatrix.daa[3][kk-1],
                    p->stress->uMatrix.daa[4][kk-1]);
          }
          goto STRESS_UPDATE_END; 
       }

       /* calculate the element nodal_load   */

        for( i = 1; i <= p->dof_per_node; i++) {
            stress[i-1][0] = p->stress->uMatrix.daa[i-1][kk-1];
        }
        nodal_load_temp = dMatrixMultRep(nodal_load_temp, B_Transpose, size, dof, stress, dof, 1);
        

        for(i = 1; i <= p->size_of_stiff; i++) {
            nodal_load[i-1][0] += nodal_load_temp[i-1][0]*jacobian*weight[ii-1]; 
        } 


STRESS_UPDATE_END: 
        ;

    }  /* end of gaussian integration */

    free((char *) x_integ_coord);
    free((char *) y_integ_coord);
    MatrixFreeIndirectDouble(shp, 3);

    MatrixFreeIndirectDouble(nodal_load_temp, p->size_of_stiff);
    MatrixFreeIndirectDouble(B_matrix, p->dof_per_node);
    MatrixFreeIndirectDouble(B_Transpose, p->size_of_stiff);
    MatrixFreeIndirectDouble(stress, p->dof_per_node);

#ifdef DEBUG 
     dMatrixPrint("nodal_load surface ", nodal_load, p->size_of_stiff, 1);
     printf(" Leaving Shell_Nodal_Load_Plane() \n");
#endif

    return(nodal_load);
}


/*
#define DEBUG 
*/

/* ======================*/
/* Shell Stiff Matrix    */
/* ======================*/

#ifdef __STDC__
void Shell_Stiff_Plane_4node(double **K, ARRAY *p, double **co_coord, double z_coord, int z_integ_pt, double EE, double nu)
#else
void Shell_Stiff_Plane_4node(K, p, co_coord, z_coord, z_integ_pt, EE, nu)
double         **K;
ARRAY           *p;
double  **co_coord;
double     z_coord;
int     z_integ_pt;      /* integration point in z-direction */
double      EE, nu;
#endif
{
int                i, j, k,n, ii, jj, kk;
static int     surface_pts, no_integ_pts;
int                            size, dof;
double                    *x_integ_coord;
double                    *y_integ_coord;
static double       weight[16], jacobian;
static double                      **shp;
static double             **stiff_matrix;
static double  **B_Transpose, **B_matrix;
static double              **temp_matrix;

static double          sum, temp1, temp2;

#ifdef DEBUG 
     printf(" Enter Shell_Stiff_Plane_4node() \n");
#endif

    shp            = MatrixAllocIndirectDouble(3,4);
    x_integ_coord  = dVectorAlloc(16);
    y_integ_coord  = dVectorAlloc(16);

    stiff_matrix   = MatrixAllocIndirectDouble(p->size_of_stiff,p->size_of_stiff);
    B_matrix       = MatrixAllocIndirectDouble(p->dof_per_node,p->size_of_stiff);
    B_Transpose    = MatrixAllocIndirectDouble(p->size_of_stiff,p->dof_per_node);
    temp_matrix    = MatrixAllocIndirectDouble(p->dof_per_node,p->size_of_stiff);

    surface_pts    = p->integ_ptr->surface_pts;
    no_integ_pts   = 0;
    size           = p->size_of_stiff;
    dof            = p->dof_per_node;

    if(surface_pts*surface_pts != no_integ_pts)
       pgauss(surface_pts, &no_integ_pts, x_integ_coord, y_integ_coord, weight);

    /* Material Matrix for elastic materials  */

     p->mater_matrix->uMatrix.daa = MATER_MAT_SHELL(p->mater_matrix->uMatrix.daa, EE, nu);   

   /* ==============================*/
   /* initialize stiffness matrix   */
   /* ==============================*/
    for ( i = 1; i <= size; i++) {
         for(j = 1; j <= size; j++) {
             K[i-1][j-1] = 0.0;
         }
    }

    for( ii = 1; ii <= no_integ_pts; ii++) { /* loop over the surface */

        elmt_shell_shape_4node(p, co_coord, x_integ_coord[ii-1],
                        y_integ_coord[ii-1],shp,&jacobian, STIFF);

        B_matrix = B_MATRIX_4Node(B_matrix, p, shp, z_coord);

#ifdef DEBUG
       printf(" jacobian = %lf \n", jacobian);
       dMatrixPrint(" shape func", shp,3,4);
       dMatrixPrint(" B_matrix", B_matrix, p->dof_per_node, p->size_of_stiff);
#endif
       for( i = 1; i <= dof; i++) {
           for(j = 1; j <= size; j++) {
               B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];
           }
       }

       /* ============================================================= */
       /* Update the Material Matrix for Elastic_Plastic Materials      */
       /* ============================================================= */

       kk = no_integ_pts*(z_integ_pt-1) +ii;
       MATER_SHELL_UPDATE(p, co_coord, EE, nu, kk);

       /* calculate the element stiffness matrix */
       temp_matrix  = dMatrixMultRep(temp_matrix,
                      p->mater_matrix->uMatrix.daa, dof, dof, B_matrix, dof, size); 
       stiff_matrix = dMatrixMultRep(stiff_matrix, B_Transpose,
                      size, dof, temp_matrix, dof, size); 

       for ( i = 1; i <= size; i++) {
            for(j = 1; j<= size; j++) {
                K[i-1][j-1] += stiff_matrix[i-1][j-1]*jacobian*weight[ii-1];
            }
       }

    }  /* end of gaussian integration */


#ifdef DEBUG
       printf(" end of integration over surface \n"); 
#endif

    free((char *) x_integ_coord);
    free((char *) y_integ_coord);
    MatrixFreeIndirectDouble(shp, 3);
    MatrixFreeIndirectDouble(B_matrix, p->dof_per_node);
    MatrixFreeIndirectDouble(B_Transpose, p->size_of_stiff);
    MatrixFreeIndirectDouble(temp_matrix, p->dof_per_node);
    MatrixFreeIndirectDouble(stiff_matrix, p->size_of_stiff);

#ifdef DEBUG 
     dMatrixPrint("stiff surface ", K, p->size_of_stiff, p->size_of_stiff);
     printf(" Leaving Shell_Stiff_Plane_4node() \n");
#endif

}



#define MASS    1

/* =================== */
/* Shell Mass Matrix    */
/* =================== */
#ifdef __STDC__
void  Shell_4Node_Mass(ARRAY *p, MATRIX *mass, double **co_coord, double density, double thickness)
#else
void  Shell_4Node_Mass(p, mass, co_coord, density, thickness)
ARRAY               *p;
MATRIX           *mass;
double      **co_coord;
double         density;
double       thickness;
#endif
{
int  i, j, k, n, ii, length, length1, length2;
int    no_integ_pts, x_integ_pts, y_integ_pts;

double             **shp, **Mass, **Mass_temp;
double         **N_Transpose, **N, **N1, **N2;

double         *x_integ_coord, *y_integ_coord;
double                   weight[16], jacobian;
double                                    sum;

double               Aera, elmt_length, width;
double            node_mass, node_Jx, node_Jy;

DIMENSIONS                      *d1, *d2, *d3;
int                   UNITS_SWITCH, UnitsType;

#ifdef DEBUG
       printf("*** Enter Shell_4Node_Mass(): \n");
#endif
     /* [a] mass initialization */

      for(i = 1; i <= p->size_of_stiff; i++) 
          for(j = 1; j <= p->size_of_stiff; j++)
                mass->uMatrix.daa[i-1][j-1] = 0.0;

    /* [b] :  mass matrix      */

#ifdef DEBUG
       printf("*** In Shell_4Node_Mass(): in main swicth() \n");
       printf("*** mtype = %d \n", p->type);
       printf("*** density   = %lf\n", density);
       printf("*** thickness = %lf\n", thickness);
#endif

    switch(p->type) {
	case LUMPED:
              
             for (i = 1; i <= p->nodes_per_elmt; i++) { 
                  elmt_length = ABS(p->coord[0][i-1].value - p->coord[0][i].value);
                  if(elmt_length != 0)
                     break;
             }
           
             width       =  p->work_section[7].value;
             Aera        =  p->work_section[10].value;
             node_mass   =  density*Aera*elmt_length/4.0;
             node_Jx     =  density*elmt_length*thickness*width*width*width/12.0/4.0;
             node_Jy     =  density*Aera*elmt_length*elmt_length*elmt_length/12.0/4.0;
#ifdef DEBUG
            printf(" width      = %lf\n", width);
            printf(" length     = %lf\n", elmt_length);
            printf(" Aera       = %lf\n", Aera);
            printf(" node_mass  = %lf\n", node_mass);
            printf(" node_Jx    = %lf\n", node_Jx);
            printf(" node_Jy    = %lf\n", node_Jy);
#endif
             
             for(i = 1; i <= p->dof_per_node; i++) {
                 for (j = 1; j <= p->nodes_per_elmt; j++) {
                      k = p->dof_per_node*(j-1) + i; 
                      if(i <= 3)
                         mass->uMatrix.daa[k-1][k-1] = node_mass;
                      else{
                         if(i == 4) 
                         mass->uMatrix.daa[k-1][k-1] = node_Jx;
                         if(i == 5) 
                         mass->uMatrix.daa[k-1][k-1] = node_Jy;
                      }
                 }
             }
             
	     break;
	case CONSISTENT:

        /* MASS : [M] = [M1]+[M3] */

             Mass = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff); 
             N              = MatrixAllocIndirectDouble(p->dof_per_node, p->size_of_stiff);
             N_Transpose    = MatrixAllocIndirectDouble(p->size_of_stiff, p->dof_per_node);
             N1             = MatrixAllocIndirectDouble(p->dof_per_node, p->size_of_stiff);
             N2             = MatrixAllocIndirectDouble(p->dof_per_node, p->size_of_stiff);
             shp            = MatrixAllocIndirectDouble(3,4);
             x_integ_coord  = dVectorAlloc(16);
             y_integ_coord  = dVectorAlloc(16);
             x_integ_pts    = p->integ_ptr->surface_pts;
             y_integ_pts    = p->integ_ptr->surface_pts;
             no_integ_pts   = 0;

             pgauss(x_integ_pts, &no_integ_pts, x_integ_coord, y_integ_coord, weight);


             for( ii = 1; ii <= no_integ_pts; ii++) {

                 elmt_shell_shape_4node(p, co_coord, x_integ_coord[ii-1],y_integ_coord[ii-1],shp,&jacobian,MASS);

                for(i = 1; i<= p->dof_per_node; i++) {
                    for(j = 1; j <= p->size_of_stiff; j++) {
                        n = (j-1)/p->dof_per_node;
                        k = j - p->dof_per_node*n;
                       if(i <= 3 || k <= 3){
                         if( i == k ) { 
                            N[i-1][j-1] = shp[2][n]; 
                         }
                         else
                            N[i-1][j-1] = 0.0;
                       }
                       else
                         N[i-1][j-1] = 0.0;
                    }
                }
            
                /* [M1] calculation  */

                for(i = 1; i <= p->dof_per_node; i++) {
                  for(j = 1; j <= p->size_of_stiff; j++) {
                     N_Transpose[j-1][i-1] = N[i-1][j-1];
                  }
                }

            Mass = dMatrixMultRep(Mass,N_Transpose,p->size_of_stiff,p->dof_per_node,N,p->dof_per_node, p->size_of_stiff);

#ifdef DEBUG
          dMatrixPrint("M1", Mass, p->size_of_stiff, p->size_of_stiff);
#endif
                for(i = 1 ; i <= p->size_of_stiff; i++) {
                    for(j = 1 ; j <= p->size_of_stiff; j++) {
                        mass->uMatrix.daa[i-1][j-1] += Mass[i-1][j-1]*jacobian*weight[ii-1]*thickness*density;
                    }
                }

                /* [M3] calculation */

                for(i = 1; i<= p->dof_per_node; i++) {
                    for(j = 1; j <= p->size_of_stiff; j++) {
                        n = (j-1)/p->dof_per_node;
                        k = j - p->dof_per_node*n;
                        if( i > 3 && k > 3 && i == k)
                            N2[i-1][j-1] = shp[2][n];
                        else
                            N2[i-1][j-1] = 0.0;
                    }
                }

                for(i = 1; i <= p->dof_per_node; i++) {
                  for(j = 1; j <= p->size_of_stiff; j++) {
                     N_Transpose[j-1][i-1] = N2[i-1][j-1];
                  }
                }

                Mass  = dMatrixMultRep(Mass,N_Transpose,p->size_of_stiff,p->dof_per_node,N2,p->dof_per_node, p->size_of_stiff);
          
#ifdef DEBUG
          dMatrixPrint("M3", Mass, p->size_of_stiff, p->size_of_stiff);
#endif

                for(i = 1 ; i <= p->size_of_stiff; i++) {
                    for(j = 1 ; j <= p->size_of_stiff; j++) {
                        mass->uMatrix.daa[i-1][j-1] += Mass[i-1][j-1]*jacobian*weight[ii-1]*thickness*thickness*thickness*density/12.0;
                    }
                }
             }
             MatrixFreeIndirectDouble(N1, p->dof_per_node); 
             MatrixFreeIndirectDouble(N2, p->dof_per_node); 
             MatrixFreeIndirectDouble(N_Transpose, p->size_of_stiff);
             MatrixFreeIndirectDouble(N, p->dof_per_node);
             MatrixFreeIndirectDouble(shp,3);
             MatrixFreeIndirectDouble(Mass, p->size_of_stiff);
             free((char *) x_integ_coord);
             free((char *) y_integ_coord);
	     break;
	default:
             FatalError("In Shell_4Node_Mass() : Type of Mass Matrix Undefined",(char *)NULL);
	     break;
    }
       
 /* ------------ MASS UNITS ---------------------------- */
 /* The units type is determined by the SetUnitsType()   */
 /* ---------------------------------------------------- */

 /* Initiation of Mass Units Buffer                      */

      UNITS_SWITCH = CheckUnits();
      UnitsType    = CheckUnitsType();
      switch( UNITS_SWITCH ) {
        case ON:
           if(UnitsType == SI || UnitsType == SI_US ) {
              d1 = DefaultUnits("Pa");
              d2 = DefaultUnits("m");
           }
           else {
              d1 = DefaultUnits("psi");
              d2 = DefaultUnits("in");
           }
           d3 = DefaultUnits("sec");

          /* node no 1 */
           UnitsMultRep( &(mass->spColUnits[0]), d1, d2 );
           UnitsCopy( &(mass->spColUnits[1]), &(mass->spColUnits[0]) );
           UnitsCopy( &(mass->spColUnits[2]), &(mass->spColUnits[0]) );
           UnitsMultRep( &(mass->spColUnits[3]), &(mass->spColUnits[0]), d2 );
           UnitsCopy( &(mass->spColUnits[4]), &(mass->spColUnits[3]) );

           UnitsPowerRep( &(mass->spRowUnits[0]), d3, 2.0, NO );
           UnitsCopy( &(mass->spRowUnits[1]), &(mass->spRowUnits[0]) );
           UnitsCopy( &(mass->spRowUnits[2]), &(mass->spRowUnits[0]) );
           UnitsMultRep( &(mass->spRowUnits[3]), d2, &(mass->spRowUnits[0]) );
           UnitsCopy( &(mass->spRowUnits[4]), &(mass->spRowUnits[3]) );

          /* node no > 1 */
           for(i = 2; i <= p->nodes_per_elmt; i++) {    
               for(j = 1; j <= p->dof_per_node; j++) {
                   k = p->dof_per_node*(i-1)+j;
                   if(j <= 3) {
                      UnitsCopy( &(mass->spColUnits[k-1]), &(mass->spColUnits[0]) );
                      UnitsCopy( &(mass->spRowUnits[k-1]), &(mass->spRowUnits[0]) );
                   }
                   else {
                      UnitsCopy( &(mass->spColUnits[k-1]), &(mass->spColUnits[3]) );
                      UnitsCopy( &(mass->spRowUnits[k-1]), &(mass->spRowUnits[3]) ); 
                   }
               }
           }
           free((char *) d1->units_name);
           free((char *) d1);
           free((char *) d2->units_name);
           free((char *) d2);
           free((char *) d3->units_name);
           free((char *) d3);
           break;
        case OFF:
           break;
        default:
           break;
     }

#ifdef DEBUG
       MatrixPrintIndirectDouble(mass);
       printf("*** leaving Shell_4Node_Mass()  \n");
#endif
}


/* ================================================== */
/* Shell_Element Shape Functions for 4-Node Element   */
/* ================================================== */

#ifdef __STDC__
void elmt_shell_shape_4node(ARRAY *p, double **coord, double ss, double tt, double **shp, double *jacobian, int MASS_FLAG)
#else
void elmt_shell_shape_4node(p, coord, ss,tt, shp,jacobian,MASS_FLAG)
ARRAY                             *p;
double                **coord, **shp;
double             ss, tt, *jacobian;
int                        MASS_FLAG;
#endif
{
int          i, j, k;
double        *s, *t;
double    **xs, **sx;
double    **shp_temp;

#ifdef DEBUG
   printf(" Enter elmt_shell_shape_4node() \n");
#endif

    s        = dVectorAlloc(4);
    t        = dVectorAlloc(4);
    xs       = MatrixAllocIndirectDouble(2,2);
    sx       = MatrixAllocIndirectDouble(2,2);
    shp_temp = MatrixAllocIndirectDouble(2,4);

    s[0] = -0.5; s[1] =  0.5;
    s[2] =  0.5; s[3] = -0.5;

    t[0] = -0.5; t[1] = -0.5;
    t[2] =  0.5; t[3] =  0.5;

  switch(p->nodes_per_elmt) { 
    case 4:

    /* form 4-node quadrilateral shape function                    */
    /* shape function:                  Ni = shape[2][i-1]         */
    /*                                  node no. i = 1, .. 4       */
    /* derivatives of shape functions:  dNi/d(ss) = shape[0][i-1]  */
    /*                                  dNi/d(tt) = shape[1][i-1]  */

    for(i = 1; i <= 4; i++){
        shp[2][i-1]      = (0.5 + s[i-1] * ss) * ( 0.5 + t[i-1] * tt); 
        shp_temp[0][i-1] = s[i-1] * (0.5 + t[i-1] * tt);                
        shp_temp[1][i-1] = t[i-1] * (0.5 + s[i-1] * ss);
    }

     /* construct jacobian matrix and its determinant */

     for(i = 1; i <=  2; i++)      /* in-plane dimen = 2 */
         for(j = 1; j <=  2; j++) {
             xs[i-1][j-1] = 0.0;
             for(k = 1; k <= p->nodes_per_elmt; k++)
                 xs[i-1][j-1] = xs[i-1][j-1] + coord[i-1][k-1]*shp_temp[j-1][k-1];
         }
     *jacobian = xs[0][0] * xs[1][1] - xs[0][1] *xs[1][0];

     if(MASS_FLAG == MASS) {
        break;
     }	

    /* compute Jacobain inverse matrix */

    sx[0][0] = xs[1][1]/ *jacobian;
    sx[1][1] = xs[0][0]/ *jacobian;
    sx[0][1] = - xs[0][1]/ *jacobian;
    sx[1][0] = - xs[1][0]/ *jacobian;


    /* form global derivatives */

    /* save dNi/dx, dNi/dy into shp[0][node] , shp[1][node]*/

    for(i = 1; i <= p->nodes_per_elmt; i++){
        shp[0][i-1] = shp_temp[0][i-1]*sx[0][0] + shp_temp[1][i-1]*sx[0][1];
        shp[1][i-1] = shp_temp[0][i-1]*sx[1][0] + shp_temp[1][i-1]*sx[1][1];
    
    }
    break;

    default:
    break;
 }
#ifdef DEBUG
    dMatrixPrint(" shp", shp, 3,4);
#endif

  free((char *) s);
  free((char *) t);
  MatrixFreeIndirectDouble(xs, 2);
  MatrixFreeIndirectDouble(sx, 2);
  MatrixFreeIndirectDouble(shp_temp, 2);
#ifdef DEBUG
   printf(" leaving elmt_shell_shape_4node() \n");
#endif

}


/* =======================================================*/
/* Calculation of Direction Matrix; Transform coordinate  */
/* into lamina coordiante system                          */
/* =======================================================*/

#ifdef __STDC__
void Lamina_Sys_Implicit(ARRAY *p, double **Direction_Matrix, double **co_coord)
#else
void Lamina_Sys_Implicit(p, Direction_Matrix, co_coord)
ARRAY                         *p;
double        **Direction_Matrix;
double                **co_coord;
#endif
{
double dof, size, temp;

double x21, x31, x42;
double y21, y31, y42;
double z21, z31, z42;

double co_x31, co_y31, co_z31;
double co_x42, co_y42, co_z42;

double co_vel_x31, co_vel_y31, co_vel_z31;
double co_vel_x42, co_vel_y42, co_vel_z42;

double co_rot_x31, co_rot_y31, co_rot_z31;
double co_rot_x42, co_rot_y42, co_rot_z42;

double Jac;          /* Determinant of Jacobian matrix */
double s_norm;

double **r21_ptr, **r31_ptr, **r42_ptr;
double **s_ptr, **e1_ptr, **e2_ptr, **e3_ptr;
double **A_matrix, **coord;

int i, j, k, ii, jj, kk, n, n1, n2, k1, k2;
    
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

    coord = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);


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
    
        MatrixFreeIndirectDouble(coord, 3);
        MatrixFreeIndirectDouble(A_matrix, 3);

#ifdef DEBUG
        dMatrixPrint("co_coord", co_coord, 3,4);
        printf(" Leaving Lamina_Sys_Implicity()\n");
#endif

}
