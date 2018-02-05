/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_shell_8n.c : Eight Node Shell Finite Element
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
 *  Written by: Xiaoguang Chen                                      December 1995
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

ARRAY *sld108( ARRAY *, int );


/* ============================================================== */
/*   Element SHELL_EIGHT_NODES                                    */
/*        Implicit code                                           */
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

ARRAY *elmt_shell_8n(p, isw)
ARRAY    *p;
int     isw;
{
#ifdef DEBUG
       printf(" enter elmt_shell_8n() \n");
#endif

   p = elmt_shell_8nodes_implicit(p, isw);

#ifdef DEBUG
       printf(" Leaving elmt_shell_8n() \n");
#endif

   return(p);
}

ARRAY *elmt_shell_8nodes_implicit(p, isw)
ARRAY    *p;
int     isw;
{
static double                                   nu; 
static QUANTITY      fy, G, E, ET, density, lambda; 
static double         Ixx, Iyy, Izz, Ixy, Ixz, Iyz;
static double        bf, tf, A, depth, h, EA, EIzz;

static DIMENSIONS *dp_length, *dp_force, *dp_moment;
static DIMENSIONS *dp_stress, *dp_degree, *dp_temperature;

double             EE, jacobian, k_bar = 0.833;     /* k_bar is shear factor */
double                   elmt_length, aspect_ratio;

double                           **co_coord = NULL;     
double                          **Direction_Matrix;
double                   **T_matrix, **T_Transpose;
double                             **e1_ptr = NULL;
double                             **e2_ptr = NULL;
double                             **e3_ptr = NULL;

static double      *integ_coord = NULL, *weight = NULL;

double                 **stress = NULL, **displ = NULL;
double   **nodal_load = NULL, **nodal_load_temp = NULL;
double                                 diff, sum, temp;
double     **Cep = NULL, **stiff = NULL, **mass = NULL; 
int             i, j, k, ii, jj, kk, n, n1, n2, k1, k2;
int          surface_pts, surf_integ_pts, no_integ_pts;
int                size, dof, length, length1, length2;
int                            UNITS_SWITCH, UnitsType;


#ifdef DEBUG
       printf("*** Enter elmt_shell_8nodes_implicit() : isw = %4d\n", isw);
#endif

     
       UNITS_SWITCH = CheckUnits();
       UnitsType    = CheckUnitsType();

       co_coord         = MatrixAllocIndirectDouble(p->no_dimen, p->nodes_per_elmt);

       for( i = 1; i <= 3; i++ ) {
           for ( j = 1; j <= p->nodes_per_elmt; j++) {
               co_coord[i-1][j-1] = p->coord[i-1][j-1].value;
           }
       }
#ifdef DEBUG
        dMatrixPrint("co_coord", co_coord, 3,8);
#endif

      no_integ_pts = p->integ_ptr->thickness_pts;
      dof          = p->dof_per_node;
      size         = p->size_of_stiff;

/*****************************************************************/
#ifdef DEBUG
    printf(" ***** In elmt_shell_8nodes_implicit(): \n");
    printf("                                      : enter main switch \n");
#endif

    switch(isw) {
        case PROPTY: 

#ifdef DEBUG
    printf(" In elmt_shell_8nodes_implicit(): \n");
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
    printf(" In elmt_shell_8nodes_implicit(): \n");
    printf("    : leaving case of PROPTY\n");
#endif
             break;
             case STIFF: /* form element stiffness */

#ifdef DEBUG
       printf("*** In elmt_shell() : start case STIFF\n");
       printf("                    : Density         = %14.4e\n", density.value);
       printf("                    : shell_thickness = %8.2f\n", h);
       printf("                    : E               = %lf\n", E.value);
       printf("                    : nu              = %lf\n", nu);
       printf("                    : no_integ_pts in z_direction = %d\n", no_integ_pts);
#endif
              /* ============================================*/
              /* initialize the stiffness matrix             */
              /* ============================================*/
                 for(i = 1; i <= p->size_of_stiff; i++) {
                    for(j = 1; j <= p->size_of_stiff; j++) {
                        p->stiff->uMatrix.daa[i-1][j-1] = 0.0;
                    }
                 }

                 stiff        = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
                 integ_coord  = dVectorAlloc(no_integ_pts + 1);
                 weight       = dVectorAlloc(no_integ_pts + 1);

                 size         = p->size_of_stiff;
                 gauss(integ_coord,weight,no_integ_pts);

                 /* Integration loop over thickness */

                 for (ii = 1; ii <= no_integ_pts; ii++) {

                      Shell_Stiff_Plane_8node(stiff, p, co_coord, integ_coord[ii], ii, E.value, nu); 	

                      for (i = 1; i <= size; i++) {
                          for (j = 1; j <= size; j++) {
                               p->stiff->uMatrix.daa[i-1][j-1] += stiff[i-1][j-1]*weight[ii];     
                          }
                      }
                 }  /* end of gaussian integration through thickness */

                 free((char *) integ_coord);
                 free((char *) weight);

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

          /* node no 1 */
           UnitsMultRep( &(p->stiff->spColUnits[0]), dp_stress, dp_length );
           UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );
           UnitsCopy( &(p->stiff->spColUnits[2]), &(p->stiff->spColUnits[0]) );
           UnitsMultRep( &(p->stiff->spColUnits[3]), &(p->stiff->spColUnits[0]), dp_length );
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
          p = sld108(p, PRESSLD); /* Equivalent Load in local_coordinate */
       }

/* ------------------------ UNITS -----------------------------------*/

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
       printf("****** In elmt_shell_8nodes_implicit() : \n");
       printf("       enter cases: STRESS_UPDATE, LOAD_MATRIX, STRESS_LOAD \n");     
       printf(" elemt_state = %d \n", p->elmt_state);
#endif
          nodal_load      = MatrixAllocIndirectDouble(p->size_of_stiff, 1); 
          nodal_load_temp = MatrixAllocIndirectDouble(p->size_of_stiff, 1); 

          h  = p->work_section[11].value;    /* thickness of the shell */

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
              = Shell_Nodal_Load_8Node(nodal_load_temp, p, co_coord,
                             integ_coord[ii], ii, E.value, nu, isw);
              for(i = 1; i<= p->size_of_stiff; i++) {
                  nodal_load[i-1][0] += nodal_load_temp[i-1][0]*weight[ii];  
              }

          } /* gaussian integration ends */

          free((char *) integ_coord);
          free((char *) weight);

      size = p->size_of_stiff;	

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
             e1_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);
             e2_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);
             e3_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);

             Lamina_Sys_8node(p, e1_ptr, e2_ptr, e3_ptr);

    Direction_Matrix = MatrixAllocIndirectDouble(6,6);
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

    T_matrix     = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
    T_Transpose  = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);

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

      MatrixFreeIndirectDouble(e1_ptr,3);
      MatrixFreeIndirectDouble(e2_ptr,3);
      MatrixFreeIndirectDouble(e3_ptr,3);
      MatrixFreeIndirectDouble(Direction_Matrix,6);
      MatrixFreeIndirectDouble(T_matrix,p->size_of_stiff);
      MatrixFreeIndirectDouble(T_Transpose,p->size_of_stiff);
    }

      MatrixFreeIndirectDouble(nodal_load,p->size_of_stiff);
      MatrixFreeIndirectDouble(nodal_load_temp,p->size_of_stiff);


   /* ------------NODAL LOAD UNITS ------------------------*/

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
           free((char *) dp_force->units_name);
           free((char *) dp_force);
           free((char *) dp_length->units_name);
           free((char *) dp_length);

         break;
      case OFF:
         break;
      default:
         break;
    }

#ifdef DEBUG
printf("*** In elmt_shell_8nodes_implicit() : end case LOAD_MATRIX, STRESS_UPDATE, STRESS_LOAD \n");
#endif
	     break;

        case MASS_MATRIX:  /* form mass matrix */

#ifdef DEBUG
       printf("*** In elmt_shell_8nodes_implicit() : start case MASS\n");
       printf("                : Density = %14.4e\n", density.value);
       printf("                : shell_thickness = %8.2f\n", h);
#endif
   /*========================== */
   /*  CALCULATING MASS MATRIX  */
   /*========================== */

      Shell_8Node_Mass(p, p->stiff, co_coord, density.value, h);

#ifdef DEBUG
       printf("*** In elmt_shell() : end case MASS\n");
#endif
             break;
        default:
             break;
    }

    MatrixFreeIndirectDouble(co_coord, p->no_dimen);

#ifdef DEBUG
       printf("*** leaving elmt_shell() \n");
#endif

    return(p);
}


/* ============================================*/
/* Equivalent Load due to distributed pressure */
/* for 8 node shell element                    */
/* ============================================*/

ARRAY *sld108(p,task)
ARRAY *p;
int task;
{
#ifdef DEBUG
      printf("**** enter sld108(): \n");
#endif
       /* add details */

#ifdef DEBUG
      printf("**** leaving sld108(): \n");
#endif
    return(p);
}


void Stress_Update_8Node(p, co_coord, B1_matrix, integ_pt)
ARRAY                   *p;
double          **co_coord;
double          **B1_matrix;
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
     printf("\n enter Stress_Update_8Node(): \n");
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
   stress           = MatrixAllocIndirectDouble(dof,1); /* stress or stress incremental  */
   stress_dev       = MatrixAllocIndirectDouble(dof,1); /* deviatric component of stress */
   stress_incr      = MatrixAllocIndirectDouble(dof,1); /* stress incremental            */
   back_stress_incr = MatrixAllocIndirectDouble(dof,1); /* stress incremental            */

#ifdef DEBUG
      printf(" In Stress_Update_8Node(): materials_name = %s \n",p->material_name);
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

     strain_incr = dMatrixMultRep(strain_incr, B1_matrix, dof,
                                  size, displ_incr, size, 1); 
     stress      = dMatrixMultRep(stress, m1, dof, dof,
                                  strain_incr, dof, 1);
#ifdef DEBUG
     dMatrixPrint("displ_incr", displ_incr, size, 1);
     dMatrixPrint("strain_incr", strain_incr, dof, 1);
     dMatrixPrint("stress", stress, dof, 1);
#endif

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
         printf("**** In Stress_Update_8Node(): elmt_no = %d\n", p->elmt_no);
         printf(" elmt_state = %d: p->elmt_state \n");
         FatalError("*****Undefine element state ()", (char *) NULL);
       break;
     }

     /* Calculate stress increment */

     strain_incr = dMatrixMultRep(strain_incr, B1_matrix, dof,
                                  size, displ_incr, size, 1); 

     stress      = dMatrixMultRep(stress, m1, dof, dof,
                                  strain_incr, dof, 1);

#ifdef DEBUG
     dMatrixPrint("displ_incr", displ_incr, size, 1);
#endif

#ifdef DEBUG
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
      printf(" element no  = %d \n", p->elmt_no);
      for (i = 1; i <= dof; i++){
         printf(" total stress[%d] = %lf \n", i, stress[i-1][0]);
         printf("incrmental stress[%d] = %lf \n", i, stress_incr[i-1][0]);
         printf("previous stress[%d] = %lf \n", i, p->stress->uMatrix.daa[i-1][kk-1]);
      }
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

       /* Step 3 Estimate number of sub-incrementations needed */

          if( ABS(p->LC_ptr->beta) < 1E-10){   /* Only for beta = 0, kinematic hardening case */
        
       /* Step 3.1 Estimate the effective plastic strain increment */

               temp = sqrt(3.0/2.0)*(A-R)/(p->LC_ptr->H[kk-1]+3.0*G);

       /* Step 3.2 Estimate the back stress increment              */
          
          /* Estimate H' */

               if(!strcmp(p->material_name, "ELASTIC_PERFECTLY_PLASTIC")){
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
          } else {
              for(i = 1; i <= dof; i++) 
                  back_stress_incr[i-1][0] = 0.0;
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
         iNo_iter_step = (int) (2.0*temp/R/Beta1) + 1 ;

#ifdef DEBUG1
    printf(" ******Plastic DEFORMATION A = %lf, R = %lf \n", A, R);
    printf(" ******Plastic DEFORMATION dA = %lf, R = %lf \n", temp, R);
    printf(" at elmt_no = %d, integ_pt = %d\n", p->elmt_no, kk);
    printf(" No of sub-Incremental steps = %d \n",iNo_iter_step);
#endif
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

               strain_incr = dMatrixMultRep(strain_incr,B1_matrix, dof, 
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
         printf(" In Stress_Update_8Node(): elmt_no \n", p->elmt_no);
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
      dMatrixPrint("p->stress in Stress_Update_8Node() leaving ",
                    p->stress->uMatrix.daa, p->dof_per_node, 12);
     printf(" Leaving Stress_Update_8Node() \n");
#endif
}

double **B_MATRIX_8Node(B_matrix, p, shp, z_coord, J_inverse)
double        **B_matrix;
ARRAY                 *p;
double             **shp;
double           z_coord;
double       **J_inverse;
{
int       i,j, k, ii, n;
double                h;
static double     **a = NULL;
static double     **b = NULL;
static double     **c = NULL;
static double     **d = NULL;
static double     **e = NULL;
static double     **g = NULL;
static double     **e1_ptr = NULL;
static double     **e2_ptr = NULL;
static double     **e3_ptr = NULL;

#ifdef DEBUG 
     printf(" enter B_MATRIX_8Node() \n");
#endif

       a = MatrixAllocIndirectDouble(p->nodes_per_elmt,1);
       b = MatrixAllocIndirectDouble(p->nodes_per_elmt,1);
       c = MatrixAllocIndirectDouble(p->nodes_per_elmt,1);
       d = MatrixAllocIndirectDouble(p->nodes_per_elmt,1);
       e = MatrixAllocIndirectDouble(p->nodes_per_elmt,1);
       g = MatrixAllocIndirectDouble(p->nodes_per_elmt,1);

       e1_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);
       e2_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);
       e3_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);

       h  = p->work_section[11].value;    /* thickness of the shell */

       for(i = 1; i <= p->nodes_per_elmt; i++) {
           a[i-1][0] = J_inverse[0][0]*shp[0][i-1] + J_inverse[0][1]*shp[1][i-1] ;  
           b[i-1][0] = J_inverse[1][0]*shp[0][i-1] + J_inverse[1][1]*shp[1][i-1] ;  
           c[i-1][0] = J_inverse[2][0]*shp[0][i-1] + J_inverse[2][1]*shp[1][i-1] ;  

           d[i-1][0] = h*0.5*(a[i-1][0]*z_coord + J_inverse[0][2]*shp[2][i-1]);  
           e[i-1][0] = h*0.5*(b[i-1][0]*z_coord + J_inverse[1][2]*shp[2][i-1]);  
           g[i-1][0] = h*0.5*(d[i-1][0]*z_coord + J_inverse[2][2]*shp[2][i-1]);  
       }

#ifdef DEBUG
       for(i = 1; i <= p->nodes_per_elmt; i++) {
           printf(" a[%d] = %lf \n", i, a[i-1][0]);
           printf(" b[%d] = %lf \n", i, b[i-1][0]);
           printf(" c[%d] = %lf \n", i, c[i-1][0]);

           printf(" d[%d] = %lf \n", i, d[i-1][0]);
           printf(" g[%d] = %lf \n", i, g[i-1][0]);
           printf(" e[%d] = %lf \n", i, e[i-1][0]);
       }
#endif
       Lamina_Sys_8node(p, e1_ptr, e2_ptr, e3_ptr);
       
      for(i = 1; i <= p->dof_per_node + 1; i++) 
          for(j = 1; j <= p->size_of_stiff; j++) 
              B_matrix[i-1][j-1] = 0.0;

      /* ------------------------------------------ */
      /* Bi'  = []6x3                               */
      /* ------------------------------------------ */

      for(j = 1; j <= p->nodes_per_elmt; j++) {
          k = p->dof_per_node*(j-1)-1;

          B_matrix[0][k+1] = a[j-1][0];
          B_matrix[0][k+2] = 0.0;
          B_matrix[0][k+3] = 0.0;

          B_matrix[1][k+1] = 0.0;
          B_matrix[1][k+2] = b[j-1][0];
          B_matrix[1][k+3] = 0.0;

          B_matrix[2][k+1] = 0.0;
          B_matrix[2][k+2] = 0.0;
          B_matrix[2][k+3] = c[j-1][0];

          B_matrix[3][k+1] = b[j-1][0]; 
          B_matrix[3][k+2] = a[j-1][0]; 
          B_matrix[3][k+3] = 0.0;

          B_matrix[4][k+1] = 0.0;
          B_matrix[4][k+2] = c[j-1][0];
          B_matrix[4][k+3] = b[j-1][0]; 

          B_matrix[5][k+1] = c[j-1][0];
          B_matrix[5][k+2] = 0.0;
          B_matrix[5][k+3] = a[j-1][0]; 


      /* ----------------------------------*/
      /*  Calculate Bi" matrix             */
      /*  Bi" = [ ] 6x2                    */
      /* ----------------------------------*/

          B_matrix[0][k+4] =  -d[j-1][0]*e2_ptr[0][j-1];
          B_matrix[0][k+5] =   d[j-1][0]*e1_ptr[0][j-1];
          
          B_matrix[1][k+4] =  -e[j-1][0]*e2_ptr[1][j-1];
          B_matrix[1][k+5] =   e[j-1][0]*e1_ptr[1][j-1];

          B_matrix[2][k+4] =  -g[j-1][0]*e2_ptr[2][j-1];
          B_matrix[2][k+5] =   g[j-1][0]*e1_ptr[2][j-1];

          B_matrix[3][k+4] =  -e[j-1][0]*e2_ptr[0][j-1] - d[j-1][0]*e2_ptr[1][j-1];
          B_matrix[3][k+5] =   e[j-1][0]*e1_ptr[0][j-1] + d[j-1][0]*e1_ptr[1][j-1];

          B_matrix[4][k+4] =  -g[j-1][0]*e2_ptr[1][j-1] - e[j-1][0]*e2_ptr[2][j-1];
          B_matrix[4][k+5] =   g[j-1][0]*e1_ptr[1][j-1] + e[j-1][0]*e1_ptr[2][j-1];

          B_matrix[5][k+4] =  -d[j-1][0]*e2_ptr[2][j-1] - g[j-1][0]*e2_ptr[0][j-1];
          B_matrix[5][k+5] =   d[j-1][0]*e1_ptr[2][j-1] + g[j-1][0]*e1_ptr[0][j-1];
      }

      MatrixFreeIndirectDouble(a, p->nodes_per_elmt);
      MatrixFreeIndirectDouble(b, p->nodes_per_elmt);
      MatrixFreeIndirectDouble(c, p->nodes_per_elmt);
      MatrixFreeIndirectDouble(d, p->nodes_per_elmt);
      MatrixFreeIndirectDouble(e, p->nodes_per_elmt);
      MatrixFreeIndirectDouble(g, p->nodes_per_elmt);

      MatrixFreeIndirectDouble(e1_ptr, 3);
      MatrixFreeIndirectDouble(e2_ptr, 3);
      MatrixFreeIndirectDouble(e3_ptr, 3);

#ifdef DEBUG 
     printf(" Leaving B_MATRIX_8Node() \n");
#endif

      return(B_matrix);
}


double **Shell_Nodal_Load_8Node(nodal_load, p, co_coord, z_coord, z_integ_pt, EE, nu, task)
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
double            **B_matrix = NULL;
double           **B1_matrix = NULL;
double            **B1_Trans = NULL;
double           **J_inverse = NULL;
double                  **TT = NULL;
double                     **stress; 
double            **nodal_load_temp;


#ifdef DEBUG 
     printf(" Enter Shell_Nodal_Load_8Node() \n");
#endif

    UNITS_SWITCH = CheckUnits();

    shp            = MatrixAllocIndirectDouble(3,8);
    x_integ_coord  = dVectorAlloc(16);
    y_integ_coord  = dVectorAlloc(16);

    TT               = MatrixAllocIndirectDouble(6, 6);
    J_inverse        = MatrixAllocIndirectDouble(3, 3);
    B_matrix         = MatrixAllocIndirectDouble(p->dof_per_node + 1, p->size_of_stiff);

    B1_matrix        = MatrixAllocIndirectDouble(p->dof_per_node, p->size_of_stiff);
    B1_Trans         = MatrixAllocIndirectDouble(p->size_of_stiff, p->dof_per_node);

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

       elmt_shell_shape_8node(p, co_coord, x_integ_coord[ii-1], y_integ_coord[ii-1],
                              z_coord, shp,&jacobian, J_inverse, TT, STIFF);

       B_matrix = B_MATRIX_8Node(B_matrix, p, shp, z_coord, J_inverse);

       /* Calculate the B_matrix in local coordinate */
       /* eliminate the 3rd row of TT matrix         */

       for(i = 1; i <= 2; i++)
           for(j = 1; j <= p->size_of_stiff; j++) {
               B1_matrix[i-1][j-1] = 0.0;
               B1_Trans[j-1][i-1]  = 0.0;
               for(k = 1; k <= p->dof_per_node + 1; k++)
                   B1_matrix[i-1][j-1] += TT[i-1][k-1]*B_matrix[k-1][j-1];
           }

       for(i = 3; i <= p->dof_per_node; i++)
           for(j = 1; j <= p->size_of_stiff; j++) {
               B1_matrix[i-1][j-1] = 0.0;
               B1_Trans[j-1][i-1]  = 0.0;
               for(k = 1; k <= p->dof_per_node + 1; k++)
                   B1_matrix[i-1][j-1] += TT[i][k-1]*B_matrix[k-1][j-1];

           }

       for( i = 1; i <= p->dof_per_node; i++) 
           for(j = 1; j <= p->size_of_stiff; j++) 
               B1_Trans[j-1][i-1] = B1_matrix[i-1][j-1];

       /* update the stress in array pointer */

       kk = no_integ_pts*(z_integ_pt-1)+ii;

#ifdef DEBUG
       printf("jj = %d, kk = %d z_integ_pt = %d, no_integ_pts = %d \n",
               ii, kk, z_integ_pt, no_integ_pts);
       printf(" In Shell_Nodal_Load_8Node(): begin Stress_Update_8Node(): \n");
#endif
    
       Stress_Update_8Node(p, co_coord, B1_matrix, kk);

#ifdef DEBUG
       printf(" end of Stress_Update_8Node() \n");
#endif

       /* IF THE TASK IS TO DO STRESS UPDATE or PRINT STRESS, */
       /* STOP HERE AND CONTINUE FOR THE NEXT LOOP            */
      
       if(task == STRESS_UPDATE)
          goto STRESS_UPDATE_END; 

       if(task == STRESS){ /* print elmement stress : local */

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

        nodal_load_temp = dMatrixMultRep(nodal_load_temp, B1_Trans, size, dof, stress, dof, 1);
        
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
    MatrixFreeIndirectDouble(B_matrix, p->dof_per_node + 1);
    MatrixFreeIndirectDouble(B1_matrix, p->dof_per_node);
    MatrixFreeIndirectDouble(B1_Trans, p->size_of_stiff);
    MatrixFreeIndirectDouble(stress, p->dof_per_node);

#ifdef DEBUG 
     dMatrixPrint("nodal_load surface ", nodal_load, p->size_of_stiff, 1);
     printf(" Leaving Shell_Nodal_Load_8Node() \n");
#endif

    return(nodal_load);
}


#ifdef __STDC__
void Lamina_Sys_8node(ARRAY *p, double **e1_ptr, double **e2_ptr, double **e3_ptr)
#else
void Lamina_Sys_8node(p, e1_ptr, e2_ptr, e3_ptr)
ARRAY                         *p;
double                  **e1_ptr;
double                  **e2_ptr;
double                  **e3_ptr;
#endif
{
double dof, size, temp;

static double x75, x86;
static double y75, y86;
static double z75, z86;

double s_norm;

static double **r75_ptr = NULL, **r86_ptr = NULL;
static double **ez_ptr  = NULL, **ey_ptr  = NULL;

static double **v1_ptr  = NULL;
static double **v2_ptr  = NULL;
static double **v3_ptr  = NULL;

double **s_ptr = NULL;

int i, j, k, ii, jj, kk, n, n1, n2, k1, k2;
    
#ifdef DEBUG
    printf(" Enter Lamina_Sys_8node() \n");
#endif

    /* =================================================== */
    /* In this case, the normal direction in each node     */
    /* is simplified as the normal direction to element    */
    /* i.e. each element is considered as a plane element  */
    /* =================================================== */

    dof  = p->dof_per_node;
    size = p->size_of_stiff;

    ez_ptr  = MatrixAllocIndirectDouble(3, 1);
    ey_ptr  = MatrixAllocIndirectDouble(3, 1);

    v1_ptr  = MatrixAllocIndirectDouble(3, 1);
    v2_ptr  = MatrixAllocIndirectDouble(3, 1);
    v3_ptr  = MatrixAllocIndirectDouble(3, 1);

    ez_ptr[0][0] = 0.0;
    ez_ptr[1][0] = 0.0;
    ez_ptr[2][0] = 1.0;

    ey_ptr[0][0] = 0.0;
    ey_ptr[1][0] = 1.0;
    ey_ptr[2][0] = 0.0;

 /* --------------------------------------------------- */
 /* [1] Orientatuions of local base vectors             */ 
 /* --------------------------------------------------- */

 /* [a] Compute e3_ptr                      */
 /*     pointer normal to the shell surface */
     
#ifdef DEBUG
    printf(" In elmt_shell_Belytschko_explicit(): \n");
    printf("    : computing e3_ptr \n");
#endif

    x75 =  p->coord[0][6].value - p->coord[0][4].value;
    y75 =  p->coord[1][6].value - p->coord[1][4].value;
    z75 =  p->coord[2][6].value - p->coord[2][4].value;

    r75_ptr = MatrixAllocIndirectDouble(3, 1);
    r75_ptr[0][0] = x75;
    r75_ptr[1][0] = y75;
    r75_ptr[2][0] = z75;

    x86 =  p->coord[0][7].value - p->coord[0][5].value;
    y86 =  p->coord[1][7].value - p->coord[1][5].value;
    z86 =  p->coord[2][7].value - p->coord[2][5].value;

    r86_ptr = MatrixAllocIndirectDouble(3, 1);
    r86_ptr[0][0] = x86;
    r86_ptr[1][0] = y86;
    r86_ptr[2][0] = z86;

#ifdef DEBUG
    dMatrixPrint("r75_ptr", r75_ptr, 3,1);
    dMatrixPrint("r86_ptr", r86_ptr, 3,1);
#endif

    s_ptr   = MatrixAllocIndirectDouble(3, 1);
    s_ptr   = dVmatrixCrossProduct(s_ptr, r75_ptr, 3, 1, r86_ptr, 3,1);

#ifdef DEBUG
    dMatrixPrint("s_ptr", s_ptr, 3,1);
#endif

    s_norm  = (double) dVmatrixL2Norm(s_ptr, 3, 1);

    v3_ptr[0][0] = s_ptr[0][0]/s_norm;
    v3_ptr[1][0] = s_ptr[1][0]/s_norm;
    v3_ptr[2][0] = s_ptr[2][0]/s_norm;

#ifdef DEBUG
    dMatrixPrint("v3_ptr", v3_ptr, 3,1);
#endif

 /*     v1_ptr = ey_ptr X v3_ptr      */

    v1_ptr = dVmatrixCrossProduct(v1_ptr, ey_ptr, 3, 1, v3_ptr, 3,1);

    s_norm  = (double) dVmatrixL2Norm(v1_ptr, 3, 1);

    /* check if ey_ptr is parallel to v3_ptr */

    if(abs(s_norm) <1E-7) { /* parallel */
       v1_ptr = dVmatrixCrossProduct(v1_ptr, ez_ptr, 3, 1, v3_ptr, 3,1);
    }

    v2_ptr = dVmatrixCrossProduct(v2_ptr, v3_ptr, 3, 1, v1_ptr, 3,1);

    /* assume that each node has the same directions in */
    /* one elements                                     */

    for (i = 1; i <= 3; i++) {
         for (k = 1; k <= p->nodes_per_elmt; k++) {
             e1_ptr[i-1][k-1] = v1_ptr[i-1][0];
             e2_ptr[i-1][k-1] = v2_ptr[i-1][0];
             e3_ptr[i-1][k-1] = v3_ptr[i-1][0];
        }
    }

    MatrixFreeIndirectDouble(s_ptr, 3);

    MatrixFreeIndirectDouble(v1_ptr, 3);
    MatrixFreeIndirectDouble(v2_ptr, 3);
    MatrixFreeIndirectDouble(v3_ptr, 3);

    MatrixFreeIndirectDouble(ez_ptr, 3);
    MatrixFreeIndirectDouble(ey_ptr, 3);

    MatrixFreeIndirectDouble(r75_ptr, 3);
    MatrixFreeIndirectDouble(r86_ptr, 3);

#ifdef DEBUG
        printf(" Leaving Lamina_Sys_8node()\n");
#endif

}

/* ======================*/
/* Shell Stiff Matrix    */
/* ======================*/

#ifdef __STDC__
void Shell_Stiff_Plane_8node(double **K, ARRAY *p, double **co_coord, double z_coord, int z_integ_pt, double EE, double nu)
#else
void Shell_Stiff_Plane_8node(K, p, co_coord, z_coord, z_integ_pt, EE, nu)
double         **K;
ARRAY           *p;
double  **co_coord;
double     z_coord;      /* gussain pt in z-direction */  
int     z_integ_pt;      /* integration point No in z-direction */
double      EE, nu;
#endif
{
int                i, j, k,n, ii, jj, kk;
static int     surface_pts, no_integ_pts;
int                            size, dof;
double                    *x_integ_coord;
double                    *y_integ_coord;
double                 **J_inverse= NULL;
static double       weight[16], jacobian;
static double               **shp = NULL;
static double      **stiff_matrix = NULL;
static double                **TT = NULL;
static double          **B_matrix = NULL;
static double         **B1_matrix = NULL; /* local B_matrix           */
static double         **B1_Trans  = NULL; /* local B_matrix Transpose */
static double       **temp_matrix = NULL;

static double    diff, sum, temp1, temp2;

#ifdef DEBUG 
     printf(" Enter Shell_Stiff_Plane_8node() \n");
#endif

    shp            = MatrixAllocIndirectDouble(3,8);
    J_inverse      = MatrixAllocIndirectDouble(3,3);
    TT             = MatrixAllocIndirectDouble(6,6);
    x_integ_coord  = dVectorAlloc(16);
    y_integ_coord  = dVectorAlloc(16);

    stiff_matrix   = MatrixAllocIndirectDouble(p->size_of_stiff,p->size_of_stiff);
    B_matrix       = MatrixAllocIndirectDouble(p->dof_per_node+ 1,p->size_of_stiff);
    B1_matrix      = MatrixAllocIndirectDouble(p->dof_per_node, p->size_of_stiff);
    B1_Trans       = MatrixAllocIndirectDouble(p->size_of_stiff,p->dof_per_node);
    temp_matrix    = MatrixAllocIndirectDouble(p->dof_per_node,p->size_of_stiff);

    surface_pts    = p->integ_ptr->surface_pts;
    no_integ_pts   = (int) 0;
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

        elmt_shell_shape_8node(p, co_coord, x_integ_coord[ii-1],
                        y_integ_coord[ii-1], z_coord, shp,&jacobian,J_inverse, TT, STIFF);

        B_matrix = B_MATRIX_8Node(B_matrix, p, shp, z_coord, J_inverse);

       /* Calculate the B_matrix in local coordinate */
       /* eliminate the 3rd row of TT matrix         */

       for(i = 1; i <= 2; i++) {
           for(j = 1; j <= size; j++) {
               B1_matrix[i-1][j-1] = 0.0;
               B1_Trans[j-1][i-1]  = 0.0;
               for(k = 1; k <= dof + 1; k++)
                   B1_matrix[i-1][j-1] += TT[i-1][k-1]*B_matrix[k-1][j-1];
           }
       }

       for(i = 3; i <= dof; i++) { 
           for(j = 1; j <= size; j++) {
               B1_matrix[i-1][j-1] = 0.0;
               B1_Trans[j-1][i-1]  = 0.0;
               for(k = 1; k <= dof + 1; k++)
                   B1_matrix[i-1][j-1] += TT[i][k-1]*B_matrix[k-1][j-1];
           }
       }

#ifdef DEBUG
       printf(" jacobian = %lf \n", jacobian);
       dMatrixPrint(" shape func", shp,3,8);
       dMatrixPrint(" B_matrix", B_matrix, p->dof_per_node + 1, p->size_of_stiff);
#endif
       for( i = 1; i <= dof; i++) {
           for(j = 1; j <= size; j++) {
               B1_Trans[j-1][i-1] = B1_matrix[i-1][j-1];
           }
       }

       /* ============================================================= */
       /* Update the Material Matrix for Elastic_Plastic Materials      */
       /* ============================================================= */

       kk = no_integ_pts*(z_integ_pt-1) +ii;
       MATER_SHELL_UPDATE(p, co_coord, EE, nu, kk);

       /* calculate the element stiffness matrix */

       temp_matrix  = dMatrixMultRep(temp_matrix,
                      p->mater_matrix->uMatrix.daa, dof, dof, B1_matrix, dof, size); 
       stiff_matrix = dMatrixMultRep(stiff_matrix, B1_Trans,
                      size, dof, temp_matrix, dof, size); 

       for ( i = 1; i <= size; i++) {
            for(j = 1; j<= size; j++) {
                K[i-1][j-1] += stiff_matrix[i-1][j-1]*jacobian*weight[ii-1];
            }
       }

    }  /* end of gaussian integration */

#ifdef DEBUG
               /* check the symmetry of the stiffness matrix */
       printf(" INSIDE CHECK \n");
                for(i = 1; i <= size; i++){
                    for(j = 1; j <= size; j++) {
                        diff = ABS(K[i-1][j-1] - K[j-1][i-1]);

                       if(diff > 1E-7 && ABS(diff/K[i-1][j-1]) > 1E-1 &&
                         ABS(K[i-1][j-1]) > 1E-7) {
                           printf("K[%d][%d] = %le \n", i, j, K[i-1][j-1]);
                           printf("K[%d][%d] = %le \n", j, i, K[j-1][i-1]);
                           printf("diff = %le (diff/K[%d][%d]) = %le\n", diff, i, j, (diff/ABS(K[i-1][j-1])));
                           printf("elmtNo = %d Stiffness matrix IS NOT SYMMETRIC \n", p->elmt_no);
                           break;
                       }
                       else {
                           if( i == size && j == size)
                              printf("elmtNo = %d Stiffness matrix IS SYMMETRIC \n", p->elmt_no);
                       }
                    }
                }
#endif



#ifdef DEBUG
       printf(" end of integration over surface \n"); 
#endif

    free((char *) x_integ_coord);
    free((char *) y_integ_coord);
    MatrixFreeIndirectDouble(shp, 3);
    MatrixFreeIndirectDouble(B_matrix, p->dof_per_node+1);
    MatrixFreeIndirectDouble(B1_matrix, p->dof_per_node);
    MatrixFreeIndirectDouble(B1_Trans, p->size_of_stiff);
    MatrixFreeIndirectDouble(temp_matrix, p->dof_per_node);
    MatrixFreeIndirectDouble(stiff_matrix, p->size_of_stiff);

#ifdef DEBUG 
     dMatrixPrint("stiff surface ", K, p->size_of_stiff, p->size_of_stiff);
     printf(" Leaving Shell_Stiff_Plane_8node() \n");
#endif

}


#define MASS    1

/* =================== */
/* Shell Mass Matrix    */
/* =================== */
#ifdef __STDC__
void  Shell_8Node_Mass(ARRAY *p, MATRIX *mass, double **co_coord, double density, double thickness)
#else
void  Shell_8Node_Mass(p, mass, co_coord, density, thickness)
ARRAY               *p;
MATRIX           *mass;
double      **co_coord;
double         density;
double       thickness;
#endif
{
int  i, j, k, n, ii,jj, length, length1, length2;
int    no_integ_pts, x_integ_pts, y_integ_pts;
int                               z_integ_pts;

double                     **J_inverse = NULL;
static double                     **TT = NULL;
static double                 **e1_ptr = NULL;
static double                 **e2_ptr = NULL;
static double                 **e3_ptr = NULL;

double             **shp, **Mass, **Mass_temp;
double                      **N1_Trans = NULL;
double                      **N2_Trans = NULL;
double               **N1 = NULL, **N2 = NULL;

double         *x_integ_coord, *y_integ_coord;
double                  *z_integ_coord = NULL;
double                       *z_weight = NULL;
double                   weight[16], jacobian;
double                                    sum;

double               Aera, elmt_length, width;
double            node_mass, node_Jx, node_Jy;

DIMENSIONS                      *d1, *d2, *d3;
int                   UNITS_SWITCH, UnitsType;

#ifdef DEBUG
       printf("*** Enter Shell_8Node_Mass(): \n");
#endif
     /* [a] mass initialization */

      for(i = 1; i <= p->size_of_stiff; i++) 
          for(j = 1; j <= p->size_of_stiff; j++)
                mass->uMatrix.daa[i-1][j-1] = 0.0;

    /* [b] :  mass matrix      */

#ifdef DEBUG
       printf("*** In Shell_8Node_Mass(): in main swicth() \n");
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
             node_mass   =  density*Aera*elmt_length/8.0;
             node_Jx     =  density*elmt_length*thickness*width*width*width/12.0/8.0;
             node_Jy     =  density*Aera*elmt_length*elmt_length*elmt_length/12.0/8.0;
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

             e1_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);
             e2_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);
             e3_ptr = MatrixAllocIndirectDouble(3,p->nodes_per_elmt);

             Mass = MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff); 
             N1             = MatrixAllocIndirectDouble(3, p->size_of_stiff);
             N2             = MatrixAllocIndirectDouble(3, p->size_of_stiff);
             N1_Trans       = MatrixAllocIndirectDouble(p->size_of_stiff, 3);
             N2_Trans       = MatrixAllocIndirectDouble(p->size_of_stiff, 3);
             shp            = MatrixAllocIndirectDouble(3,8);
             x_integ_coord  = dVectorAlloc(16);
             y_integ_coord  = dVectorAlloc(16);
             x_integ_pts    = p->integ_ptr->surface_pts;
             y_integ_pts    = p->integ_ptr->surface_pts;

             z_integ_pts    = p->integ_ptr->thickness_pts;
             z_integ_coord  = dVectorAlloc(z_integ_pts + 1);
             z_weight       = dVectorAlloc(z_integ_pts + 1);


             no_integ_pts   = 0;

             gauss(z_integ_coord, z_weight,z_integ_pts);

             pgauss(x_integ_pts, &no_integ_pts, x_integ_coord, y_integ_coord, weight);

             Lamina_Sys_8node(p, e1_ptr, e2_ptr, e3_ptr);

             /* Integration loop over thickness */

             for(jj = 1; jj <= z_integ_pts; jj++) {

             /* Integration loop over surface */

                 for( ii = 1; ii <= no_integ_pts; ii++) {

                     elmt_shell_shape_8node(p, co_coord, x_integ_coord[ii-1],y_integ_coord[ii-1],
                                           z_integ_coord[ii], shp,&jacobian, J_inverse, TT, MASS);

                 
                     for(j = 1; j <= p->nodes_per_elmt; j++) {
                         k = (j-1)*p->dof_per_node + 1;

                         N1[0][k-1] = shp[2][j-1];
                         N1[0][k]   = 0.0;
                         N1[0][k+1] = 0.0;
                         N1[0][k+2] = 0.0;
                         N1[0][k+3] = 0.0;

                         N1[1][k-1] = 0.0;
                         N1[1][k]   = shp[2][j-1];
                         N1[1][k+1] = 0.0;
                         N1[1][k+2] = 0.0;
                         N1[1][k+3] = 0.0;
  
                         N1[2][k-1] = 0.0;
                         N1[2][k]   = 0.0;
                         N1[2][k+1] = shp[2][j-1];
                         N1[2][k+2] = 0.0;
                         N1[2][k+3] = 0.0;

                         N2[0][k-1] = 0.0;
                         N2[0][k]   = 0.0;
                         N2[0][k+1] = 0.0;
                         N2[0][k+2] = -e2_ptr[0][j-1]*0.5*thickness*shp[2][j-1];
                         N2[0][k+3] =  e1_ptr[0][j-1]*0.5*thickness*shp[2][j-1];
  
                         N2[1][k-1] = 0.0;
                         N2[1][k]   = 0.0;
                         N2[1][k+1] = 0.0;
                         N2[1][k+2] = -e2_ptr[1][j-1]*0.5*thickness*shp[2][j-1];
                         N2[1][k+3] =  e1_ptr[1][j-1]*0.5*thickness*shp[2][j-1];

                         N2[2][k-1] = 0.0;
                         N2[2][k]   = 0.0;
                         N2[2][k+1] = 0.0;
                         N2[2][k+2] = -e2_ptr[2][j-1]*0.5*thickness*shp[2][j-1];
                         N2[2][k+3] =  e1_ptr[2][j-1]*0.5*thickness*shp[2][j-1];
                     }
            
                    /* [M] calculation  */

                     for(i = 1; i <= p->dof_per_node; i++) {
                         for(j = 1; j <= p->size_of_stiff; j++) {
                             N1[i-1][j-1] += z_integ_coord[jj]*N2[i-1][j-1];
                             N1_Trans[j-1][i-1] = N1[i-1][j-1];
                         }
                     }

                     Mass = dMatrixMultRep(Mass,N1_Trans,p->size_of_stiff,
                            p->dof_per_node,N1,p->dof_per_node, p->size_of_stiff);
#ifdef DEBUG
          dMatrixPrint("M", Mass, p->size_of_stiff, p->size_of_stiff);
#endif
                    for(i = 1 ; i <= p->size_of_stiff; i++) 
                        for(j = 1 ; j <= p->size_of_stiff; j++) 
                             mass->uMatrix.daa[i-1][j-1] += Mass[i-1][j-1]*jacobian*weight[ii-1]*z_weight[jj]*density;

                 } /* end of surface loop */
             } /* end of thickness loop */


             MatrixFreeIndirectDouble(e1_ptr, 3); 
             MatrixFreeIndirectDouble(e2_ptr, 3); 
             MatrixFreeIndirectDouble(e3_ptr, 3); 

             MatrixFreeIndirectDouble(N1, p->dof_per_node); 
             MatrixFreeIndirectDouble(N2, p->dof_per_node); 
             MatrixFreeIndirectDouble(N1_Trans, p->size_of_stiff);
             MatrixFreeIndirectDouble(N2_Trans, p->size_of_stiff);
             MatrixFreeIndirectDouble(shp,3);
             MatrixFreeIndirectDouble(Mass, p->size_of_stiff);
             free((char *) x_integ_coord);
             free((char *) y_integ_coord);
             free((char *) z_integ_coord);
             free((char *) z_weight);
	     break;
	default:
             FatalError("In Shell_8Node_Mass() : Type of Mass Matrix Undefined",(char *)NULL);
	     break;
    }
       
 /* ------------ MASS UNITS ---------------------------- */
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

 /* ------------====================-------------------- */

#ifdef DEBUG
       MatrixPrintIndirectDouble(mass);
       printf("*** leaving Shell_8Node_Mass()  \n");
#endif

}

/* ================================================== */
/* Shell_Element Shape Functions for 4-Node Element   */
/* ================================================== */

#ifdef __STDC__
void elmt_shell_shape_8node(ARRAY *p, double **coord, double ss, double tt, double qq, double **shp, double *jacobian, double **J_inverse, double **TT, int MASS_FLAG)
#else
void elmt_shell_shape_8node(p, coord, ss,tt, qq, shp,jacobian, J_inverse, TT, MASS_FLAG)
ARRAY                             *p;
double                   **J_inverse;
double                          **TT;
double                **coord, **shp;
double         ss, tt, qq, *jacobian;
int                        MASS_FLAG;
#endif
{
int                        i, j, k;
double        *s = NULL, *t = NULL;
double                 **xs = NULL;

static double            thickness;
static double       **s_ptr = NULL;

static double      **e1_ptr = NULL;
static double      **e2_ptr = NULL;
static double      **e3_ptr = NULL;

static double      **v1_ptr = NULL;
static double      **v2_ptr = NULL;
static double      **v3_ptr = NULL;
static double                 norm;

#ifdef DEBUG
   printf(" Enter elmt_shell_shape_8node() \n");
#endif

    s        = dVectorAlloc(8);
    t        = dVectorAlloc(8);
    xs       = MatrixAllocIndirectDouble(3,3);

    v1_ptr  = MatrixAllocIndirectDouble(3, 1);
    v2_ptr  = MatrixAllocIndirectDouble(3, 1);
    v3_ptr  = MatrixAllocIndirectDouble(3, 1);

    e1_ptr  = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);
    e2_ptr  = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);
    e3_ptr  = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);

    s_ptr   = MatrixAllocIndirectDouble(3, 1);

    thickness = p->work_section[11].value;

    s[0] = -1.0; s[1] =  1.0; s[2] =  1.0; s[3] = -1.0;
    s[4] =  0.0; s[5] =  1.0; s[6] =  0.0; s[7] = -1.0;

    t[0] = -1.0; t[1] = -1.0; t[2] =  1.0; t[3] =  1.0;
    t[4] = -1.0; t[5] =  0.0; t[6] =  1.0; t[7] =  0.0;

  switch(p->nodes_per_elmt) {
    case 8:

    /* form 8-node quadrilateral shape function                    */
    /* shape function:                  Ni = shape[2][i-1]         */
    /*                                  node no. i = 1, .. 8       */
    /* derivatives of shape functions:  dNi/d(ss) = shape[0][i-1]  */
    /*                                  dNi/d(tt) = shape[1][i-1]  */

    for(i = 1; i <= 4; i++){
        shp[2][i-1] = 0.25*(1.0+ s[i-1]*ss)*(1.0 + t[i-1] * tt)*(s[i-1] * ss + t[i-1] * tt-1.0); 
        shp[0][i-1] = s[i-1] * 0.25 * (1.0 + t[i-1] * tt)*(2.0*s[i-1] * ss + t[i-1] * tt);                
        shp[1][i-1] = t[i-1] * 0.25 * (1.0 + s[i-1] * ss)*(2.0*t[i-1] * tt + s[i-1] * ss);                
    }
    shp[2][5-1] = 0.5*(1.0 - ss * ss) * ( 1.0 + t[5-1] * tt);
    shp[2][7-1] = 0.5*(1.0 - ss * ss) * ( 1.0 + t[7-1] * tt);

    shp[0][5-1] = - ss * (1.0 + t[5-1] * tt);
    shp[1][5-1] = t[5-1] * 0.5 * (1.0 - ss * ss);

    shp[0][7-1] = - ss * (1.0 + t[7-1] * tt);
    shp[1][7-1] = t[7-1] * 0.5 * (1.0 - ss * ss);

    shp[2][6-1] = 0.5*(1.0 - tt * tt) * ( 1.0 + s[6-1] * ss);
    shp[2][8-1] = 0.5*(1.0 - tt * tt) * ( 1.0 + s[8-1] * ss);

    shp[0][6-1] =  0.5* s[6-1] * (1.0 - tt * tt);
    shp[1][6-1] = -tt *(1.0 + s[6-1] * ss);

    shp[0][8-1] =  0.5* s[8-1] * (1.0 - tt * tt);
    shp[1][8-1] = -tt *(1.0 + s[8-1] * ss);

    /* construct jacobian matrix and its determinant */

     Lamina_Sys_8node(p, e1_ptr, e2_ptr, e3_ptr);

     for(i = 1; i <=  3; i++) {
         xs[i-1][0] = 0.0;
         xs[i-1][1] = 0.0;
         xs[i-1][2] = 0.0;
         for(k = 1; k <= p->nodes_per_elmt; k++) {
             xs[i-1][0] += coord[i-1][k-1]*shp[0][k-1] +
                           thickness*0.5*qq*shp[0][k-1]*e3_ptr[i-1][k-1];

             xs[i-1][1] += coord[i-1][k-1]*shp[1][k-1] +
                           thickness*0.5*qq*shp[1][k-1]*e3_ptr[i-1][k-1];

             xs[i-1][2] += thickness*0.5*shp[2][k-1]*e3_ptr[i-1][k-1];
         }
     }

     *jacobian =  xs[0][0]*(xs[1][1] * xs[2][2] - xs[1][2] *xs[2][1]) 
                - xs[0][1]*(xs[1][0] * xs[2][2] - xs[1][2] *xs[2][0]) 
                + xs[0][2]*(xs[1][0] * xs[2][1] - xs[1][1] *xs[2][0]);

#ifdef DEBUG
        for(k = 1; k <= p->nodes_per_elmt; k++) {
           printf(" x[%d] = %lf \n", k, coord[0][k-1]);
           printf(" y[%d] = %lf \n", k, coord[1][k-1]);
           printf(" z[%d] = %lf \n", k, coord[2][k-1]);
         }
    printf(" thickness = %le, jacobian = %le \n", thickness, *jacobian);
#endif

     if(MASS_FLAG == MASS) {
        break;
     }	

    /* compute Jacobain inverse matrix */

    J_inverse[0][0] = (xs[1][1]*xs[2][2] - xs[1][2]*xs[2][1])/ *jacobian;
    J_inverse[1][1] = (xs[0][0]*xs[2][2] - xs[0][2]*xs[2][0])/ *jacobian;
    J_inverse[2][2] = (xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0])/ *jacobian;

    J_inverse[0][1] = -(xs[1][0]*xs[2][2] - xs[1][2]*xs[2][0])/ *jacobian;
    J_inverse[0][2] =  (xs[1][0]*xs[2][1] - xs[1][1]*xs[2][0])/ *jacobian;

    J_inverse[1][0] = -(xs[0][1]*xs[2][2] - xs[0][2]*xs[2][1])/ *jacobian;
    J_inverse[2][0] =  (xs[0][1]*xs[1][2] - xs[0][2]*xs[1][1])/ *jacobian;

    J_inverse[1][2] = -(xs[0][0]*xs[2][1] - xs[0][1]*xs[2][0])/ *jacobian;
    J_inverse[2][1] = -(xs[0][0]*xs[1][2] - xs[0][2]*xs[1][0])/ *jacobian;

#ifdef DEBUG
     /* Check the J-inverse */

     for(i = 1; i <= 3; i++) 
         for(j = 1; j <= 3; j++) {
             norm = 0.0;
             for(k = 1; k <= 3; k++) 
                 norm += J_inverse[i-1][k-1]*xs[k-1][j-1];
             printf("I[%d][%d] = %lf \n", i, j, norm);
         }
#endif

    /* Form Local vector */
    
    for(i = 1; i <= 3; i++) {
        v1_ptr[i-1][0] = xs[i-1][0];
        s_ptr[i-1][0]  = xs[i-1][1];
    }

    v3_ptr = dVmatrixCrossProduct(v3_ptr, v1_ptr, 3, 1, s_ptr, 3,1);
    
    norm  = (double) dVmatrixL2Norm(v1_ptr, 3, 1); 
    for(i = 1; i <= 3; i++) 
        v1_ptr[i-1][0] = v1_ptr[i-1][0]/norm;
     
    norm  = (double) dVmatrixL2Norm(v3_ptr, 3, 1); 
    for(i = 1; i <= 3; i++) 
        v3_ptr[i-1][0] = v3_ptr[i-1][0]/norm;
    
    v2_ptr = dVmatrixCrossProduct(v2_ptr, v3_ptr, 3, 1, v1_ptr, 3,1);
  
    TT[0][0] = v1_ptr[0][0]*v1_ptr[0][0];
    TT[0][1] = v1_ptr[1][0]*v1_ptr[1][0];
    TT[0][2] = v1_ptr[2][0]*v1_ptr[2][0];
    
    TT[0][3] = v1_ptr[0][0]*v1_ptr[1][0];
    TT[0][4] = v1_ptr[1][0]*v1_ptr[2][0];
    TT[0][5] = v1_ptr[2][0]*v1_ptr[0][0];

    TT[1][0] = v2_ptr[0][0]*v2_ptr[0][0];
    TT[1][1] = v2_ptr[1][0]*v2_ptr[1][0];
    TT[1][2] = v2_ptr[2][0]*v2_ptr[2][0];
   
    TT[1][3] = v2_ptr[0][0]*v2_ptr[1][0];
    TT[1][4] = v2_ptr[1][0]*v2_ptr[2][0];
    TT[1][5] = v2_ptr[2][0]*v2_ptr[0][0];

    TT[2][0] = v3_ptr[0][0]*v3_ptr[0][0];
    TT[2][1] = v3_ptr[1][0]*v3_ptr[1][0];
    TT[2][2] = v3_ptr[2][0]*v3_ptr[2][0];
    
    TT[2][3] = v3_ptr[0][0]*v3_ptr[1][0];
    TT[2][4] = v3_ptr[1][0]*v3_ptr[2][0];
    TT[2][5] = v3_ptr[2][0]*v3_ptr[0][0];

    TT[3][0] = 2.0*v1_ptr[0][0]*v2_ptr[0][0];
    TT[3][1] = 2.0*v1_ptr[1][0]*v2_ptr[1][0];
    TT[3][2] = 2.0*v1_ptr[2][0]*v2_ptr[2][0];
    
    TT[3][3] = v1_ptr[0][0]*v2_ptr[1][0] + v2_ptr[0][0]*v1_ptr[1][0];
    TT[3][4] = v1_ptr[1][0]*v2_ptr[2][0] + v2_ptr[1][0]*v1_ptr[2][0];
    TT[3][5] = v1_ptr[2][0]*v2_ptr[0][0] + v2_ptr[2][0]*v1_ptr[0][0];

    TT[4][0] = 2.0*v2_ptr[0][0]*v3_ptr[0][0];
    TT[4][1] = 2.0*v2_ptr[1][0]*v3_ptr[1][0];
    TT[4][2] = 2.0*v2_ptr[2][0]*v3_ptr[2][0];
    
    TT[4][3] = v2_ptr[0][0]*v3_ptr[1][0] + v3_ptr[0][0]*v2_ptr[1][0];
    TT[4][4] = v2_ptr[1][0]*v3_ptr[2][0] + v3_ptr[1][0]*v2_ptr[2][0];
    TT[4][5] = v2_ptr[2][0]*v3_ptr[0][0] + v3_ptr[2][0]*v2_ptr[0][0];

    TT[5][0] = 2.0*v3_ptr[0][0]*v1_ptr[0][0];
    TT[5][1] = 2.0*v3_ptr[1][0]*v1_ptr[1][0];
    TT[5][2] = 2.0*v3_ptr[2][0]*v1_ptr[2][0];
    
    TT[5][3] = v3_ptr[0][0]*v1_ptr[1][0] + v1_ptr[0][0]*v3_ptr[1][0];
    TT[5][4] = v3_ptr[1][0]*v1_ptr[2][0] + v1_ptr[1][0]*v3_ptr[2][0];
    TT[5][5] = v3_ptr[2][0]*v1_ptr[0][0] + v1_ptr[2][0]*v3_ptr[0][0];

#ifdef DEBUG
    dMatrixPrint(" J^-1", J_inverse, 3,3);
    dMatrixPrint(" shp", shp, 3,8);
#endif
  
    break;

    case 4:
    break;

    default:
    break;
 }
#ifdef DEBUG
    dMatrixPrint(" shp", shp, 3,8);
#endif

  free((char *) s);
  free((char *) t);
  MatrixFreeIndirectDouble(xs, 3);
  MatrixFreeIndirectDouble(s_ptr, 3);
  MatrixFreeIndirectDouble(v1_ptr, 3);
  MatrixFreeIndirectDouble(v2_ptr, 3);
  MatrixFreeIndirectDouble(v3_ptr, 3);
  MatrixFreeIndirectDouble(e1_ptr, 3);
  MatrixFreeIndirectDouble(e2_ptr, 3);
  MatrixFreeIndirectDouble(e3_ptr, 3);

#ifdef DEBUG
   printf(" leaving elmt_shell_shape_8node() \n");
#endif

}

/* Print SHELL_8N Element Properties */
#ifdef __STDC__
void print_property_shell_8n(EFRAME *frp, int i)
#else
void print_property_shell_8n(frp, i)
EFRAME    *frp;
int          i;                 /* elmt_attr_no */
#endif
{
int     UNITS_SWITCH;
ELEMENT_ATTR    *eap;

#ifdef DEBUG
       printf("*** Enter print_property_shell_8n()\n");
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
       printf("*** Leave print_property_shell_8n()\n");
#endif
}
