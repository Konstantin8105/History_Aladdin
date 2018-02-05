/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_psps.c : Plane Stress Plane Strain Element
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

#include <math.h>
#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "vector.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"
#include "elmt.h"

/* Function Declarations */

double **MaterialMatrixPlane( double, double, double );

/*
 *  ============================================================================ 
 *  Plane-stress/plane strain finite element         
 * 
 *  The material and section properties are: 
 * 
 *      p->work_material[0] = E;
 *      p->work_material[1] = G;
 *      p->work_material[2] = fy;
 *      p->work_material[3] = ET;
 *      p->work_material[4] = nu;
 *      p->work_material[5] = density;
 *      p->work_material[6] = fu;
 *      p->work_material[7] = alpha_thermal[0];
 *      p->work_material[8] = alpha_thermal[1];
 *
 *      p->work_section[0] = Ixx;
 *      p->work_section[1] = Iyy;
 *      p->work_section[2] = Izz;
 *      p->work_section[3] = Ixy;
 *      p->work_section[4] = Ixz;
 *      p->work_section[5] = Iyz;
 *      p->work_section[6] = weight;
 *      p->work_section[7] = bf;
 *      p->work_section[8] = tf;
 *      p->work_section[9] = depth;
 *      p->work_section[10] = area;                                 
 * 
 *  Compute (3x2) derivative matrix [B] .... 
 * 
 *       B_i[0][0] = dNi/dx, B_i[0][1] = 0  
 *       B_i[1][0] = 0,      B_i[1][1] = dNi/dy     
 *       B_i[2][0] = dNi/dy, B_i[2][2] = dNi/dx   
 *
 *       [B] = [B_1, B_2, B_3, B_4, ..., B_n] where n = no of node                             
 *
 *       The output is : shp[0][i-1] = dN_i/dx                
 *                       shp[1][i-1] = dN_i/dy                
 *       compute [B] matrix at each Gaussian integration point.
 *
 *  Note : Material properties such and "E" and "nu" must be stored as static
 *         variables.
 *  ============================================================================ 
 */

ARRAY *elmt_psps(p,isw)
ARRAY *p;
int isw;
{
static QUANTITY fy, G, E, ET, density, lambda;
static double nu;
static double temperature;
static double *alpha_thermal;
static double dFlag;
static int no_integ_pts;
static int no_stress_pts;
int ii, k, i, j,j1, k1, lint, kk;
double sg[17],tg[17],wg[17],sig[7],
       jacobian,w11,w12,w21,w22, xx, yy, dv, wd[3];
double **B_matrix;
double **B_Transpose;
double **stiff;
double **Mater_matrix;
double **body_force;
double **strain;
double **stress;
double **temp_change;
double **load;
double **displ;
double **nodal_load;
double **shp;
double **temp_matrix;

int             dof_per_elmt;
int          length1, length;
int             UNITS_SWITCH;
double         *sum_row, sum;
DIMENSIONS        *dp_stress;
DIMENSIONS        *dp_length;
DIMENSIONS     *d1, *d2, *d3;

   dof_per_elmt  = p->nodes_per_elmt*p->dof_per_node;  

    /* [a] : Allocate memory for matrices */

   strain        = MatrixAllocIndirectDouble(3, 1);
   stress        = MatrixAllocIndirectDouble(3, 1);
   displ         = MatrixAllocIndirectDouble(dof_per_elmt, 1);
   stiff         = MatrixAllocIndirectDouble(dof_per_elmt, dof_per_elmt);
   load          = MatrixAllocIndirectDouble(dof_per_elmt, 1);
   nodal_load    = MatrixAllocIndirectDouble(dof_per_elmt, 1);
   B_matrix      = MatrixAllocIndirectDouble(3, dof_per_elmt);
   shp           = MatrixAllocIndirectDouble(3, p->nodes_per_elmt);
   sum_row       = dVectorAlloc(dof_per_elmt);

   B_Transpose   = MatrixAllocIndirectDouble(dof_per_elmt, 3);
   temp_matrix   = MatrixAllocIndirectDouble(dof_per_elmt, 3);
   body_force    = MatrixAllocIndirectDouble(3, 1);

   /* [b] : deal with Individual Tasks */

   UNITS_SWITCH = CheckUnits();
   switch(isw) {
       case PROPTY: 
            alpha_thermal = dVectorAlloc(2);

            /* Input material properties */

            lint = 0;
            E.value       =  p->work_material[0].value; 
            fy.value      =  p->work_material[2].value;
            ET.value      =  p->work_material[3].value;
            nu            =  p->work_material[4].value; 
            density.value =  p->work_material[5].value;
            if( UNITS_SWITCH == ON ) {
                E.dimen       =  p->work_material[0].dimen; 
                fy.dimen      =  p->work_material[2].dimen;
                ET.dimen      =  p->work_material[3].dimen;
                density.dimen =  p->work_material[5].dimen;
            }

            /* Initialize thermal expansion coefficients in x- and y- directions */

            alpha_thermal[0] = p->work_material[7].value; 
            alpha_thermal[1] = p->work_material[8].value; 

            /* 4-node and 8-node elements */

            if(p->nodes_per_elmt == 4)
               no_integ_pts  =  2;
            if(p->nodes_per_elmt == 8)
               no_integ_pts  =  3;

            no_stress_pts = no_integ_pts; /* stress points, subject to change later */

            /* Set flag for Plane Stress/Plane Strain material matrix computations */

            dFlag = 1;
            if(strcmp("PLANE_STRAIN", p->elmt_type) == 0) {
               dFlag = 2;
            }

            break;
       case CHERROR:
            break;
       case STIFF: 

            /* Zero out the Stiffness Matrix */

            for(i = 1; i <= dof_per_elmt; i++) 
                for(j = 1; j<= dof_per_elmt; j++)
                    stiff[i-1][j-1] = 0.0;

            /* Get Gauss Integration points */

            if(no_integ_pts*no_integ_pts != lint)
               pgauss(no_integ_pts,&lint,sg,tg,wg);

            /* Start Gaussian Integration  */

            for(ii = 1; ii <= lint; ii++) { 

                /* Compute shape functions */

                shape( sg[ii-1],tg[ii-1], p->coord,shp,&jacobian, p->no_dimen,
                       p->nodes_per_elmt, p->node_connect,FALSE);

                /*
                 *  =================================================== 
                 *  Compute material matrix
                 *  
                 *  Plane Stress : 3rd argument ii = 1
                 *  Plane Strain : 3rd argument ii = 2
                 *  =================================================== 
                 */

                Mater_matrix = MaterialMatrixPlane( E.value, nu, dFlag );

                /* Form [B] matrix */
            
                for(j = 1; j <= p->nodes_per_elmt; j++) { 
                    k = 2*(j-1);
                    B_matrix[0][k]   = shp[0][j-1];
                    B_matrix[1][k+1] = shp[1][j-1];
                    B_matrix[2][k]   = shp[1][j-1];
                    B_matrix[2][k+1] = shp[0][j-1];
                }
           
                /* Multiply Jacobian determinant with weight coefficents  */
                /* and material matrix                                    */
          
                jacobian = jacobian*wg[ii-1];
                for( i = 1; i <=3 ; i++ ) 
                    for (j = 1; j <= 3; j++)  
                         Mater_matrix[i-1][j-1] *= jacobian;

                /* Transpose [B] matrix */

                for(i = 1; i <= 3; i++)
                    for(j = 1; j <= dof_per_elmt; j++) 
                        B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

                /* Compute [B]^T*[Mater] and save in [B]^T */

                temp_matrix = dMatrixMultRep(temp_matrix, B_Transpose,
                                             dof_per_elmt, 3, Mater_matrix, 3, 3);

                /* Calculate stiffness : Integral of [B]^T * [Mater] * [B] */

                stiff = dMatrixMultRep(stiff, temp_matrix, dof_per_elmt, 3,
                                       B_matrix, 3, dof_per_elmt);
 
                for (i = 1; i <= dof_per_elmt; i++) 
                for (j = 1; j <= dof_per_elmt; j++) 
                     p->stiff->uMatrix.daa[i-1][j-1]  += stiff[i-1][j-1];

                /* Free memory for material matrix with modified contents */

                MatrixFreeIndirectDouble(Mater_matrix, 3);
            }

            /* [c] : Assign units to stiffness matrix */
     
            if ( UNITS_SWITCH == ON ) {
               if( CheckUnitsType() == SI) {
                   d1 = DefaultUnits("Pa");
                   d2 = DefaultUnits("m");
               } else {
                   d1 = DefaultUnits("psi");
                   d2 = DefaultUnits("in");
               }

               /* Node 1 */

               UnitsMultRep( &(p->stiff->spColUnits[0]), d1, d2 );
               UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );

               /* Nodes 2 through 4 */

               for(i = 2; i <= p->nodes_per_elmt; i++) {
               for(j = 1; j <= p->dof_per_node; j++) {
                   k  = p->dof_per_node*(i-1) + j;
                   UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
               }
               }

               free((char *) d1->units_name);
               free((char *) d1);
               free((char *) d2->units_name);
               free((char *) d2);

               for( i=1 ; i<=p->size_of_stiff ; i++ )
                    ZeroUnits( &(p->stiff->spRowUnits[i-1]) );

            }
            break;
       case EQUIV_NODAL_LOAD: /* calculate the equivalent nodal loads */
                              /* due to distributed loadings          */
           
            /* Initialize nodal_load */
      
            for(i = 1; i <= dof_per_elmt; i++) {
                p->equiv_nodal_load->uMatrix.daa[i-1][0] = 0.0;
                nodal_load[i-1][0] = 0.0;
                load[i-1][0]       = 0.0;
            }
 
            for(i = 1; i<= no_integ_pts*no_integ_pts; i++) {
                sg[i-1] = 0.0; tg[i-1] = 0.0; wg[i-1] = 0.0;
            }

            /* Compute Gaussian Integration Points */
      
            if(no_integ_pts*no_integ_pts != lint)
               pgauss(no_integ_pts,&lint,sg,tg,wg); 

            /* Main loop for Gauss Integration */

            for(ii = 1; ii <= lint; ii++) { 

                /* Compute shape functions */

                shape(sg[ii-1],tg[ii-1],p->coord,shp,&jacobian,
                      p->no_dimen,p->nodes_per_elmt, p->node_connect,0);

                /* Form [B] matrix */
            
                for(j = 1; j <= p->nodes_per_elmt; j++) { 
                    k = 2*(j-1);
                    B_matrix[0][k]   = shp[0][j-1];
                    B_matrix[1][k+1] = shp[1][j-1];
                    B_matrix[2][k]   = shp[1][j-1];
                    B_matrix[2][k+1] = shp[0][j-1];
                }

                /* Get material matrix : see notes on dFlag above */

                Mater_matrix = MaterialMatrixPlane( E.value, nu, dFlag );
            
                /* Multiply Jacobian determinant with weight  */
                /* coefficents and material matrix            */

                jacobian *= wg[ii-1];
           
                /* [a] Calculate equivalent nodal forces due to initial strains */

                if(p->nodal_init_strain != NULL) {

                   for( i = 1; i <= 3; i++ ) 
                   for( j = 1; j <= 3; j++ )  
                        Mater_matrix[i-1][j-1] *= jacobian;

                   /* Transpose [B] matrix */

                   for(i = 1; i <= 3; i++) 
                   for(j = 1; j <= dof_per_elmt; j++)
                        B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

                   /* Calculate strain[] at gaussian points : strain = sum Ni*nodal_strain */

                   for (i = 1; i <= 3; i++)
                      strain[i-1][0]  =  0.0;
                
                   for (i = 1; i <= 3; i++) {
                   for( j = 1; j <= p->nodes_per_elmt; j++) {
                      strain[i-1][0]  += shp[2][j-1] * p->nodal_init_strain[i-1][j-1];
                   }
                   }

                   /* [Mater]_3x3 * [strain]_3x1 and save in [stress]_3x1    */

                   stress = dMatrixMultRep(stress, Mater_matrix, 3, 3, strain, 3, 1);
           
                   /* Mutiply [B]^T * [Mater]*[strain] */

                   if(nodal_load == NULL ) { 
                      nodal_load = dMatrixMultRep(nodal_load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                   } else {
                      load = dMatrixMultRep(load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                      for (i = 1; i<= dof_per_elmt; i++)  {
                          nodal_load[i-1][0] += load[i-1][0];
                      }
                   }
                }

                /* [b] Calculate equivalent nodal force due to internal stress */

                if(p->nodal_init_stress != NULL) {

                   /* Transpose [B] matrix */

                   for(i = 1; i <= 3; i++)
                   for(j = 1; j <= dof_per_elmt; j++)
                       B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

                   /* Calculate stress[] at gaussian points : stress = sum Ni*nodal_stress */
    
                   for (i = 1; i <= 3; i++)
                      stress[i-1][0]  =  0.0;
                
                   for(i = 1; i <= 3; i++) {
                      for(j = 1; j <= p->nodes_per_elmt; j++) {
                          stress[i-1][0]  += shp[2][j-1] * p->nodal_init_stress[i-1][j-1].value;
                      }
                      stress[i-1][0] *= jacobian;
                   }

                   /* Multiply [B]^T * [stress] */
               
                   if(nodal_load == NULL) {
                      nodal_load  = dMatrixMultRep(nodal_load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                   } else {
                      load        = dMatrixMultRep(load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                      for (i = 1; i<= dof_per_elmt; i++) 
                          nodal_load[i-1][0] -= load[i-1][0];
                   }
                }

                /* [c] Calculate equivalent nodal force due to body force */

                if(p->nodal_body_force != NULL) {

                   /* Calculate body_force[] at gaussian points --  */
                   /* body_force = sum of [ Ni*body_force ]         */

                   for (i = 1; i <= p->no_dimen; i++)
                      body_force[i-1][0]  =  0.0;
                
                   for(i = 1; i <= p->no_dimen; i++) {
                      for( j = 1; j <= p->nodes_per_elmt; j++) {
                          body_force[i-1][0] 
                               += shp[2][j-1] * p->nodal_body_force[i-1][j-1].value*jacobian;

                          /* Mutiply [N]^T * [body_force] */

                          k = (j-1)*p->no_dimen + i;
                          if(nodal_load == NULL) {
                             nodal_load[k-1][0] =  shp[2][j-1]*body_force[i-1][0];
                          } else
                             nodal_load[k-1][0] += shp[2][j-1]*body_force[i-1][0];
                      }
                   }
                }

                /* [d] Calculate equivalent nodal force due to temperature change */ 
  
                if(p->nodal_init_strain != NULL) {

                   for(i = 1; i <= 3 ; i++ ) 
                   for(j = 1; j <= 3; j++)  
                       Mater_matrix[i-1][j-1] *= jacobian;

                   /* Transpose [B] matrix */

                   for(i = 1; i <= 3; i++)
                   for(j = 1; j <= dof_per_elmt; j++)
                       B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

                   /*  [B]^T * [Mater]  and save in [B]^T    */

                   temp_matrix = dMatrixMultRep(temp_matrix, B_Transpose,
                                                dof_per_elmt, 3, Mater_matrix, 3, 3);
               
                   /* Calculate strain[] at gaussian points : */
                   /* strain = sum Ni*nodal_strain            */

                   for(i = 1; i <= 3; i++)
                       strain[i-1][0]  =  0.0;
                
                   for(j = 1; j <= p->nodes_per_elmt; j++) {
                       temperature  += shp[2][j-1] * p->nodal_temp[j-1].value;
                   }

                   for(i = 1; i <= 2; i++) {
                       strain[i-1][0]  = temperature *alpha_thermal[i-1];
                   }

                   strain[2][0] = 0.0;
                    
                   /* mutiply [B]^T * [Mater]* [strain] */
               
                   if(nodal_load == NULL) {
                      nodal_load  = dMatrixMultRep(nodal_load, temp_matrix,
                                                   dof_per_elmt, 3, strain, 3, 1);
                   } else {
                      load        = dMatrixMultRep(load, temp_matrix, dof_per_elmt, 3, strain, 3, 1);
                      for(i = 1; i<= dof_per_elmt; i++) 
                          nodal_load[i-1][0] += load[i-1][0];
                   }
                }

                /* Transfer nodal_load to p->equiv_nodal_load */

                for(i = 1; i <= dof_per_elmt; i++) {
                    p->equiv_nodal_load->uMatrix.daa[i-1][0] +=  nodal_load[i-1][0];
                }
            }

            /*
             *  =================================================== 
             *  Initiation of Equivalent nodal load Units Buffer   
             *  
             *  Note : If Young's Modulus E is in SI then Use SI, 
             *         otherwise, use US. 
             *  =================================================== 
             */

            if (UNITS_SWITCH == ON) {
               if( CheckUnitsType() == SI)
                   d1 = DefaultUnits("N");
               else
                   d1 = DefaultUnits("lbf");

               /* node 1 */

               UnitsCopy( &(p->equiv_nodal_load->spRowUnits[0]), d1 ); 
               UnitsCopy( &(p->equiv_nodal_load->spRowUnits[1]), d1 );

               /* node i  i > 1*/
               for(i = 2; i <= p->nodes_per_elmt; i++) {
                  for(j = 1; j <= p->dof_per_node; j++) {
                      k  = p->dof_per_node*(i-1) + j;
                      UnitsCopy( &(p->equiv_nodal_load->spRowUnits[k-1]), d1 ); 
                  }
               }

               ZeroUnits( &(p->equiv_nodal_load->spColUnits[0]) );

               free((char *) d1->units_name);
               free((char *) d1);
            }
            break;
       case STRESS_UPDATE:
            break;
       case LOAD_MATRIX:
       case STRESS:          /* compute and print element stresses */

            lint = (int ) 0;
            if(isw == STRESS)     
               ii = no_stress_pts;   /* stress pts         */
            if(isw == LOAD_MATRIX)
               ii = no_integ_pts;    /* guassian integ pts */

            /* Initilize nodal_load */
      
            for(i = 1; i <= dof_per_elmt; i++) {
               nodal_load[i-1][0] = 0.0;
               load[i-1][0]       = 0.0;
            }

            for (i = 1; i<= no_integ_pts*no_integ_pts; i++) {
               sg[i-1] = 0.0; tg[i-1] = 0.0; wg[i-1] = 0.0;
            }

            if(ii*ii != lint)
               pgauss( ii, &lint, sg, tg, wg);

            /* Get units */

            if( UNITS_SWITCH == ON ) {
               if( CheckUnitsType() == SI) {
                   dp_length = DefaultUnits("m");
                   dp_stress = DefaultUnits("Pa");
               } else {
                   dp_length = DefaultUnits("in");
                   dp_stress = DefaultUnits("psi");
               }
            }

            /* Print Element Stresses, Strains and Forces */

            if( isw == STRESS ) {

               printf("\n");
               printf("Stresses in Element No %d\n", p->elmt_no);
               printf("=======================================================================\n");
               printf("Gaussion         x          y     stress-11     stress-22     stress-12 \n");

               if(UNITS_SWITCH == ON) {
                  printf("  Points       %3s        %3s", dp_length->units_name,
                                                          dp_length->units_name );
                  printf("    %10s    %10s    %10s\n", dp_stress->units_name,
                                                       dp_stress->units_name,
                                                       dp_stress->units_name );
               } else {
                  printf("  Points \n");
               }

               printf("=======================================================================\n");
            }

            /* Compute (3x3) Constituitive Material Matrix */

            Mater_matrix = MaterialMatrixPlane( E.value, nu, dFlag );

            /* Gauss Integration */

            for( ii = 1; ii <= lint; ii++) {

                /* Get element shape function and their derivatives */

                shape( sg[ii-1], tg[ii-1], p->coord,shp, &jacobian, p->no_dimen,
                       p->nodes_per_elmt, p->node_connect,0 );

                /* Form [B] matrix */
            
                for(j = 1; j <= p->nodes_per_elmt; j++) { 
                    k = 2*(j-1);
                    B_matrix[0][k]   = shp[0][j-1];
                    B_matrix[1][k+1] = shp[1][j-1];
                    B_matrix[2][k]   = shp[1][j-1];
                    B_matrix[2][k+1] = shp[0][j-1];
                }
            
                /* Calculate strains at guassian integretion pts */

                for(i = 1;i <= 3; i++)
                    strain[i-1][0] = 0.0;

                xx = 0.0; yy = 0.0;
                for(j = 1; j <= p->nodes_per_elmt; j++) {
                    xx = xx + shp[2][j-1] * p->coord[0][j-1].value;
                    yy = yy + shp[2][j-1] * p->coord[1][j-1].value;

                    /*  converting p->displ into a array */
                    
                    for ( k = 1; k <= p->dof_per_node; k++) {
                       j1 = p->dof_per_node*(j-1) + k; 
                       displ[j1-1][0] = p->displ->uMatrix.daa[k-1][j-1];
                    }
                }

                strain = dMatrixMultRep(strain, B_matrix, 3, dof_per_elmt, displ, dof_per_elmt, 1);
                
                /* Compute Stress */

                stress = dMatrixMultRep(stress, Mater_matrix, 3, 3, strain, 3, 1);

                /* Compute equivalent nodal forces for stresses */
                /* F_equiv = integral of [B]^T stress dV        */ 

                if(isw == LOAD_MATRIX) {       
                   for(i = 1; i <= 3; i++)
                       for(j = 1; j <= dof_per_elmt; j++)
                           B_Transpose[j-1][i-1] = B_matrix[i-1][j-1];

                  dv = jacobian*wg[ii-1];

                  for (i = 1; i<= 3; i++)
                       stress[i-1][0] *= dv;

                  nodal_load = dMatrixMultRep(nodal_load, B_Transpose, dof_per_elmt, 3, stress, 3, 1);
                }

                for(i = 1; i <= dof_per_elmt; i++) 
                    p->nodal_loads[i-1].value += nodal_load[i-1][0];

                /* Print stresses  */

                if(isw == STRESS && UNITS_SWITCH == ON) {
                   printf(" %7d %10.4f %10.4f", ii,
                            xx/dp_length->scale_factor,
                            yy/dp_length->scale_factor);
                   printf(" %12.4e  %12.4e  %12.4e\n",
                            stress[0][0]/dp_stress->scale_factor,
                            stress[1][0]/dp_stress->scale_factor,
                            stress[2][0]/dp_stress->scale_factor );
                }

                if(isw == STRESS && UNITS_SWITCH == OFF) {
                   printf(" %7d %10.4f %10.4f", ii, xx, yy);
                   printf(" %12.4e  %12.4e  %12.4e\n",
                            stress[0][0], stress[1][0], stress[2][0] );
                }
            }

            /* Free memory for units used in printing stresses */

            if(p->no_dimen == 2 && isw == STRESS) {
               if(UNITS_SWITCH == ON) {
                  free( dp_length->units_name );
                  free( dp_length );
                  free( dp_stress->units_name );
                  free( dp_stress );
               }
            }

            /* 
             *  =================================================================
             *  Nodal Load Units : units type is determined by the SetUnitsType()                 
             *  =================================================================
             */ 

            if( UNITS_SWITCH == ON ) {
                if( CheckUnitsType() == SI)
                    d1 = DefaultUnits("N");
                else
                    d1 = DefaultUnits("lbf");

               /* Node no 1 */

               UnitsCopy( p->nodal_loads[0].dimen, d1 );
               UnitsCopy( p->nodal_loads[1].dimen, d1 );

               /* Node no > 1 */

               for(i = 2; i <= p->nodes_per_elmt; i++) {    
               for(j = 1; j <= p->dof_per_node; j++) {
                   k = p->dof_per_node*(i-1)+j;
                   UnitsCopy( p->nodal_loads[k-1].dimen, d1 );
               }
               }
               free((char *) d1->units_name);
               free((char *) d1);
            }
            break;
       case STRESS_MATRIX:           /* save element stresses in working array */

            lint = (int ) 0;
            if(isw == STRESS_MATRIX)     
               ii = no_stress_pts;   /* stress pts         */

            /* Initilize nodal_load */
      
            for(i = 1; i <= dof_per_elmt; i++) {
               nodal_load[i-1][0] = 0.0;
               load[i-1][0]       = 0.0;
            }

            for (i = 1; i<= no_integ_pts*no_integ_pts; i++) {
               sg[i-1] = 0.0; tg[i-1] = 0.0; wg[i-1] = 0.0;
            }

            if(ii*ii != lint)
               pgauss( ii, &lint, sg, tg, wg);

            /* Get units */

            if( UNITS_SWITCH == ON ) {
               if( CheckUnitsType() == SI) {
                   dp_length = DefaultUnits("m");
                   dp_stress = DefaultUnits("Pa");
               } else {
                   dp_length = DefaultUnits("in");
                   dp_stress = DefaultUnits("psi");
               }
            }

            /* Compute (3x3) constituitive material matrix */

            Mater_matrix = MaterialMatrixPlane( E.value, nu, dFlag );

            /* Gauss Integration */

            for( ii = 1; ii <= lint; ii++) {

                /* Get element shape function and their derivatives */

                shape( sg[ii-1], tg[ii-1], p->coord,shp, &jacobian, p->no_dimen,
                       p->nodes_per_elmt, p->node_connect,0 );

                /* Form [B] matrix */
            
                for(j = 1; j <= p->nodes_per_elmt; j++) { 
                    k = 2*(j-1);
                    B_matrix[0][k]   = shp[0][j-1];
                    B_matrix[1][k+1] = shp[1][j-1];
                    B_matrix[2][k]   = shp[1][j-1];
                    B_matrix[2][k+1] = shp[0][j-1];
                }
            
                /* Calculate strains at guassian integretion pts */

                for(i = 1;i <= 3; i++)
                    strain[i-1][0] = 0.0;

                xx = 0.0; yy = 0.0;
                for(j = 1; j <= p->nodes_per_elmt; j++) {
                    xx = xx + shp[2][j-1] * p->coord[0][j-1].value;
                    yy = yy + shp[2][j-1] * p->coord[1][j-1].value;

                    /*  converting p->displ into a array */
                    
                    for ( k = 1; k <= p->dof_per_node; k++) {
                       j1 = p->dof_per_node*(j-1) + k; 
                       displ[j1-1][0] = p->displ->uMatrix.daa[k-1][j-1];
                    }
                }

                strain = dMatrixMultRep(strain, B_matrix, 3, dof_per_elmt, displ, dof_per_elmt, 1);
                
                /* Compute Stress */

                stress = dMatrixMultRep(stress, Mater_matrix, 3, 3, strain, 3, 1);

                /* Compute equivalent nodal forces for stresses */
                /* F_equiv = integral of [B]^T stress dV        */ 

                for(i = 1; i <= dof_per_elmt; i++) 
                    p->nodal_loads[i-1].value += nodal_load[i-1][0];

                /* Save element level stresses in working array */
                /* Set column buffer units                      */

                if(ii == 1 ) {
                   UnitsCopy( &(p->stress->spColUnits[0]), dp_length ); 
                   UnitsCopy( &(p->stress->spColUnits[1]), dp_length ); 
                   UnitsCopy( &(p->stress->spColUnits[2]), dp_stress ); 
                   UnitsCopy( &(p->stress->spColUnits[3]), dp_stress ); 
                   UnitsCopy( &(p->stress->spColUnits[4]), dp_stress ); 
                }

                /* Zero out row buffer units */

                ZeroUnits( &(p->stress->spRowUnits[ii-1]) );

                /* Transfer xx and yy coordinates to to working stress matrix */

                p->stress->uMatrix.daa[ii-1][0] = xx;
                p->stress->uMatrix.daa[ii-1][1] = yy;

                /* Transfer internal stresses to working stress matrix */

                p->stress->uMatrix.daa[ii-1][2] = stress[0][0];
                p->stress->uMatrix.daa[ii-1][3] = stress[1][0];
                p->stress->uMatrix.daa[ii-1][4] = stress[2][0];
            }
            break;
       case MASS_MATRIX:
             
            /*
             *  ==================================================================
             *  Compute consistent mass matrix p->type should be -1 for consistent
             *  ==================================================================
             */

            /* Zero contents of integration points arrays */

            ii = no_integ_pts; 
            for (i = 1; i<= no_integ_pts*no_integ_pts; i++) {
               sg[i-1] = 0.0; tg[i-1] = 0.0; wg[i-1] = 0.0;
            }

            /* Get Gauss integration points */

            if(ii*ii != lint)
               pgauss(ii,&lint,sg,tg,wg);

            for(ii=1; ii <= lint; ii++) {

                /* Compute shape functions */

                shape(sg[ii-1],tg[ii-1],p->coord,shp,&jacobian,p->no_dimen,
                      p->nodes_per_elmt,p->node_connect,0);

                dv = density.value * wg[ii-1] * jacobian;

                /* For each node compute db = shape * dv  */

                j1 = 1;
                for(j = 1; j<= p->nodes_per_elmt; j++){
                    w11 = shp[2][j-1] * dv;

                    /* Compute lumped mass (store lumped mass in p->nodal_loads) */

                    p->nodal_loads[j1-1].value = p->nodal_loads[j1-1].value + w11; 

                    /* For each node compute mass matrix ( upper triangular part ) */

                    k1 = j1;
                    for(k = j; k <= p->nodes_per_elmt; k++) {
                        stiff[j1-1][k1-1] += shp[2][k-1] * w11;
                        k1 = k1 + p->dof_per_node;
                    }
                    j1 = j1 + p->dof_per_node;
                } 
                      
                for (i = 1; i <= dof_per_elmt; i++) {
                   for (j = 1; j <= dof_per_elmt; j++) {
                      p->stiff->uMatrix.daa[i-1][j-1] += stiff[i-1][j-1];
                   }
                }
            }

            /* Compute missing parts and lower part by symmetries */

            dof_per_elmt = p->nodes_per_elmt* p->dof_per_node;
            for(j = 1; j <= dof_per_elmt; j++){
                p->nodal_loads[j].value = p->nodal_loads[j-1].value;
                for(k = j; k <= p->dof_per_node; k = k + p->dof_per_node) {
                    p->stiff->uMatrix.daa[j][k]      = p->stiff->uMatrix.daa[j-1][k-1];
                    p->stiff->uMatrix.daa[k-1][j-1]  = p->stiff->uMatrix.daa[j-1][k-1];
                    p->stiff->uMatrix.daa[k][j]      = p->stiff->uMatrix.daa[j-1][k-1];
                }  
            }

            /* Units for Mass Matrix */

            if(UNITS_SWITCH == ON) {
               if( CheckUnitsType() == SI || CheckUnitsType() == SI_US ) {
                   d1 = DefaultUnits("Pa");
                   d1 = DefaultUnits("m");
               } else {
                   d1 = DefaultUnits("psi");
                   d1 = DefaultUnits("in");
               }

               d3 = DefaultUnits("sec");

               /* node no 1 */

               UnitsMultRep( &(p->stiff->spColUnits[0]), d1, d2 );
               UnitsCopy( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );

               UnitsPowerRep( &(p->stiff->spRowUnits[0]), d3, 2.0, NO );
               UnitsCopy( &(p->stiff->spRowUnits[1]), &(p->stiff->spRowUnits[0]) );

              /* node no > 1 */

              for(i = 2; i <= p->nodes_per_elmt; i++) {    
              for(j = 1; j <= p->dof_per_node; j++) {
                  k = p->dof_per_node*(i-1)+j;
                  UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                  UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[0]) );
              }
              }

              free((char *) d1->units_name);
              free((char *) d1);
              free((char *) d2->units_name);
              free((char *) d2);
              free((char *) d3->units_name);
              free((char *) d3);
           }
           break;
       default:
           break;
    }

    /* [d] : free memory and leave */
   
    free(sum_row);

    MatrixFreeIndirectDouble(strain, 3);
    MatrixFreeIndirectDouble(stress, 3);
    MatrixFreeIndirectDouble(body_force, 3);
    MatrixFreeIndirectDouble(load, dof_per_elmt);
    MatrixFreeIndirectDouble(nodal_load, dof_per_elmt);
    MatrixFreeIndirectDouble(displ, dof_per_elmt);
    MatrixFreeIndirectDouble(stiff, dof_per_elmt);
    MatrixFreeIndirectDouble(B_matrix, 3);
    MatrixFreeIndirectDouble(B_Transpose, dof_per_elmt);
    MatrixFreeIndirectDouble(temp_matrix, dof_per_elmt);
    MatrixFreeIndirectDouble(shp, 3);
    
    return(p);
}


/* 
 *  ========================================================================= 
 *  shape() : shape function                                    
 * 
 *     form 4-node quadrilateral shape function                
 *     shape function:                  Ni = shape[2][i]       
 *                                      node no. i = 1, .. 4    
 *     derivatives of shape functions:  dNi/d(ss) = shape[0][i] 
 *                                      dNi/d(tt) = shape[1][i] 
 * 
 *  Input  :  double ss        -- 
 *            double tt        -- 
 *            QUANTITY **coord -- 
 *  Output : 
 *  ========================================================================= 
 */ 

int shape(ss,tt,coord,shp,jacobian,no_dimen,nodes_per_elmt,node_connect,flg)
double  ss, tt, **shp, *jacobian, *node_connect;
QUANTITY   **coord;
int       no_dimen;
int nodes_per_elmt;
int            flg;
{
double  s[4], t[4], xs[2][2], sx[2][2], tp;
double  **shp_temp;
int i, j, k;

    /* [a] : Initialize arrays */
 
    shp_temp = MatrixAllocIndirectDouble(2, 4);
    s[0] =  0.5; s[1] = -0.5;
    s[2] = -0.5; s[3] =  0.5;
    t[0] =  0.5; t[1] =  0.5;
    t[2] = -0.5; t[3] = -0.5;

    /* [b] : form shape functions */

    switch(nodes_per_elmt) { 
        case 3: case 4:
           for(i = 1; i <= 4; i++){
               shp[2][i-1]      = (0.5 + s[i-1] * ss) * ( 0.5 + t[i-1] * tt); 
               shp_temp[0][i-1] = s[i-1] * (0.5 + t[i-1] * tt);                
               shp_temp[1][i-1] = t[i-1] * (0.5 + s[i-1] * ss);
           }

           /* Form triangle by adding third and fourth node together */

           if(nodes_per_elmt == 3) { 
               shp[2][2] = shp[2][2] + shp[2][3];
               for(i = 0; i <= 1; i++)
                   shp_temp[i][2] = shp_temp[i][2] + shp_temp[i][3];
           }

           /* Construct jacobian matrix and its determinant */

           for(i = 1; i <= no_dimen; i++)      /* no_dimen = 2 */
           for(j = 1; j <= no_dimen; j++) {
               xs[i-1][j-1] = 0.0;
               for(k = 1; k <= nodes_per_elmt; k++)
                   xs[i-1][j-1] = xs[i-1][j-1] + coord[i-1][k-1].value*shp_temp[j-1][k-1];
           }

           *jacobian = xs[0][0] * xs[1][1] - xs[0][1] *xs[1][0];

           if(flg == TRUE)
              return;

           /* Compute Jacobian inverse matrix */

           sx[0][0] = xs[1][1]/ *jacobian;
           sx[1][1] = xs[0][0]/ *jacobian;
           sx[0][1] = - xs[0][1]/ *jacobian;
           sx[1][0] = - xs[1][0]/ *jacobian;

           /* Save dNi/dx, dNi/dy into shp[2][node] */
  
           for(i = 1; i <= nodes_per_elmt; i++){
               shp[0][i-1] = shp_temp[0][i-1]*sx[0][0] + shp_temp[1][i-1]*sx[0][1];
               shp[1][i-1] = shp_temp[0][i-1]*sx[1][0] + shp_temp[1][i-1]*sx[1][1];
           }
           break;
      default:
           break;
   }

   MatrixFreeIndirectDouble(shp_temp, 2);
   return;
}


/*
 *  ============================================================================= 
 *  gauss() : Gauss Integration 
 *
 *  Input :
 *  Output :
 *  ============================================================================= 
 */

#ifdef __STDC__
gauss( double *sg, double *ag, int lt)
#else
gauss(sg, ag, lt)
double *sg, *ag;
int lt;
#endif
{
double t;
int    i;

    switch(lt) {
        case 1:
             sg[1] = 0.0;
             ag[1] = 2.0;
             break;
        case 2:
             sg[1] = -1/sqrt(3.);
             sg[2] =  -sg[1];
             ag[1] = 1.0; 
             ag[2] = 1.0; 
             break;
        case 3:
             sg[1] = -sqrt(0.6);
             sg[2] = 0.0;
             sg[3] = -sg[1];

             ag[1] = 5.0/9.0; 
             ag[2] = 8.0/9.0; 
             ag[3] = ag[1]; 
             break;
        case 4:
             t = sqrt(4.8);
             sg[1] =  sqrt((3+t)/7);
             sg[2] =  sqrt((3-t)/7);
             sg[3] = -sg[2];
             sg[4] = -sg[1];

             t  = 1/3/t;
             ag[1] = 0.5-t; 
             ag[2] = 0.5+t; 
             ag[3] = ag[2]; 
             ag[4] = ag[1]; 
             break;
        default:
             break;
    }

    return (1);
}

/*
 *  ============================================================================= 
 *  pgauss() : 
 * 
 *  Input  : int no_integ_pts --  no of gaussian integ pts in one-direction.
 *  Output : lint             --  no_integ_pts*no_integ_pts  
 *         :   sg             -- coordinates in xi direction 
 *         :   tg             -- coordinates in eta direction
 *         :   wg             -- weighting coefficients
 *  ============================================================================= 
 */

int pgauss(l, lint, r, z, w)
int l, *lint;
double r[], z[], w[];
{
int i, j, k;
double g, h;
double g4[5], h4[5], lr[10], lz[10], lw[10];

lr[1] = lr[4] =lr[8] = -1;
lr[2] = lr[3] =lr[6] =  1;
lr[5] = lr[7] =lr[9] =  0;

lz[1] = lz[2] =lz[5] = -1;
lz[3] = lz[4] =lz[7] =  1;
lz[6] = lz[8] =lz[9] =  0;

lw[3] = lw[4] = lw[1] = lw[2] = 25;
lw[5] = lw[6] = lw[7] = lw[8] = 40;
lw[9] = 64;

    if(l < 0) {
       *lint = 3;
       g = sqrt((double) 1.0/3.0);
       h = sqrt((double) 2.0/3.0);
       r[1] = -h; r[2] =  h; r[3] = 0;
       z[1] = -g; z[2] = -g; z[3] = g;
       w[1] =  1; w[2] =  1; w[3] = 2;

       return (1);
    }

    *lint = l*l;
    switch (l) {
	case 1: /*  1 x 1 integration */
	     r[0] = 0.0;
	     z[0] = 0.0;
	     w[0] = 4.0;
	     break;
        case 2: /*  2 x 2 integration */
             g = 1.0/sqrt((double) 3.0);
             for(i = 1; i<= 4; i++) {
                 r[i-1] = g * lr[i];
                 z[i-1] = g * lz[i];
                 w[i-1] = 1.0;
             }
             break;
        case 3: /* 3 x 3 integration */
             g = sqrt((double) 0.6);
             h = 1.0/81.0;
             for(i = 1; i<= 9; i++) {
                 r[i-1] = g * lr[i];
                 z[i-1] = g * lz[i];
                 w[i-1] = h * lw[i];
             }
             break;
        case 4: /* 4 x 4 integration */
             g = sqrt((double) 4.8);
             h = sqrt((double) 30.0)/36;
             g4[1] = sqrt((double) (3+g)/7.0);
             g4[4] = -g4[1];
             g4[2] = sqrt((double) (3-g)/7.0);
             g4[3] = -g4[2];
             h4[1] = 0.5 - h;
             h4[2] = 0.5 + h;
             h4[3] = 0.5 + h;
             h4[4] = 0.5 - h;
             i = 0;
             for(j = 1; j<= 4; j++) {
                 for(k = 1; k<= 4; k++) {
                     i = i +1;
                     r[i-1] = g4[k];
                     z[i-1] = g4[j];
                     w[i-1] = h4[j]* h4[k];
                 }
             }
             break;
    }

    return(1);
}

/*
 *  ===============================================================
 *  MaterialMatrixPlane() : Compute (3x3) constituitive matrix.
 *
 *  Note : Plane Stress : dFlag = 1
 *         Plane Strain : dFlag = 2
 *
 *  Input  : double     E -- Moduus of elasticity.
 *         : double    nu -- Poisson's ratio
 *         : double dFlag -- Flag for Plane Strss/Plane Strain Cases
 *  Output : double **m1  -- Pointer to constituitive matrix.
 *  ===============================================================
 */

double **MaterialMatrixPlane( double E, double nu , double dFlag ) {
double **m1;
double temp;

    m1 = MatrixAllocIndirectDouble(3, 3);
    
    temp =  E*(1 + ( 1-dFlag)*nu)/(1.0+nu)/(1.0 - dFlag*nu);

    m1[0][0] = m1[1][1] = temp;
    m1[0][1] = m1[1][0] = (nu/(1+(1-dFlag)*nu))*temp; 
    m1[2][2] = E/2.0/(1.0+nu); 

    m1[0][2] = m1[2][0] = m1[1][2] = m1[2][1] = 0.0;

    return(m1);
}


/*
 *  ==========================================================================
 *  print_property_psps() : print PLANE_STRAIN/PLANE_STRESS element properties
 *  ==========================================================================
 */

#ifdef __STDC__
void print_property_psps(EFRAME *frp, int i)
#else
void print_property_psps(frp, i)
EFRAME    *frp;
int          i;                 /* elmt_attr_no */
#endif
{
int     UNITS_SWITCH;
ELEMENT_ATTR    *eap;

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
        if( eap->work_material[0].value != 0.0 ) {
           printf("             ");
           printf("         : Young's Modulus =  E = %16.3e\n",
                            eap->work_material[0].value);
        }
        if( eap->work_material[2].value != 0.0 ) {
           printf("             ");
           printf("         : Yielding Stress = fy = %16.3e\n",
                            eap->work_material[2].value);
        }
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's ratio = nu = %16.3e   \n", eap->work_material[4].value);
        }
        if( eap->work_material[0].value != 0.0 ) {
           printf("             ");
           printf("         : Density         = %16.3e\n",
                            eap->work_material[5].value);
        }
        if( eap->work_section[2].value != 0.0 ) {
           printf("             ");
           printf("         : Inertia Izz     = %16.3e\n",
                            eap->work_section[2].value);
        }
        if( eap->work_section[10].value != 0.0 ) {
           printf("             ");
           printf("         : Area            = %16.3e\n",
                            eap->work_section[10].value);
        }
        break;
        default:
        break;
     }
}
