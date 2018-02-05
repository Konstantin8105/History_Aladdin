/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_fiber3d.c : Linear/Nonlinear 3D Fiber Element
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
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "fe_functions.h"
#include "elmt.h"

/*
#define  DEBUG
*/

#ifdef __STDC__
ARRAY *elmt_fiber_3d(ARRAY *p, int isw)
#else
ARRAY *elmt_fiber_3d(p, isw)
ARRAY *p;
int     isw;
#endif
{
static int       no_integ_pt, no_section, no_fiber, no_shear, elmt_no;
static QUANTITY  elmt_length;

int         ii, jj, kk, sec, ifib;
int         UNITS_SWITCH, UnitsType;
double      cx, cy, cz, cxz, xx, yy, zz, temp1, temp2, temp3;
double      **rot, mbar, rT;
MATRIX      *temp_m1, *temp_m2;
DIMENSIONS  *dp_length, *dp_force, *dp_stress, *dimen;

double      G, J, bf, depth;
double      *abscissas, *weights;
double      xi;
double      scale;
MATRIX      *Q, *q;
MATRIX      *F, *K, *Ke;
MATRIX      *L, *Ltrans, *R, *Rtrans;
MATRIX      *kx, *fx;
MATRIX      *bx, *bxtrans;
FIBER_ATTR  *fiber;

    UNITS_SWITCH = CheckUnits();
    UnitsType    = CheckUnitsType();

    switch( isw ) {

       case PROPTY :
          no_fiber    = p->fiber_ptr->no_fiber;
          no_shear    = p->fiber_ptr->no_shear;
          no_integ_pt = p->integ_ptr->integ_pts;
          no_section  = no_integ_pt + 2;
          elmt_no     = p->elmt_no;
          cx = p->coord[0][1].value - p->coord[0][0].value;
          cy = p->coord[1][1].value - p->coord[1][0].value;
          cz = p->coord[2][1].value - p->coord[2][0].value;
          elmt_length.value = sqrt( cx*cx + cy*cy + cz*cz );
          if( UNITS_SWITCH == ON ) {
             if( UnitsType == SI )
                dp_length = DefaultUnits("m");
             else if( UnitsType == US )
                dp_length = DefaultUnits("in");
             elmt_length.dimen = (DIMENSIONS *)MyCalloc(1, sizeof(DIMENSIONS));
             UnitsCopy( elmt_length.dimen, dp_length );
             free((char *) dp_length->units_name);
             free((char *) dp_length);
          }
          break;  /* end of case PROPTY */

       case STIFF :
          if( UNITS_SWITCH == ON )
             SetUnitsOff();

          fiber     = p->fiber_ptr->fiber;
          abscissas = (double *)MyCalloc( no_section, sizeof(double) );
          weights   = (double *)MyCalloc( no_section, sizeof(double) );
          Gauss_Lobatto( abscissas, weights, no_integ_pt );

          kx = MatrixAllocIndirect( "kx", DOUBLE_ARRAY, 5, 5 );
          bx = MatrixAllocIndirect( "bx", DOUBLE_ARRAY, 5, 5 );

          for( sec=0 ; sec < no_section; ++sec ) {
             xi = elmt_length.value/2.0*(abscissas[sec]+1);
             scale = weights[sec]*elmt_length.value/2.0;
             Force_Interpolation_Matrix_3d( bx, xi, elmt_length );
             Section_Tangent_Stiffness_3d( kx, fiber, no_fiber, no_shear, sec, elmt_no );
             fx = MatrixInverse( kx );
             bxtrans = MatrixTranspose( bx );
             temp_m1 = MatrixMult( bxtrans, fx );
             temp_m2 = MatrixMult( temp_m1, bx );
             MatrixFree( fx );
             MatrixFree( temp_m1 );
             MatrixFree( bxtrans );
             if( sec == 0 )
                F = MatrixScale( temp_m2, scale );
             else {
                temp_m1 = MatrixScale( temp_m2, scale );
                MatrixAddReplace( F, temp_m1 );
                MatrixFree( temp_m1 );
             }
             MatrixFree( temp_m2 );
          }
          K = MatrixInverse( F );

          free((char *) abscissas);
          free((char *) weights);
          MatrixFree( kx );
          MatrixFree( bx );
          MatrixFree( F );

          /* Rigid body rotation transformation */
          /* Ke = Trans(R)*K*R , element stiffness in local coordinate system */
          R = Rigid_Body_Rotation_3d( elmt_length );
          Rtrans = MatrixTranspose( R );
          temp_m1 = MatrixMult( Rtrans, K );
          Ke = MatrixMult( temp_m1, R );
          MatrixFree( temp_m1 );

          /* Calculate torsional stiffness */
          /* Assuming plane remain plane, so twisting stiffness is independent to others */
          if ( (G=p->work_material[1].value) <= 0.0 )
             G = p->work_material[0].value/2.0/(1+p->work_material[4].value);
          J = p->work_section[12].value;
          bf = p->work_section[7].value;
          depth = p->work_section[9].value;

          /* If J value not input, J calculated based on rectangular */
          /* section of size (bf x depth)                            */

          if(J == 0.0 ) {
             if(bf == 0.0 || depth == 0.0){
                printf("WARNING >> Must give 'J' or ('width' & 'depth') to calculate stiffness");
                exit(1);
             }
             /* Check bf < depth & cal J */
             if(bf < depth ) 
                J = (1.0-0.63*bf/depth)*bf*bf*bf*depth/3.0;
             else
                J = (1.0-0.63*depth/bf)*depth*depth*depth*bf/3.0;
          }

          Ke->uMatrix.daa[3][3] = G*J/elmt_length.value;
          Ke->uMatrix.daa[3][9] = -Ke->uMatrix.daa[3][3];
          Ke->uMatrix.daa[9][3] =  Ke->uMatrix.daa[3][9];
          Ke->uMatrix.daa[9][9] =  Ke->uMatrix.daa[3][3];

          /* Transform local coordinate system to global coordinate system */

          cx = (p->coord[0][1].value - p->coord[0][0].value)/elmt_length.value;
          cy = (p->coord[1][1].value - p->coord[1][0].value)/elmt_length.value;
          cz = (p->coord[2][1].value - p->coord[2][0].value)/elmt_length.value;
          cxz = sqrt( cx*cx + cz*cz );
          if( cx != 1.0 ) {
             if( cxz == 0.0 ) {
                /* Ke' = Ke*T */
                for( ii=0 ; ii < 12 ; ++ii ) {
                   for( jj=0 ; jj < 12 ; jj=jj+3 ) {
                      temp1 = Ke->uMatrix.daa[ii][jj];
                      temp2 = Ke->uMatrix.daa[ii][jj+1];
                      Ke->uMatrix.daa[ii][jj]   = -temp2*cy;
                      Ke->uMatrix.daa[ii][jj+1] =  temp1*cy;
                   }
                }
                /* Ke(global) = Trans(T)*Ke*T = Trans(T)*Ke' */
                for( ii=0 ; ii < 12 ; ii=ii+3 ) {
                   for( jj=0 ; jj < 12 ; ++jj ) {
                      temp1 = Ke->uMatrix.daa[ii][jj];
                      temp2 = Ke->uMatrix.daa[ii+1][jj];
                      Ke->uMatrix.daa[ii][jj]   = -temp2*cy;
                      Ke->uMatrix.daa[ii+1][jj] =  temp1*cy;
                   }
                }
             }
             else {
                /* Ke' = Ke*T */
                for( ii=0 ; ii < 12 ; ++ii ) {
                   for( jj=0 ; jj < 12 ; jj=jj+3 ) {
                      temp1 = Ke->uMatrix.daa[ii][jj];
                      temp2 = Ke->uMatrix.daa[ii][jj+1];
                      temp3 = Ke->uMatrix.daa[ii][jj+2];
                      Ke->uMatrix.daa[ii][jj]   = temp1*cx - temp2*cx*cy/cxz - temp3*cz/cxz;
                      Ke->uMatrix.daa[ii][jj+1] = temp1*cy + temp2*cxz;
                      Ke->uMatrix.daa[ii][jj+2] = temp1*cz - temp2*cy*cz/cxz + temp3*cx/cxz;
                   }
                }
                /* Ke(global) = Trans(T)*Ke*T = Trans(T)*Ke' */
                for( ii=0 ; ii < 12 ; ii=ii+3 ) {
                   for( jj=0 ; jj < 12 ; ++jj ) {
                      temp1 = Ke->uMatrix.daa[ii][jj];
                      temp2 = Ke->uMatrix.daa[ii+1][jj];
                      temp3 = Ke->uMatrix.daa[ii+2][jj];
                      Ke->uMatrix.daa[ii][jj]   = temp1*cx - temp2*cx*cy/cxz - temp3*cz/cxz;
                      Ke->uMatrix.daa[ii+1][jj] = temp1*cy + temp2*cxz;
                      Ke->uMatrix.daa[ii+2][jj] = temp1*cz - temp2*cy*cz/cxz + temp3*cx/cxz;
                   }
                }
             }
          }

          /* Copy stiffness matrix to p array */
          for(ii = 1; ii <= p->stiff->iNoRows; ii++)
             for(jj = 1; jj <= p->stiff->iNoColumns; jj++)
                p->stiff->uMatrix.daa[ii-1][jj-1] = Ke->uMatrix.daa[ii-1][jj-1];

          MatrixFree( K );
          MatrixFree( Rtrans );
          MatrixFree( R );
          MatrixFree( Ke );

          /* Assign units to p array stiffness */
          if( UNITS_SWITCH == ON ) {
             SetUnitsOn();
             switch( UnitsType ) {
                case SI:
                case SI_US:
                   dp_force  = DefaultUnits("N");
                   dp_length = DefaultUnits("m");
                   break;
                case US:
                   dp_force  = DefaultUnits("lbf");
                   dp_length = DefaultUnits("in");
                   break;
             }

             ZeroUnits( &(p->stiff->spRowUnits[0]) );
             ZeroUnits( &(p->stiff->spRowUnits[1]) );
             ZeroUnits( &(p->stiff->spRowUnits[2]) );
             UnitsCopy( &(p->stiff->spRowUnits[3]), dp_length );
             UnitsCopy( &(p->stiff->spRowUnits[4]), dp_length );
             UnitsCopy( &(p->stiff->spRowUnits[5]), dp_length );
             UnitsCopy( &(p->stiff->spRowUnits[6]),  &(p->stiff->spRowUnits[0]) );
             UnitsCopy( &(p->stiff->spRowUnits[7]),  &(p->stiff->spRowUnits[1]) );
             UnitsCopy( &(p->stiff->spRowUnits[8]),  &(p->stiff->spRowUnits[2]) );
             UnitsCopy( &(p->stiff->spRowUnits[9]),  &(p->stiff->spRowUnits[3]) );
             UnitsCopy( &(p->stiff->spRowUnits[10]), &(p->stiff->spRowUnits[4]) );
             UnitsCopy( &(p->stiff->spRowUnits[11]), &(p->stiff->spRowUnits[5]) );

             UnitsDivRep( &(p->stiff->spColUnits[0]), dp_force, dp_length, NO );
             UnitsCopy(   &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]) );
             UnitsCopy(   &(p->stiff->spColUnits[2]), &(p->stiff->spColUnits[0]) );
             UnitsCopy(   &(p->stiff->spColUnits[3]), dp_force );
             UnitsCopy(   &(p->stiff->spColUnits[4]), dp_force );
             UnitsCopy(   &(p->stiff->spColUnits[5]), dp_force );
             UnitsCopy( &(p->stiff->spColUnits[6]),  &(p->stiff->spColUnits[0]) );
             UnitsCopy( &(p->stiff->spColUnits[7]),  &(p->stiff->spColUnits[1]) );
             UnitsCopy( &(p->stiff->spColUnits[8]),  &(p->stiff->spColUnits[2]) );
             UnitsCopy( &(p->stiff->spColUnits[9]),  &(p->stiff->spColUnits[3]) );
             UnitsCopy( &(p->stiff->spColUnits[10]), &(p->stiff->spColUnits[4]) );
             UnitsCopy( &(p->stiff->spColUnits[11]), &(p->stiff->spColUnits[5]) );

             free((char *) dp_force->units_name);
             free((char *) dp_length->units_name);
             free((char *) dp_force);
             free((char *) dp_length);
          } /* end of units on/off for case STIFF */

          break;   /* end of case STIFF */

       case LOAD_MATRIX:
       case STRESS:

          if( UNITS_SWITCH == ON )
             SetUnitsOff();

          if( p->elmt_state == 0 ) {  /* Never pass the ElmtStateDet() */

             fiber     = p->fiber_ptr->fiber;
             abscissas = (double *)MyCalloc( no_section, sizeof(double) );
             weights   = (double *)MyCalloc( no_section, sizeof(double) );
             Gauss_Lobatto( abscissas, weights, no_integ_pt );

             kx = MatrixAllocIndirect( "kx", DOUBLE_ARRAY, 5, 5 );
             bx = MatrixAllocIndirect( "bx", DOUBLE_ARRAY, 5, 5 );

             for( sec=0 ; sec < no_section; ++sec ) {
                xi = elmt_length.value/2.0*(abscissas[sec]+1);
                scale = weights[sec]*elmt_length.value/2.0;
                Force_Interpolation_Matrix_3d( bx, xi, elmt_length );
                Section_Tangent_Stiffness_3d( kx, fiber, no_fiber, no_shear, sec, elmt_no );
                fx = MatrixInverse( kx );
                bxtrans = MatrixTranspose( bx );
                temp_m1 = MatrixMult( bxtrans, fx );
                temp_m2 = MatrixMult( temp_m1, bx );
                MatrixFree( fx );
                MatrixFree( temp_m1 );
                MatrixFree( bxtrans );
                if( sec == 0 )
                   F = MatrixScale( temp_m2, scale );
                else {
                   temp_m1 = MatrixScale( temp_m2, scale );
                   MatrixAddReplace( F, temp_m1 );
                   MatrixFree( temp_m1 );
                }
                MatrixFree( temp_m2 );
             }
             K = MatrixInverse( F );

             free((char *) abscissas);
             free((char *) weights);
             MatrixFree( kx );
             MatrixFree( bx );
             MatrixFree( F );

             /* put displacement matrix to one column form */
             /* p->displ=[px1,px2; py1,py2; pz1,pz2; rx1,rx2; ry1,ry2; rz1,rz2]6x2    */
             /* => temp_m1=pt=[px1;py1;pz1;rx1;ry1;rz1; px2;py2;pz2;rx2;ry2;rz2]12x1 */
             temp_m1  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 12, 1 );
             for( ii=0 ; ii < p->dof_per_node ; ++ii )
                for( jj=0 ; jj < p->nodes_per_elmt ; ++jj )
                   temp_m1->uMatrix.daa[ii+jj*p->dof_per_node][0] = p->displ->uMatrix.daa[ii][jj];

             L = Element_Transformation_3d( p->coord, elmt_length );
             q = MatrixMult( L, temp_m1 );
             Q = MatrixMult( K, q );
             MatrixFree( temp_m1 );
          }
          else if( p->elmt_state == 1 ) {
             L = Element_Transformation_3d( p->coord, elmt_length );
             Q = p->Q_saved;
          }

          /* Calculate torsional stiffness */
          /* Assuming plane remain plane, so twisting stiffness is independent to others */
          if ( (G=p->work_material[1].value) <= 0.0 )
             G = p->work_material[0].value/2.0/(1+p->work_material[4].value);
          J = p->work_section[12].value;
          bf = p->work_section[7].value;
          depth = p->work_section[9].value;

          /* If J value not input, J calculated based on rectangular */
          /* section of size (bf x depth)                            */

          if(J == 0.0 ) {
             if(bf == 0.0 || depth == 0.0){
                printf("WARNING >> Must give 'J' or ('width' & 'depth') to calculate stiffness");
                exit(1);
             }
             /* Check bf < depth & cal J */
             if(bf < depth ) 
                J = (1.0-0.63*bf/depth)*bf*bf*bf*depth/3.0;
             else
                J = (1.0-0.63*depth/bf)*depth*depth*depth*bf/3.0;
          }

          /* Transform local coordinate system to global coordinate system */

          cx = (p->coord[0][1].value - p->coord[0][0].value)/elmt_length.value;
          cy = (p->coord[1][1].value - p->coord[1][0].value)/elmt_length.value;
          cz = (p->coord[2][1].value - p->coord[2][0].value)/elmt_length.value;

          /* independent twisting stiffness */
          /* and element twising force in local coordinate */

          temp1 = G*J/elmt_length.value;
          temp2 = temp1*( cx*(p->displ->uMatrix.daa[3][1]-p->displ->uMatrix.daa[3][0]) + cy*(p->displ->uMatrix.daa[4][1]-p->displ->uMatrix.daa[4][0]) + cz*(p->displ->uMatrix.daa[5][1]-p->displ->uMatrix.daa[5][0]) );

          if( isw == LOAD_MATRIX ) {
             Ltrans = MatrixTranspose( L );
             temp_m1 = MatrixMult( Ltrans, Q );   /* element nodal force in global coord. */

             /* assign twisting force */
             temp_m1->uMatrix.daa[3][0]  = temp_m1->uMatrix.daa[3][0]  - cx*temp2;
             temp_m1->uMatrix.daa[4][0]  = temp_m1->uMatrix.daa[4][0]  - cy*temp2;
             temp_m1->uMatrix.daa[5][0]  = temp_m1->uMatrix.daa[5][0]  - cz*temp2;
             temp_m1->uMatrix.daa[9][0]  = temp_m1->uMatrix.daa[9][0]  + cx*temp2;
             temp_m1->uMatrix.daa[10][0] = temp_m1->uMatrix.daa[10][0] + cy*temp2;
             temp_m1->uMatrix.daa[11][0] = temp_m1->uMatrix.daa[11][0] + cz*temp2;
             MatrixFree( Ltrans );
          }
          if( isw == STRESS ) {
             R = Rigid_Body_Rotation_3d( elmt_length );
             Rtrans = MatrixTranspose( R );
             temp_m1 = MatrixMult( Rtrans, Q );   /* element nodal force in local coord. */

             /* assign twisting force */
             temp_m1->uMatrix.daa[3][0] = -temp2;
             temp_m1->uMatrix.daa[9][0] =  temp2;
             MatrixFree( R );
             MatrixFree( Rtrans );
          }
          MatrixFree( L );

          /* Assign force values */
          for( ii=0 ; ii < p->size_of_stiff ; ++ii )
             p->nodal_loads[ii].value = temp_m1->uMatrix.daa[ii][0];

          if( p->elmt_state == 0 ) {
             MatrixFree( K );
             MatrixFree( Q );
             MatrixFree( q );
          }
          MatrixFree( temp_m1 );

          /* Assign force units */
          if( UNITS_SWITCH == ON ) {
             SetUnitsOn();
             switch( UnitsType ) {
                case SI:
                case SI_US:
                   dp_force  = DefaultUnits("N");
                   dp_length = DefaultUnits("m");
                   break;
                case US:
                   dp_force  = DefaultUnits("lbf");
                   dp_length = DefaultUnits("in");
                   break;
             }
             dimen = UnitsMult( dp_force, dp_length );

             for( ii=1 ; ii <= 3 ; ++ii ) {
                UnitsCopy( p->nodal_loads[ii-1].dimen, dp_force );
                UnitsCopy( p->nodal_loads[ii-1+p->dof_per_node].dimen, dp_force );
                UnitsCopy( p->nodal_loads[ii-1+3].dimen, dimen );
                UnitsCopy( p->nodal_loads[ii-1+3+p->dof_per_node].dimen, dimen );
             }
          }

          if(isw == STRESS ) {
             xx  = 0.5*(p->coord[0][0].value + p->coord[0][1].value);   /* xx = 0.5(x1+x2) */
             yy  = 0.5*(p->coord[1][0].value + p->coord[1][1].value);   /* yy = 0.5(y1+y2) */
             zz  = 0.5*(p->coord[2][0].value + p->coord[2][1].value);   /* zz = 0.5(z1+z2) */
             printf("\n");
             printf("Elmt No %3d : ", p->elmt_no);
             switch( UNITS_SWITCH ) {
                case ON:
                   printf("Coords (X,Y,Z)= (%8.3f %s,%8.3f %s,%8.3f %s)\n\n",
                         xx/p->coord[0][0].dimen->scale_factor, p->coord[0][0].dimen->units_name,
                         yy/p->coord[1][0].dimen->scale_factor, p->coord[1][0].dimen->units_name,
                         zz/p->coord[2][0].dimen->scale_factor, p->coord[2][0].dimen->units_name);
                   /* node_i */
                   printf("            Fx1 = %13.5e %s\t Fy1 = %13.5e %s\t Fz1 = %13.5e %s\n",
                         p->nodal_loads[0].value/dp_force->scale_factor, dp_force->units_name,
                         p->nodal_loads[1].value/dp_force->scale_factor, dp_force->units_name,
                         p->nodal_loads[2].value/dp_force->scale_factor, dp_force->units_name);
                   printf("            Mx1 = %13.5e %s\t My1 = %13.5e %s\t Mz1 = %13.5e %s\n",
                         p->nodal_loads[3].value/dimen->scale_factor, dimen->units_name,
                         p->nodal_loads[4].value/dimen->scale_factor, dimen->units_name,
                         p->nodal_loads[5].value/dimen->scale_factor, dimen->units_name);
                   printf("\n");
                   /* node_j */
                   printf("            Fx2 = %13.5e %s\t Fy2 = %13.5e %s\t Fz2 = %13.5e %s\n",
                         p->nodal_loads[6].value/dp_force->scale_factor, dp_force->units_name,
                         p->nodal_loads[7].value/dp_force->scale_factor, dp_force->units_name,
                         p->nodal_loads[8].value/dp_force->scale_factor, dp_force->units_name);
                   printf("            Mx2 = %13.5e %s\t My2 = %13.5e %s\t Mz2 = %13.5e %s\n",
                         p->nodal_loads[9].value/dimen->scale_factor, dimen->units_name,
                         p->nodal_loads[10].value/dimen->scale_factor, dimen->units_name,
                         p->nodal_loads[11].value/dimen->scale_factor, dimen->units_name);
                   printf("\n");

                   printf("            Axial Force : x-direction = %13.5e %s \n",
                        -p->nodal_loads[0].value/dp_force->scale_factor, dp_force->units_name);
                   printf("            Shear Force : y-direction = %13.5e %s \n",
                         p->nodal_loads[1].value/dp_force->scale_factor, dp_force->units_name);
                   printf("                        : z-direction = %13.5e %s \n",
                         p->nodal_loads[2].value/dp_force->scale_factor, dp_force->units_name);
                   printf("\n"); 

                   free((char *) dp_force->units_name);
                   free((char *) dp_length->units_name);
                   free((char *) dimen->units_name);
                   free((char *) dp_force);
                   free((char *) dp_length);
                   free((char *) dimen);
                   break;
                case OFF:
                   printf("Coords (X,Y,Z)= (%8.3f,%8.3f,%8.3f)\n\n", xx, yy, zz);
                   /* node_i */
                   printf("            Fx1 = %13.5e\t Fy1 = %13.5e\t Fz1 = %13.5e\n",
                         p->nodal_loads[0].value, p->nodal_loads[1].value, p->nodal_loads[2].value);
                   printf("            Mx1 = %13.5e\t My1 = %13.5e\t Mz1 = %13.5e\n",
                         p->nodal_loads[3].value, p->nodal_loads[4].value, p->nodal_loads[5].value);
                   printf("\n"); 
                   /* node_j */
                   printf("            Fx2 = %13.5e\t Fy2 = %13.5e\t Fz2 = %13.5e\n",
                         p->nodal_loads[6].value, p->nodal_loads[7].value, p->nodal_loads[8].value);
                   printf("            Mx2 = %13.5e\t My2 = %13.5e\t Mz2 = %13.5e\n",
                         p->nodal_loads[9].value, p->nodal_loads[10].value, p->nodal_loads[11].value);
                   printf("\n"); 

                   printf("            Axial Force : x-direction = %13.5e \n", -p->nodal_loads[0].value);
                   printf("            Shear Force : y-direction = %13.5e \n",  p->nodal_loads[1].value);
                   printf("                        : z-direction = %13.5e \n",  p->nodal_loads[2].value);
                   printf("\n"); 
                   break;
                default:
                   break;
             }
            }
            break;  /* end of case STRESS, LOAD_MATRIX */

       case STRESS_MATRIX:

            /* save element nodal forces in working array */

            printf("--------------------------------------------------- \n");
            printf("In elmt_fiber3d.c : case STRESS_MATRIX              \n");
            printf("--------------------------------------------------- \n");
            printf("*** This block of code has not yet been implemented \n");
            printf("*** See elmt_fiber2d.c for details in 2-d case      \n");
            printf("*** Terminating Program Execution                   \n");
            exit(1);

            break;

       case MASS_MATRIX:
          if( p->work_section[6].value != 0.0 )    /* mbar = weight/gravity */
             mbar = p->work_section[6].value/9.80665;
          /* mbar = density * area */
          else if( (p->work_material[5].value > 0.0) && (p->work_section[10].value > 0.0) )
             mbar = p->work_material[5].value * p->work_section[10].value;
          else
             FatalError("\nError in input: Need density value to calculate mass matrix\n",(char *)NULL);

          cx = p->coord[0][1].value - p->coord[0][0].value;           /* Cos Term */
          cy = p->coord[1][1].value - p->coord[1][0].value;           /* Sin Term */
          cz = p->coord[2][1].value - p->coord[2][0].value;           /* Tan Term */
          p->length.value = elmt_length.value;

          /* T matrix is made here: 12*12 size */

          rot = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
          rot = (double **) tmat(rot, 6, p); 

          /* Assemble Mass Matrix */
 
	  /* Calculate radius of gyration , rT  --  m  ,  in   */
          /* original version      :  rT = p->length.value/ 1.414; */
          bf = p->work_section[7].value;
          depth = p->work_section[9].value;
          J = p->work_section[12].value;
          rT = p->work_section[13].value;

          /* If J value not input, J calculated based on rectangular */
          /* section of size (bf x depth)                            */

          if(J == 0.0 ) {
             if(bf == 0.0 || depth == 0.0){
                printf("WARNING >> Must give 'J' or ('width' & 'depth') to calculate stiffness");
                exit(1);
             }
             /* Check bf < depth & cal J */
             if(bf < depth ) 
                J = (1.0-0.63*bf/depth)*bf*bf*bf*depth/3.0;
             else
                J = (1.0-0.63*depth/bf)*depth*depth*depth*bf/3.0;
          }
          if( rT==0 && p->work_section[10].value!=0 )
             rT = sqrt( J / p->work_section[10].value );

          switch( p->type ) {  /* mass type */
             case LUMPED :
                /* use lumped mass matrix of FRAME_3D element */
                p->stiff = beamms3d(p, p->stiff,p->type, mbar, elmt_length.value, rT, rot, p->size_of_stiff, p->dof_per_node);
                MatrixFreeIndirectDouble(rot, p->size_of_stiff);
                break;

             case CONSISTENT :
             default :
                FatalError("FIBER_3D element only support LUMPED mass", (char *)NULL);
                break;
          }

          break;

       default:
          break;
    }

    return(p);
}

/*
 *  ==================================
 *  Print Fiber Element Properties
 *  ==================================
 */

#ifdef __STDC__
void print_property_fiber_3d(EFRAME *frp, int i)
#else
void print_property_fiber_3d(frp, i)
EFRAME    *frp;
int          i;                 /* elmt_attr_no */
#endif
{
int     UNITS_SWITCH;
ELEMENT_ATTR    *eap;
int             ifib;

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

     printf("             ");
     printf(" General Section and Material Properties\n");
     switch(UNITS_SWITCH) {
       case ON:
        UnitsSimplify( eap->work_material[0].dimen );
        UnitsSimplify( eap->work_material[1].dimen );
        UnitsSimplify( eap->work_material[2].dimen );
        UnitsSimplify( eap->work_material[3].dimen );
        UnitsSimplify( eap->work_material[5].dimen );
        UnitsSimplify( eap->work_material[10].dimen );
        UnitsSimplify( eap->work_material[11].dimen );
        UnitsSimplify( eap->work_section[10].dimen );
	if( eap->work_material[5].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Density                 =      %11.3e %s\n",
                           eap->work_material[5].value/eap->work_material[5].dimen->scale_factor,
                           eap->work_material[5].dimen->units_name);
	}
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's Ratio         = nu = %11.3e   \n", eap->work_material[4].value);
        }
        if( eap->work_material[0].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Young's Modulus         =  E = %11.3e %s\n",
                           eap->work_material[0].value/eap->work_material[0].dimen->scale_factor,
                           eap->work_material[0].dimen->units_name);
        }
        if( eap->work_material[3].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Young's Tangent Modulus = Et = %11.3e %s\n",
                           eap->work_material[3].value/eap->work_material[3].dimen->scale_factor,
                           eap->work_material[3].dimen->units_name);
        }
        if( eap->work_material[2].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Yielding Stress         = fy = %11.3e %s\n",
                           eap->work_material[2].value/eap->work_material[2].dimen->scale_factor,
                           eap->work_material[2].dimen->units_name);
        }
        if( eap->work_material[1].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Shear Modulus           =  G = %11.3e %s\n",
                           eap->work_material[1].value/eap->work_material[1].dimen->scale_factor,
                           eap->work_material[1].dimen->units_name);
        }
        if( eap->work_material[10].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Shear Tangent Modulus   = Gt = %11.3e %s\n",
                           eap->work_material[10].value/eap->work_material[10].dimen->scale_factor,
                           eap->work_material[10].dimen->units_name);
        }
        if( eap->work_material[11].dimen->units_name != NULL ) {
           printf("             ");
           printf("         : Shear Yielding Stress   = fv = %11.3e %s\n",
                           eap->work_material[11].value/eap->work_material[11].dimen->scale_factor,
                           eap->work_material[11].dimen->units_name);
        }
        if( eap->work_section[16].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Correction Factor = ks = %11.3e   \n", eap->work_section[16].value);
        }
	if( eap->work_section[10].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Section Area            =      %11.3e %s\n",
                           eap->work_section[10].value/eap->work_section[10].dimen->scale_factor,
                           eap->work_section[10].dimen->units_name);
	}
       break;
       case OFF:
	if( eap->work_material[5].value != 0.0 ) {
          printf("             ");
          printf("         : Density                 =      %11.3e\n", eap->work_material[5].value);
	}
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's Ratio         = nu = %11.3e\n", eap->work_material[4].value);
        }
        if( eap->work_material[0].value != 0.0 ) {
           printf("             ");
           printf("         : Young's Modulus         =  E = %11.3e\n", eap->work_material[0].value);
        }
        if( eap->work_material[3].value != 0.0 ) {
           printf("             ");
           printf("         : Young's Tangent Modulus = Et = %11.3e\n", eap->work_material[3].value);
        }
        if( eap->work_material[2].value != 0.0 ) {
           printf("             ");
           printf("         : Yielding Stress         = fy = %11.3e\n", eap->work_material[2].value);
        }
        if( eap->work_material[1].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Modulus           =  G = %11.3e\n", eap->work_material[1].value);
        }
        if( eap->work_material[10].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Tangent Modulus   = Gt = %11.3e\n", eap->work_material[10].value);
        }
        if( eap->work_material[11].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Yielding Stress   = fv = %11.3e\n", eap->work_material[11].value);
        }
        if( eap->work_section[16].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Correction Factor = ks = %11.3e   \n", eap->work_section[16].value);
        }
	if( eap->work_section[10].value != 0.0 ) {
          printf("             ");
          printf("         : Section Area            =      %11.3e\n", eap->work_section[10].value);
	}
        break;
        default:
        break;
     }

     /* Print fiber attribution */

     printf("             ");
     printf("Fiber Attribution\n");
     printf("             ");
     printf("         : No. of Fibers = %6i\n", eap->work_fiber->no_fiber );
}

/*
 *  =====================================
 *  Fiber_Elmt_State_Det_3d()
 *
 *  Input :
 *  Output : void 
 *  =====================================
 */

#ifdef  __STDC__
void Fiber_Elmt_State_Det_3d( ARRAY *p, HISTORY_DATA *hp, int *flag )
#else
void Fiber_Elmt_State_Det_3d( p , hp, flag)
ARRAY *p;
HISTORY_DATA *hp;
int *flag;
#endif
{
int         no_integ_pt, no_section, no_fiber, no_shear, elmt_no;
int         total_fiber;
QUANTITY    length;

int         ii, jj, kk, sec, ifib;
int         UNITS_SWITCH, UnitsType;
double      cs, sn, tn;
MATRIX      *temp_m1, *temp_m2;
DIMENSIONS  *dp_length, *dp_force, *dp_stress, *dimen;

double      *abscissas, *weights;
double      xi;
double      scale;
MATRIX      *Q, *q;
MATRIX      *dpe, *dQ, *dq, *s;
MATRIX      *Dx, *dx, *ex, *rx;
MATRIX      *dDx, *ddx, *dex;
MATRIX      *DRx, *DUx;
MATRIX      *F, *K;
MATRIX      *L, *Ltrans, *R, *Rtrans;
MATRIX      *kx, *fx;
MATRIX      *lx, *bx, *bxtrans;
MATRIX      *stress, *tangent;
FIBER_ATTR  *fiber;
SECTION_DATA *saved_data;

    UNITS_SWITCH = CheckUnits();
    UnitsType    = CheckUnitsType();

    /* Assign Element Property */

    fiber       = p->fiber_ptr->fiber;
    no_fiber    = p->fiber_ptr->no_fiber;
    no_shear    = p->fiber_ptr->no_shear;
    total_fiber = no_fiber + (p->no_dimen-1)*no_shear;
    no_integ_pt = p->integ_ptr->integ_pts;
    no_section  = no_integ_pt + 2;
    elmt_no     = p->elmt_no;
    cs = p->coord[0][1].value - p->coord[0][0].value;
    sn = p->coord[1][1].value - p->coord[1][0].value;
    tn = p->coord[2][1].value - p->coord[2][0].value;
    length.value = sqrt( cs*cs + sn*sn + tn*tn );
    if( UNITS_SWITCH == ON ) {
       if( UnitsType == SI )
          dp_length = DefaultUnits("m");
       else if( UnitsType == US )
          dp_length = DefaultUnits("in");
       length.dimen = (DIMENSIONS *)MyCalloc(1, sizeof(DIMENSIONS));
       UnitsCopy( length.dimen, dp_length );
       free((char *) dp_length->units_name);
       free((char *) dp_length);

       SetUnitsOff();
    }

    /* Gauss-Lobatto Integration */

    abscissas = (double *)MyCalloc( no_section, sizeof(double) );
    weights   = (double *)MyCalloc( no_section, sizeof(double) );
    Gauss_Lobatto( abscissas, weights, no_integ_pt );

    /* Put Displacement Matrix To One Column Form */
    /* p->displ=[dpx1,dpx2; dpy1,dpy2; dpz1,dpz2; drx1,drx2; dry1,dry2; drz1,drz2]6x2 */
    /*  =>  dpe=[dpx1;dpy1;dpz1;drx1;dry1;drz1; dpx2;dpy2;dpz2;drx2;dry2;drz2]12x1   */

    dpe = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 12, 1 );
    for( ii=0 ; ii < p->dof_per_node ; ++ii ) {
       for( jj=0 ; jj < p->nodes_per_elmt ; ++jj ) {
          kk = ii + jj*p->dof_per_node;
          dpe->uMatrix.daa[kk][0] = p->displ->uMatrix.daa[ii][jj];
       }
    }

    /* Use the same address in array, therefore  */
    /* frame->element[elmt_no-1]->rp->Q_saved, q_saved will update automatically. */
    Q = p->Q_saved;
    q = p->q_saved;
    ex = hp->strain;   /* ex[no_section][total_fiber] */
    stress = hp->stress;
    tangent = hp->tangent;
#ifdef DEBUG
printf("\n Q_saved\n");
MatrixPrintVar( Q, (MATRIX *)NULL );
printf("\n q_saved\n");
MatrixPrintVar( q, (MATRIX *)NULL );
printf("\n stress\n");
MatrixPrintVar( stress, (MATRIX *)NULL );
printf("\n tangent\n");
MatrixPrintVar( tangent, (MATRIX *)NULL );
#endif

    /* Calculate The Initial Tangent Stiffness, And Each Section Related Matrix */

    saved_data = (SECTION_DATA *)MyCalloc( no_section, sizeof(SECTION_DATA) );
    kx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 5, 5 );
    bx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 5, 5 );
    lx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, total_fiber, 5 );

    for( sec=0 ; sec < no_section; ++sec ) {
       saved_data[sec].xi = length.value/2.0*(abscissas[sec]+1);
       saved_data[sec].wi = weights[sec]*length.value/2.0;

       Force_Interpolation_Matrix_3d( bx, saved_data[sec].xi, length );
       Section_Tangent_Stiffness_3d( kx, fiber, no_fiber, no_shear, sec, elmt_no );

       saved_data[sec].fx = MatrixInverse( kx );
       saved_data[sec].Dx = MatrixMult( bx, Q );
       dx = MatrixMult( saved_data[sec].fx, saved_data[sec].Dx );
       saved_data[sec].rx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 5, 1 );
       MatrixFree( dx );

       bxtrans = MatrixTranspose( bx );
       temp_m1 = MatrixMult( bxtrans, saved_data[sec].fx );
       temp_m2 = MatrixMult( temp_m1, bx );
       MatrixFree( temp_m1 );
       MatrixFree( bxtrans );
       if( sec == 0 )
          F = MatrixScale( temp_m2, saved_data[sec].wi );
       else {
          temp_m1 = MatrixScale( temp_m2, saved_data[sec].wi );
          MatrixAddReplace( F, temp_m1 );
          MatrixFree( temp_m1 );
       } 
       MatrixFree( temp_m2 );
    }
    K = MatrixInverse( F );

    /* Element Deformation Increments */

    L = Element_Transformation_3d( p->coord, length );
    dq = MatrixMult( L, dpe );
    MatrixAddReplace( q, dq );
#ifdef DEBUG
printf("\n dpe\n");
MatrixPrintVar( dpe, (MATRIX *)NULL );
printf("\n dq\n");
MatrixPrintVar( dq, (MATRIX *)NULL );
#endif

    /* Initial Residual Deformation Matrix, Resisting Force Matrix */

    s   = MatrixAllocIndirect(   "s",  DOUBLE_ARRAY, 5, 1 );
    DRx = MatrixAllocIndirect( "DRx",  DOUBLE_ARRAY, 5, 1 );

    /* Element Converge, j ; At Least Do Once */

    do
    {
       /* Reset Flexibility Matrix F, And Residual Matrix s */

       for( ii=0 ; ii < F->iNoRows ; ++ii )
       for( jj=0 ; jj < F->iNoColumns ; ++jj )
          F->uMatrix.daa[ii][jj] = 0.0;

       for( ii=0 ; ii < s->iNoRows ; ++ii ) {
          s->uMatrix.daa[ii][0] = 0.0;
       }

       /* Section State Determination */

       dQ = MatrixMult( K, dq );
       MatrixAddReplace( Q, dQ );
#ifdef DEBUG
printf("\n K in the do loop\n");
MatrixPrintVar( K, (MATRIX *)NULL );
printf("\n dq\n");
MatrixPrintVar( dq, (MATRIX *)NULL );
printf("\n dQ\n");
MatrixPrintVar( dQ, (MATRIX *)NULL );
#endif

       for( sec=0 ; sec < no_section; ++sec ) {
          xi = saved_data[sec].xi;
          scale = saved_data[sec].wi;
          fx = saved_data[sec].fx;
          Dx = saved_data[sec].Dx;
          rx = saved_data[sec].rx;
#ifdef DEBUG
printf("\n sec = %i\n",sec); 
printf("\n fx at the beginnig of do loop\n");
MatrixPrintVar( fx, (MATRIX *)NULL );
printf("\n Dx\n");
MatrixPrintVar( Dx, (MATRIX *)NULL );
#endif

          Force_Interpolation_Matrix_3d( bx, xi, length );
          Linear_Geometric_Matrix_3d( lx, fiber, no_fiber, no_shear );

          dDx = MatrixMult( bx, dQ );
          MatrixAddReplace( Dx, dDx );

          temp_m1 = MatrixMult( fx, dDx );
          ddx = MatrixAdd( rx, temp_m1 );
          dex = MatrixMult( lx, ddx );
          MatrixFree( temp_m1 );
          for( ifib=0 ; ifib < total_fiber ; ++ifib )
             ex->uMatrix.daa[sec][ifib] += dex->uMatrix.daa[ifib][0];

#ifdef DEBUG
printf("\n dDx\n");
MatrixPrintVar( dDx, (MATRIX *)NULL );
printf("\n ddx\n");
MatrixPrintVar( ddx, (MATRIX *)NULL );
printf("\n dex\n");
MatrixPrintVar( dex, (MATRIX *)NULL );
#endif

          Stress_Strain_Relationship( hp, p, fiber, ex, dex, no_fiber, no_shear, sec, flag[sec] );
          flag[sec]++;
#ifdef DEBUG
printf("\n stress\n");
MatrixPrintVar( stress, (MATRIX *)NULL );
printf("\n tangent\n");
MatrixPrintVar( tangent, (MATRIX *)NULL );
#endif

          MatrixFree( dDx );
          MatrixFree( ddx );
          MatrixFree( dex );

          Section_Tangent_Stiffness_3d( kx, fiber, no_fiber, no_shear, sec, elmt_no );
          MatrixFree( fx );
          fx = MatrixInverse( kx );
          saved_data[sec].fx = fx;

          Section_Resisting_Force_3d( DRx, stress, fiber, no_fiber, no_shear, sec );
          DUx = MatrixSub( Dx, DRx );
          MatrixFree( rx );
          rx = MatrixMult( fx, DUx );
          saved_data[sec].rx = rx;
          MatrixFree( DUx );
#ifdef DEBUG
printf("\n fx after Tangent\n");
MatrixPrintVar( fx, (MATRIX *)NULL );
printf("\n rx\n");
MatrixPrintVar( rx, (MATRIX *)NULL );
#endif

          bxtrans = MatrixTranspose( bx );
          temp_m1 = MatrixMult( bxtrans, fx );
          temp_m2 = MatrixMult( temp_m1, bx );
          MatrixFree( temp_m1 );
          temp_m1 = MatrixScale( temp_m2, scale );
          MatrixAddReplace( F, temp_m1 );
          MatrixFree( temp_m1 );
          MatrixFree( temp_m2 );

          temp_m1 = MatrixMult( bxtrans, rx );
          temp_m2 = MatrixScale( temp_m1, scale );
          MatrixAddReplace( s, temp_m2 );
          MatrixFree( temp_m1 );
          MatrixFree( temp_m2 );
          MatrixFree( bxtrans );
       }
       MatrixFree( K );
       K = MatrixInverse( F );

       for( ii=0 ; ii < dq->iNoRows ; ++ii )
          dq->uMatrix.daa[ii][0] = -s->uMatrix.daa[ii][0];

       MatrixFree( dQ );

    }  while ( dVmatrixL2Norm(s->uMatrix.daa, s->iNoRows, s->iNoColumns) > 0.001 );

    for( sec=0 ; sec < no_section ; ++sec ) {
       MatrixFree( saved_data[sec].fx );
       MatrixFree( saved_data[sec].Dx );
       MatrixFree( saved_data[sec].rx );
    }
    free((char *) saved_data);
    free((char *) abscissas);
    free((char *) weights);
    MatrixFree( dpe );
    MatrixFree( kx );
    MatrixFree( bx );
    MatrixFree( lx );
    MatrixFree( F );
    MatrixFree( K );
    MatrixFree( L );
    MatrixFree( s );
    MatrixFree( dq );
    MatrixFree( DRx );

    if( UNITS_SWITCH == ON )
       SetUnitsOn();
}


#ifdef  __STDC__
void Force_Interpolation_Matrix_3d( MATRIX *bx, double x, QUANTITY L )
#else
void Force_Interpolation_Matrix_3d( bx, x, L )
MATRIX *bx;
double   x;
QUANTITY L;
#endif
{
    bx->uMatrix.daa[0][0] = 1.0;
    bx->uMatrix.daa[0][1] = 0.0;
    bx->uMatrix.daa[0][2] = 0.0;
    bx->uMatrix.daa[0][3] = 0.0;
    bx->uMatrix.daa[0][4] = 0.0;

    bx->uMatrix.daa[1][0] = 0.0;
    bx->uMatrix.daa[1][1] = x/L.value-1.0;
    bx->uMatrix.daa[1][2] = x/L.value;
    bx->uMatrix.daa[1][3] = 0.0;
    bx->uMatrix.daa[1][4] = 0.0;

    bx->uMatrix.daa[2][0] = 0.0;
    bx->uMatrix.daa[2][1] = 0.0;
    bx->uMatrix.daa[2][2] = 0.0;
    bx->uMatrix.daa[2][3] = x/L.value-1.0;
    bx->uMatrix.daa[2][4] = x/L.value;

    bx->uMatrix.daa[3][0] = 0.0;
    bx->uMatrix.daa[3][1] = 0.0;
    bx->uMatrix.daa[3][2] = 0.0;
    bx->uMatrix.daa[3][3] = 1.0/L.value;
    bx->uMatrix.daa[3][4] = 1.0/L.value;

    bx->uMatrix.daa[4][0] = 0.0;
    bx->uMatrix.daa[4][1] = -1.0/L.value;
    bx->uMatrix.daa[4][2] = -1.0/L.value;
    bx->uMatrix.daa[4][3] = 0.0;
    bx->uMatrix.daa[4][4] = 0.0;
}

#ifdef  __STDC__
void Linear_Geometric_Matrix_3d( MATRIX *lx, FIBER_ATTR *fiber, int no_fiber, int no_shear )
#else
void Linear_Geometric_Matrix_3d( lx, fiber, no_fiber, no_shear )
MATRIX	*lx;
FIBER_ATTR  *fiber;
int  no_fiber, no_shear;
#endif
{
int	ifib;

    for( ifib=0 ; ifib < no_fiber ; ++ifib )
    {
       lx->uMatrix.daa[ifib][0] =  1.0;
       lx->uMatrix.daa[ifib][1] =  fiber[ifib].z.value;
       lx->uMatrix.daa[ifib][2] = -fiber[ifib].y.value;
       lx->uMatrix.daa[ifib][3] =  0.0;
       lx->uMatrix.daa[ifib][4] =  0.0;
    }
    for( ifib=no_fiber ; ifib < (no_fiber+no_shear) ; ++ifib )
    {
       lx->uMatrix.daa[ifib][0] =  0.0;
       lx->uMatrix.daa[ifib][1] =  0.0;
       lx->uMatrix.daa[ifib][2] =  0.0;
       lx->uMatrix.daa[ifib][3] =  1.0;
       lx->uMatrix.daa[ifib][4] =  0.0;
    }
    for( ifib=(no_fiber+no_shear) ; ifib < (no_fiber+2*no_shear) ; ++ifib )
    {
       lx->uMatrix.daa[ifib][0] =  0.0;
       lx->uMatrix.daa[ifib][1] =  0.0;
       lx->uMatrix.daa[ifib][2] =  0.0;
       lx->uMatrix.daa[ifib][3] =  0.0;
       lx->uMatrix.daa[ifib][4] =  1.0;
    }
}

#ifdef  __STDC__
void Section_Tangent_Stiffness_3d( MATRIX *kx, FIBER_ATTR *fiber, int no_fiber, int no_shear, int sec, int elmt_no )
#else
void Section_Tangent_Stiffness_3d( kx, fiber, no_fiber, no_shear, sec, elmt_no )
MATRIX *kx;
FIBER_ATTR *fiber;
int no_fiber, no_shear, sec;
int elmt_no;
#endif
{
int	  ifib, i, j;
double	y, z, A, Exi;
HISTORY_DATA     *hp;
MATRIX      *tangent;

    for( i=0 ; i < kx->iNoRows ; ++i )
       for( j=0 ; j < kx->iNoColumns ; ++j )
          kx->uMatrix.daa[i][j] = 0.0;

    /* according to strain at xi, get E(xi) for the fiber */
    hp = FiberElmtHistory( elmt_no );
    tangent = hp->tangent;

    for( ifib=0 ; ifib < no_fiber ; ++ifib )
    {
       y = fiber[ifib].y.value;
       z = fiber[ifib].z.value;
       A = fiber[ifib].area.value;
       Exi = tangent->uMatrix.daa[sec][ifib];

       kx->uMatrix.daa[0][0] = kx->uMatrix.daa[0][0] + Exi*A;
       kx->uMatrix.daa[0][1] = kx->uMatrix.daa[0][1] + Exi*A*z;
       kx->uMatrix.daa[0][2] = kx->uMatrix.daa[0][2] - Exi*A*y;
       kx->uMatrix.daa[1][1] = kx->uMatrix.daa[1][1] + Exi*A*z*z;
       kx->uMatrix.daa[1][2] = kx->uMatrix.daa[1][2] - Exi*A*y*z;
       kx->uMatrix.daa[2][2] = kx->uMatrix.daa[2][2] + Exi*A*y*y;
    }
    kx->uMatrix.daa[1][0] = kx->uMatrix.daa[0][1];
    kx->uMatrix.daa[2][0] = kx->uMatrix.daa[0][2];
    kx->uMatrix.daa[2][1] = kx->uMatrix.daa[1][2];

    if( no_shear == 1 ) {
       kx->uMatrix.daa[3][3] = tangent->uMatrix.daa[sec][no_fiber]*fiber[no_fiber].area.value;
       kx->uMatrix.daa[4][4] = tangent->uMatrix.daa[sec][no_fiber+1]*fiber[no_fiber+1].area.value;
    }
    else {
       for( ifib=no_fiber ; ifib < (no_fiber+no_shear) ; ++ifib )
          kx->uMatrix.daa[3][3] = kx->uMatrix.daa[3][3] + tangent->uMatrix.daa[sec][ifib]*fiber[ifib-no_fiber].area.value;

       for( ifib=(no_fiber+no_shear) ; ifib < (no_fiber+2*no_shear) ; ++ifib )
          kx->uMatrix.daa[4][4] = kx->uMatrix.daa[4][4] + tangent->uMatrix.daa[sec][ifib]*fiber[ifib-no_fiber-no_shear].area.value;
    }
}

#ifdef  __STDC__
void Section_Resisting_Force_3d( MATRIX *DRx, MATRIX *stress, FIBER_ATTR *fiber, int no_fiber, int no_shear, int sec )
#else
void Section_Resisting_Force_3d( DRx, stress, fiber, no_fiber, no_shear, sec )
MATRIX *DRx;
MATRIX *stress;
FIBER_ATTR *fiber;
int  no_fiber, no_shear, sec;
#endif
{
int  ifib, i, j;

    for( i=0 ; i < DRx->iNoRows ; ++i )
       DRx->uMatrix.daa[i][0] = 0.0;

    for( ifib=0 ; ifib < no_fiber ; ++ifib )
    {
       DRx->uMatrix.daa[0][0] = DRx->uMatrix.daa[0][0] + stress->uMatrix.daa[sec][ifib]*fiber[ifib].area.value;
       DRx->uMatrix.daa[1][0] = DRx->uMatrix.daa[1][0] + stress->uMatrix.daa[sec][ifib]*fiber[ifib].area.value*fiber[ifib].z.value;
       DRx->uMatrix.daa[2][0] = DRx->uMatrix.daa[2][0] - stress->uMatrix.daa[sec][ifib]*fiber[ifib].area.value*fiber[ifib].y.value;
    }

    if( no_shear == 1 ) {
       DRx->uMatrix.daa[3][0] = stress->uMatrix.daa[sec][no_fiber]*fiber[no_fiber].area.value;
       DRx->uMatrix.daa[4][0] = stress->uMatrix.daa[sec][no_fiber+1]*fiber[no_fiber+1].area.value;
    }
    else {
       for( ifib=no_fiber ; ifib < (no_fiber+no_shear) ; ++ifib )
          DRx->uMatrix.daa[3][0] = DRx->uMatrix.daa[3][0] + stress->uMatrix.daa[sec][ifib]*fiber[ifib-no_fiber].area.value;

       for( ifib=(no_fiber+no_shear) ; ifib < (no_fiber+2*no_shear) ; ++ifib )
          DRx->uMatrix.daa[4][0] = DRx->uMatrix.daa[4][0] + stress->uMatrix.daa[sec][ifib]*fiber[ifib-no_fiber-no_shear].area.value;
    }
}

#ifdef  __STDC__
MATRIX *Rigid_Body_Rotation_3d( QUANTITY length )
#else
MATRIX *Rigid_Body_Rotation_3d( length )
QUANTITY length;
#endif
{
int  i, j;
MATRIX *R;

    R  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 5, 12 );

    for( i=0 ; i < R->iNoRows ; ++i ) {
       for( j=0 ; j < R->iNoColumns ; ++j ) {
          R->uMatrix.daa[i][j] = 0.0;
       }
    }
    R->uMatrix.daa[0][0]  = -1.0;
    R->uMatrix.daa[0][6]  =  1.0;
    R->uMatrix.daa[1][2]  = -1.0/length.value;
    R->uMatrix.daa[1][4]  =  1.0;
    R->uMatrix.daa[1][8]  =  1.0/length.value;
    R->uMatrix.daa[2][2]  = -1.0/length.value;
    R->uMatrix.daa[2][8]  =  1.0/length.value;
    R->uMatrix.daa[2][10] =  1.0;
    R->uMatrix.daa[3][1]  =  1.0/length.value;
    R->uMatrix.daa[3][5]  =  1.0;
    R->uMatrix.daa[3][7]  = -1.0/length.value;
    R->uMatrix.daa[4][1]  =  1.0/length.value;
    R->uMatrix.daa[4][7]  = -1.0/length.value;
    R->uMatrix.daa[4][11] =  1.0;

    return( R );
}

#ifdef  __STDC__
MATRIX *Element_Transformation_3d( QUANTITY **coord, QUANTITY length )
#else
MATRIX *Element_Transformation_3d( coord, length )
QUANTITY **coord;
QUANTITY length;
#endif
{
MATRIX *Lele;
double cx, cy, cz, cxz;
double temp1, temp2, temp3;
int    i, j, k;

    /* Lele5x12 = Rigid5x12 * Trans12x12 */

    cx = (coord[0][1].value - coord[0][0].value)/length.value;
    cy = (coord[1][1].value - coord[1][0].value)/length.value;
    cz = (coord[2][1].value - coord[2][0].value)/length.value;
    cxz = sqrt( cx*cx + cz*cz );
    Lele = Rigid_Body_Rotation_3d( length );

    if( cx == 1.0 )
       return( Lele );

    /* rotation of vertical axes */
    if( cxz == 0.0 ) {
       for( i=0 ; i < 5 ; ++i ) {  /* [Q]5x1 */
          for( j=0 ; j < 12 ; j=j+3 ) {
             temp1 = Lele->uMatrix.daa[i][j];
             temp2 = Lele->uMatrix.daa[i][j+1];
             Lele->uMatrix.daa[i][j]   = -temp2*cy;
             Lele->uMatrix.daa[i][j+1] =  temp1*cy;
          }
       }
       return ( Lele );
    }

    for( i=0 ; i < 5 ; ++i ) {  /* [Q]5x1 */
       for( j=0 ; j < 12 ; j=j+3 ) {
          temp1 = Lele->uMatrix.daa[i][j];
          temp2 = Lele->uMatrix.daa[i][j+1];
          temp3 = Lele->uMatrix.daa[i][j+2];
          Lele->uMatrix.daa[i][j]   = temp1*cx - temp2*cx*cy/cxz - temp3*cz/cxz;
          Lele->uMatrix.daa[i][j+1] = temp1*cy + temp2*cxz;
          Lele->uMatrix.daa[i][j+2] = temp1*cz - temp2*cy*cz/cxz + temp3*cx/cxz;
       }
    }
    return ( Lele );
}
