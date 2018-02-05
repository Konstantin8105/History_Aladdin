/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_fiber2d.c : Linear/Nonlinear 2D Fiber Element
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
 *  --------------------------------------------------------------------
 *  Convention for Nodal Forces
 *             +ve M     -  anticlockwise(RHT Rule)
 *             +ve X,Y,Z -  along +ve axis
 *  Convention for Member End  Forces
 *             +ve M     -  Sagging Moments
 *             +ve SF    -  LHS upwards
 *             +ve AF    -  Tension(LHS outwards)
 *  --------------------------------------------------------------------
 *  Written by: Wane-Jang Lin                                   May 1996
 *  Modified by: Mark Austin                                  March 2000
 *  ====================================================================
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

/* 
 *  ============================================================== 
 *  Element FIBER_2D                                            
 *
 *  Input Properties:                                      
 *
 *       p->work_material[0] = E;
 *       p->work_material[1] = G;
 *       p->work_material[2] = fy;
 *       p->work_material[3] = ET;
 *       p->work_material[4] = nu;
 *       p->work_material[5] = density;
 *       p->work_material[6] = fu;
 *
 *       p->work_section[0] = Ixx;
 *       p->work_section[1] = Iyy;
 *       p->work_section[2] = Izz;
 *       p->work_section[3] = Ixy;
 *       p->work_section[4] = Ixz;
 *       p->work_section[5] = Iyz;
 *       p->work_section[6] = weight;
 *       p->work_section[7] = bf;
 *       p->work_section[8] = tf;
 *       p->work_section[9] = depth;                                  
 *       p->work_section[10] = area;
 *       p->work_section[11] = plate_thickness;
 *       p->work_section[12] = J;
 *       p->work_section[13] = rT;
 *       p->work_section[14] = width;
 *       p->work_section[15] = tw;                                  
 * 
 *  Input  :  ARRAY *p  -- pointer to working ARRAY data structure
 *         :  int isw   -- flag for task to be computed.
 *  Output :  ARRAY *p  -- pointer to working ARRAY data structure
 *  ============================================================== 
 */

#ifdef __STDC__
ARRAY *elmt_fiber_2d(ARRAY *p, int isw)
#else
ARRAY *elmt_fiber_2d(p, isw)
ARRAY *p;
int     isw;
#endif
{
static int       no_integ_pt, no_section, no_fiber, no_shear, elmt_no;
static QUANTITY  elmt_length;

int         ii, jj, kk, sec, ifib;
int         UNITS_SWITCH, UnitsType;
double      cs, sn, tn, xx, yy;
double      mbar;
MATRIX      *temp_m1, *temp_m2;
DIMENSIONS  *dp_length, *dp_force, *dp_stress, *dimen;

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

#ifdef DEBUG
       printf(" Enter elmt_fiber_2d() : isw = %d\n", isw );
#endif

    switch( isw ) {
       case PROPTY :
            no_fiber    = p->fiber_ptr->no_fiber;
            no_shear    = p->fiber_ptr->no_shear;
            no_integ_pt = p->integ_ptr->integ_pts;
            no_section  = no_integ_pt + 2;
            elmt_no     = p->elmt_no;
            cs = p->coord[0][1].value - p->coord[0][0].value;
            sn = p->coord[1][1].value - p->coord[1][0].value;
            elmt_length.value = sqrt( cs*cs + sn*sn );

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

            /* Allocate working arrays for element flexibility */

            kx = MatrixAllocIndirect( "kx", DOUBLE_ARRAY, 3, 3 );
            bx = MatrixAllocIndirect( "bx", DOUBLE_ARRAY, 3, 3 );

            /* Compute element flexibility matrix */

            for( sec=0 ; sec < no_section; ++sec ) {
                xi    = elmt_length.value/2.0*(abscissas[sec]+1);
                scale = weights[sec]*elmt_length.value/2.0;

                Force_Interpolation_Matrix_2d( bx, xi, elmt_length );

                Section_Tangent_Stiffness_2d( kx, fiber, no_fiber, no_shear, sec, elmt_no );

                fx      = MatrixInverse( kx );
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

            /* Invert element flexibility to get stiffness matrix */

            K = MatrixInverse( F );

            free((char *) abscissas);
            free((char *) weights);
            MatrixFree( kx );
            MatrixFree( bx );
            MatrixFree( F );

            /* Rigid body rotation and transform local coordinate to global */

            L       = Element_Transformation_2d( p->coord, elmt_length );
            Ltrans  = MatrixTranspose( L );
            temp_m1 = MatrixMult( Ltrans, K );
            Ke      = MatrixMult( temp_m1, L );

            /* Copy stiffness matrix to p array */

            for(ii = 1; ii <= p->stiff->iNoRows; ii++)
            for(jj = 1; jj <= p->stiff->iNoColumns; jj++)
                p->stiff->uMatrix.daa[ii-1][jj-1] = Ke->uMatrix.daa[ii-1][jj-1];

            MatrixFree( temp_m1 );
            MatrixFree( K );
            MatrixFree( L );
            MatrixFree( Ltrans );
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
                UnitsCopy( &(p->stiff->spRowUnits[2]), dp_length );
                UnitsCopy( &(p->stiff->spRowUnits[3]),
                           &(p->stiff->spRowUnits[0]) );
                UnitsCopy( &(p->stiff->spRowUnits[4]),
                           &(p->stiff->spRowUnits[1]) );
                UnitsCopy( &(p->stiff->spRowUnits[5]),
                           &(p->stiff->spRowUnits[2]) );

                UnitsDivRep( &(p->stiff->spColUnits[0]),
                             dp_force, dp_length, NO );
                UnitsCopy( &(p->stiff->spColUnits[1]),
                           &(p->stiff->spColUnits[0]) );
                UnitsCopy( &(p->stiff->spColUnits[2]), dp_force );
                UnitsCopy( &(p->stiff->spColUnits[3]),
                           &(p->stiff->spColUnits[0]) );
                UnitsCopy( &(p->stiff->spColUnits[4]),
                           &(p->stiff->spColUnits[1]) );
                UnitsCopy( &(p->stiff->spColUnits[5]),
                           &(p->stiff->spColUnits[2]) );

                free((char *) dp_force->units_name);
                free((char *) dp_length->units_name);
                free((char *) dp_force);
                free((char *) dp_length);
            } /* end of units on/off for case STIFF */

            break;   /* end of case STIFF */

       case LOAD_MATRIX:  /* compute internal load vector       */
       case STRESS:       /* compute and print element stresses */

            if( UNITS_SWITCH == ON )
               SetUnitsOff();

            /* --------------------------------------------------------- */
            /* Two cases to deal with:                                   */
            /*                                                           */
            /* p->elmt_state == 0 : ElmtStateDet() has not been called.  */
            /*    Therefore, element stresses and loads must be computed */
            /*    from scratch.                                          */
            /* p->elmt_state == 1 : ElmtStateDet() has been called.      */
            /*    Therefore, element stresses and loads are saveed and   */
            /*    only need to be retrieved from Aladdin database.       */
            /* --------------------------------------------------------- */

            if( p->elmt_state == 0 ) { 

            fiber     = p->fiber_ptr->fiber;
            abscissas = (double *)MyCalloc( no_section, sizeof(double) );
            weights   = (double *)MyCalloc( no_section, sizeof(double) );
            Gauss_Lobatto( abscissas, weights, no_integ_pt );

            /* Allocate memory for element force/stiffness */

            kx = MatrixAllocIndirect( "kx", DOUBLE_ARRAY, 3, 3 );
            bx = MatrixAllocIndirect( "bx", DOUBLE_ARRAY, 3, 3 );

            /* Construct element flexibility matrix by looping over sections */

            for( sec=0 ; sec < no_section; ++sec ) {
                xi    = elmt_length.value/2.0*(abscissas[sec]+1);
                scale = weights[sec]*elmt_length.value/2.0;

                Force_Interpolation_Matrix_2d( bx, xi, elmt_length );
                Section_Tangent_Stiffness_2d( kx, fiber, no_fiber, no_shear, sec, elmt_no );

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

            /* Invert element flexibility to get element-level stiffness */

            K = MatrixInverse( F );

            /* Clean up ..... */

            free((char *) abscissas);
            free((char *) weights);
            MatrixFree( kx );
            MatrixFree( bx );
            MatrixFree( F );

            /* --------------------------------------------------------- */
            /* Put displacement matrix into one column format i.e.,      */
            /* Transfer                                                  */
            /*                                                           */
            /*  p->displ=[ px1, px2;   to temp_m1 = pe = [ px1;          */
            /*             py1, py2;                       py1;          */
            /*             pz1, pz2 ]                      pz1;          */
            /*                                             px2;          */
            /*                                             py2;          */
            /*                                             pz2 ]         */
            /* --------------------------------------------------------- */

            temp_m1  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 6, 1 );

            for( ii=0 ; ii < p->dof_per_node ; ++ii )
            for( jj=0 ; jj < p->nodes_per_elmt ; ++jj )
                 temp_m1->uMatrix.daa[ii+jj*p->dof_per_node][0] =
                    p->displ->uMatrix.daa[ii][jj];

            L = Element_Transformation_2d( p->coord, elmt_length );
            q = MatrixMult( L, temp_m1 );
            Q = MatrixMult( K, q );
            MatrixFree( temp_m1 );

            }

            if( p->elmt_state == 1 ) {
               L = Element_Transformation_2d( p->coord, elmt_length );
               Q = p->Q_saved;
            }

            if( isw == LOAD_MATRIX ) {
               Ltrans  = MatrixTranspose( L );
               temp_m1 = MatrixMult( Ltrans, Q );   /* element nodal force in global coord. */
               MatrixFree( Ltrans );
            }

            if( isw == STRESS ) {
               R = Rigid_Body_Rotation_2d( elmt_length );
               Rtrans = MatrixTranspose( R );
               temp_m1 = MatrixMult( Rtrans, Q );   /* element nodal force in local coord. */
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

                /* node 1 */

                UnitsCopy( p->nodal_loads[0].dimen, dp_force );
                UnitsCopy( p->nodal_loads[1].dimen, dp_force );
                UnitsMultRep( p->nodal_loads[2].dimen, dp_force, dp_length );

                /* node 2 */

                UnitsCopy( p->nodal_loads[3].dimen,
                           p->nodal_loads[0].dimen );
                UnitsCopy( p->nodal_loads[4].dimen,
                           p->nodal_loads[1].dimen );
                UnitsCopy( p->nodal_loads[5].dimen,
                           p->nodal_loads[2].dimen );

                free((char *) dp_force->units_name);
                free((char *) dp_length->units_name);
                free((char *) dp_force);
                free((char *) dp_length);
            }

            if(isw == STRESS ) {
               xx  = 0.5*(p->coord[0][0].value + p->coord[0][1].value);   /* xx = 0.5(x1+x2) */
               yy  = 0.5*(p->coord[1][0].value + p->coord[1][1].value);   /* yy = 0.5(y1+y2) */
               printf("\n");
               printf("Elmt No %3d : \n", p->elmt_no);
               switch( UNITS_SWITCH ) {
                  case ON :
                       printf("Coords (X,Y) = (%10.3f %s, %10.3f %s)\n", 
                          xx/elmt_length.dimen->scale_factor,elmt_length.dimen->units_name,
                          yy/elmt_length.dimen->scale_factor,elmt_length.dimen->units_name);
                       printf("\n");

                       /* node i */

                       printf(" Fx1 = %13.5e %s  Fy1 = %13.5e %s  Mz1 = %13.5e %s\n",
                          p->nodal_loads[0].value/p->nodal_loads[0].dimen->scale_factor,
                          p->nodal_loads[0].dimen->units_name,
                          p->nodal_loads[1].value/p->nodal_loads[1].dimen->scale_factor,
                          p->nodal_loads[1].dimen->units_name,
                          p->nodal_loads[2].value/p->nodal_loads[2].dimen->scale_factor,
                          p->nodal_loads[2].dimen->units_name);

                       /* node j */

                       printf(" Fx2 = %13.5e %s  Fy2 = %13.5e %s  Mz2 = %13.5e %s\n",
                          p->nodal_loads[3].value/p->nodal_loads[3].dimen->scale_factor,
                          p->nodal_loads[3].dimen->units_name,
                          p->nodal_loads[4].value/p->nodal_loads[4].dimen->scale_factor,
                          p->nodal_loads[4].dimen->units_name,
                          p->nodal_loads[5].value/p->nodal_loads[5].dimen->scale_factor,
                          p->nodal_loads[5].dimen->units_name);
                          printf("\n");

                      /* Member Forces */

                      printf(" Axial Force : x-direction = %13.5e %s\n",
                         -p->nodal_loads[0].value/p->nodal_loads[0].dimen->scale_factor,
                          p->nodal_loads[0].dimen->units_name);
                      printf(" Shear Force : y-direction = %13.5e %s\n",
                          p->nodal_loads[1].value/p->nodal_loads[1].dimen->scale_factor,
                          p->nodal_loads[1].dimen->units_name);
                      printf("\n");
                      break;
                 case OFF :
                      printf("Coords (X,Y) = (%10.3f , %10.3f )\n", xx, yy);
                      printf("\n");

                      /* node i */

                      printf(" Fx1 = %13.5e   Fy1 = %13.5e   Mz1 = %13.5e \n",
                         p->nodal_loads[0].value, 
                         p->nodal_loads[1].value, 
                         p->nodal_loads[2].value);

                      /* node j */

                      printf(" Fx2 = %13.5e   Fy2 = %13.5e   Mz2 = %13.5e \n",
                         p->nodal_loads[3].value,
                         p->nodal_loads[4].value,
                         p->nodal_loads[5].value);
                      printf("\n");

                      /* Member Forces */

                      printf(" Axial Force : x-direction = %13.5e \n", -p->nodal_loads[0].value);
                      printf(" Shear Force : y-direction = %13.5e \n",  p->nodal_loads[1].value);
                      printf("\n");
                      break;
                  default:
                      break;
               }
            }
            break;  /* end of case STRESS, LOAD_MATRIX */

       case STRESS_MATRIX:

            /* save element nodal forces in working array */

            if( UNITS_SWITCH == ON )
                SetUnitsOff();

            /* --------------------------------------------------------- */
            /* Two cases to deal with:                                   */
            /*                                                           */
            /* p->elmt_state == 0 : ElmtStateDet() has not been called.  */
            /*    Therefore, element stresses and loads must be computed */
            /*    from scratch.                                          */
            /* p->elmt_state == 1 : ElmtStateDet() has been called.      */
            /*    Therefore, element stresses and loads are saveed and   */
            /*    only need to be retrieved from Aladdin database.       */
            /* --------------------------------------------------------- */

            if( p->elmt_state == 1 ) {
                L = Element_Transformation_2d( p->coord, elmt_length );
                Q = p->Q_saved;
            }

            if( p->elmt_state == 0 ) { 

                fiber     = p->fiber_ptr->fiber;
                abscissas = (double *)MyCalloc( no_section, sizeof(double) );
                weights   = (double *)MyCalloc( no_section, sizeof(double) );
                Gauss_Lobatto( abscissas, weights, no_integ_pt );

                /* Allocate memory for element force/stiffness */

                kx = MatrixAllocIndirect( "kx", DOUBLE_ARRAY, 3, 3 );
                bx = MatrixAllocIndirect( "bx", DOUBLE_ARRAY, 3, 3 );

                /* Construct element flexibility matrix by looping over sections */

                for( sec=0 ; sec < no_section; ++sec ) {
                   xi    = elmt_length.value/2.0*(abscissas[sec]+1);
                   scale = weights[sec]*elmt_length.value/2.0;

                   Force_Interpolation_Matrix_2d( bx, xi, elmt_length );
                   Section_Tangent_Stiffness_2d( kx, fiber, no_fiber, no_shear, sec, elmt_no );

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

                /* Invert element flexibility to get element-level stiffness */

                K = MatrixInverse( F );

                /* Clean up ..... */

                free((char *) abscissas);
                free((char *) weights);
                MatrixFree( kx );
                MatrixFree( bx );
                MatrixFree( F );

                /* --------------------------------------------------------- */
                /* Put displacement matrix into one column format i.e.,      */
                /* Transfer                                                  */
                /*                                                           */
                /*  p->displ=[ px1, px2;   to temp_m1 = pe = [ px1;          */
                /*             py1, py2;                       py1;          */
                /*             pz1, pz2 ]                      pz1;          */
                /*                                             px2;          */
                /*                                             py2;          */
                /*                                             pz2 ]         */
                /* --------------------------------------------------------- */

                temp_m1  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 6, 1 );
                for( ii=0 ; ii < p->dof_per_node ; ++ii )
                for( jj=0 ; jj < p->nodes_per_elmt ; ++jj )
                   temp_m1->uMatrix.daa[ii+jj*p->dof_per_node][0] = p->displ->uMatrix.daa[ii][jj];

                L = Element_Transformation_2d( p->coord, elmt_length );
                q = MatrixMult( L, temp_m1 );
                Q = MatrixMult( K, q );
                MatrixFree( temp_m1 );
            }

            /* Compute nodal forces in local element coordinates */

            R       = Rigid_Body_Rotation_2d( elmt_length );
            Rtrans  = MatrixTranspose( R );
            temp_m1 = MatrixMult( Rtrans, Q );

            /* Assign force values */

            for( ii=0 ; ii < p->size_of_stiff ; ++ii )
                 p->nodal_loads[ii].value = temp_m1->uMatrix.daa[ii][0];

            /* Cleanup memory from working arrays */

            if( p->elmt_state == 0 ) {
                MatrixFree( K );
                MatrixFree( Q );
                MatrixFree( q );
            }

            MatrixFree( R );
            MatrixFree( Rtrans );
            MatrixFree( L );
            MatrixFree( temp_m1 );

            /* Assign force units */

            if( UNITS_SWITCH == ON ) {
               SetUnitsOn();
               if( CheckUnitsType() == SI ) {
                   dp_length = DefaultUnits("m");
                   dp_force  = DefaultUnits("N");
               } else if (CheckUnitsType() == US ) {
                   dp_length = DefaultUnits("in");
                   dp_force  = DefaultUnits("lbf");
               } 

               /* node 1 */

               UnitsCopy( p->nodal_loads[0].dimen, dp_force );
               UnitsCopy( p->nodal_loads[1].dimen, dp_force );
               UnitsMultRep( p->nodal_loads[2].dimen, dp_force, dp_length );

               /* node 2 */

               UnitsCopy( p->nodal_loads[3].dimen, p->nodal_loads[0].dimen );
               UnitsCopy( p->nodal_loads[4].dimen, p->nodal_loads[1].dimen );
               UnitsCopy( p->nodal_loads[5].dimen, p->nodal_loads[2].dimen );

               /* Transfer nodal coordinates and forces/moment to working array */

               UnitsCopy( &(p->stress->spColUnits[0]), dp_length );
               UnitsCopy( &(p->stress->spColUnits[1]), dp_length );

               UnitsCopy( &(p->stress->spColUnits[2]), dp_force  );
               UnitsCopy( &(p->stress->spColUnits[3]), dp_force  );
               UnitsMultRep( &(p->stress->spColUnits[4]), dp_force, dp_length );

               /* Zero out units buffer */

               ZeroUnits( &(p->stress->spRowUnits[0]) );
               ZeroUnits( &(p->stress->spRowUnits[1]) );

               /* Release working memory */

               free((char *) dp_force->units_name);
               free((char *) dp_length->units_name);
               free((char *) dp_force);
               free((char *) dp_length);
            }

            /* Transfer coordinates to working stress matrix */

            p->stress->uMatrix.daa[0][0] = p->coord[0][0].value;
            p->stress->uMatrix.daa[0][1] = p->coord[1][0].value;
            p->stress->uMatrix.daa[1][0] = p->coord[0][1].value;
            p->stress->uMatrix.daa[1][1] = p->coord[1][1].value;

            /* Transfer internal loads to working stress matrix */

            p->stress->uMatrix.daa[0][2] = p->nodal_loads[0].value;
            p->stress->uMatrix.daa[0][3] = p->nodal_loads[1].value;
            p->stress->uMatrix.daa[0][4] = p->nodal_loads[2].value;
            p->stress->uMatrix.daa[1][2] = p->nodal_loads[3].value;
            p->stress->uMatrix.daa[1][3] = p->nodal_loads[4].value;
            p->stress->uMatrix.daa[1][4] = p->nodal_loads[5].value;


#ifdef DEBUG
       printf("*** Finished case STRESS_MATRIX : isw = %d\n", isw );

       for (ii = 1; ii <= 2; ii = ii + 1 ) {
            for (jj = 1; jj <= 5; jj = jj + 1 ) {
                 printf(" %10.4f ", p->stress->uMatrix.daa[ ii-1 ][ jj-1 ]);
            }
            printf("\n");
       }

       /* MatrixPrintVar ( p->stress ); */
       printf("*** Leaving STRESS_MATRIX case\n" );
#endif

            break; 
       case MASS_MATRIX:   /* compute element mass matrix */

          /* mbar = density * area */

          if( p->work_section[6].value != 0.0 )    /* mbar = weight/gravity */
             mbar = p->work_section[6].value/9.80665;
          else if( (p->work_material[5].value > 0.0) && (p->work_section[10].value > 0.0) )
             mbar = p->work_material[5].value * p->work_section[10].value;
          else
             FatalError("\nError in input: Need density value to calculate mass matrix\n",(char *)NULL);

          cs = p->coord[0][1].value - p->coord[0][0].value;
          sn = p->coord[1][1].value - p->coord[1][0].value;
          cs = cs / elmt_length.value;
          sn = sn / elmt_length.value;

          switch( p->type ) {  /* mass type */
             case LUMPED :
                /* use lumped mass matrix of FRAME_2D element */
                p->stiff = beamms( p, p->stiff, p->type, mbar, elmt_length.value, cs, sn,
                                   p->size_of_stiff, p->dof_per_node );
                break;

             case CONSISTENT :
             default :
                FatalError("FIBER_2D element only support LUMPED mass", (char *)NULL);
                break;
          }

          break;

       default:
          printf("*** ERROR in elmt_fiber2d() : isw = %d not defined\n", isw );
          break;
    }

    return(p);
}


/* ========================================================== */
/* print_property_fiber_2d() : Print Fiber Element Properties */
/*                                                            */
/* Input:                                                     */
/*                                                            */
/*   EFRAME *frp :                                            */
/*   int i       :                                            */
/*                                                            */
/* Output: void.                                              */
/* ========================================================== */

#ifdef __STDC__
void print_property_fiber_2d(EFRAME *frp, int i)
#else
void print_property_fiber_2d(frp, i)
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
           printf("         : Shear Correction Factor = ks = %11.3e\n",
                 eap->work_section[16].value);
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
           printf("         : Density                 =      %11.3e\n",
                 eap->work_material[5].value);
	}
        if( eap->work_material[4].value != 0.0 ) {
           printf("             ");
           printf("         : Poisson's Ratio         = nu = %11.3e\n",
                 eap->work_material[4].value);
        }
        if( eap->work_material[0].value != 0.0 ) {
           printf("             ");
           printf("         : Young's Modulus         =  E = %11.3e\n",
                 eap->work_material[0].value);
        }
        if( eap->work_material[3].value != 0.0 ) {
           printf("             ");
           printf("         : Young's Tangent Modulus = Et = %11.3e\n",
                 eap->work_material[3].value);
        }
        if( eap->work_material[2].value != 0.0 ) {
           printf("             ");
           printf("         : Yielding Stress         = fy = %11.3e\n",
                 eap->work_material[2].value);
        }
        if( eap->work_material[1].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Modulus           =  G = %11.3e\n",
                 eap->work_material[1].value);
        }
        if( eap->work_material[10].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Tangent Modulus   = Gt = %11.3e\n",
                 eap->work_material[10].value);
        }
        if( eap->work_material[11].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Yielding Stress   = fv = %11.3e\n",
                 eap->work_material[11].value);
        }
        if( eap->work_section[16].value != 0.0 ) {
           printf("             ");
           printf("         : Shear Correction Factor = ks = %11.3e\n",
                 eap->work_section[16].value);
        }
	if( eap->work_section[10].value != 0.0 ) {
          printf("             ");
          printf("         : Section Area            =      %11.3e\n",
                 eap->work_section[10].value);
	}
        break;
        default:
        break;
     }

     /* Print Fiber Attributes1 */

     printf("             ");
     printf("Fiber Attributes\n");
     printf("             ");
     printf("         : No. of Fibers = %6i\n", eap->work_fiber->no_fiber );
}


/*
 *  =================================================================
 *  State Determination for 2D fiber element
 *
 *  Input  : Array        *p  -- pointer to working array.
 *         : HISTORY_DATA *pp -- pointer to history info.
 *         : int *flag        -- flag.
 *  Output : void
 *  =================================================================
 */

/* #define DEBUG2 */

#ifdef  __STDC__
void Fiber_Elmt_State_Det_2d( ARRAY *p, HISTORY_DATA *hp, int *flag )
#else
void Fiber_Elmt_State_Det_2d( p , hp, flag)
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
double      cs, sn, tn, xx, yy;
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

#ifdef DEBUG
       printf("*** ===============================\n");
       printf("*** ENTER Fiber_Elmt_State_Det_2d()\n");
       printf("*** ===============================\n");
#endif

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
    length.value = sqrt( cs*cs + sn*sn );
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

    /* Put Displacement Matrix To One Column Form    */
    /* p->displ=[dpx1,dpx2; dpy1,dpy2; dpz1,dpz2]3x2 */
    /*  =>  dpe=[dpx1;dpy1;dpz1;dpx2;dpy2;dpz2]6x1   */

    dpe = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 6, 1 );
    for( ii=0 ; ii < p->dof_per_node ; ++ii ) {
       for( jj=0 ; jj < p->nodes_per_elmt ; ++jj ) {
          kk = ii + jj*p->dof_per_node;
          dpe->uMatrix.daa[kk][0] = p->displ->uMatrix.daa[ii][jj];
       }
    }

    /* Use the same address in array, therefore                                   */
    /* frame->element[elmt_no-1]->rp->Q_saved, q_saved will update automatically. */

    Q = p->Q_saved;
    q = p->q_saved;
    ex = hp->strain;   /* ex[no_section][total_fiber] */
    stress = hp->stress;
    tangent = hp->tangent;

#ifdef DEBUG2
printf("\n Q_saved\n");
MatrixPrintVar( Q, (MATRIX *)NULL );
printf("\n q_saved\n");
MatrixPrintVar( q, (MATRIX *)NULL );
printf("\n stress\n");
MatrixPrintVar( stress, (MATRIX *)NULL );
printf("\n tangent\n");
MatrixPrintVar( tangent, (MATRIX *)NULL );
#endif

    /* Calculate the initial tangent stiffness and related section matrix */

    saved_data = (SECTION_DATA *)MyCalloc( no_section, sizeof(SECTION_DATA) );
    kx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 3, 3 );
    bx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 3, 3 );
    lx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, total_fiber, 3 );

    for( sec=0 ; sec < no_section; ++sec ) {
       saved_data[sec].xi = length.value/2.0*(abscissas[sec]+1);
       saved_data[sec].wi = weights[sec]*length.value/2.0;

       Force_Interpolation_Matrix_2d( bx, saved_data[sec].xi, length );
       Section_Tangent_Stiffness_2d( kx, fiber, no_fiber, no_shear, sec, elmt_no );

       saved_data[sec].fx = MatrixInverse( kx );
       saved_data[sec].Dx = MatrixMult( bx, Q );
       dx = MatrixMult( saved_data[sec].fx, saved_data[sec].Dx );
       saved_data[sec].rx = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 3, 1 );
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

    /* Element deformation increments */

    L = Element_Transformation_2d( p->coord, length );
    dq = MatrixMult( L, dpe );
    MatrixAddReplace( q, dq );

#ifdef DEBUG2
printf("\n dpe\n");
MatrixPrintVar( dpe, (MATRIX *)NULL );
printf("\n dq\n");
MatrixPrintVar( dq, (MATRIX *)NULL );
#endif

    /* Initial Residual Deformation Matrix, Resisting Force Matrix */

    s   = MatrixAllocIndirect(   "s",  DOUBLE_ARRAY, 3, 1 );
    DRx = MatrixAllocIndirect( "DRx",  DOUBLE_ARRAY, 3, 1 );

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

#ifdef DEBUG2
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

#ifdef DEBUG2
printf("\n sec = %i\n",sec); 
printf("\n fx at the beginnig of do loop\n");
MatrixPrintVar( fx, (MATRIX *)NULL );
printf("\n Dx\n");
MatrixPrintVar( Dx, (MATRIX *)NULL );
#endif

          Force_Interpolation_Matrix_2d( bx, xi, length );
          Linear_Geometric_Matrix_2d( lx, fiber, no_fiber, no_shear );

          dDx = MatrixMult( bx, dQ );
          MatrixAddReplace( Dx, dDx );

          temp_m1 = MatrixMult( fx, dDx );
          ddx = MatrixAdd( rx, temp_m1 );
          dex = MatrixMult( lx, ddx );
          MatrixFree( temp_m1 );
          for( ifib=0 ; ifib < total_fiber ; ++ifib )
               ex->uMatrix.daa[sec][ifib] += dex->uMatrix.daa[ifib][0];

#ifdef DEBUG2
printf("\n dDx\n");
MatrixPrintVar( dDx, (MATRIX *)NULL );
printf("\n ddx\n");
MatrixPrintVar( ddx, (MATRIX *)NULL );
printf("\n dex\n");
MatrixPrintVar( dex, (MATRIX *)NULL );
#endif

          Stress_Strain_Relationship( hp, p, fiber, ex, dex, no_fiber, no_shear, sec, flag[sec] );
          flag[sec]++;

#ifdef DEBUG2
printf("\n stress\n");
MatrixPrintVar( stress, (MATRIX *)NULL );
printf("\n tangent\n");
MatrixPrintVar( tangent, (MATRIX *)NULL );
#endif

          MatrixFree( dDx );
          MatrixFree( ddx );
          MatrixFree( dex );

          Section_Tangent_Stiffness_2d( kx, fiber, no_fiber, no_shear, sec, elmt_no );
          MatrixFree( fx );
          fx = MatrixInverse( kx );
          saved_data[sec].fx = fx;

          Section_Resisting_Force_2d( DRx, stress, fiber, no_fiber, no_shear, sec );
          DUx = MatrixSub( Dx, DRx );
          MatrixFree( rx );
          rx = MatrixMult( fx, DUx );
          saved_data[sec].rx = rx;
          MatrixFree( DUx );

#ifdef DEBUG2
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

#ifdef DEBUG
    printf("*** ===============================\n");
    printf("*** LEAVE Fiber_Elmt_State_Det_2d()\n");
    printf("*** ===============================\n");
#endif

}



/*
 *  =================================================================
 *  Force Interpolation Matrix (2d case)
 *
 *  Input  : Array  *bx -- pointer to .....
 *         : double x   -- ....
 *         : QUANTITY L -- 
 *  Output : void
 *  =================================================================
 */



#ifdef  __STDC__
void Force_Interpolation_Matrix_2d( MATRIX *bx, double x, QUANTITY L )
#else
void Force_Interpolation_Matrix_2d( bx, x, L )
MATRIX *bx;
double   x;
QUANTITY L;
#endif
{
    bx->uMatrix.daa[0][0] = 1.0;
    bx->uMatrix.daa[0][1] = 0.0;
    bx->uMatrix.daa[0][2] = 0.0;
    bx->uMatrix.daa[1][0] = 0.0;
    bx->uMatrix.daa[1][1] = x/L.value-1.0;
    bx->uMatrix.daa[1][2] = x/L.value;
    bx->uMatrix.daa[2][0] = 0.0;
    bx->uMatrix.daa[2][1] = 1.0/L.value;
    bx->uMatrix.daa[2][2] = 1.0/L.value;
}


/*
 *  =================================================================
 *  Linear Geometric Matrix (2d case).
 *
 *  Input  : MATRIX *lx        --
 *         : FIBER_ATTR *fiber --
 *         : int no_fiber      --
 *           int no_shear      --
 *  Output : void
 *  =================================================================
 */

#ifdef  __STDC__
void Linear_Geometric_Matrix_2d( MATRIX *lx, FIBER_ATTR *fiber, int no_fiber, int no_shear )
#else
void Linear_Geometric_Matrix_2d( lx, fiber, no_fiber, no_shear )
MATRIX	*lx;
FIBER_ATTR  *fiber;
int  no_fiber, no_shear;
#endif
{
int	ifib, total_fiber;

    total_fiber = no_fiber + no_shear;
    for( ifib=0 ; ifib < no_fiber ; ++ifib )
    {
       lx->uMatrix.daa[ifib][0] =  1.0;
       lx->uMatrix.daa[ifib][1] = -fiber[ifib].y.value;
       lx->uMatrix.daa[ifib][2] =  0.0;
    }
    for( ifib=no_fiber ; ifib < total_fiber ; ++ifib )
    {
       lx->uMatrix.daa[ifib][0] =  0.0;
       lx->uMatrix.daa[ifib][1] =  0.0;
       lx->uMatrix.daa[ifib][2] =  1.0;
    }
}


/*
 *  =================================================================
 *  Section Tangent Stiffness Matrix
 *
 *  Input  : MATRIX *kx        --
 *         : FIBER_ATTR *fiber --
 *         : int no_fiber      --
 *         : int no_shear      --
 *         : int sec           -- 
 *         : int elmt_no       --
 *  Output : void
 *  =================================================================
 */


#ifdef  __STDC__
void
Section_Tangent_Stiffness_2d( MATRIX *kx, FIBER_ATTR *fiber, int no_fiber,
                              int no_shear, int sec, int elmt_no )
#else
void
Section_Tangent_Stiffness_2d( kx, fiber, no_fiber, no_shear, sec, elmt_no )
MATRIX *kx;
FIBER_ATTR *fiber;
int no_fiber, no_shear, sec;
int elmt_no;
#endif
{
int	  ifib, i, j;
int      total_fiber;
double	y, z, A, Exi;
HISTORY_DATA     *hp;
MATRIX      *tangent;

    /* Zero tangent stiffness matrix */

    for( i=0 ; i < kx->iNoRows ; ++i )
       for( j=0 ; j < kx->iNoColumns ; ++j )
          kx->uMatrix.daa[i][j] = 0.0;

    /* For strain at xi, get E(xi) for the fiber */

    hp = FiberElmtHistory( elmt_no );
    tangent = hp->tangent;

    for( ifib=0 ; ifib < no_fiber ; ++ifib ) {

       y   = fiber[ifib].y.value;
       A   = fiber[ifib].area.value;
       Exi = tangent->uMatrix.daa[sec][ifib];

       kx->uMatrix.daa[0][0] = kx->uMatrix.daa[0][0] + Exi*A;
       kx->uMatrix.daa[0][1] = kx->uMatrix.daa[0][1] - Exi*A*y;
       kx->uMatrix.daa[1][0] = kx->uMatrix.daa[0][1];
       kx->uMatrix.daa[1][1] = kx->uMatrix.daa[1][1] + Exi*A*y*y;
    }

    if( no_shear == 1 ) {
       kx->uMatrix.daa[2][2] = tangent->uMatrix.daa[sec][no_fiber]*fiber[no_fiber].area.value;
    } else {
       for( ifib=no_fiber ; ifib < total_fiber ; ++ifib )
          kx->uMatrix.daa[2][2] = kx->uMatrix.daa[2][2] +
                                  tangent->uMatrix.daa[sec][ifib]*fiber[ifib-no_fiber].area.value;
    }
}


/*
 *  =================================================================
 *  Section Resisting Force
 *
 *  Input  : MATRIX *DRx       --
 *         : MATRIX *stress    --
 *         : FIBER_ATTR *fiber --
 *         : int no_fiber      --
 *         : int no_shear      --
 *         : int sec           --
 *  Output : void
 *  =================================================================
 */

#ifdef  __STDC__
void
Section_Resisting_Force_2d( MATRIX *DRx, MATRIX *stress, FIBER_ATTR *fiber,
                            int no_fiber, int no_shear, int sec )
#else
void
Section_Resisting_Force_2d( DRx, stress, fiber, no_fiber, no_shear, sec )
MATRIX *DRx;
MATRIX *stress;
FIBER_ATTR *fiber;
int  no_fiber, no_shear, sec;
#endif
{
int  ifib, i, j;
int  total_fiber;

    total_fiber = no_fiber + no_shear;
    for( i=0 ; i < DRx->iNoRows ; ++i )
       DRx->uMatrix.daa[i][0] = 0.0;

    for( ifib=0 ; ifib < no_fiber ; ++ifib )
    {
       DRx->uMatrix.daa[0][0] = DRx->uMatrix.daa[0][0] +
            stress->uMatrix.daa[sec][ifib]*fiber[ifib].area.value;
       DRx->uMatrix.daa[1][0] = DRx->uMatrix.daa[1][0] -
            stress->uMatrix.daa[sec][ifib]*fiber[ifib].area.value*fiber[ifib].y.value;
    }

    if( no_shear == 1 ) {
       DRx->uMatrix.daa[2][0] = stress->uMatrix.daa[sec][no_fiber]*fiber[no_fiber].area.value;
    }
    else {
    for( ifib=no_fiber ; ifib < total_fiber ; ++ifib )
       DRx->uMatrix.daa[2][0] = DRx->uMatrix.daa[2][0] +
            stress->uMatrix.daa[sec][ifib]*fiber[ifib-no_fiber].area.value;
    }
}


/*
 *  =================================================================
 *  Rigid Body Rotation
 *
 *  Input  : QUANTITY length  --
 *  Output : MATRIX           --
 *  =================================================================
 */

#ifdef  __STDC__
MATRIX *Rigid_Body_Rotation_2d( QUANTITY length )
#else
MATRIX *Rigid_Body_Rotation_2d( length )
QUANTITY length;
#endif
{
int  i, j;
MATRIX *R;

    R  = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 3, 6 );

    for( i=0 ; i < R->iNoRows ; ++i ) {
       for( j=0 ; j < R->iNoColumns ; ++j ) {
          R->uMatrix.daa[i][j] = 0.0;
       }
    }
    R->uMatrix.daa[0][0] = -1.0;
    R->uMatrix.daa[0][3] =  1.0;
    R->uMatrix.daa[1][1] =  1.0/length.value;
    R->uMatrix.daa[1][4] = -1.0/length.value;
    R->uMatrix.daa[1][2] =  1.0;
    R->uMatrix.daa[2][1] =  1.0/length.value;
    R->uMatrix.daa[2][4] = -1.0/length.value;
    R->uMatrix.daa[2][5] =  1.0;

    return( R );
}


/*
 *  =================================================================
 *  Element Transformation (2d case)
 *
 *  Input  : QUANTITY **coord  --
 *         : QUANTITY length   --
 *  Output : MATRIX            --
 *  =================================================================
 */

#ifdef  __STDC__
MATRIX *Element_Transformation_2d( QUANTITY **coord, QUANTITY length )
#else
MATRIX *Element_Transformation_2d( coord, length )
QUANTITY **coord;
QUANTITY length;
#endif
{
MATRIX *Lele;
double cs, sn;
double temp;
int    i, j, k;

    /* Lele3x6 = Rigid3x6 * Trans6x6 */

    cs = (coord[0][1].value - coord[0][0].value)/length.value;
    sn = (coord[1][1].value - coord[1][0].value)/length.value;
    Lele = Rigid_Body_Rotation_2d( length );

    if( cs == 1.0 )
       return( Lele );

    for( i=0 ; i < 3 ; ++i ) {    /* [Q]3x1 */
       for( j=0 ; j < 6 ; j=j+3 ) {
          k = j + 1;
          temp = cs*Lele->uMatrix.daa[i][j] - sn*Lele->uMatrix.daa[i][k];
          Lele->uMatrix.daa[i][k] = sn*Lele->uMatrix.daa[i][j] + cs*Lele->uMatrix.daa[i][k];
          Lele->uMatrix.daa[i][j] = temp;
       }
    }

    return ( Lele );
}


/* 
 *  =======================================================================================
 *  Stress_Strain_Relationship() : Find fiber stress and tangent values according to strain
 *
 *  Input  : HISTORY_DATA  *hp --
 *         : ARRAY      *array --
 *         : FIBER_ATTR *fiber --
 *         : MATRIX        *ex --
 *         : MATRIX       *dex --
 *         : int no_fiber      --
 *         : int no_shear      --
 *         : int sec           --
 *         : int indx          --
 *  Output : void
 *  =======================================================================================
 */ 

#ifdef  __STDC__
void
Stress_Strain_Relationship( HISTORY_DATA *hp, ARRAY *array, FIBER_ATTR *fiber,
                            MATRIX *ex, MATRIX *dex, int no_fiber,
                            int no_shear, int sec, int indx )
#else
void Stress_Strain_Relationship( hp, array, fiber, ex, dex, no_fiber, no_shear, sec, indx )
HISTORY_DATA *hp;
ARRAY *array;
FIBER_ATTR *fiber;
MATRIX  *ex, *dex;
int  no_fiber, no_shear;
int sec, indx;
#endif
{
int  ifib, i, j;
int  total_fiber;
double  fy, ey, ks, kt;
double  k_x, s_x, e_x;

    total_fiber = no_fiber + (array->no_dimen-1)*no_shear;
    for( ifib=0 ; ifib < total_fiber ; ++ifib ) {
       hp->yielding[sec][ifib]  = array->yielding_saved[sec][ifib];
       hp->pre_range[sec][ifib] = array->pre_range_saved[sec][ifib];
       hp->pre_load[sec][ifib]  = array->pre_load_saved[sec][ifib];

       if( indx == 0 ) {
          if( dex->uMatrix.daa[ifib][0] > 0.0 )
             hp->loading[sec][ifib] = 1;
          else if( dex->uMatrix.daa[ifib][0] < 0.0 )
             hp->loading[sec][ifib] = -1;
          else
             hp->loading[sec][ifib] = hp->pre_load[sec][ifib];
       }

       hp->sr->uMatrix.daa[sec][ifib] = array->sr_saved->uMatrix.daa[sec][ifib];
       hp->er->uMatrix.daa[sec][ifib] = array->er_saved->uMatrix.daa[sec][ifib];
       hp->s0->uMatrix.daa[sec][ifib] = array->s0_saved->uMatrix.daa[sec][ifib];
       hp->e0->uMatrix.daa[sec][ifib] = array->e0_saved->uMatrix.daa[sec][ifib];
    }

    for( ifib=0 ; ifib < total_fiber ; ++ifib ) {
       if( ifib < no_fiber )
       {
          ks = fiber[ifib].Es.value;
          kt = fiber[ifib].Et.value;
          fy = fiber[ifib].fy.value;
          ey = fy/ks;
       }
       else if( ifib < 2*no_fiber )
       {
          if( no_shear==1 ) {
             ks = fiber[ifib].Es.value;
             kt = fiber[ifib].Et.value;
             fy = fiber[ifib].fy.value;
             ey = fy/ks;
          }
          else {
             ks = fiber[ifib-no_fiber].Gs.value;
             kt = fiber[ifib-no_fiber].Gt.value;
             fy = fiber[ifib-no_fiber].fv.value;
             ey = fy/ks;
          }
       }
       else
       {
          ks = fiber[ifib-no_fiber-no_shear].Gs.value;
          kt = fiber[ifib-no_fiber-no_shear].Gt.value;
          fy = fiber[ifib-no_fiber-no_shear].fv.value;
          ey = fy/ks;
       }

       if( hp->yielding[sec][ifib] == 0 )  /* first elastic or plastic */
       {
          if( ABS(ex->uMatrix.daa[sec][ifib]) <= ey ) /* first elastic */
          {
             k_x = ks;
             s_x = 0.0;
             e_x = 0.0;
             hp->yielding[sec][ifib] = 0;
             hp->pre_range[sec][ifib] = 0;
          } /* end of first elastic */
          else  /* first from elastic -> plastic, first start yielding */
          {
             k_x = kt;
             hp->s0->uMatrix.daa[sec][ifib] = fy*hp->loading[sec][ifib];
             hp->e0->uMatrix.daa[sec][ifib] = ey*hp->loading[sec][ifib];
             s_x = hp->s0->uMatrix.daa[sec][ifib];
             e_x = hp->e0->uMatrix.daa[sec][ifib];
             hp->yielding[sec][ifib] = 1;
             hp->pre_range[sec][ifib] = 1;
             array->elmt_state = 1;
          } /* end of first start yielding */
       } /* end of yielding[sec][ifib]==0, first elastic or plastic */

       else  /* yielding==1, plastic residual will occurs */
       {
          if( hp->pre_load[sec][ifib] != hp->loading[sec][ifib] )  /* load reversed */
          {
             if( hp->pre_range[sec][ifib] == 1 )
             {
                hp->sr->uMatrix.daa[sec][ifib] = array->sx_saved->uMatrix.daa[sec][ifib];
                hp->er->uMatrix.daa[sec][ifib] = array->ex_saved->uMatrix.daa[sec][ifib];
                k_x = ks;
                s_x = hp->sr->uMatrix.daa[sec][ifib];
                e_x = hp->er->uMatrix.daa[sec][ifib];
                hp->pre_range[sec][ifib] = 0;
             } /* end of pre_range==1 */
             else  /* pre_range==0 */
             {
                if( hp->loading[sec][ifib]*ex->uMatrix.daa[sec][ifib]
                <=  hp->loading[sec][ifib]*hp->er->uMatrix.daa[sec][ifib] )
                {
                   k_x = ks;
                   s_x = hp->sr->uMatrix.daa[sec][ifib];
                   e_x = hp->er->uMatrix.daa[sec][ifib];
                   hp->pre_range[sec][ifib] = 0;
                }
                else
                {
                   if( hp->loading[sec][ifib]*array->ex_saved->uMatrix.daa[sec][ifib]
                   >=  hp->loading[sec][ifib]*hp->er->uMatrix.daa[sec][ifib] )
                   {
                      k_x = ks;
                      s_x = hp->sr->uMatrix.daa[sec][ifib];
                      e_x = hp->er->uMatrix.daa[sec][ifib];
                      hp->pre_range[sec][ifib] = 0;
                   }
                   else
                   {
                      k_x = kt;
                      s_x = hp->sr->uMatrix.daa[sec][ifib];
                      e_x = hp->er->uMatrix.daa[sec][ifib];
                      hp->pre_range[sec][ifib] = 1;
                   }
                }
             } /* end of pre_range==0 */
          }  /* end of pre_load != loading, load reversed */

          else  /* pre_load==loading, add load in same direction */
          {
             if( hp->pre_range[sec][ifib] == 1 )
             {
                k_x = kt;
                s_x = hp->s0->uMatrix.daa[sec][ifib];
                e_x = hp->e0->uMatrix.daa[sec][ifib];
                hp->pre_range[sec][ifib] = 1;
             }
             else  /* pre_range=0 */
             {
                if( hp->loading[sec][ifib]*array->ex_saved->uMatrix.daa[sec][ifib]
                <=  hp->loading[sec][ifib]*hp->er->uMatrix.daa[sec][ifib] )
                {
                   if( hp->loading[sec][ifib]*ex->uMatrix.daa[sec][ifib]
                   <=  hp->loading[sec][ifib]*hp->er->uMatrix.daa[sec][ifib] )
                   {
                      k_x = ks;
                      s_x = hp->sr->uMatrix.daa[sec][ifib];
                      e_x = hp->er->uMatrix.daa[sec][ifib];
                      hp->pre_range[sec][ifib] = 0;
                   }
                   else
                   {
                      k_x = kt;
                      s_x = hp->sr->uMatrix.daa[sec][ifib];
                      e_x = hp->er->uMatrix.daa[sec][ifib];
                      hp->pre_range[sec][ifib] = 1;
                   }
                }
                else
                {
                   if( ABS(ex->uMatrix.daa[sec][ifib] - hp->er->uMatrix.daa[sec][ifib]) <= 2*ey )
                   {
                      k_x = ks;
                      s_x = hp->sr->uMatrix.daa[sec][ifib];
                      e_x = hp->er->uMatrix.daa[sec][ifib];
                      hp->pre_range[sec][ifib] = 0;
                   }
                   else
                   {
                      hp->s0->uMatrix.daa[sec][ifib] =
                         hp->sr->uMatrix.daa[sec][ifib] + hp->loading[sec][ifib]*2*fy;
                      hp->e0->uMatrix.daa[sec][ifib] =
                         hp->er->uMatrix.daa[sec][ifib] + hp->loading[sec][ifib]*2*ey;
                      k_x = kt;
                      s_x = hp->s0->uMatrix.daa[sec][ifib];
                      e_x = hp->e0->uMatrix.daa[sec][ifib];
                      hp->pre_range[sec][ifib] = 1;
                   }
                }
             } /* end of pre_range=0 */
          }  /* end of pre_load == loading */

       } /* end of yielding==1 */

       hp->tangent->uMatrix.daa[sec][ifib] = k_x;
       hp->stress->uMatrix.daa[sec][ifib]  = s_x + k_x*(ex->uMatrix.daa[sec][ifib]-e_x);
    }
}


/*
 *  =========================================================
 *  Get Gauss-Lobatto abscissas and weights
 * 
 *  Input : abscissas = section position in one element
 *        : weights   = array of weighting coefficients.
 *        : n = integration points = number of section in one
 *              element
 *  Output : void.
 *  =========================================================
 */

#ifdef  __STDC__
void Gauss_Lobatto( double *abscissas, double *weights , int n )
#else
void Gauss_Lobatto( abscissas, weights , n )
double *abscissas, *weights;
int n;
#endif
{
	/* Gauss-Lobatto integration */
	/* integral(-1,1)f(x)dx = A*f(-1) + sum(1,n)Ai*f(xi) + A*f(1) */

	switch(n)
	{
	   case 2:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.4472135954;
	      abscissas[2] =  0.4472135954;
	      abscissas[3] =  1.0;

	      weights[0] = 0.1666666666;
	      weights[1] = 0.8333333333;
	      weights[2] = 0.8333333333;
	      weights[3] = 0.1666666666;
	      break;

	   case 3:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.6546536707;
	      abscissas[2] =  0.0;
	      abscissas[3] =  0.6546536707;
	      abscissas[4] =  1.0;

	      weights[0] = 0.1000000000;
	      weights[1] = 0.5444444444;
	      weights[2] = 0.7111111111;
	      weights[3] = 0.5444444444;
	      weights[4] = 0.1000000000;
	      break;

	   case 4:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.7650553239;
	      abscissas[2] = -0.2852315164;
	      abscissas[3] =  0.2852315164;
	      abscissas[4] =  0.7650553239;
	      abscissas[5] =  1.0;

	      weights[0] = 0.0666666666;
	      weights[1] = 0.3784749562;
	      weights[2] = 0.5548583770;
	      weights[3] = 0.5548583770;
	      weights[4] = 0.3784749562;
	      weights[5] = 0.0666666666;
	      break;

	   case 5:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.8302238962;
	      abscissas[2] = -0.4688487934;
	      abscissas[3] =  0.0;
	      abscissas[4] =  0.4688487934;
	      abscissas[5] =  0.8302238962;
	      abscissas[6] =  1.0;

	      weights[0] = 0.0476190476;
	      weights[1] = 0.2768260473;
	      weights[2] = 0.4317453812;
	      weights[3] = 0.4876190476;
	      weights[4] = 0.4317453812;
	      weights[5] = 0.2768260473;
	      weights[6] = 0.0476190476;
	      break;

	   case 6:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.8717401485;
	      abscissas[2] = -0.5917001814;
	      abscissas[3] = -0.2092992179;
	      abscissas[4] =  0.2092992179;
	      abscissas[5] =  0.5917001814;
	      abscissas[6] =  0.8717401485;
	      abscissas[7] =  1.0;

	      weights[0] = 0.0357142857;
	      weights[1] = 0.2107042271;
	      weights[2] = 0.3411226924;
	      weights[3] = 0.4124587946;
	      weights[4] = 0.4124587946;
	      weights[5] = 0.3411226924;
	      weights[6] = 0.2107042271;
	      weights[7] = 0.0357142857;
	      break;

	   case 7:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.8997579954;
	      abscissas[2] = -0.6771862795;
	      abscissas[3] = -0.3631174638;
	      abscissas[4] =  0.0;
	      abscissas[5] =  0.3631174638;
	      abscissas[6] =  0.6771862795;
	      abscissas[7] =  0.8997579954;
	      abscissas[8] =  1.0;

	      weights[0] = 0.0277777777;
	      weights[1] = 0.1654953615;
	      weights[2] = 0.2745387125;
	      weights[3] = 0.3464285109;
	      weights[4] = 0.3715192743;
	      weights[5] = 0.3464285109;
	      weights[6] = 0.2745387125;
	      weights[7] = 0.1654953615;
	      weights[8] = 0.0277777777;
	      break;

	   case 8:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.9195339081;
	      abscissas[2] = -0.7387738651;
	      abscissas[3] = -0.4779249498;
	      abscissas[4] = -0.1652789576;
	      abscissas[5] =  0.1652789576;
	      abscissas[6] =  0.4779249498;
	      abscissas[7] =  0.7387738651;
	      abscissas[8] =  0.9195339081;
	      abscissas[9] =  1.0;

	      weights[0] = 0.0222222222;
	      weights[1] = 0.1333059908;
	      weights[2] = 0.2248893420;
	      weights[3] = 0.2920426836;
	      weights[4] = 0.3275397611;
	      weights[5] = 0.3275397611;
	      weights[6] = 0.2920426836;
	      weights[7] = 0.2248893420;
	      weights[8] = 0.1333059908;
	      weights[9] = 0.0222222222;
	      break;

	   case 9:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.9340014304;
	      abscissas[2] = -0.7844834736;
	      abscissas[3] = -0.5652353269;
	      abscissas[4] = -0.2957581355;
	      abscissas[5] =  0.0;
	      abscissas[6] =  0.2957581355;
	      abscissas[7] =  0.5652353269;
	      abscissas[8] =  0.7844834736;
	      abscissas[9] =  0.9340014304;
	      abscissas[10]=  1.0;

	      weights[0] = 0.0181818181;
	      weights[1] = 0.1096122732;
	      weights[2] = 0.1871698817;
	      weights[3] = 0.2480481042;
	      weights[4] = 0.2868791247;
	      weights[5] = 0.3002175954;
	      weights[6] = 0.2868791247;
	      weights[7] = 0.2480481042;
	      weights[8] = 0.1871698817;
	      weights[9] = 0.1096122732;
	      weights[10]= 0.0181818181;
	      break;

	   case 10:
	      abscissas[0] = -1.0;
	      abscissas[1] = -0.9448992722;
	      abscissas[2] = -0.8192793216;
	      abscissas[3] = -0.6328761530;
	      abscissas[4] = -0.3995309409;
	      abscissas[5] = -0.1365529328;
	      abscissas[6] =  0.1365529328;
	      abscissas[7] =  0.3995309409;
	      abscissas[8] =  0.6328761530;
	      abscissas[9] =  0.8192793216;
	      abscissas[10]=  0.9448992722;
	      abscissas[11]=  1.0;

	      weights[0] = 0.0151515151;
	      weights[1] = 0.0916845174;
	      weights[2] = 0.1579747055;
	      weights[3] = 0.2125084177;
	      weights[4] = 0.2512756031;
	      weights[5] = 0.2714052409;
	      weights[6] = 0.2714052409;
	      weights[7] = 0.2512756031;
	      weights[8] = 0.2125084177;
	      weights[9] = 0.1579747055;
	      weights[10]= 0.0916845174;
	      weights[11]= 0.0151515151;
	      break;

	   default:
	      printf("\nThis program doesn's support n=%d integral points in Gauss-Lobatto integration.\n");
	      printf("Please modify n (n=2,3,4,5,6,7,8,9,10) before continue\n");
	      exit(1);
	      break;
	}
}
