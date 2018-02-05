/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_frame3d.c : Three Dimensional Frame Element
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
 *  Three Dimensional Frame Element                                    
 *                                                                    
 *  Convention for Nodal Forces                                        
 *             +ve M     -  anticlockwise(RHT Rule)                    
 *             +ve X,Y,Z -  along +ve axis                             
 *  Convention for Member End  Forces                                  
 *             +ve M     -  Sagging Moments                           
 *             +ve SF    -  LHS upwards                               
 *             +ve AF    -  Tension(LHS outwards)                     
 *  ------------------------------------------------------------------- 
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
/*
#define DEBUG
*/

/* function declarations */

ARRAY *sld05( ARRAY * , int );


/* ============================================================== */
/*   Element FRAME_3D                                             */
/*   3D   Frame Element                                           */
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
      p->work_section[11] = plate_thickness;
      p->work_section[12] = J;
      p->work_section[13] = rT;
      p->work_section[14] = width;
      p->work_section[15] = tw;                                   */
/* ============================================================== */

/* macro for J constant */

#define J_Const(x,y) ((1 - 0.63 * x/y)*( x * x* x* y/3))

#ifdef __STDC__
ARRAY *elmt_frame_3d(ARRAY *p, int isw)
#else
ARRAY *elmt_frame_3d(p, isw)
ARRAY *p;
int   isw;
#endif
{
static double  nu;
static QUANTITY  fy, G, E, ET, density;
static double  Ixx, Iyy, Izz, Ixy, Ixz, Iyz, bf, tf, A, depth, weight, EA, EIzz, J, rT;

double d1;
double cs, sn,tn,xl,xx,yy,zz,vv,xn,xm,mbar;
double **rot, **trot, **fr, **dlocal;
int    i,j,k,l;

DIMENSIONS  *dp_length, *dp_force, *dp_moment;
DIMENSIONS  *dp_stress, *dp_degree, *dp_temperature;
int          UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter elmt_frame_3d() : isw = %4d\n", isw);
#endif

   H_Print = 0;
   UNITS_SWITCH = CheckUnits();

   switch(isw) {
       case PROPTY:  /* MAT PROPS */

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

          /* (1)   check  poi_ratio value */

          if( nu == 0.0 || nu > 0.5 ) {
              printf("WARNING >> ... In 3d beam element() nu value = %9.4f reset to 0.3!\n", nu);
              nu = 0.3;    /* default poi_ratio value */
          }

          /* (2)   calculate  G value */
            
/*        if(E.value/((1.0 - 2.0*nu)) != p->work_material[1].value) {
              printf(" elmt_frame_3d(): WARNING: G is not equal to E/(1-2nu), check G for homogeneous material \n");
              printf("                : ignore this message for non-homogeneous materials \n");
          }
*/
          G.value = p->work_material[1].value = E.value/(1.0 - 2.0*nu) ;
          if(UNITS_SWITCH == ON)  G.dimen = E.dimen;

          Ixx    = p->work_section[0].value;
          Iyy    = p->work_section[1].value;
          Izz    = p->work_section[2].value;
          Ixy    = p->work_section[3].value;
          Ixz    = p->work_section[4].value;
          Iyz    = p->work_section[5].value;
          weight = p->work_section[6].value;
          bf     = p->work_section[7].value;
          tf     = p->work_section[8].value;
          depth  = p->work_section[9].value;
          A      = p->work_section[10].value;
          J      = p->work_section[12].value;
          rT     = p->work_section[13].value;

          EA     = E.value*A;
          EIzz   = E.value*Izz;

          /* (3) If J value not input, J calculated based on rectangular */
          /*     section of size (bf x depth)                            */

          if(J == 0.0 ) {
             if(bf == 0.0 || depth == 0.0){
                printf("WARNING >> Must give 'J' or ('width' & 'depth') to calculate stiffness");
                exit(1);
             }

             /* Check bf < depth & cal J */

             if(bf < depth ) 
                 J = p->work_section[12].value = J_Const(bf,depth); 
             else
                 J = p->work_section[12].value = J_Const(depth,bf); 
          }
          break;
       case CHERROR:
            break;
       case STRESS_UPDATE:
            break;
       case PRESSLD:
            break;

       case STIFF:

            cs = p->coord[0][1].value - p->coord[0][0].value;           /* Cos Term */
            sn = p->coord[1][1].value - p->coord[1][0].value;           /* Sin Term */
            tn = p->coord[2][1].value - p->coord[2][0].value;           /* Tan Term */
            xl = sqrt(cs * cs + sn * sn + tn * tn); /* Calculate  Length */
            p->length.value = xl;

            /* T matrix is made here: 12*12 size */

            rot = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
            rot = (double **) tmat(rot, 6, p); 

            if( UNITS_SWITCH==ON )   p->stiff->spColUnits[0].units_type = CheckUnitsType();

            p->stiff = beamst3d(p, p->stiff, EA, EIzz, E.value*Iyy, G.value*J ,xl, rot,
                                p->size_of_stiff, p->dof_per_node); 

            MatrixFreeIndirectDouble(rot, p->size_of_stiff);
            break;

       case MASS_MATRIX:

            cs = p->coord[0][1].value - p->coord[0][0].value;           /* Cos Term */
            sn = p->coord[1][1].value - p->coord[1][0].value;           /* Sin Term */
            tn = p->coord[2][1].value - p->coord[2][0].value;           /* Tan Term */
            xl = sqrt(cs * cs + sn * sn + tn * tn); /* Calculate  Length */
            p->length.value = xl;

            /* T matrix is made here: 12*12 size */

            rot = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
            rot = (double **) tmat(rot, 6, p); 

            /* Assemble Mass Matrix */
 
            /* Calculate mbar =  mass/length                 */
            /* in units of (kg/m) or (lbf*sec^2/in/in)       */
            /* if no units, assume gravity g=9.80665 m/sec^2 */

             if( weight != 0.0 )
               mbar = weight/9.80665;
             else
               if( density.value > 0 )  mbar = A * density.value ;
             else {
                printf("\nError in input: Need density value to calculate mass matrix\n");
                exit(1);
             }

	    /* Calculate radius of gyration , rT  --  m  ,  in   */
            /* original version      :  rT = p->length.value/ 1.414; */

            if( rT==0 && A!=0 )    rT = sqrt( J / A );

            if( UNITS_SWITCH == ON ) 
                p->stiff->spColUnits[0].units_type = CheckUnitsType();

            p->stiff = beamms3d(p, p->stiff,p->type, mbar, xl, rT, rot, p->size_of_stiff, p->dof_per_node);

            MatrixFreeIndirectDouble(rot, p->size_of_stiff);
            break;

       case STRESS:
       case LOAD_MATRIX:

            cs = p->coord[0][1].value - p->coord[0][0].value;           /* Cos Term */
            sn = p->coord[1][1].value - p->coord[1][0].value;           /* Sin Term */
            tn = p->coord[2][1].value - p->coord[2][0].value;           /* Tan Term */
            xl = sqrt(cs * cs + sn * sn + tn * tn); /* Calculate  Length */
            p->length.value = xl;

            /* T matrix is made here: 12*12 size */

            rot = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
            rot = (double **) tmat(rot, 6, p); 

            /* --------------------------- */
            /* Output Stresses and Strains */
            /* --------------------------- */

            cs = rot[0][0];
            rot[0][0] = 1.0;

            if (UNITS_SWITCH==ON)  p->stiff->spColUnits[0].units_type = CheckUnitsType();
            p->stiff = beamst3d(p, p->stiff, EA, EIzz, E.value*Iyy,
                       G.value*J, xl, rot, p->size_of_stiff, p->dof_per_node); 

            rot[0][0] = cs;

            fr = MatrixAllocIndirectDouble(p->size_of_stiff,1); /* Here size_of_stiff = 12; i.e 12x1 mat */

            for(l = 1; l<= p->nodes_per_elmt; l++) {
                for(k = 1; k <= p->dof_per_node; k++) {
                    j = (l-1)*p->dof_per_node + k;
                    fr[j-1][0] = p->displ->uMatrix.daa[k-1][l-1];
                }
            }

            dlocal = (double **)  dMatrixMult(rot, p->size_of_stiff, p->size_of_stiff, fr, p->size_of_stiff, 1);

            fr = (double **)  dMatrixMultRep(fr, p->stiff->uMatrix.daa, p->size_of_stiff, p->size_of_stiff, dlocal, p->size_of_stiff, 1);

            xx = 0.5 *(p->coord[0][0].value +p->coord[0][1].value);
            yy = 0.5 *(p->coord[1][0].value +p->coord[1][1].value);
            zz = 0.5 *(p->coord[2][0].value +p->coord[2][1].value);

            if(H_Print == YES && isw == STRESS){
               printf( "\n");
               printf( "3D Frame Element   Stresses\n");
               printf( "---------------------------\n");
               H_Print = NO;
            }

            if(isw == STRESS ) {
               printf("\n");
               printf("Elmt# %3d : ", p->elmt_no);
               if(UNITS_SWITCH == ON)
                  printf(": Coords (X,Y,Z)= (%8.3f %s,%8.3f %s,%8.3f %s)\n",
                            xx/p->coord[0][0].dimen->scale_factor, p->coord[0][0].dimen->units_name,
                            yy/p->coord[1][0].dimen->scale_factor, p->coord[1][0].dimen->units_name,
                            zz/p->coord[2][0].dimen->scale_factor, p->coord[2][0].dimen->units_name);
	       else
                  printf(": Coords (X,Y,Z)= (%8.3f,%8.3f,%8.3f)\n", xx, yy, zz);

               printf("\n");
	    }

            /*--------------------------------------------------------*/
            /* nodal forces   & member end forces                     */
            /*--------------------------------------------------------*/

            if(p->elmt_load_ptr != NULL ) { /* calculate FEF */
               printf("Fixed End Loads; \n");
               p = sld05(p, STRESS);

               /* Add FEF to joint p->[  ]  orces */
                for(i = 1; i <= 12; i++)
                  fr[i-1][0] = fr[i-1][0] - p->nodal_loads[i-1].value;
            }

            /* -------------------- */
            /* Print Element forces */
            /* -------------------- */

            for(j = 1; j <= p->size_of_stiff; j++) 
                p->nodal_loads[j-1].value = fr[j-1][0];

            /* Assign element forces's units */

            if( UNITS_SWITCH == ON ) {
                if( CheckUnitsType() == SI) {
                    dp_length = DefaultUnits("m");
                    dp_force  = DefaultUnits("N");
                }
                if( CheckUnitsType() == US) {
                    dp_length = DefaultUnits("in");
                    dp_force  = DefaultUnits("lbf");
                }
                dp_moment = UnitsMult( dp_force, dp_length );

                for(j=1;j<=3;j++) {
		    UnitsCopy( p->nodal_loads[j-1].dimen, dp_force );
		    UnitsCopy( p->nodal_loads[j-1+p->dof_per_node].dimen, dp_force );
		    UnitsCopy( p->nodal_loads[j-1+3].dimen, dp_moment );
		    UnitsCopy( p->nodal_loads[j-1+3+p->dof_per_node].dimen, dp_moment );
		}
            }

            if(isw == STRESS ) {
	       switch(UNITS_SWITCH) {
		 case ON:
                  /* node_i */
                  printf("            Fx1 = %13.5e %s\t Fy1 = %13.5e %s\t Fz1 = %13.5e %s\n",
                         p->nodal_loads[0].value/dp_force->scale_factor, dp_force->units_name,
                         p->nodal_loads[1].value/dp_force->scale_factor, dp_force->units_name,
                         p->nodal_loads[2].value/dp_force->scale_factor, dp_force->units_name);
                  printf("            Mx1 = %13.5e %s\t My1 = %13.5e %s\t Mz1 = %13.5e %s\n",
                         p->nodal_loads[3].value/dp_moment->scale_factor, dp_moment->units_name,
                         p->nodal_loads[4].value/dp_moment->scale_factor, dp_moment->units_name,
                         p->nodal_loads[5].value/dp_moment->scale_factor, dp_moment->units_name);
                  printf("\n");
                  /* node_j */
                  printf("            Fx2 = %13.5e %s\t Fy2 = %13.5e %s\t Fz2 = %13.5e %s\n",
                        p->nodal_loads[6].value/dp_force->scale_factor, dp_force->units_name,
                        p->nodal_loads[7].value/dp_force->scale_factor, dp_force->units_name,
                        p->nodal_loads[8].value/dp_force->scale_factor, dp_force->units_name);
                  printf("            Mx2 = %13.5e %s\t My2 = %13.5e %s\t Mz2 = %13.5e %s\n",
                        p->nodal_loads[9].value/dp_moment->scale_factor, dp_moment->units_name,
                        p->nodal_loads[10].value/dp_moment->scale_factor, dp_moment->units_name,
                        p->nodal_loads[11].value/dp_moment->scale_factor, dp_moment->units_name);
                  printf("\n");

                  printf("            Axial Force : x-direction = %13.5e %s \n",
                            -p->nodal_loads[0].value/dp_force->scale_factor,
                             dp_force->units_name);
                  printf("            Shear Force : y-direction = %13.5e %s \n",
                             p->nodal_loads[1].value/dp_force->scale_factor,
                             dp_force->units_name);
                  printf("                        : z-direction = %13.5e %s \n",
                             p->nodal_loads[2].value/dp_force->scale_factor,
                             dp_force->units_name);
                  printf("\n"); 
                 break;
                 case OFF:
                  /* node_i */
                  printf("            Fx1 = %13.5e\t Fy1 = %13.5e\t Fz1 = %13.5e\n",
                           p->nodal_loads[0].value,
                           p->nodal_loads[1].value,
                           p->nodal_loads[2].value); 
                  printf("            Mx1 = %13.5e\t My1 = %13.5e\t Mz1 = %13.5e\n",
                           p->nodal_loads[3].value,
                           p->nodal_loads[4].value,
                           p->nodal_loads[5].value); 
                  printf("\n"); 
                  /* node_j */
                  printf("            Fx2 = %13.5e\t Fy2 = %13.5e\t Fz2 = %13.5e\n",
                           p->nodal_loads[6].value,
                           p->nodal_loads[7].value,
                           p->nodal_loads[8].value); 
                  printf("            Mx2 = %13.5e\t My2 = %13.5e\t Mz2 = %13.5e\n",
                           p->nodal_loads[9].value,
                           p->nodal_loads[10].value,
                           p->nodal_loads[11].value); 
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

            if(isw == LOAD_MATRIX ) {
                for(j = 1; j <= p->size_of_stiff; j++) 
                    fr[j-1][0] = 0.0;
                for( i=1 ; i<=p->size_of_stiff; i++ ) {
                    for( j=1 ; j<=p->size_of_stiff; j++ ) {
                        fr[i-1][0] += rot[j-1][i-1]*p->nodal_loads[j-1].value;
                    }
                }
                for(j = 1; j <= p->size_of_stiff; j++) 
                    p->nodal_loads[j-1].value = fr[j-1][0];
            }

            if( UNITS_SWITCH==ON ) {
                free((char *) dp_length->units_name);
                free((char *) dp_length);
                free((char *) dp_force->units_name);
                free((char *) dp_force);
                free((char *) dp_moment->units_name);
                free((char *) dp_moment);
            }

            MatrixFreeIndirectDouble(fr, p->size_of_stiff);
            MatrixFreeIndirectDouble(dlocal, p->size_of_stiff);
            MatrixFreeIndirectDouble(rot, p->size_of_stiff);
            break;
      default:
            break;
   }

   return(p);
}


/* =============================================================== */
/* 3D Frame Stiffness: This function returns a 12x12 rotated       */
/* stiffness matrix for 3D frame elmts.                            */
/*                                                                 */
/* Input:                                                          */
/*   ea->EA ; s->matrix pointer; eiz->EI_zaxis;  eiy->EI_yaxis;    */
/*   gj->G.J; xl->length       ; rot->rotation_matrix(12x12)       */
/*   nst->num_of_deg_of free per element; ndf->dof's per node;     */
/* =============================================================== */

#ifdef __STDC__
MATRIX *beamst3d(ARRAY *p, MATRIX *s, double ea, double eiz, double eiy, double gj, double xl, double **rot, int nst, int no_dof)
#else
MATRIX *beamst3d(p,s, ea, eiz, eiy, gj, xl, rot, nst, no_dof)
ARRAY  *p;
MATRIX *s;
double **rot;
double ea, eiz,eiy,gj, xl;
int    nst, no_dof;
#endif
{
int    i, j, k,ii,jj,kk;
double t;
DIMENSIONS     *d1, *d2;

    i  = no_dof + 1;
    j  = no_dof + 2; 
    k  = no_dof + 3;
    ii = no_dof + 4;
    jj = no_dof + 5;
    kk = no_dof + 6;

    t = ea/xl; 
#ifdef DEBUG

    printf(" in beamst3d() : \n");
    printf("               : EA = %lf, xl = %lf \n", ea, xl);
    printf("               : t = %lf\n", t);
#endif

    s->uMatrix.daa[0][0] =  t;
    s->uMatrix.daa[i-1][i-1] =  t;
    s->uMatrix.daa[0][i-1] = -t;
    s->uMatrix.daa[i-1][0] = -t;

    t = 12 * eiz/xl/xl/xl ;

    s->uMatrix.daa[1][1] = t;
    s->uMatrix.daa[j-1][j-1] = t;
    s->uMatrix.daa[1][j-1] = -t;
    s->uMatrix.daa[j-1][1] = -t;

    t = 12 * eiy/xl/xl/xl ;

    s->uMatrix.daa[2][2] = t;
    s->uMatrix.daa[k-1][k-1] = t;
    s->uMatrix.daa[2][k-1] = -t;
    s->uMatrix.daa[k-1][2] = -t;

    t = gj / xl;

    s->uMatrix.daa[3][3] = s->uMatrix.daa[ii-1][ii-1] =  t;
    s->uMatrix.daa[3][ii-1] = s->uMatrix.daa[ii-1][3] = -t;

    t = (eiy+ eiy) / xl;

    s->uMatrix.daa[4][4] = s->uMatrix.daa[jj-1][jj-1] = t+t;
    s->uMatrix.daa[4][jj-1] = s->uMatrix.daa[jj-1][4] = t;

    t = 6 * eiy/xl/xl;

    s->uMatrix.daa[2][4] = s->uMatrix.daa[4][2] = -t;
    s->uMatrix.daa[2][jj-1] = s->uMatrix.daa[jj-1][2] = -t;
    s->uMatrix.daa[k-1][4] = s->uMatrix.daa[4][k-1] =  t;
    s->uMatrix.daa[jj-1][k-1] = s->uMatrix.daa[k-1][jj-1] =  t;

    t = (eiz+ eiz) / xl;

    s->uMatrix.daa[5][5] = s->uMatrix.daa[kk-1][kk-1] = t+t;
    s->uMatrix.daa[5][kk-1] = s->uMatrix.daa[kk-1][5] = t;

    t = 6 * eiz/xl/xl;

    s->uMatrix.daa[1][5] = s->uMatrix.daa[5][1] = t;
    s->uMatrix.daa[1][kk-1] = s->uMatrix.daa[kk-1][1] = t;
    s->uMatrix.daa[5][j-1] = s->uMatrix.daa[j-1][5] = -t;
    s->uMatrix.daa[j-1][kk-1] = s->uMatrix.daa[kk-1][j-1] = -t;

    /* ==================== */
    /* Initial units buffer */
    /* ==================== */

    if( CheckUnits() == ON ) {
       if( CheckUnitsType() == SI) {
           d1 = DefaultUnits("Pa");
           d2 = DefaultUnits("m");
       }
       else {
           d1 = DefaultUnits("psi");
           d2 = DefaultUnits("in");
       }

       UnitsMultRep( &(s->spColUnits[0]), d1, d2 );
       UnitsCopy( &(s->spColUnits[1]), &(s->spColUnits[0]) );
       UnitsCopy( &(s->spColUnits[2]), &(s->spColUnits[0]) );
       UnitsMultRep( &(s->spColUnits[3]), &(s->spColUnits[0]), d2 );
       UnitsCopy( &(s->spColUnits[4]), &(s->spColUnits[3]) );
       UnitsCopy( &(s->spColUnits[5]), &(s->spColUnits[3]) );

       ZeroUnits( &(s->spRowUnits[0]) );
       ZeroUnits( &(s->spRowUnits[1]) );
       ZeroUnits( &(s->spRowUnits[2]) );
       UnitsCopy( &(s->spRowUnits[3]), d2 ); 
       UnitsCopy( &(s->spRowUnits[4]), d2 ); 
       UnitsCopy( &(s->spRowUnits[5]), d2 ); 

       UnitsCopy( &(s->spColUnits[6]), &(s->spColUnits[0]) ); 
       UnitsCopy( &(s->spColUnits[7]), &(s->spColUnits[1]) ); 
       UnitsCopy( &(s->spColUnits[8]), &(s->spColUnits[2]) ); 
       UnitsCopy( &(s->spColUnits[9]), &(s->spColUnits[3]) ); 
       UnitsCopy( &(s->spColUnits[10]), &(s->spColUnits[4]) ); 
       UnitsCopy( &(s->spColUnits[11]), &(s->spColUnits[5]) ); 

       UnitsCopy( &(s->spRowUnits[6]), &(s->spRowUnits[0]) ); 
       UnitsCopy( &(s->spRowUnits[7]), &(s->spRowUnits[1]) ); 
       UnitsCopy( &(s->spRowUnits[8]), &(s->spRowUnits[2]) ); 
       UnitsCopy( &(s->spRowUnits[9]), &(s->spRowUnits[3]) ); 
       UnitsCopy( &(s->spRowUnits[10]), &(s->spRowUnits[4]) ); 
       UnitsCopy( &(s->spRowUnits[11]), &(s->spRowUnits[5]) ); 

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
    }

    /* ===================== */
    /* Rotate stiff matrix   */
    /* ===================== */

    s->uMatrix.daa = (double **) rotate3d(s->uMatrix.daa, rot, no_dof);

#ifdef DEBUG
  printf("flag: in beamst3d() : STIFF after rotation \n");
  MatrixPrintIndirectDouble( s );
#endif
    
     return(s);
}


/* =============================================================== */
/* function   beamms3d()                                           */
/* 3D Frame Mass Matrix: This function returns a 12x12 rotated     */
/* mass matrix for 3D frame elmts.                                 */
/*                                                                 */
/* Input:                                                          */
/*   mbar->Density * X_area ;                                      */ 
/*   s->matrix pointer;                                            */
/*   xl->length       ; rot->rotation_matrix(12x12)                */
/*   nst->num_of_deg_of free per element; ndf->dof's per node;     */
/* =============================================================== */

#ifdef __STDC__
MATRIX *beamms3d(ARRAY *p, MATRIX *s, int mtype, double mbar, double xl, double rg , double **rot, int nst, int no_dof)
#else
MATRIX *beamms3d(p, s, mtype, mbar, xl, rg , rot, nst, no_dof)
ARRAY  *p;
MATRIX *s;
double **rot;
double mbar, rg , xl;
int    nst, no_dof, mtype;
#endif
{
int    n,m, i, j, k;
double t;
DIMENSIONS *d1,*d2,*d3;

#ifdef DEBUG
       printf("INFO >> In beamms3d() : no_dof= %d   \n",no_dof);
#endif

    switch(mtype) {
        case LUMPED:
	     t = mbar * xl/2;
             s->uMatrix.daa[0][0] = s->uMatrix.daa[1][1] = s->uMatrix.daa[2][2]  =  t;
	     s->uMatrix.daa[3][3] = t*rg*rg;
	     s->uMatrix.daa[4][4] = s->uMatrix.daa[5][5] = t*xl*xl/12.0;

             s->uMatrix.daa[no_dof][no_dof] = s->uMatrix.daa[no_dof+1][no_dof+1] = s->uMatrix.daa[no_dof+2][no_dof+2]  =  t;
	     s->uMatrix.daa[no_dof+3][no_dof+3] = t*rg*rg;
	     s->uMatrix.daa[no_dof+4][no_dof+4] = s->uMatrix.daa[no_dof+5][no_dof+5] = t*xl*xl/12.0;

             break;

        case CONSISTENT:

	    /* -------------------------------------- */
	    /* M terms =|  Mjj   Mjk |                */
	    /*          |  Mkj   Mkk |                */
	    /* -------------------------------------- */
	
	    /* -------------------------------------- */
	    /* Mjj  terms                             */
	    /* -------------------------------------- */

	    t = mbar*xl/420; 

	    s->uMatrix.daa[0][0]           =  t* 140;
	    s->uMatrix.daa[1][1] = s->uMatrix.daa[2][2] =  t* 156;
	    s->uMatrix.daa[3][3]           =  t* 140* rg * rg ;
	    s->uMatrix.daa[4][4] = s->uMatrix.daa[5][5] =  t* 4 * xl * xl;

	    s->uMatrix.daa[4][2] = s->uMatrix.daa[2][4] = -  t* 22 * xl;
	    s->uMatrix.daa[5][1] = s->uMatrix.daa[1][5] =    t* 22 * xl;

	    for(m=1; m <= no_dof; m++){
	        for(n = 1; n <= no_dof; n++){
	            if(n > m)                        /* sym terms */
	               s->uMatrix.daa[m-1][n-1] = s->uMatrix.daa[n-1][m-1] ;         
	            s->uMatrix.daa[n+no_dof-1][m+no_dof-1] = s->uMatrix.daa[n-1][m-1]; /* Mkk */         
	        }
	    }

	    /* ------------------------------------ */
	    /* off diagonal terms & sign adjustment */
	    /* ------------------------------------ */

	    s->uMatrix.daa[4+no_dof][2+no_dof] = s->uMatrix.daa[2+no_dof][4+no_dof] =   t* 22 * xl;
	    s->uMatrix.daa[5+no_dof][1+no_dof] = s->uMatrix.daa[1+no_dof][5+no_dof] = - t* 22 * xl;

	    /*-------------------------------------------*/
	    /* Mkj, Mjk terms - off main diagonal terms  */
	    /*-------------------------------------------*/

            s->uMatrix.daa[no_dof][0]                          =  t* 70;
            s->uMatrix.daa[no_dof+ 1][1] = s->uMatrix.daa[no_dof+2][2]   =  t* 54;
            s->uMatrix.daa[no_dof+ 3][3]                       =  t*  70 * rg * rg;
            s->uMatrix.daa[no_dof+ 4][4] = s->uMatrix.daa[no_dof+ 5][5]  = - t* 3 * xl * xl;
            s->uMatrix.daa[no_dof+ 4][2] = s->uMatrix.daa[no_dof+ 1][5]  =   t* 13 * xl;
            s->uMatrix.daa[no_dof+ 5][1] = s->uMatrix.daa[no_dof+ 2][4]  = - t* 13 * xl;

            /* Symmmetric terms */

            for(m = 1; m <= no_dof; m++)
                for(n = 1; n <= no_dof; n++)
                    s->uMatrix.daa[n-1][m+no_dof-1] = s->uMatrix.daa[n+no_dof-1][m-1];

	    break;

        default:
             FatalError("In elmt_frame3d() : beamms() : Type of Mass Matrix Undefined",(char *)NULL);
             break;
    }


#ifdef DEBUG
    printf("elmt mass[12][12]");
    MatrixPrintIndirectDouble(s);
#endif

    /* ---------------------- */
    /* initial unit of matrix */
    /* ---------------------- */

    if( CheckUnits() == ON ) {
       if( CheckUnitsType() == SI) {
           d1 = DefaultUnits("Pa");
           d2 = DefaultUnits("m");
       }
       else {
           d1 = DefaultUnits("psi");
           d2 = DefaultUnits("in");
       }
       d3 = DefaultUnits("sec");

       UnitsMultRep( &(s->spColUnits[0]), d1, d2 );
       UnitsCopy( &(s->spColUnits[1]), &(s->spColUnits[0]) );
       UnitsCopy( &(s->spColUnits[2]), &(s->spColUnits[0]) );
       UnitsMultRep( &(s->spColUnits[3]), &(s->spColUnits[0]), d2 );
       UnitsCopy( &(s->spColUnits[4]), &(s->spColUnits[3]) );
       UnitsCopy( &(s->spColUnits[5]), &(s->spColUnits[3]) );

       UnitsPowerRep( &(s->spRowUnits[0]), d3, 2.0, NO );
       UnitsCopy( &(s->spRowUnits[1]), &(s->spRowUnits[0]) );
       UnitsCopy( &(s->spRowUnits[2]), &(s->spRowUnits[0]) );
       UnitsMultRep( &(s->spRowUnits[3]), d2, &(s->spRowUnits[0]) );
       UnitsCopy( &(s->spRowUnits[4]), &(s->spRowUnits[3]) );
       UnitsCopy( &(s->spRowUnits[5]), &(s->spRowUnits[3]) );

       UnitsCopy( &(s->spColUnits[6]), &(s->spColUnits[0]) );
       UnitsCopy( &(s->spColUnits[7]), &(s->spColUnits[1]) );
       UnitsCopy( &(s->spColUnits[8]), &(s->spColUnits[2]) );
       UnitsCopy( &(s->spColUnits[9]), &(s->spColUnits[3]) );
       UnitsCopy( &(s->spColUnits[10]), &(s->spColUnits[4]) );
       UnitsCopy( &(s->spColUnits[11]), &(s->spColUnits[5]) );

       UnitsCopy( &(s->spRowUnits[6]), &(s->spRowUnits[0]) );
       UnitsCopy( &(s->spRowUnits[7]), &(s->spRowUnits[1]) );
       UnitsCopy( &(s->spRowUnits[8]), &(s->spRowUnits[2]) );
       UnitsCopy( &(s->spRowUnits[9]), &(s->spRowUnits[3]) );
       UnitsCopy( &(s->spRowUnits[10]), &(s->spRowUnits[4]) );
       UnitsCopy( &(s->spRowUnits[11]), &(s->spRowUnits[5]) );

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) d3->units_name);
       free((char *) d3);
    }

    /* ---------------------- */
    /* rotate to global axis  */
    /* ---------------------- */

    s->uMatrix.daa = (double **) rotate3d(s->uMatrix.daa, rot, no_dof);

    return(s);
}


/* ====================================================== */
/* function   tmat(rot, type, p)                          */
/*    Transformation matrix for 2,3d element              */
/*    Generates a transformation matrix                   */ 
/* Input:                                                 */
/*        rot : rotation matrix (12 x 12)                 */
/*        type: elmt type                                 */
/*           p: ARRAY * storing elmt info                 */
/* ====================================================== */

#ifdef __STDC__
double **tmat(double **rot, int type, ARRAY *p)
#else
double **tmat(rot, type, p)
double **rot;
int    type;
ARRAY  *p;
#endif
{
int    i,j,n,sr,ii,jj,kk,nf, ndff;
double t[4][4];
double alpha, cx,cy,cz,csa,sna,elone, el;

    alpha = p->eangle.value;
    sr    = p->size_of_stiff;

    for(i = 0; i <= 2; i++)
        for(j = 0; j <= 2; j++)
            t[i][j] = 0.0;

    for(i = 1; i <= sr; i++)
        for(n =1; n<=sr;n++)
            rot[i-1][n-1] = 0.0;

    cx = p->coord[0][1].value - p->coord[0][0].value;
    cy = p->coord[1][1].value - p->coord[1][0].value;

    switch(type) {
        case 1: /* plane truss */
	case 2:
	case 3:
	case 4: /* plane frame, with axial force */
             el = sqrt(cx * cx + cy * cy );
             cx = cx/el; 
             cy = cy/el;
             if(type == 4) {
                t[0][0] = cx;
                t[0][1] = cy;
                t[1][0] = -cy;
                t[1][1] = cx;
                if(type == 1) break;
                   t[2][2] = 1;
             }
             else {
                t[0][0] = cx;
                t[0][2] = cy;
                t[1][1] = 1;
                t[2][0] = -cy;
                t[2][2] = cx;
             } 
             break;
        case 5: /* space truss */
        case 6: /* space frame */
             cz = p->coord[2][1].value - p->coord[2][0].value;
             el = sqrt(cx * cx + cy * cy + cz * cz);
             cx = cx/el; 
             cy = cy/el;
             cz = cz/el;
             if(type == 5){
                t[0][0] = cx*cx;
                t[0][1] = cx*cy;
                t[0][2] = cx*cz;
                t[1][1] = cy*cy;
                t[1][2] = cy*cz;
                t[2][2] = cz*cz;
                for(i = 2; i <= 3; i++){
                    ii = i -1;
                    for(j = 1; j <= ii; j++) 
                        t[i-1][j-1] = t[j-1][i-1];
                }
             }
             else {
                csa = cos(alpha);
                sna = sin(alpha);
                t[0][1] = cy;
                if(cx == 0.0 && cz == 0.0){  /* vertical space frame */
                   t[1][0] = -cy*csa;
                   t[1][2] = sna;
                   t[2][0] = cy * sna;
                   t[2][2] = csa;
                }
		else {
                   elone = sqrt(cx * cx + cz * cz);
                   t[0][0] = cx;
                   t[0][2] = cz;
                   t[1][0] = (-cx*cy*csa - cz*sna)/elone;
                   t[1][2] = (-cy*cz*csa + cx*sna)/elone;
                   t[1][1] = csa * elone;
                   t[2][0] = ( cx*cy*sna - cz*csa)/elone;
                   t[2][1] = -sna * elone;
                   t[2][2] = ( cy*cz*sna + cx*csa)/elone;
               }
            }
            break;
        default:
            break;
    }

    if(type == 6){
       nf = 4;
       ndff = 3;
    }
    else{
       nf   = 2;
       ndff = p->dof_per_node;
    }

    for(i=1;i<=nf;i++) {
        jj = (i-1) * ndff; 
        kk = (i-1) * ndff;
        for(j=1;j<=ndff;j++){
            jj = jj+ 1;
            for(n=1;n<=ndff;n++){
                kk = kk  + 1;
                rot[jj-1][kk-1] = t[j-1][n-1];
            }
            kk = (i-1) * ndff;
        }
    }

    return(rot);
}


/* ===================================================================== */
/* MATRIX *rotate3d(s, r ,ndf)                                           */
/* Rotation matrix  for 3d frame element                                 */
/* Multiplies  the stiffness matrix "s" by the transformation matrix "r" */
/* (ie. Kg = Tt*Kl*T)  Return Kg                                         */
/* ===================================================================== */

#ifdef __STDC__
double **rotate3d(double **s, double **r , int no_dof)
#else
double **rotate3d(s, r ,no_dof)
double **s, **r;
int    no_dof;
#endif
{
int    i,j,n;
double t;
double **rt, **asj, **sj, **sij, **kj;

     if(r[0][0] == 1.0) return(s);

  /* transpose of rotation matrix  r */

     rt = MatrixAllocIndirectDouble(no_dof, no_dof);
     for(i = 1; i <= no_dof; i++)
       for(n =1; n<= no_dof; n++)
           rt[n-1][i-1] = r[i-1][n-1];

  /* sii terms of transformed stiffness matrix s */

     kj = MatrixAllocIndirectDouble(no_dof, no_dof);
     for(i = 1; i <= no_dof; i++)
         for(n =1; n<= no_dof; n++)
             kj[i-1][n-1] = s[i-1][n-1];

     asj = (double **) dMatrixMult( kj, no_dof, no_dof, r, no_dof, no_dof);
     sj  = (double **) dMatrixMult( rt, no_dof, no_dof, asj, no_dof, no_dof);

     for(i=1;i<= no_dof; i++)
         for(n=1;n<=no_dof; n++)
             s[i-1][n-1] = sj[i-1][n-1];

     MatrixFreeIndirectDouble(asj, no_dof);
     MatrixFreeIndirectDouble( sj, no_dof);

  /* sjj terms of transformed stiffness matrix s */

     for(i = 1; i <= no_dof; i++)
         for(n =1; n<=no_dof;n++)
             kj[i-1][n-1] = s[i+no_dof-1][n+no_dof-1];

     asj = (double **) dMatrixMult( kj, no_dof, no_dof, r, no_dof, no_dof);
     sj  = (double **) dMatrixMult( rt, no_dof, no_dof, asj, no_dof, no_dof);

     for(i=1;i<= no_dof; i++)
         for(n=1;n<=no_dof; n++)
             s[i+no_dof-1][n+no_dof-1] = sj[i-1][n-1];

     MatrixFreeIndirectDouble(asj, no_dof);
     MatrixFreeIndirectDouble( sj, no_dof);

  /* sij terms of transformed stiffness matrix s */

     for(i = 1; i <= no_dof; i++)
         for(n =1; n<= no_dof; n++)
             kj[i-1][n-1] = s[i-1][n+no_dof-1];

     asj = (double **) dMatrixMult( kj, no_dof, no_dof, r, no_dof, no_dof);
     sj  = (double **) dMatrixMult( rt, no_dof, no_dof, asj, no_dof, no_dof);

     for(i=1;i<= no_dof; i++)
         for(n=1;n<= no_dof ;n++)
             s[i-1][n+no_dof-1] = s[n+no_dof-1][i-1] = sj[i-1][n-1];

     MatrixFreeIndirectDouble(asj, no_dof);
     MatrixFreeIndirectDouble( sj, no_dof);
     MatrixFreeIndirectDouble( kj, no_dof);
     MatrixFreeIndirectDouble( rt, no_dof);

     return(s);
}


/* 
 *  =============================================================
 *  print_property_frame_3d() : Print FRAME_3D Element Properties
 *  =============================================================
 */ 

#ifdef __STDC__
void print_property_frame_3d(EFRAME *frp, int i)
#else
void print_property_frame_3d(frp, i)
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
        UnitsSimplify( eap->work_section[1].dimen );
        UnitsSimplify( eap->work_section[2].dimen );
        UnitsSimplify( eap->work_section[12].dimen );
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
	if( eap->work_section[1].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Inertia Iyy     = %16.3e %s\n",
                           eap->work_section[1].value/eap->work_section[1].dimen->scale_factor,
                           eap->work_section[1].dimen->units_name);
	}
	if( eap->work_section[2].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Inertia Izz     = %16.3e %s\n",
                           eap->work_section[2].value/eap->work_section[2].dimen->scale_factor,
                           eap->work_section[2].dimen->units_name);
	}
	if( eap->work_section[12].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Torsional Constant J = %16.3e %s\n",
                           eap->work_section[12].value/eap->work_section[12].dimen->scale_factor,
                           eap->work_section[12].dimen->units_name);
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
        if( eap->work_section[1].value != 0.0 ) {
           printf("             ");
           printf("         : Inertia Iyy     = %16.3e\n",
                            eap->work_section[1].value);
        }
        if( eap->work_section[2].value != 0.0 ) {
           printf("             ");
           printf("         : Inertia Izz     = %16.3e\n",
                            eap->work_section[2].value);
        }
        if( eap->work_section[12].value != 0.0 ) {
           printf("             ");
           printf("         : Torsional Constant J = %16.3e\n",
                            eap->work_section[12].value);
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


/* ================================================== */
/* Equivalent Loading Procedure for 3 D frame element */
/* ================================================== */

#ifdef __STDC__
ARRAY *sld05(ARRAY *p, int task)
#else
ARRAY *sld05(p,task)
ARRAY *p;
int   task;
#endif
{
ELEMENT_LOADS   *elsptr;
ELOAD_LIB          *elp;
double  P, a ,b ;
double  L, load[15];
double  px,py,pz, mx,my,mz, bx,by,bz, ze,ze2,ze3;
double  X1,Y1,Z1,MX1,MY1,MZ1,
        X2,Y2,Z2,MX2,MY2,MZ2, 
        MCZ, MCY, MCZT, MCYT;
double  f1,f2,f3,f4,f5,f6,
        df1x, df3x, df5x,df2x, df6x;
double  **rmat, **rt;
int     *da;
int     inc,i,j;

/* printf(">>>>In sld05 : 3d elmt    \n");  */

       /* Initialize total load */

       for(inc=1; inc <= 12; inc++)
          load[inc-1] = 0.0;
       MCZT = 0.0;
       MCYT = 0.0;

    switch(task){
        case PRESSLD: case STRESS:
             L = p->length.value;  elsptr  =  p->elmt_load_ptr;

             for(j=1; j<= elsptr->no_loads_faces; j++) {
                 elp  =  &elsptr->elib_ptr[j-1];

                 P =  elp->P.value;
                 a =  elp->a.value;
                 b =  elp->b.value;

                 /* Pt loads */

                 px =  elp->px.value;
                 py =  elp->py.value;
                 pz =  elp->pz.value;

                 /* moments */

                 mx =  elp->mx.value;
                 my =  elp->my.value;
                 mz =  elp->mz.value;

                 /* distributed loading */

                 bx =  elp->bx.value;
                 by =  elp->by.value;
                 bz =  elp->bz.value;

                 if(a > L) /* error message */
                    printf(">>ERROR in sld; Elmt Load dist. 'a' > Elmt Length; El_no= %d\n",p->elmt_no);  

                 if (elp->type == -1) { /* Distributed loading Condition */
                     inc = 0;

                 /* set default values */
                 if(b == 0.0) b = L; /* dist loading acts on entire length */

                 /* first calc f(b) */
                    ze = b/L;    
SHP_START:
                    ze2 = ze * ze; ze3 = ze2 * ze;
                    f1 =    1   -  ze2/2;
                    f2 =    ze3*ze/2  - ze3 + ze ;
                    f3 =    (ze3*ze/4 - 2*ze3/3 + ze2/2) * L;
                    f4 =    ze2/2 ;
                    f5 =    -ze *ze3/2 +  ze3;
                    f6 =    (ze3*ze/4  - ze3/3) * L;
                    inc++;
                 
                    if(inc == 1) {
                       /* temp hold of values f(b) */
                       X1 = f1; Y1 = f2; Z1 = f3; 
                       X2 = f4; Y2 = f5; Z2 = f6; 
                       ze = a/L;
                       goto SHP_START;
                    }
                    else{
                       /* f() = f(b) - f(a)  */
                       f1 = X1 - f1; f2 =  Y1 - f2; f3 = Z1 - f3; 
                       f4 = X2 - f4; f5 =  Y2 - f5; f6 = Z2 - f6; 
                    }

                    X1 = bx * f1 * L;
                    Y1 = by * f2  * L;
                    Z1 = bz * f2 * L;
                    MX1 = 0.0 ;
                    MY1 =-bz * f3 * L;
                    MZ1 = by * f3 * L;

                    X2 = bx * f4 * L;
                    Y2 = by * f5 * L;
                    Z2 = bz * f5 * L;
                    MX2 = 0.0 ;
                    MY2 =-bz * f6 * L;
                    MZ2 = by * f6 * L;

                    if (task == STRESS){
                        /* +ve  simply support moment at center */
                        if(b==L && a== 0.0) {/* udl acting on entire length */
                           MCZ = -by * (L * L)/8;   
                           MCY = -bz * (L * L)/8;   
                        }
                    else { /* approximate mom at center */
                        MCZ = -by *(b-a)* (L - (a+b)/2)/2;   
                        MCY = -bz *(b-a)* (L - (a+b)/2)/2;   
                    }
                 }
               }      /* end of dist loading */

               else {                 /* Concentrated Loading Condition */
                   /* shape functions */
                  ze = a/L;      ze2 = ze * ze;        ze3 = ze2 * ze;

                  f1 =     1  -  ze;
                  f2 =     2 * ze3  -  3 * ze2 + 1;
                  f3 =    (ze3 - 2 * ze2 + ze) * L;
                  f4 =     ze ;
                  f5 =    -2 * ze3  +  3 * ze2;
                  f6 =    (ze3 - ze2 ) * L;

                  if(my != 0.0 || mz != 0.0){
                     /* derivatives of shape function */  
                     df2x =  6 *( ze2  -  ze) / L;
                     df3x =  3 *ze2  - 4*ze  + 1;
                     df5x = - df2x; 
                     df6x =  3 *ze2  - 2*ze;
                  }

                  X1 = px * f1 + 0;
                  Y1 = py * f2 + mz *  df2x;
                  Z1 = pz * f2 - my *  df2x;
                  MX1 =       mx * f1 ;
                  MY1 =-pz * f3 + my *  df3x;
                  MZ1 = py * f3 + mz *  df3x;

                  X2 = px * f4 + 0;
                  Y2 = py * f5 + mz *  df5x;
                  Z2 = pz * f5 - my *  df5x;
                  MX2 =       mx * f4 ;
                  MY2 =-pz * f6 + my *  df6x;
                  MZ2 = py * f6 + mz *  df6x;

                  if(task == STRESS) { /* +ve  simply support moment at center */
                     MCY = -pz * (L -a)/2 ;   
                     MCZ = -py * (L -a)/2 ;   
                  }
              }

              /* Add Contributation to Total Equivalent Load */

                 load[0] = load[0] + X1;
                 load[1] = load[1] + Y1;
                 load[2] = load[2] + Z1;
                 load[3] = load[3] + MX1;
                 load[4] = load[4] + MY1;
                 load[5] = load[5] + MZ1;
                 load[6] = load[6] + X2;
                 load[7] = load[7] + Y2;
                 load[8] = load[8] + Z2;
                 load[9] = load[9] + MX2;
                 load[10] = load[10] + MY2;
                 load[11] = load[11] + MZ2;

                 if(task == STRESS) { /* mid pt moment */
                    MCYT = MCYT+ MCY;
                    MCZT = MCZT+ MCZ;
                 }
              }
              break;
         default:
              break;
    }
  
    /* Rotate to Global */
 
    rmat = (double **) MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
    rmat = (double **) tmat(rmat,6,p);
    rt   = (double **)MatrixAllocIndirectDouble(p->size_of_stiff, p->size_of_stiff);
    for( i=1 ; i<=p->size_of_stiff ; i++ )
       for( j=1 ; j<=p->size_of_stiff ; j++ )
          rt[j-1][i-1] = rmat[i-1][j-1];

    for (inc=1; inc <= p->size_of_stiff; inc++){
         p->nodal_loads[inc-1].value = 0.0;
         for (j=1; j <= p->size_of_stiff; j++)
              p->nodal_loads[inc-1].value = p->nodal_loads[inc-1].value + rt[inc-1][j-1]* (double) load[j-1];
    }

    MatrixFreeIndirectDouble(rmat, p->size_of_stiff);
    MatrixFreeIndirectDouble(rt, p->size_of_stiff);

    return(p);
} 

