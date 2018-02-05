/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_plate.c : Quadrilateral DKT Plate Bending Element
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
 *  Element 8 - Quadrilateral DKT Plate Bending Element               
 *            - Written by B.K. Voon (1990)                           
 *                                                                    
 *  Note : 2-dimensional plate element from Berkeley.               
 *       : Only works for elements lying in x-z plane. Here we      
 *         assume that y-axis corresponds to vertical direction.   
 *  ============================================================================= 
 */

#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "fe_functions.h"
#include "elmt.h"

/*
#define DEBUG 
*/

/* ================================================== */
/* elmt_plate() : Quadrilateral DKT Plate Bending Element */
/* ================================================== */
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
/* ---------------------------------------------------------------*/
/*    p->work_section[i] for i>11 is used to fit each elmts needs.*/
/*    Therefore, p->work_section[i] may be different for differnt */
/*    elmts when i > 11. For DKQ(T) plate element :               */ 
/* ---------------------------------------------------------------*/
/*    p->work_section[12] = D = (E/(1-nu^2))*t^3/12;              */
/*    p->work_section[13] = D*nu;                                 */
/*    p->work_section[14] = D*(1-nu)/2;                           */

/*    p->work_section[16] = flag   ?                              */
/* ============================================================== */
ARRAY *elmt_plate(p,isw)
ARRAY *p;
int isw;
{
static double nu;
static QUANTITY fy, G, E, ET, density;
static double  Ixx, Iyy, Izz, Ixy, Ixz, Iyz, bf, tf, A, depth, weight, EA, thickness;

double ia[3][27],eps[4],sig[4],aa[5],bb[5],cc[5],dd[5],ee[5];
double sg[10],tg[10],wg[10];
double d1,d2,d3,d4,d5,d6,thick,xsj,dn1,dn2,dn3;
double        xx,yy;
double        **bmu;
double        **shp;
double         **bm;
int i,j,l,lint,i1,k;
int    length1, length2, unit_length;
int          UNITS_SWITCH, UnitsType;
int           ii, jj, kk, k1, k2, k3;
DIMENSIONS  *dimen, *dimen1, *dimen2;

   H_Print = 0;
   UNITS_SWITCH = CheckUnits();
   UnitsType    = CheckUnitsType();

#ifdef DEBUG
       printf("In elmt_plate() : starting with isw = %10d\n", isw);
#endif

   if(isw != PROPTY) {
      shp = MatrixAllocIndirectDouble(3,8);
      bm  = MatrixAllocIndirectDouble(3,12);  
   }

   switch(isw) {
       case PROPTY:   /* Input The Material Properties */
            if(p->work_material[16].value == 0.0 ) {
              /* d[10] has been changed to work_material[16].value    */
              /* what is d[10] any way ??? a flag ? array start at 1  */
              /* p->d[1];   Modulus    : p->d[2];   Poisson   */
              /* p->d[3];   Thickness  : p->d[4];   q-loading */

             E.value       =  p->work_material[0].value;
             fy.value      =  p->work_material[2].value;
             ET.value      =  p->work_material[3].value;
             nu            =  p->work_material[4].value;
             density.value =  p->work_material[5].value;
             thickness     =  p->work_section[11].value;
             p->work_section[12].value
                           = E.value/(1.0-nu*nu)*thickness*thickness*thickness/12.0; 
             p->work_section[13].value 
                           = nu*p->work_section[12].value; 
             p->work_section[14].value 
                           = 0.5*(1-nu)*p->work_section[12].value; 
             p->work_section[16].value = 1.0;

             switch(UNITS_SWITCH) {
                case ON:
                   E.dimen       =  p->work_material[0].dimen;
                   fy.dimen      =  p->work_material[2].dimen;
                   ET.dimen      =  p->work_material[3].dimen;
                   density.dimen =  p->work_material[5].dimen;

                   dimen            = UnitsPower( p->work_section[11].dimen, 3.0, NO );
                   UnitsMultRep( p->work_section[12].dimen, E.dimen, dimen );
                   free((char *) dimen->units_name);
                   free((char *) dimen);
                   UnitsCopy( p->work_section[13].dimen, p->work_section[12].dimen );
                   UnitsCopy( p->work_section[14].dimen, p->work_section[12].dimen );
                break;
                case OFF:
                break;
                default:
                break;
             }
     

/* ---------------------------------------------------------*/
/* In fortran FEAP version                                  */
/* construct rotation parameters: u_x = 1, u_y = 2          */
/*                      ia[1][iel] = 1,  iel = elmt_type no */
/*                      ia[2][iel] = 2                      */
/* construct rotation parameters: theta_x = 4, theta_y = 5  */
/*                      ir[1][iel] = 4,  iel = elmt_type no */
/*                      ir[2][iel] = 5                      */
/* ---------------------------------------------------------*/
/*-----
               ia[0][p->eiel-1] = 2;         
               ia[1][p->eiel-1] = 3;        
-------*/

#ifdef DEBUG
               printf("In elmt_plate() : modulus   = %10.4f\n",E.value);
               printf("            : poisson   = %10.4f\n",nu);
               printf("            : thickness = %10.4f\n",thickness);
               printf("            : D         = %10.4f\n", p->work_section[12].value);
               printf("            : D*v       = %10.4f\n", p->work_section[13].value);
               printf("            : D(1-v)/2  = %10.4f\n", p->work_section[14].value);
               printf("            : p->work_section[16] = %10.4f\n\n", p->work_section[16].value);
#endif

            }
            break;
       case CHERROR:
       case MASS_MATRIX:
            break;
       case STIFF:   /* Compute The Element Tangent Array */

#ifdef DEBUG
       printf("In elmt_plate() : Start to compute Element Tangent Array\n");
#endif

            jacqud(p->coord,aa,bb,cc,dd,ee);
            l = 3;
            pgauss(l,&lint,sg,tg,wg);

            for(l=1; l <= lint; l++) {

#ifdef DEBUG
       printf("\n\n *** STARTING INT POINT : LINT %3d\n\n", l);
#endif

                shp = qushp8(sg[l-1],tg[l-1],shp,p->coord,&xsj);
                xsj = xsj*wg[l-1];
                bm  = dktqbm(shp,bm,aa,bb,cc,dd,ee);

                /* Compute Weighted Jacobian And D-matrix Constants */

                d1 = p->work_section[12].value * xsj;
                d2 = p->work_section[13].value * xsj;
                d3 = p->work_section[14].value * xsj;

                /* Compute The Element Load Vector */

/* nodal loads 
                for(i = 1; i <= p->nodes_per_elmt; i++) {
                    j = (i-1)*3+1; 
                    p->nodal_loads[j-1].value += (p->d[8])*xsj*(1.-sg[l-1])*(1.-tg[l-1]);
                }
*/
               /* Compute Contribution To Element Stiffness For This Point */

                for(i=1;i<=12;i++) {
                    dn1 = d1*bm[0][i-1]+d2*bm[1][i-1];
                    dn2 = d2*bm[0][i-1]+d1*bm[1][i-1];
                    dn3 = d3*bm[2][i-1];
                    for(j=i;j<=12;j++) 
                        p->stiff->uMatrix.daa[i-1][j-1] += dn1*bm[0][j-1]+dn2*bm[1][j-1]+dn3*bm[2][j-1];
                }
            }

            /* Make Stiffness Symmetric */

            for(i= 2; i<=12; i++) {
                i1 = i-1;
                for(j=1;j<=i1;j++)
                    p->stiff->uMatrix.daa[i-1][j-1] = p->stiff->uMatrix.daa[j-1][i-1];
            }

            /* ========================== feature not included in program yet ============
               Modify Load Vector for Non-Zero Displacement Boundary Conditions : p->u[12][1] 
            for(j=1;j<=p->size_of_stiff; j++) {
                if(p->displ->uMatrix.daa[j-1][0] != 0)
                   for(i=1;i<=p->size_of_stiff;i++)
                       p->nodal_loads[i-1].value -= p->stiff->uMatrix.daa[i-1][j-1]*p->displ->uMatrix.daa[j-1][0];

            } 
              ========================== end of feature =================================== */
           
            /**************************************************/
            /* Assign Units to Stiffness Matrix               */
            /**************************************************/

            /* Initiation of Stiffness Units Buffer           */

            if( CheckUnits() == ON ) {
               if( UnitsType == SI) {
                  dimen1 = DefaultUnits("Pa");
                  dimen2 = DefaultUnits("m");
               }
               else {
                  dimen1 = DefaultUnits("psi");
                  dimen2 = DefaultUnits("in");
               }

               /* node 1*/

               UnitsMultRep( &(p->stiff->spColUnits[0]), dimen1, dimen2 );
               UnitsMultRep( &(p->stiff->spColUnits[1]), &(p->stiff->spColUnits[0]), dimen2 );
               UnitsCopy( &(p->stiff->spColUnits[2]), &(p->stiff->spColUnits[1]) );

               ZeroUnits( &(p->stiff->spRowUnits[0]) ); 
               UnitsCopy( &(p->stiff->spRowUnits[1]), dimen2 ); 
               UnitsCopy( &(p->stiff->spRowUnits[2]), dimen2 ); 

               /* node i  i > 1*/

               for(i = 2; i <= p->nodes_per_elmt; i++) {
                   kk = p->dof_per_node*(i-1) + 1;
                   for(j = 1; j <= p->dof_per_node; j++) {
                       k  = p->dof_per_node*(i-1) + j;
                       if( k <= kk) {
                           UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[0]) );
                           UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[0]) );
                       }
                       if(k > kk) {
                           UnitsCopy( &(p->stiff->spColUnits[k-1]), &(p->stiff->spColUnits[1]) );
                           UnitsCopy( &(p->stiff->spRowUnits[k-1]), &(p->stiff->spRowUnits[1]) );
                       }
                   }
               }

               free((char *) dimen1->units_name);
               free((char *) dimen1);
               free((char *) dimen2->units_name);
               free((char *) dimen2);
            }
     
            break;
       case STRESS:   /* Compute And Output The Element Variables */

            jacqud(p->coord,aa,bb,cc,dd,ee);
            shp = qushp8(0.,0.,shp,p->coord,&xsj);
            bm  = dktqbm(shp,bm,aa,bb,cc,dd,ee);

            for(i=1; i<=3; i++) {
                eps[i-1] = 0.0;
                for(j=1; j <= 4; j++) {
                    for(k=1; k <= 3; k++) 
			eps[i-1] += bm[i-1][(int) (3*(j-1)+k-1)]*p->displ->uMatrix.daa[k-1][j-1];
                }
            }

            sig[0] = (p->work_section[12].value*eps[0])+(p->work_section[13].value*eps[1]);
            sig[1] = (p->work_section[13].value*eps[0])+(p->work_section[12].value*eps[1]);
            sig[2] =  p->work_section[14].value*eps[2]/2.0;

            xx = 0.25*(p->coord[0][0].value+p->coord[0][1].value+p->coord[0][2].value+p->coord[0][3].value);
            yy = 0.25*(p->coord[2][0].value+p->coord[2][1].value+p->coord[2][2].value+p->coord[2][3].value);

            H_Print = H_Print-1;
            if(H_Print <= 0) {
               H_Print = 50;
               printf("---------------------------------------------------------------------------\n");
               printf("DKQ Plate Bending Element\n");
               printf("elmt   mat   x-coord   z-coord      Mxx/length     Mxz/length    Mzz/legnth\n");
               printf("---------------------------------------------------------------------------\n");
               if(UNITS_SWITCH == ON) {   
                  printf("Units");
	          switch( UnitsType ) {
	            case US:
                       dimen = DefaultUnits("lbf");
                       break;
	            case SI:
	            case SI_US:
	            default:
                       dimen = DefaultUnits("N");
                       break;
                  }
                  printf("           %s        %s           %s            %s             %s",
                         p->coord[0][0].dimen->units_name,
                         p->coord[2][0].dimen->units_name,
                         dimen->units_name, dimen->units_name, dimen->units_name);
                  printf("\n");
               }
            }

            printf("%4d %s %9.3f %9.3f", p->elmt_no, p->material_name,
                   xx/p->coord[0][0].dimen->scale_factor,
                   yy/p->coord[2][0].dimen->scale_factor);

            printf("%14.5e %14.5e %14.5e\n",
                   sig[0]/dimen->scale_factor,
                   sig[1]/dimen->scale_factor,
                   sig[2]/dimen->scale_factor);
            printf("\n");
            free((char *) dimen->units_name);
            free((char *) dimen);

            break;
       case LOAD_MATRIX:   /* Compute The Element Residual Vector */
            break;
    }

    if(isw != PROPTY){
       MatrixFreeIndirectDouble(shp, 3);
       MatrixFreeIndirectDouble(bm, 3);
    }

    return(p);
}


/* ============= */
/* Jacqud(x,ndm) */
/* ============= */

int jacqud(x,aa,bb,cc,dd,ee)
QUANTITY **x;
double aa[5],bb[5],cc[5],dd[5],ee[5];
{
int i,k;
double b,c,sql;

   for(i=1 ; i<=4; i++) {
       k = i%4 + 1;

       b = x[2][k-1].value-x[2][i-1].value;
       c = x[0][i-1].value-x[0][k-1].value;
       sql = b*b+c*c;

       aa[i-1] = 1.5*c/sql;
       bb[i-1] = 0.75*b*c/sql;
       cc[i-1] = (0.25*c*c-0.5*b*b)/sql;
       dd[i-1] = -1.5*b/sql;
       ee[i-1] = (0.25*b*b-0.5*c*c)/sql;
   }

   return;
}


/* ============================= */
/* dktqbm(shm,bm,aa,bb,cc,dd,ee) */
/* ============================= */

double **dktqbm(shm,bm,aa,bb,cc,dd,ee)
double **shm,**bm;
double aa[5],bb[5],cc[5],dd[5],ee[5];
{
int i,j,i1,i2,i3;

   i1 = 1;
   for(i=1; i<=4; i++) {
      j  = (i+2)%4 +1;
      i2 = i1+1;
      i3 = i2+1;
      
      bm[0][i1-1] =  aa[i-1]*shm[0][i+4-1] - aa[j-1]*shm[0][j+4-1];
      bm[0][i2-1] =  bb[i-1]*shm[0][i+4-1] + bb[j-1]*shm[0][j+4-1];
      bm[0][i3-1] =  cc[i-1]*shm[0][i+4-1] + cc[j-1]*shm[0][j+4-1] - shm[0][i-1];
      bm[1][i1-1] =  dd[i-1]*shm[1][i+4-1] - dd[j-1]*shm[1][j+4-1];
      bm[1][i2-1] = -ee[i-1]*shm[1][i+4-1] - ee[j-1]*shm[1][j+4-1] + shm[1][i-1];
      bm[1][i3-1] = -bb[i-1]*shm[1][i+4-1] - bb[j-1]*shm[1][j+4-1];
      bm[2][i1-1] =  aa[i-1]*shm[1][i+4-1] - aa[j-1]*shm[1][j+4-1]
                    +dd[i-1]*shm[0][i+4-1] - dd[j-1]*shm[0][j+4-1];
      bm[2][i2-1] = -ee[i-1]*shm[0][i+4-1] - ee[j-1]*shm[0][j+4-1] + shm[0][i-1] - bm[1][i3-1];
      bm[2][i3-1] =  cc[i-1]*shm[1][i+4-1] + cc[j-1]*shm[1][j+4-1] - shm[1][i-1] - bm[0][i2-1];
      i1 = i1+3;
   }

   return(bm);
} 


/* ========================== */
/* qushp8(s,t,shp,x,xsj)      */
/* ========================== */

double **qushp8(s,t,shp,x,xsj) 
double           s,t;
double         **shp;
QUANTITY    **x;
double          *xsj;
{
int i;
double xs,xt,ys,yt,ss,tt,sn,tn,si[5],ti[5]; 

#ifdef DEBUG
       printf("*** In qushp8() : s   = %10.5f t   = %10.5f\n", s,t);
#endif

   si[0] = -1.0;
   si[1] =  1.0;
   si[2] =  1.0;
   si[3] = -1.0;

   ti[0] = -1.0;
   ti[1] = -1.0;
   ti[2] =  1.0;
   ti[3] =  1.0;

   xs = 0.0; xt = 0.0;
   ys = 0.0; yt = 0.0;

   for(i=1;i<=4;i++) {
       ss = si[i-1]*s;
       tt = ti[i-1]*t;
       sn = si[i-1]*(1.+tt);
       tn = ti[i-1]*(1.+ss);
       xs = xs+sn*x[0][i-1].value;
       xt = xt+tn*x[0][i-1].value;
       ys = ys+sn*x[2][i-1].value;
       yt = yt+tn*x[2][i-1].value;

       shp[0][i-1] = 0.25*sn*(ss+ss+tt);
       shp[1][i-1] = 0.25*tn*(ss+tt+tt);
       shp[2][i-1] = 0.25*(1.+ss)*(1.+tt)*(-1.+ss+tt);
   }

   *xsj = (xs*yt-ys*xt)/4.0;

   xs = xs/ *xsj;
   xt = xt/ *xsj;
   ys = ys/ *xsj;
   yt = yt/ *xsj;
   *xsj = *xsj/4.0;

   for(i=5;i<=7;i=i+2) {
       ss = si[i-4-1]*s;
       tt = ti[i-4-1]*t;
       shp[0][i-1] =     -s*(1.+tt);
       shp[1][i-1] =    0.5*ti[i-4-1]*(1.-s*s);
       shp[2][i-1] =    0.5*(1.-s*s)*(1.+tt);
       shp[0][i] = -0.5*si[i-4-1]*(1.-t*t);
       shp[1][i] =   -t*(1.-ss);
       shp[2][i] =  0.5*(1.-ss)*(1.-t*t);
   }

   for(i=1; i<= 8; i++) {
      sn = yt*shp[0][i-1]-ys*shp[1][i-1];
      shp[1][i-1] = xs*shp[1][i-1]-xt*shp[0][i-1];
      shp[0][i-1] = sn;
   }

#ifdef DEBUG
       printf("\n");
       printf("In Qushp8 () : yt = %12.5f ys = %12.5f\n", yt, ys);
       printf("In Qushp8 () : xt = %12.5f xs = %12.5f\n", xt, xs);
       dMatrixPrint("shape func", shp, 3, 8);
#endif

   return(shp);
}

/*
 *  ===============================================
 *  Print DKT_PLATE Element Properties
 *  ===============================================
 */

#ifdef __STDC__
void print_property_plate(EFRAME *frp, int i)
#else
void print_property_plate(frp, i)
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
	if( eap->work_section[11].dimen->units_name != NULL ) {
          printf("             ");
          printf("         : Plate Thickness = %16.3e %s\n",
                           eap->work_section[11].value/eap->work_section[11].dimen->scale_factor,
                           eap->work_section[11].dimen->units_name);
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
        if( eap->work_section[11].value != 0.0 ) {
           printf("             ");
           printf("         : Plate Thickness = %16.3e\n",
                            eap->work_section[11].value);
        }
        break;
        default:
        break;
     }
}

ARRAY *sld08(p, isw)
ARRAY *p;
int isw;
{
    printf("ERROR >> ***In sld08() : elmt no =%3d : isw= %3d\n",p->elmt_no, isw);
    return(p);
}
