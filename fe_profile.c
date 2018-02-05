/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_profile.c : Major fuctions for FE Profiler
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
 *  Written by: Mark Austin, Xiaoguang Chen, and Wane-Jang Lin           May 1997
 *  ============================================================================= 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "vector.h"
#include "elmt.h"

/*
#define DEBUG 
*/

static int max_size_of_stiff;
static int max_no_dof;
static int max_nodes_per_elmt;
static int max_no_dimen;

/* ================================== */
/* profile() :  Input: EFRAME *frp    */
/* ================================== */

#ifdef __STDC__
EFRAME *profile(EFRAME *frame, int prt, int iop)
#else
EFRAME *profile(frame, prt, iop)
EFRAME *frame;
int    prt, iop;
#endif
{
int n, i,j,ii,jj,m,neq,nad;
int idl_vector[30];
int nmin;

#ifdef DEBUG
       printf("*** Enter profile() : iop        = %4d \n", iop);
       printf("                    : frame->no_dof = %4d \n", frame->no_dof);
       printf("                    : frame->no_nodes = %4d \n", frame->no_nodes);
#endif

   prt = PRINT_PROFILE; /* PRINT_PROFILE = ON  has been set in */
                        /* fe_setflags.c prt = PRINT_PROFILE;  */
                        /* option allows profile to be printed */
                        /* all the time                        */

   switch(iop){
       case EQNOS:
            nmin = 0; nad = 0;

            for(i = 1; i <= frame->no_nodes; i++) {
                for(j = 1; j <= frame->no_dof; j++) {
	            if(frame->node[i-1].bound_id[j-1] < nmin)
	               nmin = frame->node[i-1].bound_id[j-1];
             }
            }
         
            neq = -nmin;
            for(n = 1 ; n <= frame->no_nodes; n++) {
                for(i = 1; i <= frame->no_dof; i++)  {
                    j =  frame->node[n-1].bound_id[i-1];
                    if(j < 0) {    
                       frame->node[n-1].bound_id[i-1] 
                       = -frame->node[n-1].bound_id[i-1]; 
                    }             
                    else if(j == 0) {
                       neq = neq + 1;
                       frame->node[n-1].bound_id[i-1] = neq; 

                       /* unconstrained nodes have positive bound_id */

                       CheckJdiagSpace(frame, neq);
                       frame->jdiag[neq-1]      = 0; 
                    }
                    else {
                       nad = nad - 1;
                       frame->node[n-1].bound_id[i-1]  = nad; 
                        /* change restrained node ID from 1 to -1 */
                    }
                }
            }
            frame->no_eq = neq;

            /* compute diagonal pointers for profile */

            nad = 1;  
            frame->jdiag[0] = 1;
            if(frame->no_eq == 1)
               return(frame);

            CheckJdiagSpace(frame, frame->no_eq);
            for(n = 2; n <= frame->no_eq; n++)
                frame->jdiag[n-1] = frame->jdiag[n-1] + frame->jdiag[n-2] +1;

            nad = frame->jdiag[frame->no_eq-1];

            /**********    frame->nad = nad;           *******/
            break;
        default:
            break;
   } 

   /* 
    *  ===============
    *  Print profile()
    *  ===============
    */ 

   if(prt != ON)
      return (frame);

   printf("\n");
   printf("===================================\n");
   printf("           PROFILE                 \n");  
   printf("===================================\n");
   printf("\n");

   printf("node");
   for(i=1; i <= frame->no_dof; i++) {
       switch(i) {
           case 1:
                printf("    dof1");
                break;
           case 2:
                printf("    dof2");
                break;
           case 3:
                printf("    dof3");
                break;
           case 4:
                printf("    dof4");
                break;
           case 5:
                printf("    dof5");
                break;
           case 6:
                printf("    dof6");
                break;
           default:
                break;
         }
   }

   printf("\n");
   printf("----------------------------------------------------\n");

   for(i = 1 ; i <= frame->no_nodes; i++) {
       printf("\n%4d ",i);
       for(j = 1; j <= frame->no_dof; j++)
           printf("\t%4d", frame->node[i-1].bound_id[j-1]);
   }
   printf("\n");
   printf("----------------------------------------------------\n\n");

   return(frame);
}


/* ===========================================================*/
/*  rplink                                                    */
/*  Input: EFRAME *frp                                        */
/*  objective: This subroutine is used to link rigid body and */
/*  flexible columns connected by nodes                       */ 
/* ===========================================================*/

#ifdef __STDC__
EFRAME *rplink(EFRAME *frp, MATRIX *dp, int *idl)
#else
EFRAME *rplink(frp,dp,idl)
EFRAME *frp;
MATRIX *dp;
int    *idl;
#endif
{
MATRIX *sp;
int j,ii,i,nmax,ntemp,nelim, ip, dofl;
int prt = PRINT_PLINK;
int n1;

#ifdef DEBUG
       printf("Enter rplink() : idl[]\n");
       printf("               : ");
       for(i = 1; i<= frp->no_dof;i++)
           printf("%4d ", idl[i-1]);
       printf("\n");
#endif

    /* Find max -ve index in Nodal bound_id[] array */

    nmax = 0;
    for(i = 1; i <= frp->no_nodes; i++) {
        for(j = 1; j <= frp->no_dof; j++)
            if(frp->node[i-1].bound_id[j] < nmax)
               nmax = frp->node[i-1].bound_id[j];
    }

    return(frp);
}


/* =================================================== */
/*  plink                                              */
/*  Input: EFRAME *frp                                 */
/* =================================================== */

#ifdef __STDC__
EFRAME *plink(EFRAME *frp, MATRIX *dofp, int *idl)
#else
EFRAME *plink(frp, dofp, idl)
EFRAME *frp;
MATRIX *dofp;
int    *idl;
#endif
{
MATRIX *sp;
int j,ii,i,nmax,ntemp,nelim, ip, dofl;
int prt = PRINT_PLINK;
int n1;

#ifdef DEBUG
       printf("Enter plink() : idl[]\n");
       printf("              : ");
       for(i = 1; i <= frp->no_dof; i++)
           printf("%4d ", idl[i-1]);
       printf("\n");
#endif

    /* Find max -ve index in Nodal id[fcp->no_dof] array */

       nmax = 0;
       for(i = 1; i <= frp->no_nodes; i++) {
           for(j = 1; j <= frp->no_dof; j++)
	       if(frp->node[i-1].bound_id[j] < nmax)
	          nmax = frp->node[i-1].bound_id[j];
       }

    return(frp);
}


/* ===================================================== */
/*  pload  ;form load vector in compact form             */
/*  Input: EFRAME *frp,                                  */
/*         load vector f(1,TDOF)                         */
/*         multiplier factor p                           */
/*         nneq = no_nodes * ndf = TDOF                  */
/* ===================================================== */

#ifdef __STDC__
QUANTITY *pload(EFRAME *frp, QUANTITY *f, QUANTITY *b, int nneq, QUANTITY p)
#else
QUANTITY *pload(frp,f,b,nneq,p)
EFRAME   *frp;
QUANTITY *f, *b, p;
int      nneq;
#endif
{
int n, i,j,jj;

    for(n = 1; n <= nneq; n++)
        b[n].value = 0.0;
        for(n =1; n <=frp->no_nodes; n++){
        for(i =1; i <=frp->no_dof; i++){
            j = frp->node[n-1].bound_id[i-1];
            if(j >0){
               jj = (n-1) * frp->no_dof+ i;
               b[j].value = f[jj].value * p.value + b[j].value;
            }
        }
        }

    return(b);

}


/* 
 *  ----------------------------------------------------------------- 
 *  Function : ARRAY  *Assign_p_Array()                           
 *                                                                 
 *  Transfers info from EFRAME *frp to working array ARRAY *p
 *  for a particular element elmt_no.                          
 *  Depending on task, elmt arrays in ARRAY *p              
 *  are checked to suit for element arrays &                    
 *  
 *  Input:    EFRAME *frame                                     
 *            int    elmt_no                                 
 *            ARRAY  *p                                        
 *            int        task                                     
 *  
 *  Output:   ARRAY  *p                                         
 *  ----------------------------------------------------------------- 
 */ 

#ifdef __STDC__
ARRAY *Assign_p_Array(EFRAME *frp, int elmt_no, ARRAY *p, int task)
#else
ARRAY *Assign_p_Array(frp, elmt_no, p, task)
EFRAME *frp;
ARRAY  *p;
int    elmt_no, task;
#endif
{
DIMENSIONS  *dp_length_temp, *dp_force_temp;
NODE                                    *np;
ELEMENT                                 *el;
ELEMENT_ATTR                           *eap;
MATRIX                                   *m;
int      ma,nen,j,i,node_no, jj,num_node, k;
int                         n, elmt_attr_no;
int                            dof_per_elmt;
int                                  length;
int                            UNITS_SWITCH;

SYMBOL                                 *slp;
int                        iInPlaneIntegPts;
int                      iThicknessIntegPts;
int                           iNO_INTEG_pts;

#ifdef DEBUG
       printf("Enter Assign_p_Array() : elmt_no    = %4d : task = %4d\n", elmt_no, task);
       printf("Enter Assign_p_Array() : p->elmt_no = %4d\n", p->elmt_no );
#endif

   UNITS_SWITCH = CheckUnits();

   slp = lookup("InPlaneIntegPts");   /* number of integration pts in plane/surface      */
   if(slp == NULL ){
      iInPlaneIntegPts = UNIT_IN_PLANE_INTEG_PTS;        /* 2x2 as default */
   }
   else {
      iInPlaneIntegPts = (int) slp->u.q->value;
   }

   slp = lookup("ThicknessIntegPts"); /* number of integration pts in thickness direction*/

   if(slp == NULL){
      iThicknessIntegPts = UNIT_INTEG_PTS; /* 2 as default */
   }
   else {
      iThicknessIntegPts = (int) slp->u.q->value;
   }

   iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;

   el              = &frp->element[elmt_no-1];  /* element ptr               */
   elmt_attr_no    =  el->elmt_attr_no;
   eap             = &frp->eattr[elmt_attr_no -1]; 
   p->elmt_state   =  el->esp->state;

   for (i = 0; i < NO_ELEMENTS_IN_LIBRARY; i++) { 
      if(!strcmp(elmt_library[i].name, eap->elmt_type)) { 
           n = i;
           break;
      }
   }

   if(task != PROPTY) {
      if(p->elmt_no != elmt_no) { /* check if already values assigned */
         p->elmt_no        =  elmt_no;
         p->elmt_attr_no   =  elmt_attr_no;
         p->dof_per_node   =  elmt_library[n].no_dof;
         p->no_dimen       =  elmt_library[n].no_dimen;
         p->elmt_type      =  elmt_library[n].name;
         p->nodes_per_elmt = (int) MIN(frp->no_nodes_per_elmt,
                                       elmt_library[n].no_node_per_elmt);
         p->size_of_stiff  =  p->nodes_per_elmt* p->dof_per_node;

         free((char *) p->material_name);
         p->material_name = SaveString( frp->eattr[elmt_attr_no-1].material );

         /* p->work_section[] = working arrary of section properties for element */
         /* p->work_material[] = working arrary of material properties for element */

         for(i = 1; i <= UNIT_MATERIAL_ATTR; i++) {  
             if( eap->work_material[i-1].value != 0 ) {
                 p->work_material[i-1].value = eap->work_material[i-1].value;
                 if( UNITS_SWITCH == ON )
                    UnitsCopy( p->work_material[i-1].dimen, eap->work_material[i-1].dimen );
             }
         }
         for(i = 1; i <= UNIT_SECTION_ATTR; i++) {  
             if( eap->work_section[i-1].value != 0 ) {
                 p->work_section[i-1].value = eap->work_section[i-1].value;
                 if( UNITS_SWITCH == ON )
                    UnitsCopy( p->work_section[i-1].dimen, eap->work_section[i-1].dimen );
             }
         }
         if(el->LC_ptr->name != (char *)NULL){
            free((char *) p->LC_ptr->name);
            p->LC_ptr->name  = SaveString(el->LC_ptr->name);
            p->LC_ptr->n     = el->LC_ptr->n;
            p->LC_ptr->alpha = el->LC_ptr->alpha;
            p->LC_ptr->beta  = el->LC_ptr->beta;
         }
         else{
            free((char *) p->LC_ptr->name);
            p->LC_ptr->name  = (char *)NULL;
            p->LC_ptr->n     = 0.0;
            p->LC_ptr->alpha = 0.0;
            p->LC_ptr->beta  = 0.0;
         }
         p->LC_ptr->ialph   = el->LC_ptr->ialph;
         p->LC_ptr->pen     = el->LC_ptr->pen;
         p->LC_ptr->load[0] = el->LC_ptr->load[0];
         p->LC_ptr->load[1] = el->LC_ptr->load[1];
         p->LC_ptr->load[2] = el->LC_ptr->load[2];
         p->LC_ptr->load[3] = el->LC_ptr->load[3];
         p->LC_ptr->load[4] = el->LC_ptr->load[4];
         p->LC_ptr->load[5] = el->LC_ptr->load[5];

         for(j = 1; j <= iNO_INTEG_pts; j++) {
             p->LC_ptr->R[j-1] = el->LC_ptr->R[j-1];
             p->LC_ptr->H[j-1] = el->LC_ptr->H[j-1];
             for(i = 1; i <= 6; i++)
                 p->LC_ptr->back_stress[i-1][j-1] 
                 = el->LC_ptr->back_stress[i-1][j-1];
         }

	 /* Assign work fiber attribution in frame to p array */
	 if( !(strcmp(p->elmt_type, "FIBER_2D"))  || !(strcmp(p->elmt_type, "FIBER_3D"))
	  || !(strcmp(p->elmt_type, "FIBER_2DS")) || !(strcmp(p->elmt_type, "FIBER_3DS")) ) {
            p->integ_ptr->integ_pts = frp->no_integ_pt;
	    /* Use the data stored in frame DIRECTLY , use carefully, do not free them */
	    p->fiber_ptr = eap->work_fiber;

	    p->Q_saved  = el->rp->Q_saved;
	    p->q_saved  = el->rp->q_saved;
	    p->sr_saved = el->rp->sr_saved;
	    p->er_saved = el->rp->er_saved;
	    p->s0_saved = el->rp->s0_saved;
	    p->e0_saved = el->rp->e0_saved;
	    p->sx_saved = el->rp->sx_saved;
	    p->ex_saved = el->rp->ex_saved;

	    p->yielding_saved  = el->esp->yielding_saved;
	    p->pre_range_saved = el->esp->pre_range_saved;
	    p->pre_load_saved  = el->esp->pre_load_saved;
	 }
      }
   }

/*
 *  ==================================================================
 *  Check Allocation of elmt arrays & load values :  depends on 'task'
 *  ==================================================================
 */

    if(p->size_of_stiff > max_size_of_stiff) {

       /* Stiffness Matrix */

       MatrixFreeIndirect(p->stiff);
       p->stiff = MatrixAllocIndirect("stiff/mass",DOUBLE_ARRAY,
                  p->size_of_stiff, p->size_of_stiff);

       if(UNITS_SWITCH == ON ) {

          /* Applied Nodal Loads */

          for( i=1 ; i <= max_size_of_stiff ; i++ ) {
               free((char *) p->nodal_loads[i-1].dimen->units_name);
	       free((char *) p->nodal_loads[i-1].dimen);
          }
       }

       free((char *) p->nodal_loads);
       p->nodal_loads  = (QUANTITY *) MyCalloc(p->size_of_stiff, sizeof(QUANTITY));

       if( UNITS_SWITCH == ON ) {
           for (i = 1; i <= p->size_of_stiff; i++)
                p->nodal_loads[i-1].dimen =(DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS)); 
       }

       MatrixFreeIndirect( p->equiv_nodal_load );
       p->equiv_nodal_load = MatrixAllocIndirect("equiv_nodal_loads",
                             DOUBLE_ARRAY,p->size_of_stiff,1); 

       max_size_of_stiff = p->size_of_stiff;
   }

   if( p->nodes_per_elmt > max_nodes_per_elmt ) {
       if(p->no_dimen > max_no_dimen ) {
	           for(i=1 ; i<=max_no_dimen ; i++ ) {
                       if(UNITS_SWITCH == ON ) {
	                  for(j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                      free((char *) p->coord[i-1][j-1].dimen->units_name);
	                      free((char *) p->coord[i-1][j-1].dimen);
	                      free((char *) p->nodal_traction[i-1][j-1].dimen->units_name);
	                      free((char *) p->nodal_traction[i-1][j-1].dimen);
	                  }
                       }
	               free((char *) p->coord[i-1]);
	               free((char *) p->nodal_traction[i-1]);
	           }
	           free((char *) p->coord);
	           free((char *) p->nodal_traction);

                   p->coord = (QUANTITY **) MyCalloc(p->no_dimen, sizeof(QUANTITY *));
                   p->nodal_traction =(QUANTITY **) MyCalloc(p->no_dimen,sizeof(QUANTITY *));
                   for(i = 1; i <= p->no_dimen; i++) {
                       p->coord[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt,
                                                             sizeof(QUANTITY));
                       p->nodal_traction[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt,
                                                                      sizeof(QUANTITY));  
                       if(UNITS_SWITCH == ON ) {
                          for(j = 1; j <= p->nodes_per_elmt; j++) {
                              p->coord[i-1][j-1].dimen 
                              = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                              p->nodal_traction[i-1][j-1].dimen 
                              = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));  
	                   }
                        }
                    }
                    max_no_dimen = p->no_dimen;
	        }
	        else { /* p->no_dimen <= max_no_dimen &&         */
                       /* p->nodes_per_elmt > max_nodes_per_elmt */

	            for( i=1 ; i<=max_no_dimen ; i++ ) {
                       if(UNITS_SWITCH == ON ) {
	                  for(j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                      free((char *) p->coord[i-1][j-1].dimen->units_name);
	                      free((char *) p->coord[i-1][j-1].dimen);
	                      free((char *) p->nodal_traction[i-1][j-1].dimen->units_name);
	                      free((char *) p->nodal_traction[i-1][j-1].dimen);
	                  }
                       }
	               free((char *) p->coord[i-1]);
	               free((char *) p->nodal_traction[i-1]);
	            }

                    for(i = 1; i <= max_no_dimen; i++) {
                        p->coord[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt,
                                                              sizeof(QUANTITY));
                        p->nodal_traction[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt,
                                                              sizeof(QUANTITY));  
                        if( UNITS_SWITCH == ON ) {
                           for (j = 1; j <= p->nodes_per_elmt; j++) {
                               p->coord[i-1][j-1].dimen  
                               = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                               p->nodal_traction[i-1][j-1].dimen
                               = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));  
	                   }
                        }
                    }
	        }

	        if( p->dof_per_node > max_no_dof ) {

                    MatrixFreeIndirect( p->displ );
                    p->displ = MatrixAllocIndirect("displacement", DOUBLE_ARRAY,
                               p->dof_per_node, p->nodes_per_elmt);

	            for(i=1 ; i<=max_no_dof ; i++ ) {
                        if(UNITS_SWITCH == ON ) {
	                   for(j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                       free((char *) p->nodal_body_force[i-1][j-1].dimen->units_name);
	                       free((char *) p->nodal_body_force[i-1][j-1].dimen);
	                   }
                        }
	                free((char *) p->nodal_body_force[i-1]);
	            }
	            free((char *) p->nodal_body_force);
                    p->nodal_body_force = (QUANTITY **) MyCalloc(p->dof_per_node, sizeof(QUANTITY *));  

                    for(i = 1; i <= p->dof_per_node; i++) {
                        p->nodal_body_force[i-1] 
                        = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  
                        if( UNITS_SWITCH == ON ) {
                           for(j = 1; j <= p->nodes_per_elmt; j++)
                               p->nodal_body_force[i-1][j-1].dimen 
                               = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                        }
                    }

	            MatrixFreeIndirectDouble( p->nodal_init_strain, max_no_dof );
                    p->nodal_init_strain = 
                    (double **) MatrixAllocIndirectDouble(p->dof_per_node, p->nodes_per_elmt);

	            for( i=1 ; i<=max_no_dof ; i++ ) {
                        if(UNITS_SWITCH == ON ) {
	                   for(j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                       free((char *)p->nodal_init_stress[i-1][j-1].dimen->units_name);
	                       free((char *)p->nodal_init_stress[i-1][j-1].dimen);
	                   }
                        }
	                free((char *) p->nodal_init_stress[i-1]);
	            }
	            free((char *) p->nodal_init_stress);
                    p->nodal_init_stress 
                    = (QUANTITY **) MyCalloc(p->dof_per_node, sizeof(QUANTITY *));  
                    for(i = 1; i <= p->dof_per_node; i++) {
                        p->nodal_init_stress[i-1]
                        = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  
                        if( UNITS_SWITCH == ON ) {
                           for(j = 1; j <= p->nodes_per_elmt; j++) 
                               p->nodal_init_stress[i-1][j-1].dimen 
                               = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                        }
                    }
                    if( UNITS_SWITCH == ON ) {
	               for( j=1 ; j<= max_no_dof ; j++ ) {
	                   free((char *) p->nodal_temp[j-1].dimen->units_name);
	                   free((char *) p->nodal_temp[j-1].dimen);
	               }
                    }
	            free((char *) p->nodal_temp);
                    p->nodal_temp = (QUANTITY *) MyCalloc(p->dof_per_node, sizeof(QUANTITY));  
                    if( UNITS_SWITCH == ON ) {
                       for (i = 1; i <= p->dof_per_node; i++)
                           p->nodal_temp[i-1].dimen 
                           = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                    }
                    max_no_dof         = p->dof_per_node;
	        }
	        else { /* p->dof_per_node <= max_no_dof &&       */
                       /* p->nodes_per_elmt > max_nodes_per_elmt */

                    MatrixFreeIndirect( p->displ );
                    p->displ = MatrixAllocIndirect("displacement", 
                               DOUBLE_ARRAY, max_no_dof, p->nodes_per_elmt);

	            for( i=1 ; i <= max_no_dof ; i++ ) {
                        if( UNITS_SWITCH == ON ) {
	                   for( j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                          free((char *)p->nodal_body_force[i-1][j-1].dimen->units_name);
	                          free((char *)p->nodal_body_force[i-1][j-1].dimen);

	                   }
                        }
	                free((char *) p->nodal_body_force[i-1]);
	            }
                    for (i = 1; i <= max_no_dof; i++) {
                        p->nodal_body_force[i-1] 
                        = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  
                        if( UNITS_SWITCH == ON ) {
                           for(j = 1; j <= p->nodes_per_elmt; j++)
                               p->nodal_body_force[i-1][j-1].dimen 
                               = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                        }
                    }

	            for( i=1; i <= max_no_dof; i++ ) {
	                free((char *) p->nodal_init_strain[i-1] );
	                p->nodal_init_strain[i-1] 
                        = (double *)MyCalloc(p->nodes_per_elmt, sizeof(double));
	            }

	            for( i=1 ; i<=max_no_dof ; i++ ) {
                        if( UNITS_SWITCH == ON ) {
	                   for( j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                          free((char *)p->nodal_init_stress[i-1][j-1].dimen->units_name);
     	                          free((char *)p->nodal_init_stress[i-1][j-1].dimen);
	                   }
                        }
	                free((char *) p->nodal_init_stress[i-1]);
	            }
                    for(i = 1; i <= max_no_dof; i++) {
                        p->nodal_init_stress[i-1]
                        = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  
                        if( UNITS_SWITCH == ON ) {
                           for(j = 1; j <= p->nodes_per_elmt; j++) 
                               p->nodal_init_stress[i-1][j-1].dimen
                               = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                        }
                    }
	        }
                free((char *) p->node_connect);
                p->node_connect = iVectorAlloc(p->nodes_per_elmt);

                max_nodes_per_elmt = p->nodes_per_elmt;
    } else { 

      /* p->nodes_per_elmt <= max_nodes_per_elmt */

      if( p->no_dimen > max_no_dimen ) {
          for( i=1 ; i<=max_no_dimen ; i++ ) {
               if( UNITS_SWITCH == ON ) {
	           for( j=1 ; j<= max_nodes_per_elmt ; j++ ) {
	                free((char *) p->coord[i-1][j-1].dimen->units_name);
	                free((char *) p->coord[i-1][j-1].dimen);
	                free((char *) p->nodal_traction[i-1][j-1].dimen->units_name);
	                free((char *) p->nodal_traction[i-1][j-1].dimen);
	           }
               }
	       free((char *) p->coord[i-1]);
	       free((char *) p->nodal_traction[i-1]);
	  }

	  free((char *) p->coord);
	  free((char *) p->nodal_traction);

          p->coord = (QUANTITY **) MyCalloc(p->no_dimen, sizeof(QUANTITY *));
          p->nodal_traction = (QUANTITY **)MyCalloc(p->no_dimen,sizeof(QUANTITY *));

          for( i = 1; i <= p->no_dimen; i++) {
               p->coord[i-1]  = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));
               p->nodal_traction[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  

               if( UNITS_SWITCH == ON ) {
                   for(j = 1; j <= p->nodes_per_elmt; j++) {
                       p->coord[i-1][j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                       p->nodal_traction[i-1][j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));  
	           }
               }
          }
          max_no_dimen = p->no_dimen;
       }

       if( p->dof_per_node > max_no_dof ) {
           MatrixFreeIndirect( p->displ );
           p->displ = MatrixAllocIndirect("displacement", DOUBLE_ARRAY,
                                          p->dof_per_node, p->nodes_per_elmt );

           for( i=1 ; i<=max_no_dof ; i++ ) {
                if( UNITS_SWITCH == ON ) {
                    for( j=1 ; j<= max_nodes_per_elmt ; j++ ) {
                         free((char *)p->nodal_body_force[i-1][j-1].dimen->units_name);
                         free((char *)p->nodal_body_force[i-1][j-1].dimen);
                    }
                }
                free((char *) p->nodal_body_force[i-1]);
           }

           free((char *) p->nodal_body_force);
           p->nodal_body_force = (QUANTITY **) MyCalloc(p->dof_per_node, sizeof(QUANTITY *));  

           for(i = 1; i <= p->dof_per_node; i++) {
               p->nodal_body_force[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  
               if(UNITS_SWITCH == ON ) {
                  for(j = 1; j <= p->nodes_per_elmt; j++)
                      p->nodal_body_force[i-1][j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
               }
           }

	   MatrixFreeIndirectDouble( p->nodal_init_strain, max_no_dof );
           p->nodal_init_strain = (double **)MatrixAllocIndirectDouble(p->dof_per_node, p->nodes_per_elmt); 

           for( i=1 ; i<=max_no_dof ; i++ ) {
                if(UNITS_SWITCH == ON ) {
	           for(j=1 ; j<= max_nodes_per_elmt ; j++ ) {
                       if(p->nodal_init_stress[i-1][j-1].dimen != NULL &&
                          p->nodal_init_stress[i-1][j-1].dimen->units_name != NULL)
	                  free((char *)p->nodal_init_stress[i-1][j-1].dimen->units_name);
	                  free((char *)p->nodal_init_stress[i-1][j-1].dimen);
	               }
                }
	        free((char *) p->nodal_init_stress[i-1]);
	   }

	   free((char *) p->nodal_init_stress);
           p->nodal_init_stress = (QUANTITY **) MyCalloc(p->dof_per_node, sizeof(QUANTITY *));  

           for(i = 1; i <= p->dof_per_node; i++) {
               p->nodal_init_stress[i-1] = (QUANTITY *) MyCalloc(p->nodes_per_elmt, sizeof(QUANTITY));  
               if( UNITS_SWITCH == ON ) {
                   for(j = 1; j <= p->nodes_per_elmt; j++) 
                       p->nodal_init_stress[i-1][j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
                   }
           }

           if( UNITS_SWITCH == ON ) {
	       for( j=1 ; j<= max_no_dof ; j++ ) {
	            free((char *) p->nodal_temp[j-1].dimen->units_name);
	            free((char *) p->nodal_temp[j-1].dimen);
	       }
           }

	   free((char *) p->nodal_temp);
           p->nodal_temp = (QUANTITY *) MyCalloc(p->dof_per_node, sizeof(QUANTITY));  

           if(UNITS_SWITCH == ON ) {
              for (i = 1; i <= p->dof_per_node; i++)
                   p->nodal_temp[i-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           }
           max_no_dof = p->dof_per_node;
        }
   }

   switch (task) {
      case PROPTY:
           for(i = 1; i <= UNIT_SECTION_ATTR; i++) {
               p->work_section[i-1].value 
               = frp->eattr[elmt_attr_no-1].work_section[i-1].value;
           }

           if(UNITS_SWITCH == ON ){
              for(i = 1; i <= UNIT_SECTION_ATTR; i++) 
                  UnitsCopy( p->work_section[i-1].dimen,
                  frp->eattr[elmt_attr_no-1].work_section[i-1].dimen );
           }

           for(i = 1; i <= UNIT_MATERIAL_ATTR; i++) {
               p->work_material[i-1].value 
               = frp->eattr[elmt_attr_no-1].work_material[i-1].value;
              if(UNITS_SWITCH == ON )
                 UnitsCopy( p->work_material[i-1].dimen,
                 frp->eattr[elmt_attr_no-1].work_material[i-1].dimen );
            }

            if(el->LC_ptr->name != (char *)NULL) {
               free((char *) p->LC_ptr->name);
               p->LC_ptr->name  = SaveString(el->LC_ptr->name);
               p->LC_ptr->n     = el->LC_ptr->n;
               p->LC_ptr->alpha = el->LC_ptr->alpha;
               p->LC_ptr->beta  = el->LC_ptr->beta;
            }else {
               free((char *) p->LC_ptr->name);
               p->LC_ptr->name  = (char *)NULL;
               p->LC_ptr->n     = 0.0;
               p->LC_ptr->alpha = 0.0;
               p->LC_ptr->beta  = 0.0;
            }
            p->LC_ptr->ialph   = el->LC_ptr->ialph;
            p->LC_ptr->pen     = el->LC_ptr->pen;
            p->LC_ptr->load[0] = el->LC_ptr->load[0];
            p->LC_ptr->load[1] = el->LC_ptr->load[1];
            p->LC_ptr->load[2] = el->LC_ptr->load[2];
            p->LC_ptr->load[3] = el->LC_ptr->load[3];
            p->LC_ptr->load[4] = el->LC_ptr->load[4];
            p->LC_ptr->load[5] = el->LC_ptr->load[5];

            for(j = 1; j <= iNO_INTEG_pts; j++) {
                p->LC_ptr->R[j-1] = el->LC_ptr->R[j-1];
                p->LC_ptr->H[j-1] = el->LC_ptr->H[j-1];
                for(i = 1; i <= 6; i++)
                    p->LC_ptr->back_stress[i-1][j-1] = el->LC_ptr->back_stress[i-1][j-1];
            }

	    /* Assign work fiber attribution in frame to p array */

	    if( !(strcmp(p->elmt_type, "FIBER_2D"))  ||
                !(strcmp(p->elmt_type, "FIBER_3D"))  ||
                !(strcmp(p->elmt_type, "FIBER_2DS")) ||
                !(strcmp(p->elmt_type, "FIBER_3DS")) ) {

                p->integ_ptr->integ_pts = frp->no_integ_pt;

                /* Use the data stored in frame DIRECTLY , use carefully, do not free them */

                p->fiber_ptr = eap->work_fiber;

	        p->Q_saved  = el->rp->Q_saved;
	        p->q_saved  = el->rp->q_saved;
	        p->sr_saved = el->rp->sr_saved;
	        p->er_saved = el->rp->er_saved;
	        p->s0_saved = el->rp->s0_saved;
	        p->e0_saved = el->rp->e0_saved;
	        p->sx_saved = el->rp->sx_saved;
	        p->ex_saved = el->rp->ex_saved;

	        p->yielding_saved  = el->esp->yielding_saved;
	        p->pre_range_saved = el->esp->pre_range_saved;
	        p->pre_load_saved  = el->esp->pre_load_saved;
	    }
            break;

       case STIFF:
       case MASS_MATRIX:

            if(task == MASS_MATRIX)
               for(i = 1;i <= p->size_of_stiff; i++) {
                   p->nodal_loads[i-1].value = 0.0;
                   if( UNITS_SWITCH == ON )
                      ZeroUnits( p->nodal_loads[i-1].dimen );
               }

            /* p->node_connect[] = list of nodal connectivities */

            for(j = 1; j <= p->nodes_per_elmt;j++)
                p->node_connect[j-1] = el->node_connect[j-1];

            for(j = 1; j <= iNO_INTEG_pts; j++) {
                p->effect_pl_strain[j-1]   = el->rp->effect_pl_strain[j-1];
                p->eff_pl_strain_incr[j-1] = 0.0;
                for(i = 1; i <= 9; i++){
                    p->stress->uMatrix.daa[i-1][j-1] 
                    = el->rp->stress->uMatrix.daa[i-1][j-1];
                    p->strain_pl->uMatrix.daa[i-1][j-1] 
                    = el->rp->strain_pl->uMatrix.daa[i-1][j-1];
                    p->strain_pl_incr->uMatrix.daa[i-1][j-1] 
                    = 0.0;
                }
            }

            for(j = 1; j <= frp->no_nodes_per_elmt; j++) {
                node_no = el->node_connect[j-1];
                if(node_no != 0) {
                   for(i=1;i <= p->no_dimen; i++) {
                       p->coord[i-1][j-1].value 
                       = frp->node[node_no -1].coord[i-1].value;
                       if(UNITS_SWITCH == ON )
                          UnitsCopy( p->coord[i-1][j-1].dimen,
                          frp->node[node_no -1].coord[i-1].dimen );
                   }
                }
            }

            /*  initialize and Transfer Element Displacements */

            if( UNITS_SWITCH == ON ) {
               for(i = 1; i <= p->dof_per_node; i++) {
                   ZeroUnits( &(p->displ->spRowUnits[i-1]) );
               }
               for(i = 1; i <= p->nodes_per_elmt; i++) {
                   ZeroUnits( &(p->displ->spColUnits[i-1]) );
               }
            }

            for(i=1; i <= p->dof_per_node; i++) {
                for(j=1; j <= p->nodes_per_elmt; j++) {
                    p->displ->uMatrix.daa[i-1][j-1]      = el->rp->displ->uMatrix.daa[i-1][j-1];
                    p->displ_incr->uMatrix.daa[i-1][j-1] = 0.0;
                }
            }
            break;

       case LOAD_MATRIX:
       case EQUIV_NODAL_LOAD:

            for(i = 1; i <= p->size_of_stiff; i++) 
	        p->nodal_loads[i-1].value   = 0.0;

            if(UNITS_SWITCH == ON ) {
               for(i = 1; i <= p->size_of_stiff; i++) {
                   ZeroUnits( &(p->equiv_nodal_load->spRowUnits[i-1]) );
                   ZeroUnits( p->nodal_loads[i-1].dimen );
               }
               ZeroUnits( &(p->equiv_nodal_load->spColUnits[0]) );
            }

            dof_per_elmt = p->nodes_per_elmt*p->dof_per_node;
            for(j = 1; j <= dof_per_elmt; j++) {
                p->equiv_nodal_load->uMatrix.daa[j-1][0] = 0.0;
            }
    
            for(j = 1; j <= p->nodes_per_elmt; j++) {
               p->nodal_temp[j-1].value = frp->eforces->elib_ptr[j-1].temp_change.value; 
               if(UNITS_SWITCH == ON )
                  UnitsCopy( p->nodal_temp[j-1].dimen,
                  frp->eforces->elib_ptr[j-1].temp_change.dimen );
               for( i = 1; i <= p->no_dimen; i++) {
                   /*== subjected to change later ==*/
                    p->nodal_traction[i-1][j-1].value = 0.0;
                    if(UNITS_SWITCH == ON )
                       ZeroUnits( p->nodal_traction[i-1][j-1].dimen );
                    if(frp->eforces->elib_ptr[j-1].body_force != NULL) {
                       p->nodal_body_force[i-1][j-1].value 
                       = frp->eforces->elib_ptr[j-1].body_force[i-1].value; 
                       if( UNITS_SWITCH == ON )
                          UnitsCopy( p->nodal_body_force[i-1][j-1].dimen,
                          frp->eforces->elib_ptr[j-1].body_force[i-1].dimen );
                    }
                    else {
                       p->nodal_body_force[i-1][j-1].value  = 0.0;
                       if(UNITS_SWITCH == ON )
                          ZeroUnits( p->nodal_body_force[i-1][j-1].dimen );
                    }
               }
            }

            for(j = 1; j <= p->nodes_per_elmt; j++) {
               for( i = 1; i <= p->dof_per_node; i++) { 
                 if(frp->eforces->elib_ptr[j-1].init_strain != NULL)
                    p->nodal_init_strain[i-1][j-1] 
                    = frp->eforces->elib_ptr[j-1].init_strain[i-1];
                 else
                    p->nodal_init_strain[i-1][j-1] = 0.0;
                 if(frp->eforces->elib_ptr[j-1].init_stress != NULL) {
                    p->nodal_init_stress[i-1][j-1].value 
                    = frp->eforces->elib_ptr[j-1].init_stress[i-1].value; 
                    if( UNITS_SWITCH == ON )
                       UnitsCopy( p->nodal_init_stress[i-1][j-1].dimen,
                       frp->eforces->elib_ptr[j-1].init_stress[i-1].dimen );
                 }
                 else {
                   p->nodal_init_stress[i-1][j-1].value  = 0.0;
                   if(UNITS_SWITCH == ON )
                      ZeroUnits( p->nodal_init_stress[i-1][j-1].dimen );
                 }
               }
            }

            for(j=1;j <= p->nodes_per_elmt;j++) {
               node_no = el->node_connect[j-1];
               if(node_no != 0){
                  for(i=1;i<=p->no_dimen;i++) {
                      p->coord[i-1][j-1].value 
                      = frp->node[node_no -1].coord[i-1].value;    
                      if( UNITS_SWITCH == ON )
                         UnitsCopy( p->coord[i-1][j-1].dimen,
                         frp->node[node_no -1].coord[i-1].dimen );
                  }
               }
            }

            /* Assign the nodal load and stress in previous step to p->nodal Load */

            for(k = 1; k <= p->size_of_stiff; k++)
                    p->nodal_loads[k-1].value = 0.0;

            for(j = 1; j <= iNO_INTEG_pts; j++) {
                p->effect_pl_strain[j-1]   = el->rp->effect_pl_strain[j-1];
                p->eff_pl_strain_incr[j-1] = 0.0;
                for(i = 1; i <= 9; i++){
                    p->stress->uMatrix.daa[i-1][j-1]
                    = el->rp->stress->uMatrix.daa[i-1][j-1];
                    p->strain_pl->uMatrix.daa[i-1][j-1] 
                    = el->rp->strain_pl->uMatrix.daa[i-1][j-1];
                    p->strain_pl_incr->uMatrix.daa[i-1][j-1] 
                    = 0.0;
                }
            }

            for(i=1; i <= p->dof_per_node; i++) {
            for(j=1; j <= p->nodes_per_elmt; j++) {
                p->displ->uMatrix.daa[i-1][j-1]      = el->rp->displ->uMatrix.daa[i-1][j-1];
                p->displ_incr->uMatrix.daa[i-1][j-1] = 0.0;
            }
            }
            break;

       case PRESSLD:
       case STRESS_LOAD:

            for(i=1;i<=p->size_of_stiff ;i++){
               p->nodal_loads[i-1].value  = 0.0;
               if( UNITS_SWITCH == ON )
                  ZeroUnits( p->nodal_loads[i-1].dimen );
            }

            for(j = 1; j <= iNO_INTEG_pts; j++){
                p->effect_pl_strain[j-1]   = el->rp->effect_pl_strain[j-1];
                p->eff_pl_strain_incr[j-1] = 0.0;
                for(i = 1; i <= 9; i++){
                    p->stress->uMatrix.daa[i-1][j-1] = el->rp->stress->uMatrix.daa[i-1][j-1];
                    p->strain_pl->uMatrix.daa[i-1][j-1]      = el->rp->strain_pl->uMatrix.daa[i-1][j-1];
                    p->strain_pl_incr->uMatrix.daa[i-1][j-1] = 0.0;
                }
            }

            for(i=1; i <= p->dof_per_node; i++) {
            for(j=1; j <= p->nodes_per_elmt; j++) {
                p->displ->uMatrix.daa[i-1][j-1] = el->rp->displ->uMatrix.daa[i-1][j-1];
                p->displ_incr->uMatrix.daa[i-1][j-1] = 0.0;
            }
            }

            for(j=1;j<= p->nodes_per_elmt ;j++) {
               node_no = el->node_connect[j-1];
               if(node_no != 0){
                  for(i=1;i<=p->no_dimen;i++) {
                      p->coord[i-1][j-1].value = frp->node[node_no -1].coord[i-1].value;
                      if(UNITS_SWITCH == ON )
                         UnitsCopy( p->coord[i-1][j-1].dimen,
                         frp->node[node_no -1].coord[i-1].dimen );
                  }
               }
            }

            if(task == STRESS_LOAD) {
               /* Get_elmt_forces_and_strains(p,frp,elmt_no); */ ;
            } else{  /* PRESSLD */
               i = 0;
               for(j=1;j <= frp->no_element_loads; j++) {
                   if(elmt_no == frp->eforces[i].elmt_no ){ 
                      p->elmt_load_ptr =  &frp->eforces[i-1] ;
                      goto PRESS_END;
                   }
                   else p->elmt_load_ptr = (ELEMENT_LOADS *) NULL;
                   i = i+1;
               }
            }

            PRESS_END:
            break;
       case STRESS:

            /* For nonlinear problems, get element forces and strains */

            for(i = 1; i <= p->size_of_stiff; i++) {
                p->nodal_loads[i-1].value = 0.0;
                if( UNITS_SWITCH == ON )
                    ZeroUnits( p->nodal_loads[i-1].dimen );
            }

	    /* Element Displacements */

            if( UNITS_SWITCH == ON ) {
                for(i = 1; i <= p->dof_per_node; i++)
                    ZeroUnits( &(p->displ->spRowUnits[i-1]) );
                for(i = 1; i <= p->nodes_per_elmt; i++)
                    ZeroUnits( &(p->displ->spColUnits[i-1]) );
            }

            /* Zero and Transfer Element Displacements */

            for(i=1; i <= p->dof_per_node; i++) {
            for(j=1; j <= p->nodes_per_elmt; j++) {
                p->displ->uMatrix.daa[i-1][j-1]      = el->rp->displ->uMatrix.daa[i-1][j-1];
                p->displ_incr->uMatrix.daa[i-1][j-1] = 0.0;
            }
            }

            /* Element Coordinates */

            for(j = 1; j <= p->nodes_per_elmt; j++) {
                node_no = el->node_connect[j-1];
                if(node_no != 0){
                   for(i=1; i <= p->no_dimen; i++) {
                       p->coord[i-1][j-1].value 
                       = frp->node[node_no -1].coord[i-1].value;
                       if(UNITS_SWITCH == ON ) {
                          UnitsCopy( p->coord[i-1][j-1].dimen,
                          frp->node[node_no -1].coord[i-1].dimen );
                       }
                   }
                }
            }

            /* [g] : Get Elmt Load Array */

            i = 0;
            for(j = 1; j <= frp->no_element_loads; j++) {
                if(elmt_no == frp->eforces[i-1].elmt_no) { 
                   p->elmt_load_ptr = &frp->eforces[i-1];
                   goto  TRESS_END;
                }
                else 
                   p->elmt_load_ptr = (ELEMENT_LOADS *) NULL;
                i = i+1;
            }

            /* [h] : Element stress & strain  */

            for(j = 1; j <= iNO_INTEG_pts; j++) {

            p->effect_pl_strain[j-1]   = el->rp->effect_pl_strain[j-1];
            p->eff_pl_strain_incr[j-1] = 0.0;

            for(i = 1; i <= 9; i++) {
                p->stress->uMatrix.daa[i-1][j-1]         = el->rp->stress->uMatrix.daa[i-1][j-1];
                p->strain_pl->uMatrix.daa[i-1][j-1]      = el->rp->strain_pl->uMatrix.daa[i-1][j-1];
                p->strain_pl_incr->uMatrix.daa[i-1][j-1] = 0.0;
            }
            }
            break;

            TRESS_END:
            break;
       default:
            break; 
    }
    return(p);

#ifdef DEBUG
       printf("*** Leave Assign_p_Array()\n");
#endif
}


/* 
 *  ===================================================
 *  Alloc_p_Array() : Allocate memory for working array
 *  ===================================================
 */ 

ARRAY *Alloc_p_Array()
{
ARRAY               *pp;
int   i,j, no_integ_pts;
int                temp;
int        UNITS_SWITCH;

SYMBOL             *slp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;

   slp = lookup("InPlaneIntegPts");   /* number of integration pts in plane/surface      */
   if(slp == NULL) {
      iInPlaneIntegPts = UNIT_IN_PLANE_INTEG_PTS;        /* 2x2 as default */
   }
   else{
      iInPlaneIntegPts = (int) slp->u.q->value;
   }

   slp = lookup("ThicknessIntegPts"); /* number of integration pts in thickness direction*/
   if(slp == NULL)
      iThicknessIntegPts = 2;        /* 2 as default */
   else
      iThicknessIntegPts = (int) slp->u.q->value;

   iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;

   slp = lookup("NDimension");         /* No of dimensions : 2 or 3 */
   max_no_dimen = (int) slp->u.q->value;

   slp = lookup("NDofPerNode");        /* No gdof per node */
   max_no_dof = (int) slp->u.q->value;

   slp = lookup("MaxNodesPerElement"); /* Max no nodes per element */

   max_nodes_per_elmt = (int) slp->u.q->value;
   max_size_of_stiff  = max_nodes_per_elmt*max_no_dof;
   no_integ_pts       = iThicknessIntegPts;

   UNITS_SWITCH = CheckUnits();
   pp = (ARRAY *) MyCalloc(1,sizeof(ARRAY));

   pp->integ_ptr                = (INTEG_PTS *) MyCalloc(1,sizeof(INTEG_PTS));
   pp->integ_ptr->surface_pts   = (int) sqrt((double)iInPlaneIntegPts);  
   pp->integ_ptr->thickness_pts =  iThicknessIntegPts;
   pp->integ_ptr->integ_pts     =  no_integ_pts;

   pp->elmt_type      = (char *) NULL;
   pp->material_name  = (char *) NULL;
   pp->displ          = MatrixAllocIndirect("displ",DOUBLE_ARRAY,
                        max_no_dof,max_nodes_per_elmt);
   pp->work_section   = (QUANTITY *) MyCalloc(UNIT_SECTION_ATTR,sizeof(QUANTITY)); 
   pp->work_material  = (QUANTITY *) MyCalloc(UNIT_MATERIAL_ATTR, sizeof(QUANTITY));
   pp->d_array        = iVectorAlloc(max_size_of_stiff); 
   pp->coord          = (QUANTITY **) MyCalloc(max_no_dimen, sizeof(QUANTITY *));

   for (i = 1; i <= max_no_dimen; i++)
       pp->coord[i-1] = (QUANTITY *) MyCalloc(max_nodes_per_elmt, sizeof(QUANTITY));

   pp->node_connect   = iVectorAlloc(max_nodes_per_elmt);  
   pp->nodal_loads    = (QUANTITY *) MyCalloc(max_size_of_stiff, sizeof(QUANTITY));
   pp->nodal_temp     = (QUANTITY *) MyCalloc(max_no_dof, sizeof(QUANTITY)); 
   pp->nodal_body_force  = (QUANTITY **) MyCalloc(max_no_dof, sizeof(QUANTITY *));   

   for( i = 1; i <= max_no_dof; i++)
       pp->nodal_body_force[i-1] = (QUANTITY *) MyCalloc(max_nodes_per_elmt, sizeof(QUANTITY));

   pp->nodal_init_strain = (double **)MatrixAllocIndirectDouble(max_no_dof, max_nodes_per_elmt); 
   pp->nodal_init_stress = (QUANTITY **) MyCalloc(max_no_dof, sizeof(QUANTITY *));   

   for( i = 1; i <= max_no_dof; i++)
       pp->nodal_init_stress[i-1] = (QUANTITY *) MyCalloc(max_nodes_per_elmt, sizeof(QUANTITY));

   pp->nodal_traction    = (QUANTITY **) MyCalloc(max_no_dimen, sizeof(QUANTITY *));   

   for( i = 1; i <= max_no_dimen; i++)
        pp->nodal_traction [i-1] = (QUANTITY *)
                                   MyCalloc(max_nodes_per_elmt, sizeof(QUANTITY));

   pp->stiff            = MatrixAllocIndirect("stiff/mass",
                          DOUBLE_ARRAY, max_size_of_stiff, max_size_of_stiff);
   pp->equiv_nodal_load = MatrixAllocIndirect("equivalent_load",
                          DOUBLE_ARRAY,max_size_of_stiff, 1);
   pp->mater_matrix     = MatrixAllocIndirect("material_matrix",
                          DOUBLE_ARRAY,max_no_dof, max_no_dof);
   pp->displ_incr       = MatrixAllocIndirect("incremt_displ",
                          DOUBLE_ARRAY, max_no_dof, max_nodes_per_elmt);
   pp->strain_rate      = MatrixAllocIndirect("strain_rate", DOUBLE_ARRAY, 9, 1);
   pp->stress_rate      = MatrixAllocIndirect("stress_rate", DOUBLE_ARRAY, 9, 1);

   pp->effect_pl_strain   = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
   pp->eff_pl_strain_incr = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
   pp->direc_cos        = (double *) MyCalloc(2, sizeof(double));
   pp->direc_cos[0]     = 1.0;
   pp->direc_cos[1]     = 0.0;

   pp->LC_ptr           = (MATER_LOAD_CURVE *) MyCalloc(1,sizeof(MATER_LOAD_CURVE));
   pp->LC_ptr->name     = (char *)NULL;
   pp->LC_ptr->R        = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
   pp->LC_ptr->H        = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
   pp->LC_ptr->back_stress = MatrixAllocIndirectDouble(6, iNO_INTEG_pts);

   pp->stress           = MatrixAllocIndirect("stress", DOUBLE_ARRAY, 9, iNO_INTEG_pts);
   pp->strain_pl        = MatrixAllocIndirect("strain", DOUBLE_ARRAY, 9, iNO_INTEG_pts);
   pp->strain_pl_incr   = MatrixAllocIndirect("strain", DOUBLE_ARRAY, 9, iNO_INTEG_pts);

/* element load pointer is not allocated until the load is present */

   if( UNITS_SWITCH == ON ) {

      for(i = 1; i <= UNIT_SECTION_ATTR; i++) {
          pp->work_section[i-1].dimen 
          = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS)); 
      }

      for(i = 1; i <= UNIT_MATERIAL_ATTR; i++) {
          pp->work_material[i-1].dimen
          = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS)); 
      }

      for (i = 1; i <= max_no_dimen; i++) {
          for (j = 1; j <= max_nodes_per_elmt; j++) {
             pp->coord[i-1][j-1].dimen 
             = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
      }

      for (i = 1; i <= max_size_of_stiff; i++) {
         pp->nodal_loads[i-1].dimen 
         = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
      }

      for (i = 1; i <= max_no_dof; i++) {
         pp->nodal_temp[i-1].dimen 
         = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
      }

      for( i = 1; i <= max_no_dof; i++) {
           for (j = 1; j <= max_nodes_per_elmt; j++) {
              pp->nodal_body_force[i-1][j-1].dimen 
              = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           }
      }

      for( i = 1; i <= max_no_dof; i++) {
          for (j = 1; j <= max_nodes_per_elmt; j++) {
             pp->nodal_init_stress[i-1][j-1].dimen 
             = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
      }

      for( i = 1; i <= max_no_dimen; i++) {
          for (j = 1; j <= max_nodes_per_elmt; j++) {
             pp->nodal_traction[i-1][j-1].dimen 
             = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
      }

      pp->eangle.dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
      pp->ealpha.dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
      pp->length.dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
   }
   return(pp);

#ifdef DEBUG
       printf("*** Leave fe_profile.c : Alloc_p_Array()\n\n");
#endif

}


/* ==================================================== */
/* Get_elmt_forces_and_strains(p,frp,elmt_no)           */
/* ==================================================== */

#ifdef __STDC__
Get_elmt_forces_and_strains(ARRAY *p, EFRAME *frp, int elmt_no)
#else
Get_elmt_forces_and_strains(p,frp,elmt_no)
ARRAY  *p;
EFRAME *frp;
int    elmt_no; 
#endif
{
 
int i, j, k, no_dof, no_nodes_per_elmt;
char                   *elmt_type_name;
ELEMENT_ATTR                      *eap;
int                       elmt_attr_no;
int                             length;
int                       UNITS_SWITCH;

    UNITS_SWITCH = CheckUnits();
    p->elmt_state = frp->element[elmt_no -1].esp->state;
/*
    p->eep        = frp->element[elmt_no -1].ep;
    p->ealpha     = frp->element[elmt_no -1].alpha;
*/
    
    elmt_attr_no = frp->element[elmt_no-1].elmt_attr_no;
    eap = &frp->eattr[elmt_attr_no-1];

    /* get forces from EFRAME *    */
    
    for (i = 1; i <= NO_ELEMENTS_IN_LIBRARY; i++) {
        if(!strcmp(elmt_library[i-1].name, frp->eattr[elmt_attr_no-1].elmt_type)) {
        elmt_type_name    = elmt_library[i-1].name;
        no_dof            = elmt_library[i-1].no_dof;
        no_nodes_per_elmt = elmt_library[i-1].no_node_per_elmt;
        }
    }
    
    for(i = 1; i <= no_dof+1; i++){ 
       for(j = 1; j <= no_nodes_per_elmt; j++) {

          k = no_dof*(j-1)+i;
          p->nodal_loads[k-1].value   = frp->element[elmt_no -1].rp->Forces->uMatrix.daa[i-1][j-1]; 
          if( UNITS_SWITCH == ON ) {
              UnitsCopy( p->nodal_loads[k-1].dimen,
              &(frp->element[elmt_no -1].rp->Forces->spRowUnits[i-1]) );
          }
       }
    }
}
