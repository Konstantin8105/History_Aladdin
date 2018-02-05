/*
 *  ============================================================================= 
 *  ALADDIN Version 2.0 :
 *                                                                     
 *  fe_checkmemory.c : Check and Reallocate Memory for Frame Data Structure
 *                                                                     
 *  Copyright (C) 1995-1997 by Mark Austin, Xiaoguang Chen, and Wane-Jang Lin
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
 *  3. Altered versions must be plainly marked as such, and must not
 *     be misrepresented as being the original software.
 *  4. This notice is to remain intact.
 *                                                                    
 *  Written by: Mark Austin, Xiaoguang Chen, and Wane-Jang Lin           May 1997
 *  ============================================================================= 
 */

#include <stdio.h>
#include <ctype.h>

#include "defs.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "vector.h"
#include "fe_functions.h"

/*
#define DEBUG
*/


/* ================ */
/* Check Node Space */
/* ================ */

#ifdef __STDC__
EFRAME *CheckNodeSpace(EFRAME *frp, int node_no)
#else
EFRAME *CheckNodeSpace(frp, node_no)
EFRAME *frp;
int    node_no;
#endif
{
NODE         *np;
int         i, j;  
int UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter CheckNodeSpace()");
#endif

    UNITS_SWITCH = CheckUnits();
    if(node_no > max_no_nodes) {

       frp->no_nodes = node_no;
       frp->node  = (NODE *) realloc(frp->node,(max_no_nodes+node_no)*sizeof(NODE));
       frp->jdiag = (int *) realloc(frp->jdiag,(((max_no_nodes+node_no)*UNIT_NDF*sizeof(int))));

       for(i = max_no_nodes + 1; i <= max_no_nodes + node_no; i++) {

           np = &frp->node[i-1];
           np->coord    = (QUANTITY *) MyCalloc(frp->no_dimen, sizeof(QUANTITY));
           np->bound_id =      (int *) MyCalloc(frp->no_dof, sizeof(int));
           np->disp     = (QUANTITY *) MyCalloc(frp->no_dof, sizeof(QUANTITY));
           np->rb_num   = (int) 0;

           if( UNITS_SWITCH == ON ) {
              for(j = 1; j <= frp->no_dimen; j++) {
                 np->coord[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              }
              for (j = 1; j <= frp->no_dof; j++) {
                 np->disp[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              }
           }

           np->TrT  = (MATRIX *) MyCalloc(1,sizeof(MATRIX));
       }

       max_no_nodes = max_no_nodes + node_no;
    }

#ifdef DEBUG
       printf("*** Leave CheckNodeSpace()\n");
#endif

    return (frp);
}


/* ================== */
/* Check jdiag Space  */
/* ================== */

#ifdef __STDC__
EFRAME *CheckJdiagSpace(EFRAME *frp, int iEqNo)
#else
EFRAME *CheckJdiagSpace(frp, iEqNo)
EFRAME      *frp;
int        iEqNo;
#endif
{
int  i, j;  

#ifdef DEBUG
       printf("*** Enter CheckJdiagSpace()");
       printf("  frp->no_dimen = %4d \n", frp->no_dimen);
       printf(" : iEqNo = %4d : max_no_eqs = %4d\n", iEqNo, max_no_eqs);
#endif
   
    if(iEqNo > max_no_eqs) {
       frp->jdiag = (int *) realloc(frp->jdiag, (max_no_eqs+iEqNo)*sizeof(int));

#ifdef DEBUG
       printf("*** In CheckJdiagSpace() : REALLOCATE SPACE()\n");
#endif

       for(i = max_no_eqs + 1; i <= max_no_eqs + iEqNo; i++) {
           frp->jdiag[i-1] = (int) 0;
       }

       max_no_eqs += iEqNo;
    }

#ifdef DEBUG
       printf("*** Leave CheckJdiagSpace()\n");
#endif

    return (frp);
}

/* =================== */
/* Check Element Space */
/* =================== */
 
#ifdef __STDC__
EFRAME *CheckElementSpace(EFRAME *frp, int element_no)
#else
EFRAME *CheckElementSpace(frp, element_no)
EFRAME *frp;
int    element_no;
#endif
{
ELEMENT             *el;
int         i, j, k, ii;
int        UNITS_SWITCH;

SYMBOL             *slp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;

#ifdef DEBUG
       printf("*** Enter CheckElementSpace()");
       printf(" : element_no = %4d : max_no_elements = %4d\n", element_no, max_no_elements);
#endif

    slp = lookup("InPlaneIntegPts");   /* number of integration pts in plane/surface      */
    if(slp == NULL)
       iInPlaneIntegPts = 2*2;        /* 2x2 as default */
    else
       iInPlaneIntegPts = (int) slp->u.q->value;

    slp = lookup("ThicknessIntegPts"); /* number of integration pts in thickness direction*/
    if(slp == NULL)
       iThicknessIntegPts = 2;        /* 2 as default */
    else
       iThicknessIntegPts = (int) slp->u.q->value;

    iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;

    UNITS_SWITCH = CheckUnits();
    if(element_no > max_no_elements) {
       frp->no_elements = element_no;
       frp->element = (ELEMENT *) realloc(frp->element,
                      (max_no_elements + element_no)*sizeof(ELEMENT));
 
#ifdef DEBUG
       printf("*** In CheckElementSpace() : REALLOCATE SPACE()\n");
#endif

       for(i=max_no_elements + 1 ; i <= max_no_elements + element_no; i++) {
#ifdef DEBUG
       printf("*** In CheckElementSpace() : i = %d\n", i);
#endif
           el = &frp->element[i-1];
           el->elmt_attr_name = (char *)NULL;
           el->elmt_attr_no   = (int) 0;
           el->node_connect   = iVectorAlloc(frp->no_nodes_per_elmt + 1);
           el->node_connect   = (int *) MyCalloc((frp->no_nodes_per_elmt + 1), sizeof(int));
           el->d_array        = (int *) MyCalloc(frp->no_nodes_per_elmt, sizeof(int));
           el->rp             = (RESPONSE *) MyCalloc(1,sizeof(RESPONSE));
           el->esp            = (ELEMENT_STATE *) MyCalloc(1,sizeof(ELEMENT_STATE));
           el->esp->state     = (int) 0;

           el->rp->Forces     = MatrixAllocIndirect("nodal_force",
                                DOUBLE_ARRAY, frp->no_dof, UNIT_NODES);
           el->LC_ptr         = (MATER_LOAD_CURVE *) MyCalloc(1,sizeof(MATER_LOAD_CURVE));
           el->LC_ptr->name   = (char *)NULL;
           el->LC_ptr->R      = (double *) MyCalloc(iNO_INTEG_pts,sizeof(double));
           el->LC_ptr->H      = (double *) MyCalloc(iNO_INTEG_pts,sizeof(double));
           el->LC_ptr->back_stress = MatrixAllocIndirectDouble(6, iNO_INTEG_pts);
           el->rp->stress     = MatrixAllocIndirect("nodalforce",DOUBLE_ARRAY,9,iNO_INTEG_pts);
           el->rp->displ      = MatrixAllocIndirect("displ",DOUBLE_ARRAY,frp->no_dof,UNIT_NODES);
           el->rp->strain_pl  = MatrixAllocIndirect("strain_pl",DOUBLE_ARRAY,9,iNO_INTEG_pts);
           el->rp->strain_pl_incr  = MatrixAllocIndirect("strain_pl_incr",DOUBLE_ARRAY,9,iNO_INTEG_pts);
           el->rp->effect_pl_strain      = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
           el->rp->eff_pl_strain_incr = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));

           if(UNITS_SWITCH == ON) {
             el->rp->min_moment.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             el->rp->max_moment.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             el->rp->min_shear.dimen   = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             el->rp->max_shear.dimen   = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             el->rp->Mzc.dimen         = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           }
       }  

       max_no_elements = max_no_elements + element_no;
    }

#ifdef DEBUG
       printf("*** Leave CheckElementSpace()\n");
#endif

    return(frp);
}


/* ======================================= */
/* Check the space for Rigid Body Elements */
/* ======================================= */

#ifdef __STDC__
EFRAME *CheckRigidSpace(EFRAME *frp, int rigid_no)
#else
EFRAME *CheckRigidSpace(frp, rigid_no)
EFRAME *frp;
int    rigid_no;
#endif
{
RIGID  *rig;
int    i;

#ifdef DEBUG
       printf("*** Enter CheckRigidSpace()");
       printf(" : rigid_no = %4d : max_no_rigid = %4d\n", rigid_no, max_no_rigid);
#endif

    if(rigid_no > max_no_rigid ) {
       frp->no_rigid  = rigid_no; /* Allocate space till the 'rigid_no' input */
       frp->rigid = (RIGID *) realloc(frp->rigid,(max_no_rigid + rigid_no)*sizeof(RIGID));
 
       for(i=max_no_rigid  + 1 ; i<= (max_no_rigid + rigid_no); i++) {
           rig           = &frp->rigid[i-1];
           rig->rest_dof = iVectorAlloc(6);
       }  
       max_no_rigid = max_no_rigid + rigid_no;
    }

#ifdef DEBUG
       printf("*** Leave CheckRigidSpace()\n");
#endif

    return(frp);
}

/* ======================= */
/* Check load_nforce space */
/* ======================= */
 
#ifdef __STDC__
EFRAME *CheckNforcesSpace(EFRAME *frp, int nforce_no)
#else
EFRAME *CheckNforcesSpace(frp, nforce_no)
EFRAME *frp;
int    nforce_no;
#endif
{
int    i, j, k;
int  UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter CheckNforecesSpace()");
       printf(" : nforce_no = %4d : max_no_nforces = %4d\n", nforce_no, max_no_nforces);
#endif

     UNITS_SWITCH = CheckUnits();
    if(nforce_no > max_no_nforces) {
       frp->no_node_loads = nforce_no;
       frp->nforces = (NODE_LOADS *) realloc(frp->nforces, 
                      (max_no_nforces + nforce_no)*sizeof(NODE_LOADS));

       for(i = max_no_nforces + 1; i <= (max_no_nforces + nforce_no); i++){
           frp->nforces[i-1].fn = (QUANTITY *) MyCalloc(frp->no_dof, sizeof(QUANTITY));
           if( UNITS_SWITCH == ON ) {
              for(j = 1; j <= frp->no_dof; j++) {
                  frp->nforces[i-1].fn[j-1].dimen  = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
              }
           }
       }

       max_no_nforces += nforce_no;
    }

#ifdef DEBUG
       printf("*** Leave CheckNforecesSpace()\n");
#endif
    return(frp);
}

/* ======================= */
/* Check load_eforce space */
/* ======================= */

#ifdef __STDC__
EFRAME *CheckEforcesSpace(EFRAME *frp, int loaded_elmts)
#else
EFRAME *CheckEforcesSpace(frp, loaded_elmts)
EFRAME *frp;
int    loaded_elmts;
#endif
{
ELEMENT_LOADS   *elsp;
int                 i;
int      UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter CheckEforcesSpace()");
       printf(" : loaded_elmts = %4d : max_no_loaded_elmts = %4d\n", loaded_elmts, max_no_loaded_elmts);
#endif

  UNITS_SWITCH = CheckUnits();
  if(loaded_elmts > max_no_loaded_elmts) {
     frp->eforces = (ELEMENT_LOADS *) realloc(frp->eforces,
                    (loaded_elmts+max_no_loaded_elmts)*sizeof(ELEMENT_LOADS));

     for(i= max_no_loaded_elmts + 1 ; i <= loaded_elmts+max_no_loaded_elmts; i++){
        elsp = &frp->eforces[i -1];
        elsp->no_loads_faces  = 1;                        
        elsp->elib_ptr = (ELOAD_LIB *)calloc(elsp->no_loads_faces,sizeof(ELOAD_LIB));
        elsp->face_direc           = (double *) MyCalloc(frp->no_dimen, sizeof(double));
        elsp->elib_ptr             = (ELOAD_LIB *) MyCalloc(UNIT_NODES, sizeof(ELOAD_LIB));
        elsp->elib_ptr->name       =  NULL;
        elsp->elib_ptr->nopl       = (int *) MyCalloc(4, sizeof(int));
        elsp->elib_ptr->pr         = (MATRIX *) MyCalloc(1,sizeof(MATRIX));
        elsp->elib_ptr->pr->iNoRows      = 4;
        elsp->elib_ptr->pr->iNoColumns   = 3;
        elsp->elib_ptr->pr->cpMatrixName = (char *)NULL;
        if( UNITS_SWITCH == ON ) {
           elsp->elib_ptr->pr->spColUnits   = BufferInit(3);
           elsp->elib_ptr->pr->spRowUnits   = BufferInit(4);
        }
     }
     max_no_loaded_elmts += loaded_elmts;
  }

#ifdef DEBUG
       printf("*** Leave CheckEforcesSpace()\n");
#endif
 
  return(frp);
}
