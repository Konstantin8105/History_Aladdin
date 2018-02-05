/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_allocate.c : Allocation Functions for Finite Elements
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

#include <stdio.h>
#include <ctype.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"

/*
#define DEBUG
*/

static int num_elmt_attrs = 0;

int  max_no_nodes;
int  max_no_eqs;
int  max_no_elements;
int  max_no_rigid;
int  max_no_nforces;
int  max_no_loaded_nodes;
int  max_no_loaded_elmts;


/* =========================================================== */
/* Allocate space for Element, Section and Material Attributes */
/* =========================================================== */

ELEMENT_ATTR *Alloc_Element_Attr_Item() {
ELEMENT_ATTR *eap;
int             i;

#ifdef DEBUG
       printf("*** In Alloc_Element_Attr_Item ()\n");
#endif

      eap                   = (ELEMENT_ATTR *) MyCalloc(1,sizeof(ELEMENT_ATTR));
      eap->name             = (char *)NULL;
      eap->elmt_type        = (char *)NULL;
      eap->section          = (char *)NULL;
      eap->material         = (char *)NULL;
      eap->fiber_attr_name  = (char *)NULL;
      eap->map_ldof_to_gdof = (int  *) MyCalloc(6, sizeof(int));
      eap->work_material    = (QUANTITY *) MyCalloc(UNIT_MATERIAL_ATTR, sizeof(QUANTITY));
      eap->work_section     = (QUANTITY *) MyCalloc(UNIT_SECTION_ATTR, sizeof(QUANTITY));

#ifdef DEBUG
       printf("*** leaving Alloc_Element_Attr_Item ()\n");
#endif
      return(eap);
}

SECTION_ATTR *Alloc_Section_Attr_Item() {
SECTION_ATTR *sp;

      sp = (SECTION_ATTR *) MyCalloc(1, sizeof(SECTION_ATTR));
      return(sp);
}

MATERIAL_ATTR *Alloc_Material_Attr_Item() {
MATERIAL_ATTR       *mp;
SYMBOL             *slp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;

      slp = lookup("InPlaneIntegPts");   /* number of integration pts in plane/surface      */
      if(slp == NULL)
         iInPlaneIntegPts = 2*2;         /* 2x2 as default */
      else
         iInPlaneIntegPts = (int) slp->u.q->value;

      slp = lookup("ThicknessIntegPts"); /* number of integration pts in thickness direction*/
      if(slp == NULL)
         iThicknessIntegPts = 2;        /* 2 as default */
      else
         iThicknessIntegPts = (int) slp->u.q->value;

      iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;

      mp                = (MATERIAL_ATTR *) MyCalloc(1, sizeof(MATERIAL_ATTR));
      mp->alpha_thermal = (QUANTITY *) MyCalloc(3, sizeof(QUANTITY));
      mp->LC_ptr        = (MATER_LOAD_CURVE *) MyCalloc(1,sizeof(MATER_LOAD_CURVE));
      mp->LC_ptr->back_stress = MatrixAllocIndirectDouble(6, iNO_INTEG_pts);
      mp->LC_ptr->R     = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
      mp->LC_ptr->H     = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
      mp->LC_ptr->ialph = 0;
      mp->LC_ptr->pen   = 0;
      mp->LC_ptr->load[0] = 0;
      mp->LC_ptr->load[1] = 0;
      mp->LC_ptr->load[2] = 0;
      mp->LC_ptr->load[3] = 0;
      mp->LC_ptr->load[4] = 0;
      mp->LC_ptr->load[5] = 0;

      return(mp);
}

#ifdef __STDC__
FIBER_ELMT *Alloc_Fiber_Elmt_Attr_Item( int no_fiber ) 
#else
FIBER_ELMT *Alloc_Fiber_Elmt_Attr_Item( no_fiber ) 
int	no_fiber;
#endif
{
FIBER_ELMT *fep;

#ifdef DEBUG
       printf("*** In Alloc_Fiber_Elmt_Attr_Item ()\n");
#endif

      fep            = (FIBER_ELMT *) MyCalloc(1,sizeof(FIBER_ELMT));
      fep->no_fiber  = no_fiber;
      fep->fiber     = (FIBER_ATTR *) MyCalloc(no_fiber, sizeof(FIBER_ATTR));

#ifdef DEBUG
       printf("*** Leaving Alloc_Fiber_Elmt_Attr_Item ()\n");
#endif
      return(fep);
}

/*
 *  --------------------------------
 *  Mapping from global to local dof
 *  --------------------------------
 */

#ifdef __STDC__
void Ldof_to_gdof(ELEMENT_ATTR *eap, MATRIX  *m_local, MATRIX *m_global)
#else
void Ldof_to_gdof(eap, m_local, m_global)
ELEMENT_ATTR                 *eap;
MATRIX        *m_local, *m_global;
#endif
{
int    i, j, k, n;
int ldof, gdof;

#ifdef DEBUG
   printf(" #### Enter Ldof_to_gdof() \n");
#endif

  if(m_global->iNoRows == 1) {
     n = MAX(m_global->iNoRows, m_global->iNoColumns);
  }
  if( n  < 6) {
      if(m_global->iNoRows == 1) {
         for(j = 1; j <= m_global->iNoColumns; j++)
             eap->map_ldof_to_gdof[(int)(m_local->uMatrix.daa[0][j-1]-1)] 
             = (int) m_global->uMatrix.daa[0][j-1];
      }
      if(m_global->iNoColumns == 1) {
         for(j = 1; j <= m_global->iNoRows; j++)
             eap->map_ldof_to_gdof[(int)(m_local->uMatrix.daa[j-1][0]-1)] 
             = (int) m_global->uMatrix.daa[j-1][0];
      }
  }
  else
     if(n > 6)
        FatalError("number maping dofs can not more than 6",(char *)NULL);
#ifdef DEBUG
   printf(" #### Leaving Ldof_to_gdof() \n");
#endif
}


/* ===================================================== */
/* Alloc Frame                                           */
/* function to allocate space for datastructure EFRAME   */
/* for  a fixed size ---as defined by  parameters in     */
/* defs.h file                                           */
/*                  UNIT_NDM,UNIT_NDF,UNIT_NEN,UNIT_NAD  */
/*                  UNIT_ELEMENTS                        */
/*                  UNIT_NODES,UNIT_NFORCES,UNIT_EFORCES */
/*                  UNIT_SECTION_ATTR,UNIT_MATERIAL_ATTR */
/* ====================================================  */

EFRAME *FrameAlloc()
{
EFRAME            *frp;
NODE               *np;
ELEMENT            *ep;
ELEMENT_ATTR      *eap;
ELEMENT_LOADS    *elsp;
ELOAD_LIB         *elp;
NODE_LOADS        *nlp;
RIGID             *rig;
SYMBOL             *sp;
int  i,j,k, el_attr_no;

int       UNITS_SWITCH;
int      iNO_INTEG_pts;
int   iInPlaneIntegPts;
int iThicknessIntegPts;

#ifdef DEBUG
       printf("*** Enter FrameAlloc()\n");
#endif

     UNITS_SWITCH = CheckUnits();

    /* --------------------------------------------- */
    /* Allocate and initialize EFRAME data structure */
    /* --------------------------------------------- */

    frp = (EFRAME *) MyCalloc(1,sizeof(EFRAME));

    sp = lookup("NDimension");         /* No of dimensions : 2 or 3 */
    frp->no_dimen = (int) sp->u.q->value;

    sp = lookup("NDofPerNode");        /* No gdof per node */
    frp->no_dof = (int) sp->u.q->value;

    sp = lookup("MaxNodesPerElement"); /* Max no nodes per element */
    frp->no_nodes_per_elmt = (int) sp->u.q->value;

    sp = lookup("GaussIntegPts"); /* section no.(integration pts) along element */
    if(sp == NULL)
       frp->no_integ_pt = UNIT_INTEG_PTS;        /* 2 as default */
    else {
       if( (int)sp->u.q->value < 2 ) {
          printf("\nWARNING : The value of GaussIntegPts can not less than 2!");
          printf("\n          Use 2 Gauss Integration points.\n");
          sp->u.q->value = (double) 2;
       }
       else if( (int)sp->u.q->value > 10 ) {
          printf("\nWARNING : ALADDIN does not support GaussIntegPts more than 10!");
          printf("\n          Use 10 Gauss Integration points.\n");
          sp->u.q->value = (double) 10;
       }
       frp->no_integ_pt = (int) sp->u.q->value;
    }

    sp = lookup("InPlaneIntegPts");   /* number of integration pts in plane/surface      */
    if(sp == NULL) 
       iInPlaneIntegPts = 2*2;        /* 2x2 as default */
    else
       iInPlaneIntegPts = (int) sp->u.q->value;

    sp = lookup("ThicknessIntegPts"); /* number of integration pts in thickness direction*/
    if(sp == NULL) 
       iThicknessIntegPts = 2;        /* 2 as default */
    else
       iThicknessIntegPts = (int) sp->u.q->value;

    /* ---------------------------------------------- */
    /* Allocate Pointers to component data structures */
    /* ---------------------------------------------- */

    frp->node      = (NODE *)    MyCalloc( UNIT_NODES,    sizeof(NODE));
    frp->element   = (ELEMENT *) MyCalloc( UNIT_ELEMENTS, sizeof(ELEMENT));
    frp->eattr     = (ELEMENT_ATTR *) MyCalloc(UNIT_ELEMENT_ATTR, sizeof(ELEMENT_ATTR));
    frp->rigid     = (RIGID *) MyCalloc( UNIT_RIGIDS,   sizeof(RIGID));
    frp->jdiag     = (int *)   MyCalloc((UNIT_NODES * UNIT_NDF), sizeof(int));
    frp->nforces   = (NODE_LOADS *) MyCalloc( UNIT_NFORCES, sizeof(NODE_LOADS));
    frp->eforces   = (ELEMENT_LOADS *) MyCalloc( UNIT_ELEMENTS, sizeof(ELEMENT_LOADS));

    /* ----------------------------------------------------- */
    /* Allocation of Arrays within Component Data Structures */ 
    /* ----------------------------------------------------- */

    for(i=1; i <= UNIT_NODES; i++) {

        np = &frp->node[i-1];
        np->coord = (QUANTITY *) MyCalloc(frp->no_dimen, sizeof(QUANTITY));
        if( UNITS_SWITCH == ON ) {
           for(j = 1; j <= frp->no_dimen; j++){
              np->coord[j-1].dimen    =(DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           }
        }

        np->bound_id = (int *) MyCalloc(frp->no_dof, sizeof(int));
        np->rb_num   = (int ) 0;

        np->disp     = (QUANTITY *) MyCalloc(frp->no_dof, sizeof(QUANTITY));
        if( UNITS_SWITCH == ON ) {
           for (j = 1; j <= frp->no_dof; j++) {
              np->disp[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           }
        }

        np->TrT = (MATRIX *) MyCalloc(1,sizeof(MATRIX));
    }

    iNO_INTEG_pts = iThicknessIntegPts*iInPlaneIntegPts;

    for(i=1; i<= UNIT_ELEMENTS; i++) {

        ep                 = &frp->element[i-1];
        ep->elmt_attr_name = (char *) NULL;
        ep->elmt_attr_no   = (int) 0;
        ep->node_connect   = (int  *) MyCalloc((frp->no_nodes_per_elmt+1),sizeof(int));
        ep->d_array        = (int  *) MyCalloc(frp->no_nodes_per_elmt, sizeof(int));
        ep->rp             = (RESPONSE *) MyCalloc(1,sizeof(RESPONSE));
        ep->esp            = (ELEMENT_STATE *) MyCalloc(1,sizeof(ELEMENT_STATE));
        ep->esp->state     = (int) 0;

        ep->LC_ptr         = (MATER_LOAD_CURVE *) MyCalloc(1,sizeof(MATER_LOAD_CURVE));
        ep->LC_ptr->name   = (char *) NULL;
        ep->LC_ptr->R      = (double *) MyCalloc(iNO_INTEG_pts,sizeof(double));
        ep->LC_ptr->H      = (double *) MyCalloc(iNO_INTEG_pts,sizeof(double));
        ep->LC_ptr->back_stress = MatrixAllocIndirectDouble(6, iNO_INTEG_pts);

        ep->rp->Forces     = MatrixAllocIndirect("nodal_force",DOUBLE_ARRAY,frp->no_dof,UNIT_NODES);

        ep->rp->displ     = MatrixAllocIndirect("displ",    DOUBLE_ARRAY, frp->no_dof, UNIT_NODES);
        ep->rp->stress    = MatrixAllocIndirect("stress",   DOUBLE_ARRAY, 9, iNO_INTEG_pts );
        ep->rp->strain_pl = MatrixAllocIndirect("strain_pl",DOUBLE_ARRAY, 9, iNO_INTEG_pts );
        ep->rp->strain_pl_incr     = MatrixAllocIndirect("strain_pl_incr",DOUBLE_ARRAY,9,iNO_INTEG_pts);
        ep->rp->effect_pl_strain   = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));
        ep->rp->eff_pl_strain_incr = (double *) MyCalloc(iNO_INTEG_pts, sizeof(double));

        if(UNITS_SWITCH == ON) {
           ep->rp->min_moment.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           ep->rp->max_moment.dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           ep->rp->min_shear.dimen   = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           ep->rp->max_shear.dimen   = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           ep->rp->Mzc.dimen         = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
        }
    }

    for(i=1; i<= UNIT_ELEMENT_ATTR; i++) {

        eap                   = &frp->eattr[i-1];
        eap->name             = (char *) NULL;
        eap->elmt_type        = (char *) NULL;
        eap->map_ldof_to_gdof = (int  *) MyCalloc(6,sizeof(int));
        eap->material         = (char *) NULL;
        eap->section          = (char *) NULL;
        eap->work_material    = (QUANTITY *) MyCalloc(UNIT_MATERIAL_ATTR,sizeof(QUANTITY)); 
        eap->work_section     = (QUANTITY *) MyCalloc(UNIT_SECTION_ATTR,sizeof(QUANTITY));

        if( UNITS_SWITCH == ON ) {
           for(j = 1; j <= UNIT_MATERIAL_ATTR; j++) {
              eap->work_material[j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              ZeroUnits(eap->work_material[j-1].dimen);
           }
           for(j = 1; j <= UNIT_SECTION_ATTR; j++) {
              eap->work_section[j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
              ZeroUnits(eap->work_section[j-1].dimen);
           }
        }
    }

   /* In the following: 10 is arbitrary and is subjected to change*/

   for(i=1; i <= UNIT_RIGIDS;i++) {
        rig                  = &frp->rigid[i-1];
        rig->rbody_attr_name = (char *) NULL;
        rig->in              = (int *) MyCalloc(10, sizeof(int)); 
        rig->rest_dof        = (int *) MyCalloc(6, sizeof(int));
        rig->per_nodes       = (int *) MyCalloc(10, sizeof(int));
        rig->prop            = (QUANTITY *) MyCalloc(10, sizeof(QUANTITY));  
        if( UNITS_SWITCH == ON ) {
           for(j = 1; j <= 10; j++) 
              rig->prop[j-1].dimen = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
        }
   }

   for(i=1; i <= UNIT_ELEMENTS; i++) {
        elsp                       = &frp->eforces[i-1];
        elsp->no_loads_faces       = 1;                        
        elsp->face_direc           = (double *) MyCalloc(frp->no_dimen, sizeof(double));
        elsp->elib_ptr             = (ELOAD_LIB *) MyCalloc(UNIT_NODES, sizeof(ELOAD_LIB));
        elsp->elib_ptr->name       = (char *) NULL;
        elsp->elib_ptr->nopl       = (int *) MyCalloc(4, sizeof(int));
        elsp->elib_ptr->pr         = (MATRIX *) MyCalloc(1,sizeof(MATRIX));
        elsp->elib_ptr->face_no           = (int) 0;
        elsp->elib_ptr->numnp_face        = (int) 0;
        elsp->elib_ptr->pr->iNoRows       = 4;
        elsp->elib_ptr->pr->iNoColumns    = 3;
        elsp->elib_ptr->pr->cpMatrixName  = (char *)NULL;
        if( UNITS_SWITCH == ON ) {
            elsp->elib_ptr->pr->spRowUnits = BufferInit(4);
            elsp->elib_ptr->pr->spColUnits = BufferInit(3);
        }
    }

    for(i=1; i <= UNIT_NODES; i++) {
       elp                       = &elsp->elib_ptr[i-1];
       elp->body_force           = (QUANTITY *) MyCalloc(frp->no_dimen, sizeof(QUANTITY));
       elp->init_stress          = (QUANTITY *) MyCalloc(frp->no_dimen, sizeof (QUANTITY));
       elp->init_strain          = (double *) MyCalloc(frp->no_dimen, sizeof (double));
       elp->traction             = (QUANTITY *) MyCalloc(frp->no_dimen, sizeof(QUANTITY));
       if( UNITS_SWITCH == ON ) {
          for(j = 1; j <= frp->no_dimen; j++) {
              elp->body_force[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
          for(j = 1; j <= frp->no_dimen; j++) {
              elp->init_stress[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
          for(j = 1; j <= frp->no_dimen; j++) {
              elp->traction[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
       }
    }

    for (i=1; i <= UNIT_NFORCES; i++) {
       nlp         = &frp->nforces[i-1];
       nlp->fn     = (QUANTITY *) MyCalloc(frp->no_dof, sizeof(QUANTITY));
       nlp->node_f = (int) 0; 
       if( UNITS_SWITCH == ON ) {
          for(j = 1; j <= frp->no_dof; j++) {
             nlp->fn[j-1].dimen  = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          }
       }
    }

    /*------------------------------------------------------*/
    /* Assignment of max_no_ parameters for CheckSpace()    */
    /*------------------------------------------------------*/

    max_no_nodes        = UNIT_NODES;
    max_no_elements     = UNIT_ELEMENTS;
    max_no_nforces      = UNIT_NFORCES;
    max_no_loaded_elmts = UNIT_EFORCES;
    max_no_rigid        = UNIT_RIGIDS;
    max_no_eqs          = UNIT_NODES * UNIT_NDF;

#ifdef DEBUG
       printf("*** Leave FrameAlloc()\n\n");
#endif

    return(frp);
}

