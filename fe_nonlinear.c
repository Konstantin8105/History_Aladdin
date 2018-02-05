/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_nonlinear.c : Functions Needed for Nonlinear Finite Element Analysis
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
 *  Modified by: Mark Austin                                           March 2000
 *  ============================================================================= 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "fe_functions.h"
#include "elmt.h"

extern ARRAY     *array;
extern EFRAME    *frame;

/* #define DEBUG */

/* ======================================================================= */
/* Setup pointers to working arrays : notes static flag limits             */
/* scope of arrays to fe_nonlinear.c (i.e., don't break this file apart!). */
/* ======================================================================= */

static RESPONSE                *RespondBuff;
static int                   *ElmtStateBuff;
static MATER_LOAD_CURVE          *LoadCurve;
static FIBER_RESPONSE   *FiberRespondBuffer;
static int                  **FiberLoadFlag;

/*
 * =======================================================================
 * Element State Determination
 *
 * Input:  -- MATRIX *d_displ : matrix of global displacements.
 * Output: -- void.
 * =======================================================================
 */

#ifdef  __STDC__
void Elmt_State_Det( MATRIX *d_displ )
#else
void Elmt_State_Det( d_displ )
MATRIX *d_displ;
#endif
{
MATRIX  *dp;
ELEMENT  *el;
ELEMENT_ATTR  *eap;
int  dof_per_elmt;
int  node_no, elmt_no;
int  i, j, k, ii, jj; 
HISTORY_DATA  *hp;
int  UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter Elmt_State_Det()\n");
#endif

    UNITS_SWITCH = CheckUnits();

    for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) { 
       el  = &frame->element[elmt_no-1];
       eap = &(frame->eattr[el->elmt_attr_no-1]);

       /* If it is a fiber element, start element state determation */
       /* This block of code should be extended later               */

       if( !(strcmp(eap->elmt_type, "FIBER_2D"))  ||
           !(strcmp(eap->elmt_type, "FIBER_3D"))  ||
           !(strcmp(eap->elmt_type, "FIBER_2DS")) ||
           !(strcmp(eap->elmt_type, "FIBER_3DS")) ) {

          /* Assign element propty */

          for (i = 0; i < NO_ELEMENTS_IN_LIBRARY; i++) { 
             if(!strcmp(elmt_library[i].name, eap->elmt_type)) { 
                k = i;
                break;
             }
          }

          if(array->elmt_no != elmt_no) { /* check if already values assigned */
             array->elmt_no        =  elmt_no;
             array->elmt_attr_no   =  el->elmt_attr_no;
             array->dof_per_node   =  elmt_library[k].no_dof;
             array->no_dimen       =  elmt_library[k].no_dimen;
             array->elmt_type      =  elmt_library[k].name;
             array->nodes_per_elmt = (int) MIN(frame->no_nodes_per_elmt,
                                           elmt_library[k].no_node_per_elmt);
             array->size_of_stiff = array->nodes_per_elmt* array->dof_per_node;

             free((char *) array->material_name);
             array->material_name = SaveString( eap->material );
          }
          array = Assign_p_Array(frame, elmt_no, array, PROPTY);

          /* Assign element coordinates */

          for(j = 1; j <= array->nodes_per_elmt; j++) {
             node_no = el->node_connect[j-1];
             if(node_no != 0) {
                for(i=1;i <= array->no_dimen; i++) {
                   array->coord[i-1][j-1].value = frame->node[node_no -1].coord[i-1].value;
                   if(UNITS_SWITCH == ON )
                      UnitsCopy( array->coord[i-1][j-1].dimen,
                                 frame->node[node_no -1].coord[i-1].dimen );
                }
             }
          }

          /* Transfer fixed displacements */

          for(i=1; i <= array->nodes_per_elmt; i++) {
             k = 1;
             node_no = el->node_connect[i-1];
             for(j = 1; j <= frame->no_dof; j++) {
                jj = frame->node[node_no - 1].bound_id[j-1];  
                if(jj > 0) {
                   array->displ->uMatrix.daa[j-1][i-1] = d_displ->uMatrix.daa[jj-1][0];
                   if( UNITS_SWITCH == ON ) {
                      UnitsCopy(&(array->displ->spRowUnits[j-1]),
                                &(d_displ->spRowUnits[jj-1]));
                      UnitsCopy(&(array->displ->spColUnits[i-1]),
                                &(d_displ->spColUnits[0]));
                   }
                } else {
                   array->displ->uMatrix.daa[j-1][i-1] = frame->node[node_no -1].disp[j-1].value;
                   if( UNITS_SWITCH == ON ) {
                      UnitsCopy(&(array->displ->spRowUnits[j-1]),
                                frame->node[node_no -1].disp[j-1].dimen);
                      UnitsCopy(&(array->displ->spColUnits[i-1]),
                                &(d_displ->spColUnits[0]));
                   }
                }
             }
          }

          /* Retrieve Previous History From The RespondBuffer */

          if( FiberRespondBuffer == (FIBER_RESPONSE *)NULL )
              FatalError("Must have fiber element to use ElmtStateDet()!",
                        (char *) NULL);

          for( ii=0 ; ii < FiberRespondBuffer->total_fiber_elmt ; ++ii ) {
               if( elmt_no == FiberRespondBuffer->history[ii].elmt_no )
                   break;
          }

          hp = &FiberRespondBuffer->history[ii];

          if( !(strcmp(eap->elmt_type, "FIBER_2D")) ||
              !(strcmp(eap->elmt_type, "FIBER_2DS")) )
             Fiber_Elmt_State_Det_2d( array, hp, FiberLoadFlag[ii] );
          else
             Fiber_Elmt_State_Det_3d( array, hp, FiberLoadFlag[ii] );

          el->esp->state = array->elmt_state;
       }
    }

#ifdef DEBUG
       printf("*** Leave Elmt_State_Det()\n");
#endif

}


/*
 * =======================================================================
 * UpdateResponse() : Save element response to frame data structures.
 *
 * Input:  -- void.
 * Output: -- void.
 * =======================================================================
 */

void UpdateResponse() {
ELEMENT             *ep;
ELEMENT_ATTR       *eap;
int             elmt_no;
int        elmt_attr_no;
int                i, j;
int    total_shell_elmt;

SYMBOL             *slp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;

#ifdef DEBUG
       printf("*** Enter UpdateResponse() \n");
#endif

  total_shell_elmt = 0;

  for(i=1; i<=frame->no_elements; i++) {

    elmt_attr_no = frame->element[i-1].elmt_attr_no;
    eap = &(frame->eattr[elmt_attr_no-1]);

    /*
     * -------------------------------------------------------------------
     * Fiber elements exist, update the load history.
     * Including the flags, stresses and strains for each fiber element.
     * The updated load history will be assign to frame.
     * The related details are in elmt_fiber.c and be done in elmt_fiber.c
     * -------------------------------------------------------------------
     */

    if( !(strcmp(eap->elmt_type, "FIBER_2D"))  ||
        !(strcmp(eap->elmt_type, "FIBER_3D"))  ||
        !(strcmp(eap->elmt_type, "FIBER_2DS")) ||
        !(strcmp(eap->elmt_type, "FIBER_3DS")) ) {
          ep = &(frame->element[i-1]);
          j  = eap->work_fiber->no_fiber +
               (frame->no_dimen-1)*eap->work_fiber->no_shear;
          SaveFiberRespondBuffer( ep, i, frame->no_integ_pt+2, j );
    }

    if( !(strcmp(eap->elmt_type, "SHELL_4N")) ||
        !(strcmp(eap->elmt_type, "SHELL_8N")) )
        total_shell_elmt++;
  }

  /* Update load history for SHELL_4N or SHELL_8N elements */

  if( total_shell_elmt != 0 ) {

      /* no of integration pts in plane/surface */

      slp = lookup("InPlaneIntegPts");  
      if(slp == NULL)
         iInPlaneIntegPts = 2*2;        /* 2x2 as default */
      else
         iInPlaneIntegPts = (int) slp->u.q->value;

      /* no integration pts in thickness direction */

      slp = lookup("ThicknessIntegPts"); 
      if(slp == NULL)
         iThicknessIntegPts = 2;        /* 2 as default */
      else
         iThicknessIntegPts = (int) slp->u.q->value;

      iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;

      for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {
          ep = &frame->element[elmt_no-1];
          save_action(array, ep, &frame->eattr[ep->elmt_attr_no -1],
                      elmt_no, iNO_INTEG_pts);
      }
  }

#ifdef DEBUG
       printf("*** Leaving UpdateResponse() \n");
#endif

}


/*
 * =======================================================================
 * Save Action ()
 *
 * Input:  -- ARRAY    *p       --
 *         -- ELEMENT *ep       --
 *         -- ELEMENT_ATTR *eap --
 *         -- int elmt_no       --
 *         -- int iNO_INTEG_pts --
 * Output: -- void.
 * =======================================================================
 */

#ifdef __STDC__
void save_action(ARRAY *p, ELEMENT *ep, ELEMENT_ATTR *eap,
                 int elmt_no, int iNO_INTEG_pts )
#else
void save_action(p, ep, eap, elmt_no, iNO_INTEG_pts)
ELEMENT          *ep;
ELEMENT_ATTR    *eap;
ARRAY             *p;
int          elmt_no;
int    iNO_INTEG_pts;

#endif
{
char    *name;
int   i, j, k;
DIMENSIONS *d;

#ifdef DEBUG
       printf("*** Enter save_actions() \n");
#endif

  name = eap->elmt_type;  
  for(i = 1; i <= frame->no_dof; i++ ) {
  for(j = 1; j <= frame->no_nodes_per_elmt; j++) {
       k = frame->no_dof*(j-1) + i;
       ep->rp->Forces->uMatrix.daa[i-1][j-1]
       = RespondBuff[elmt_no-1].Forces->uMatrix.daa[i-1][j-1];
       ep->rp->displ->uMatrix.daa[i-1][j-1] 
       = RespondBuff[elmt_no-1].displ->uMatrix.daa[i-1][j-1];
  }
  }

  ep->rp->max_moment.value = RespondBuff[elmt_no-1].max_moment.value;
  ep->esp->state           = ElmtStateBuff[elmt_no-1];

  /* state, strains parameters */

  switch(ep->esp->state) {
     case 0:   /* ELASTIC state    */
          for(j = 1; j <= iNO_INTEG_pts; j++) {
          ep->rp->effect_pl_strain[j-1] = 0.0;
          for(i = 1; i <= 9; i++) 
              ep->rp->stress->uMatrix.daa[i-1][j-1] 
              = RespondBuff[elmt_no-1].stress->uMatrix.daa[i-1][j-1];
          }
          break;
    case 1:   /* Perfectly plastic or  */
              /* elastic plastic state */
         if(LoadCurve[elmt_no-1].name != NULL) 
            ep->LC_ptr->name = SaveString(LoadCurve[elmt_no-1].name);
         else 
            ep->LC_ptr->name = (char *)NULL;

         for(j = 1; j <= iNO_INTEG_pts; j++) {
            ep->rp->effect_pl_strain[j-1]
               += RespondBuff[elmt_no-1].eff_pl_strain_incr[j-1];

            for(i = 1; i <= 9; i++) {
                ep->rp->stress->uMatrix.daa[i-1][j-1] 
                = RespondBuff[elmt_no-1].stress->uMatrix.daa[i-1][j-1];
                ep->rp->strain_pl->uMatrix.daa[i-1][j-1] 
                += RespondBuff[elmt_no-1].strain_pl_incr->uMatrix.daa[i-1][j-1];
            }
            ep->LC_ptr->R[j-1] = LoadCurve[elmt_no-1].R[j-1];
            ep->LC_ptr->H[j-1] = LoadCurve[elmt_no-1].H[j-1];

            for(i = 1; i <= 6; i++) 
                ep->LC_ptr->back_stress[i-1][j-1] =
                    LoadCurve[elmt_no-1].back_stress[i-1][j-1];
         }
         break;
    default:
         break;
  }

  if(CheckUnits() == ON) {
       switch(CheckUnitsType()) {
         case SI:
              d = DefaultUnits("Pa");
              break;
         case US:
              d = DefaultUnits("psi");
              break;
       }
       for(i = 1; i <= 9; i++) 
           UnitsCopy( &(ep->rp->stress->spRowUnits[i-1]), d );
       free((char *) d->units_name);
       free((char *) d);
  }

#ifdef DEBUG
        printf("******Leaving save_action(): \n");
#endif

}


/*
 *  ===================================================
 *  SetUpRespondBuffer() : Create response buffer......
 *
 *  Input:  -- void.
 *  Output: -- void.
 *  ===================================================
 */

#ifdef __STDC__
void SetUpRespondBuffer()
#else
void SetUpRespondBuffer()
#endif
{
int                i, j;
int         iNoIntegPts; 

SYMBOL             *slp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;

ELEMENT_ATTR   *eap;
FIBER_ELMT     *fep;

int	total_fiber_elmt;
int	total_shell_elmt;

#ifdef DEBUG
   printf("*** Enter SetUpRespondBuffer() \n");
   printf("*** No_elements = %d           \n", frame->no_elements);
#endif

  total_fiber_elmt = 0;
  total_shell_elmt = 0;

  for(i=1; i<=frame->no_elements; i++) {

     eap = lookup(frame->element[i-1].elmt_attr_name)->u.eap;
     if(eap == NULL)
        FatalError("Elmt_Attribute name not found",(char *)NULL);

     if( !(strcmp(eap->elmt_type, "FIBER_2D")) ||
         !(strcmp(eap->elmt_type, "FIBER_3D")) ||
         !(strcmp(eap->elmt_type, "FIBER_2DS")) ||
         !(strcmp(eap->elmt_type, "FIBER_3DS")) )
           total_fiber_elmt++;

     if( !(strcmp(eap->elmt_type, "SHELL_4N")) ||
         !(strcmp(eap->elmt_type, "SHELL_8N")) )
           total_shell_elmt++;
  }

  /*
   * ----------------------------------------------------------------------
   * Fiber elements exist, setup the needed space for storing load history.
   * Including the flags, stresses and strains for each fiber element.
   * The storage space will be static and assigned to frame until updated.
   * The related details are in elmt_fiber.c and be done in elmt_fiber.c
   * ----------------------------------------------------------------------
   */

  if( total_fiber_elmt != 0 )
      SetUpFiberRespondBuffer( total_fiber_elmt, frame );

  /* SHELL_4N or SHELL_8N elements exist, setup the needed space for storing load history. */

  if( total_shell_elmt != 0 ) {

    /* Number of integration pts in plane/surface */

    slp = lookup("InPlaneIntegPts"); 
    if(slp == NULL)
       iInPlaneIntegPts = 2*2;    /* 2x2 as default */
    else
       iInPlaneIntegPts = (int) slp->u.q->value;

    /* Number of integration pts in thickness direction*/

    slp = lookup("ThicknessIntegPts");
    if(slp == NULL)
       iThicknessIntegPts = 2;    /* 2 as default */
    else
       iThicknessIntegPts = (int) slp->u.q->value;

    iNoIntegPts = iInPlaneIntegPts*iThicknessIntegPts;


    RespondBuff = (RESPONSE *) MyCalloc(frame->no_elements, sizeof(RESPONSE)); 
    ElmtStateBuff = (int *)    MyCalloc(frame->no_elements, sizeof(int)); 
    LoadCurve     = (MATER_LOAD_CURVE *) MyCalloc(frame->no_elements,
                                           sizeof(MATER_LOAD_CURVE)); 

    for(i = 1; i <= frame->no_elements; i++) {
        RespondBuff[i-1].Forces 
          = MatrixAllocIndirect("NodalForce",DOUBLE_ARRAY,6,frame->no_nodes);
        RespondBuff[i-1].displ
          = MatrixAllocIndirect("NodalDispl",DOUBLE_ARRAY,6,frame->no_nodes);
        RespondBuff[i-1].stress
          = MatrixAllocIndirect("Stress",DOUBLE_ARRAY,9,iNoIntegPts);
        RespondBuff[i-1].strain_pl
          = MatrixAllocIndirect("PlasticStrain",DOUBLE_ARRAY,9,iNoIntegPts);
        RespondBuff[i-1].strain_pl_incr
          = MatrixAllocIndirect("PlasticStrainIncr",DOUBLE_ARRAY,9,iNoIntegPts);
        RespondBuff[i-1].effect_pl_strain
          = (double *) MyCalloc(iNoIntegPts, sizeof(double));
        RespondBuff[i-1].eff_pl_strain_incr 
          = (double *) MyCalloc(iNoIntegPts, sizeof(double));

        LoadCurve[i-1].name = (char *) MyCalloc(1,sizeof(char));
        LoadCurve[i-1].R    = (double *) MyCalloc(iNoIntegPts, sizeof(double));
        LoadCurve[i-1].H    = (double *) MyCalloc(iNoIntegPts, sizeof(double));
        LoadCurve[i-1].back_stress = MatrixAllocIndirectDouble(6, iNoIntegPts);
    }
  }  /* end of total_shell_elmt != 0, SHELL_4N or SHELL_8N elements exist */

#ifdef DEBUG
   printf("*** Leaving SetUpRespondBuffer(): \n");
#endif

}

/*
 * =======================================================================
 * SaveRespondBuffer() : Save response buffer..... 
 *
 * Input:  -- ARRAy *p     --
 *         -- int integ_pt --
 * Output: -- void.
 * =======================================================================
 */

#ifdef __STDC__
void SaveRespondBuffer(ARRAY *p, int integ_pt) 
#else
void SaveRespondBuffer(p, integ_pt) 
ARRAY     *p;
int integ_pt;
#endif
{
int i, j, k, kk;
      
#ifdef DEBUG
       printf("*** Enter SaveRespondBuffer(): \n");
#endif

    kk = integ_pt;
    ElmtStateBuff[p->elmt_no-1] = p->elmt_state;

    if(p->LC_ptr->name != NULL) 
       LoadCurve[p->elmt_no-1].name = SaveString(p->LC_ptr->name);
     else
       LoadCurve[p->elmt_no-1].name = NULL;
    
    LoadCurve[p->elmt_no-1].R[kk-1] = p->LC_ptr->R[kk-1];
    LoadCurve[p->elmt_no-1].H[kk-1] = p->LC_ptr->H[kk-1];

    for(i = 1; i <= 6; i++)
        LoadCurve[p->elmt_no-1].back_stress[i-1][kk-1]
        = p->LC_ptr->back_stress[i-1][kk-1];

    RespondBuff[p->elmt_no-1].effect_pl_strain[kk-1]
    = p->effect_pl_strain[kk-1];

    RespondBuff[p->elmt_no-1].eff_pl_strain_incr[kk-1]
    = p->eff_pl_strain_incr[kk-1];

    for(i = 1; i <= p->dof_per_node; i++) {
        RespondBuff[p->elmt_no-1].stress->uMatrix.daa[i-1][kk-1]
        = p->stress->uMatrix.daa[i-1][kk-1];
        RespondBuff[p->elmt_no-1].strain_pl->uMatrix.daa[i-1][kk-1]
        = p->strain_pl->uMatrix.daa[i-1][kk-1];
        RespondBuff[p->elmt_no-1].strain_pl_incr->uMatrix.daa[i-1][kk-1]
        = p->strain_pl_incr->uMatrix.daa[i-1][kk-1];
    }

    for(i = 1; i <= p->dof_per_node; i++ ) {
    for(j = 1; j <= p->nodes_per_elmt; j++) {
        k = p->dof_per_node*(j-1) + i;
        RespondBuff[p->elmt_no-1].Forces->uMatrix.daa[i-1][j-1] =
                    p->nodal_loads[k-1].value;
        RespondBuff[p->elmt_no-1].displ->uMatrix.daa[i-1][j-1] =
                    p->displ->uMatrix.daa[i-1][j-1];
    }
    }

    /* What's this ??? */

    if((p->elmt_type != NULL) && !strcmp(p->elmt_type, "FRAME_2D"))
        RespondBuff[p->elmt_no-1].max_moment.value 
           = MAX( p->nodal_loads[2].value, p->nodal_loads[5].value);
    if((p->elmt_type != NULL) && !strcmp(p->elmt_type, "FRAME_3D"))
        RespondBuff[p->elmt_no-1].max_moment.value
           = MAX( p->nodal_loads[5].value, p->nodal_loads[11].value);

#ifdef DEBUG
        printf("******Leaving SaveRespondBuffer(): \n");
#endif

}

/*
 * =======================================================================
 * SaveFiberRespondBuffer() : Save the static flags, stresses and strains
 *                            of each fiber element to frame.
 *
 * Input:  -- ELEMENT *el,
 *         -- int elmt_no,
 *         -- int no_section,
 *         -- int total_no_fiber
 * Output: -- void.
 * =======================================================================
 */

#ifdef  __STDC__
void SaveFiberRespondBuffer( ELEMENT *el, int elmt_no,
                             int no_section, int total_no_fiber )
#else
void SaveFiberRespondBuffer( el, elmt_no, no_section, total_no_fiber )
ELEMENT  *el;
int  elmt_no, no_section, total_no_fiber;
#endif
{
HISTORY_DATA      *hp;
int     ii, ifib, sec;

    for( ii=0 ; ii < FiberRespondBuffer->total_fiber_elmt ; ++ii ) {
       if( elmt_no == FiberRespondBuffer->history[ii].elmt_no )
          break;
    }

    hp = &FiberRespondBuffer->history[ii];

    for( sec=0 ; sec < no_section ; ++sec ) {
    for( ifib=0 ; ifib < total_no_fiber ; ++ifib ) {
       hp->pre_load[sec][ifib] = hp->loading[sec][ifib];
       el->esp->yielding_saved[sec][ifib]  = hp->yielding[sec][ifib];
       el->esp->pre_range_saved[sec][ifib] = hp->pre_range[sec][ifib];
       el->esp->pre_load_saved[sec][ifib]  = hp->pre_load[sec][ifib];

       el->rp->sr_saved->uMatrix.daa[sec][ifib] = hp->sr->uMatrix.daa[sec][ifib];
       el->rp->er_saved->uMatrix.daa[sec][ifib] = hp->er->uMatrix.daa[sec][ifib];
       el->rp->s0_saved->uMatrix.daa[sec][ifib] = hp->s0->uMatrix.daa[sec][ifib];
       el->rp->e0_saved->uMatrix.daa[sec][ifib] = hp->e0->uMatrix.daa[sec][ifib];
       el->rp->sx_saved->uMatrix.daa[sec][ifib] = hp->stress->uMatrix.daa[sec][ifib];
       el->rp->ex_saved->uMatrix.daa[sec][ifib] = hp->strain->uMatrix.daa[sec][ifib];
    }
    }

    /* Initial flag matrix, for each load step */
    for( sec=0 ; sec < no_section; ++sec )
       FiberLoadFlag[ii][sec] = 0;
}

/*
 * ========================================================================
 * SetUpFiberRespondBuffer() : Setup the static flags, stresses and strains
 *                             to store load history. 
 *
 * Input:  -- int total_fiber_elmt 
 *         -- EFRAME *frp 
 * Output: -- void.
 * ========================================================================
 */

#ifdef  __STDC__
void SetUpFiberRespondBuffer( int total_fiber_elmt, EFRAME *frp )
#else
void SetUpFiberRespondBuffer( total_fiber_elmt, frp )
int  total_fiber_elmt;
EFRAME *frp;
#endif
{
HISTORY_DATA    *hp;
int      ii, jj, kk;
int       ifib, sec;
int         elmt_no;
int      no_section;
int        no_fiber;
int        no_shear;
int  total_no_fiber;
int    elmt_attr_no;
int    UNITS_SWITCH;
QUANTITY         Es;
DIMENSIONS   *dimen;

    UNITS_SWITCH = CheckUnits();
    if( UNITS_SWITCH == ON )
       SetUnitsOff();

    FiberRespondBuffer = (FIBER_RESPONSE *)MyCalloc(1, sizeof(FIBER_RESPONSE));
    FiberRespondBuffer->total_fiber_elmt = total_fiber_elmt;
    FiberRespondBuffer->history = (HISTORY_DATA *) MyCalloc( total_fiber_elmt,
                                                             sizeof(HISTORY_DATA));

    no_section = frp->no_integ_pt + 2;  /* include 2 end sections */

    /* Initial flag matrix, for each load step */

    FiberLoadFlag = (int **)MyCalloc( total_fiber_elmt, sizeof(int *) );
    for( ii=0 ; ii < total_fiber_elmt; ++ii ) {
       FiberLoadFlag[ii] = (int *)MyCalloc( no_section, sizeof(int) );
       for( sec=0 ; sec < no_section; ++sec )
          FiberLoadFlag[ii][sec] = 0;
    }

    jj = 0;
    for( ii=1 ; ii <= frp->no_elements; ++ii ) {

       elmt_attr_no = frp->element[ii-1].elmt_attr_no;
       if( !(strcmp(frp->eattr[elmt_attr_no-1].elmt_type, "FIBER_2D"))
        || !(strcmp(frp->eattr[elmt_attr_no-1].elmt_type, "FIBER_3D"))
        || !(strcmp(frp->eattr[elmt_attr_no-1].elmt_type, "FIBER_2DS"))
        || !(strcmp(frp->eattr[elmt_attr_no-1].elmt_type, "FIBER_3DS")) )
       { 
	  /* total no. of fiber includes shear fibers */
	  no_fiber = frp->eattr[elmt_attr_no-1].work_fiber->no_fiber;
	  no_shear = frp->eattr[elmt_attr_no-1].work_fiber->no_shear;
	  total_no_fiber = no_fiber + (frp->no_dimen-1)*no_shear;
          elmt_no = ii;
          jj++;
          hp = &FiberRespondBuffer->history[jj-1];
          hp->elmt_no = elmt_no;
          hp->sr = (MATRIX *) MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                  no_section, total_no_fiber);
          hp->er = (MATRIX *) MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                  no_section, total_no_fiber);
          hp->s0 = (MATRIX *) MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                  no_section, total_no_fiber);
          hp->e0 = (MATRIX *) MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                  no_section, total_no_fiber);
          hp->stress = (MATRIX *) MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                       no_section, total_no_fiber);
          hp->strain = (MATRIX *) MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                      no_section, total_no_fiber);
          hp->tangent = (MATRIX *)MatrixAllocIndirect((char *)NULL, DOUBLE_ARRAY,
                                                      no_section, total_no_fiber);

          /* Initial the fiber tangent value into kx */

          for( ifib=0 ; ifib < no_fiber ; ++ifib ) {
             Es = frp->eattr[elmt_attr_no-1].work_fiber->fiber[ifib].Es;
             for( sec=0 ; sec < no_section ; ++sec )
                hp->tangent->uMatrix.daa[sec][ifib] = Es.value;
          }

          if( !(strcmp(frp->eattr[elmt_attr_no-1].elmt_type, "FIBER_2D")) 
           || !(strcmp(frp->eattr[elmt_attr_no-1].elmt_type, "FIBER_3D")) ) {

             for( ifib=no_fiber ; ifib < total_no_fiber ; ++ifib ) {
                Es = frp->eattr[elmt_attr_no-1].work_fiber->fiber[ifib].Es;
                for( sec=0 ; sec < no_section ; ++sec )
                   hp->tangent->uMatrix.daa[sec][ifib] = Es.value;
             }

          } else {

             for( ifib=no_fiber ; ifib < (no_fiber+no_shear) ; ++ifib ) {
                Es = frp->eattr[elmt_attr_no-1].work_fiber->fiber[ifib-no_fiber].Gs;
                for( sec=0 ; sec < no_section ; ++sec )
                   hp->tangent->uMatrix.daa[sec][ifib] = Es.value;
             }
             for( ifib=(no_fiber+no_shear) ; ifib < total_no_fiber ; ++ifib ) {
                Es = frp->eattr[elmt_attr_no-1].work_fiber->fiber[ifib-no_fiber-no_shear].Gs;
                for( sec=0 ; sec < no_section ; ++sec )
                   hp->tangent->uMatrix.daa[sec][ifib] = Es.value;
             }

          }

          hp->yielding  = (int **)MyCalloc( no_section, sizeof(int *) );
          hp->pre_load  = (int **)MyCalloc( no_section, sizeof(int *) );
          hp->pre_range = (int **)MyCalloc( no_section, sizeof(int *) );
          hp->loading   = (int **)MyCalloc( no_section, sizeof(int *) );

          for( kk=0 ; kk < no_section ; ++kk ) {

	     hp->yielding[kk]  = (int *)MyCalloc( total_no_fiber, sizeof(int) );
	     hp->pre_load[kk]  = (int *)MyCalloc( total_no_fiber, sizeof(int) );
	     hp->pre_range[kk] = (int *)MyCalloc( total_no_fiber, sizeof(int) );
	     hp->loading[kk]   = (int *)MyCalloc( total_no_fiber, sizeof(int) );

          }
       } /* end of if loop for fiber element */
    } /* end of for loop for all element */

    if( UNITS_SWITCH == ON )
       SetUnitsOn();
}

/*
 * =======================================================================
 * FiberElmtHistory() : Fiber element history .....
 *
 * Input:  -- int elmt_no :
 * Output: -- void.
 * =======================================================================
 */

#ifdef  __STDC__
HISTORY_DATA *FiberElmtHistory( int elmt_no )
#else
HISTORY_DATA *FiberElmtHistory( elmt_no )
int elmt_no;
#endif
{
int  ii;
HISTORY_DATA *hp;

    for( ii=0 ; ii < FiberRespondBuffer->total_fiber_elmt ; ++ii ) {
       if( elmt_no == FiberRespondBuffer->history[ii].elmt_no )
          break;
    }

    hp = &FiberRespondBuffer->history[ii];

    return( hp );
}
