/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_library.c : Utility Functions for Element Library
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
#include <string.h>
#include <stdlib.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "elmt.h"

#define Streq(s1, s2) (strcmp(s1, s2) == 0)

extern ARRAY     *array;
extern EFRAME    *frame;

#ifdef __STDC__
ARRAY *elmlib(ARRAY *p, int isw)
#else
ARRAY *elmlib(p, isw)
ARRAY    *p;
int     isw;
#endif
{
int i,l,m,k;
static int FLAG;

#ifdef DEBUG
       printf(" Enter elmtlib() \n");
#endif

    /* For selected element tasks, zero stiffness matrix */

    switch ( isw ) {
        case STIFF:
        case STRESS:
        case MASS_MATRIX:
        case LOAD_MATRIX:
             k = p->size_of_stiff;
             if(isw >= STIFF) { 
                for(i=1; i <= k; i++) {
                    for(m=1; m <= k; m++)
                        p->stiff->uMatrix.daa[i-1][m-1] = 0.0;
                }
                if( CheckUnits()==ON ) {
                    for(i = 1; i <= k; i++) {
                        ZeroUnits( &(p->stiff->spRowUnits[i-1]) );
                        ZeroUnits( &(p->stiff->spColUnits[i-1]) );
                        p->stiff->spRowUnits[i-1].units_type = SI_US;
                        p->stiff->spColUnits[i-1].units_type = SI_US;
                    }
                }
             }
             break;
        default:
             break;
    }

    /* Switch to element tasks() */

    FLAG = FALSE;
    for(i = 0; i < NO_ELEMENTS_IN_LIBRARY; i++) { 
       if(elmt_library[i].name != NULL && Streq(p->elmt_type, elmt_library[i].name)) {
          p =  (*(ARRAY * (*) ()) (elmt_library[i].elmt_lib_func))(p, isw);
          FLAG = TRUE;
          break;
       }
    }

    if(FLAG == FALSE) {
      printf("FATAL ERROR >> In elmtlib(): \n");
      printf("FATAL ERROR >> elmt_type = %s is not defined in element library: elmt.h \n", p->elmt_type);
      exit(1);
    }


#ifdef DEBUG
       printf(" Leaving elmtlib() \n");
#endif

   return(p);
}

/* ============================================================== */
/*   Assign Properties                                            */
/*   Print material and element information                       */
/* ============================================================== */

#ifdef __STDC__
EFRAME *assign_properties(EFRAME *frame, ELEMENT_ATTR *eap, MATERIAL_ATTR *map,
                          SECTION_ATTR *sap, FIBER_ELMT *fep, int elmt_no, int n)
#else
EFRAME *assign_properties(frame, eap, map, sap, fep, elmt_no, n)
EFRAME       *frame;
ELEMENT_ATTR   *eap;
MATERIAL_ATTR  *map;
SECTION_ATTR   *sap;
FIBER_ELMT     *fep;
int         elmt_no;
int               n;                     /* elmt attr no. */
#endif
{
int        no_dof, i, j;
int              length;
int        UNITS_SWITCH;
static int         FLAG;

SYMBOL             *slp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;
int         total_fiber;

#ifdef DEBUG
       printf(" enter assign_properties\n");
#endif

      UNITS_SWITCH = CheckUnits();

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

         /* material property */

         if(map->LC_ptr->name != (char *)NULL) {
            frame->element[elmt_no-1].LC_ptr->name  = SaveString(map->LC_ptr->name);
            frame->element[elmt_no-1].LC_ptr->alpha = map->LC_ptr->alpha;
            frame->element[elmt_no-1].LC_ptr->n     = map->LC_ptr->n;
            frame->element[elmt_no-1].LC_ptr->beta  = map->LC_ptr->beta;
         } else {
            frame->element[elmt_no-1].LC_ptr->name  = (char *)NULL;
            frame->element[elmt_no-1].LC_ptr->alpha = 0.0;
            frame->element[elmt_no-1].LC_ptr->n     = 0.0;
            frame->element[elmt_no-1].LC_ptr->beta  = 0.0;
         }
         frame->element[elmt_no-1].LC_ptr->ialph   = map->LC_ptr->ialph;
         frame->element[elmt_no-1].LC_ptr->pen     = map->LC_ptr->pen;
         frame->element[elmt_no-1].LC_ptr->load[0] = map->LC_ptr->load[0];
         frame->element[elmt_no-1].LC_ptr->load[1] = map->LC_ptr->load[1];
         frame->element[elmt_no-1].LC_ptr->load[2] = map->LC_ptr->load[2];
         frame->element[elmt_no-1].LC_ptr->load[3] = map->LC_ptr->load[3];
         frame->element[elmt_no-1].LC_ptr->load[4] = map->LC_ptr->load[4];
         frame->element[elmt_no-1].LC_ptr->load[5] = map->LC_ptr->load[5];

         for(j = 1; j <= iNO_INTEG_pts; j++) {
             frame->element[elmt_no-1].LC_ptr->H[j-1] = map->LC_ptr->H[j-1];
             frame->element[elmt_no-1].LC_ptr->R[j-1] = map->LC_ptr->R[j-1];
             for(i = 1; i <= 6; i++) {
                 frame->element[elmt_no-1].LC_ptr->back_stress[i-1][j-1]
                 = map->LC_ptr->back_stress[i-1][j-1];
             }
         }

         frame->eattr[n-1].work_material[0].value  = map->E.value;
         frame->eattr[n-1].work_material[1].value  = map->G.value;
         frame->eattr[n-1].work_material[2].value  = map->fy.value;
         frame->eattr[n-1].work_material[3].value  = map->ET.value;
         frame->eattr[n-1].work_material[4].value  = map->nu;
         frame->eattr[n-1].work_material[5].value  = map->density.value;
         frame->eattr[n-1].work_material[6].value  = map->fu.value;
         frame->eattr[n-1].work_material[7].value  = map->alpha_thermal[0].value;
         frame->eattr[n-1].work_material[8].value  = map->alpha_thermal[1].value;
         frame->eattr[n-1].work_material[9].value  = map->alpha_thermal[2].value;
         frame->eattr[n-1].work_material[10].value  = map->Gt.value;
         frame->eattr[n-1].work_material[11].value  = map->fv.value;

    switch(UNITS_SWITCH) {
      case ON:
        if(map->E.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[0].dimen, map->E.dimen );
        if(map->G.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[1].dimen, map->G.dimen );
        if(map->fy.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[2].dimen, map->fy.dimen );
        if(map->ET.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[3].dimen, map->ET.dimen );
        ZeroUnits( frame->eattr[n-1].work_material[4].dimen );
        if(map->density.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[5].dimen, map->density.dimen );
        if(map->fu.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[6].dimen, map->fu.dimen );
        if(map->alpha_thermal[0].dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[7].dimen, map->alpha_thermal[0].dimen );
        if(map->alpha_thermal[1].dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[8].dimen, map->alpha_thermal[1].dimen );
        if(map->alpha_thermal[2].dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[9].dimen, map->alpha_thermal[2].dimen );
        if(map->Gt.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[10].dimen, map->Gt.dimen );
        if(map->fv.dimen != NULL)
           UnitsCopy( frame->eattr[n-1].work_material[11].dimen, map->fv.dimen );

      break;
      case OFF:
      break;
      default:
      break;
    }
         /* section property  */

         frame->eattr[n-1].work_section[0].value  = sap->Ixx.value;
         frame->eattr[n-1].work_section[1].value  = sap->Iyy.value;
         frame->eattr[n-1].work_section[2].value  = sap->Izz.value;
         frame->eattr[n-1].work_section[3].value  = sap->Ixy.value;
         frame->eattr[n-1].work_section[4].value  = sap->Ixz.value;
         frame->eattr[n-1].work_section[5].value  = sap->Iyz.value;
         frame->eattr[n-1].work_section[6].value  = sap->weight.value;
         frame->eattr[n-1].work_section[7].value  = sap->bf.value;
         frame->eattr[n-1].work_section[8].value  = sap->tf.value;
         frame->eattr[n-1].work_section[9].value  = sap->depth.value;
         if(sap->area.value != 0.0)
            frame->eattr[n-1].work_section[10].value = sap->area.value;
         else
            if( sap->bf.value != 0 )
                frame->eattr[n-1].work_section[10].value  = sap->bf.value*sap->depth.value;
            else
                frame->eattr[n-1].work_section[10].value  = sap->width.value*sap->depth.value;
         frame->eattr[n-1].work_section[11].value  = sap->plate_thickness.value;
         frame->eattr[n-1].work_section[12].value  = sap->tor_const.value;
         frame->eattr[n-1].work_section[13].value  = sap->rT.value;
         frame->eattr[n-1].work_section[14].value  = sap->width.value;
         frame->eattr[n-1].work_section[15].value  = sap->tw.value;
         frame->eattr[n-1].work_section[16].value  = sap->ks;

    switch(UNITS_SWITCH) {
      case ON:
      if(sap->Ixx.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[0].dimen, sap->Ixx.dimen );
      if(sap->Iyy.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[1].dimen, sap->Iyy.dimen );
      if(sap->Izz.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[2].dimen, sap->Izz.dimen );
      if(sap->Ixy.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[3].dimen, sap->Ixy.dimen );
      if(sap->Ixz.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[4].dimen, sap->Ixz.dimen );
      if(sap->Iyz.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[5].dimen, sap->Iyz.dimen );
      if(sap->weight.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[6].dimen, sap->weight.dimen );
      if(sap->bf.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[7].dimen, sap->bf.dimen );
      if(sap->tf.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[8].dimen, sap->tf.dimen );
      if(sap->depth.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[9].dimen, sap->depth.dimen );

      if(sap->area.value != 0.0) {
         if(sap->area.dimen != NULL)
            UnitsCopy( frame->eattr[n-1].work_section[10].dimen, sap->area.dimen );
      }
      else {
        if(sap->depth.dimen != NULL && sap->bf.dimen !=NULL)
           UnitsMultRep( frame->eattr[n-1].work_section[10].dimen,
                         sap->bf.dimen , sap->depth.dimen );
        else if(sap->depth.dimen != NULL && sap->width.dimen !=NULL)
           UnitsMultRep( frame->eattr[n-1].work_section[10].dimen,
                         sap->width.dimen , sap->depth.dimen );
      }

      if(sap->plate_thickness.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[11].dimen, sap->plate_thickness.dimen );
      if(sap->tor_const.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[12].dimen, sap->tor_const.dimen );
      if(sap->rT.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[13].dimen, sap->rT.dimen );
      if(sap->width.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[14].dimen, sap->width.dimen );
      if(sap->tw.dimen != NULL)
         UnitsCopy( frame->eattr[n-1].work_section[15].dimen, sap->tw.dimen );
      ZeroUnits( frame->eattr[n-1].work_section[16].dimen );

      break;
      case OFF:
      break;
      default:
      break;
    }

    /* store hash table data to frame */
    if( !(strcmp(eap->elmt_type, "FIBER_2D"))  || !(strcmp(eap->elmt_type, "FIBER_3D")) 
    ||  !(strcmp(eap->elmt_type, "FIBER_2DS")) || !(strcmp(eap->elmt_type, "FIBER_3DS")) ) {

       if( !(strcmp(eap->elmt_type, "FIBER_2D")) || !(strcmp(eap->elmt_type, "FIBER_3D")) ) {
          fep->no_shear = 1;  /* global shear spring for each element in one direction */
          total_fiber = fep->no_fiber + (frame->no_dimen-1)*fep->no_shear;

          frame->eattr[n-1].work_fiber = (FIBER_ELMT *)MyCalloc(1,sizeof(FIBER_ELMT));
          frame->eattr[n-1].work_fiber->no_fiber = fep->no_fiber;
          frame->eattr[n-1].work_fiber->no_shear = fep->no_shear;
          frame->eattr[n-1].work_fiber->fiber = 
             (FIBER_ATTR *)MyCalloc(total_fiber, sizeof(FIBER_ATTR));
          for( i=0 ; i < total_fiber ; ++i ) {
             if( i < fep->no_fiber ) {  /* flexual fibers */
                frame->eattr[n-1].work_fiber->fiber[i].y.value = fep->fiber[i].y.value;
                frame->eattr[n-1].work_fiber->fiber[i].z.value = fep->fiber[i].z.value;
                frame->eattr[n-1].work_fiber->fiber[i].area.value = fep->fiber[i].area.value;
                frame->eattr[n-1].work_fiber->fiber[i].Es.value = fep->fiber[i].Es.value;
                frame->eattr[n-1].work_fiber->fiber[i].Et.value = fep->fiber[i].Et.value;
                frame->eattr[n-1].work_fiber->fiber[i].fy.value = fep->fiber[i].fy.value;
                frame->eattr[n-1].work_fiber->fiber[i].Gs.value = 0.0;
                frame->eattr[n-1].work_fiber->fiber[i].Gt.value = 0.0;
                frame->eattr[n-1].work_fiber->fiber[i].fv.value = 0.0;
             }
             else {  /* shear spring */
                frame->eattr[n-1].work_fiber->fiber[i].y.value  = 0.0;
                frame->eattr[n-1].work_fiber->fiber[i].z.value  = 0.0;
                frame->eattr[n-1].work_fiber->fiber[i].Gs.value = 0.0;
                frame->eattr[n-1].work_fiber->fiber[i].Gt.value = 0.0;
                frame->eattr[n-1].work_fiber->fiber[i].fv.value = 0.0;

                if( sap->ks <= 0.0 )   /* default shear correction factor ks = 1.2 */
                   frame->eattr[n-1].work_fiber->fiber[i].area.value
                    = frame->eattr[n-1].work_section[10].value / 1.2;
                else
                   frame->eattr[n-1].work_fiber->fiber[i].area.value
                    = frame->eattr[n-1].work_section[10].value / sap->ks;

                if( map->nu <= 0.0 )   /* default nu = 0.3 */
                   frame->eattr[n-1].work_material[4].value  = 0.3;

                if( map->G.value <= 0.0 )   /* G = E/2/(1+nu) */
                   if( map->E.value <= 0.0 )
                      frame->eattr[n-1].work_fiber->fiber[i].Es.value
                       = frame->eattr[n-1].work_fiber->fiber[0].Es.value/2.0/(1.0+frame->eattr[n-1].work_material[4].value);
                   else
                      frame->eattr[n-1].work_fiber->fiber[i].Es.value
                       = map->E.value/2.0/(1.0+frame->eattr[n-1].work_material[4].value);
                else
                   frame->eattr[n-1].work_fiber->fiber[i].Es.value = map->G.value;

                if( map->Gt.value <= 0.0 ) {
                   if( map->ET.value <= 0.0 )   /* Et not exist, assume always linear, Et=Es */
                      frame->eattr[n-1].work_fiber->fiber[i].Et.value
                       = frame->eattr[n-1].work_fiber->fiber[i].Es.value;
                   else
                      frame->eattr[n-1].work_fiber->fiber[i].Et.value
                       = map->ET.value/2.0/(1.0+frame->eattr[n-1].work_material[4].value);
                }
                else
                   frame->eattr[n-1].work_fiber->fiber[i].Et.value = map->Gt.value;

                if( map->fv.value <= 0.0 )   /* use fy as shear yielding stress */
                   if( map->fy.value <= 0.0 )
                      frame->eattr[n-1].work_fiber->fiber[i].fy.value
                       = frame->eattr[n-1].work_fiber->fiber[0].fy.value;
                   else
                      frame->eattr[n-1].work_fiber->fiber[i].fy.value = map->fy.value;
                else
                   frame->eattr[n-1].work_fiber->fiber[i].fy.value = map->fv.value;
             }

             if(UNITS_SWITCH == ON) {
                frame->eattr[n-1].work_fiber->fiber[i].y.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].z.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].area.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Es.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Et.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].fy.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Gs.dimen = (DIMENSIONS *) NULL;
                frame->eattr[n-1].work_fiber->fiber[i].Gt.dimen = (DIMENSIONS *) NULL;
                frame->eattr[n-1].work_fiber->fiber[i].fv.dimen = (DIMENSIONS *) NULL;

                if( i < fep->no_fiber ) {  /* flexual fibers */
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].y.dimen, fep->fiber[i].y.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].z.dimen, fep->fiber[i].z.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].area.dimen, fep->fiber[i].area.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Es.dimen, fep->fiber[i].Es.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Et.dimen, fep->fiber[i].Et.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].fy.dimen, fep->fiber[i].fy.dimen );
                }
                else {  /* shear spring */
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].y.dimen, fep->fiber[0].y.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].z.dimen, fep->fiber[0].z.dimen );
                   UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].area.dimen, frame->eattr[n-1].work_section[10].dimen );

                   if( map->G.value <= 0.0 )
                      if( map->E.value <= 0.0 )
                         UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Es.dimen,
                                    frame->eattr[n-1].work_fiber->fiber[0].Es.dimen );
                      else
                         UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Es.dimen, map->E.dimen );
                   else
                      UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Es.dimen, map->G.dimen );

                   if( map->Gt.value <= 0.0 ) {
                      if( map->ET.value <= 0.0 )
                         UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Et.dimen,
                                    frame->eattr[n-1].work_fiber->fiber[i].Es.dimen );
                      else
                         UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Et.dimen, map->ET.dimen );
                   }
                   else
                      UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Et.dimen, map->Gt.dimen );

                   if( map->fv.value <= 0.0 )
                      if( map->fy.value <= 0.0 )
                         UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].fy.dimen,
                                    frame->eattr[n-1].work_fiber->fiber[0].fy.dimen );
                      else
                         UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].fy.dimen, map->fy.dimen );
                   else
                      UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].fy.dimen, map->fv.dimen );
                }
             } /* UNIT_SWITCH == ON */
          } /* i=0, i<total_fiber */
       } /* FIBER_2D and FIBER_3D */

       else {  /* FIBER_2DS and FIBER_3DS */
          fep->no_shear = fep->no_fiber;  /* individual shear spring for each fiber in one direction */
          frame->eattr[n-1].work_fiber = (FIBER_ELMT *)MyCalloc(1,sizeof(FIBER_ELMT));
          frame->eattr[n-1].work_fiber->no_fiber = fep->no_fiber;
          frame->eattr[n-1].work_fiber->no_shear = fep->no_shear;
          frame->eattr[n-1].work_fiber->fiber = 
             (FIBER_ATTR *)MyCalloc(fep->no_fiber,sizeof(FIBER_ATTR));
          for( i=0 ; i < fep->no_fiber ; ++i ) {
             frame->eattr[n-1].work_fiber->fiber[i].y.value = fep->fiber[i].y.value;
             frame->eattr[n-1].work_fiber->fiber[i].z.value = fep->fiber[i].z.value;
             frame->eattr[n-1].work_fiber->fiber[i].area.value = fep->fiber[i].area.value;
             frame->eattr[n-1].work_fiber->fiber[i].Es.value = fep->fiber[i].Es.value;
             frame->eattr[n-1].work_fiber->fiber[i].Et.value = fep->fiber[i].Et.value;
             frame->eattr[n-1].work_fiber->fiber[i].fy.value = fep->fiber[i].fy.value;
             frame->eattr[n-1].work_fiber->fiber[i].Gs.value = fep->fiber[i].Gs.value;
             frame->eattr[n-1].work_fiber->fiber[i].Gt.value = fep->fiber[i].Gt.value;
             frame->eattr[n-1].work_fiber->fiber[i].fv.value = fep->fiber[i].fv.value;

             if(UNITS_SWITCH == ON) {
                frame->eattr[n-1].work_fiber->fiber[i].y.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].z.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].area.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Es.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Et.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].fy.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Gs.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].Gt.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );
                frame->eattr[n-1].work_fiber->fiber[i].fv.dimen
                 = (DIMENSIONS *)MyCalloc( 1, sizeof(DIMENSIONS) );

                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].y.dimen, fep->fiber[i].y.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].z.dimen, fep->fiber[i].z.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].area.dimen, fep->fiber[i].area.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Es.dimen, fep->fiber[i].Es.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Et.dimen, fep->fiber[i].Et.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].fy.dimen, fep->fiber[i].fy.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Gs.dimen, fep->fiber[i].Gs.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].Gt.dimen, fep->fiber[i].Gt.dimen );
                UnitsCopy( frame->eattr[n-1].work_fiber->fiber[i].fv.dimen, fep->fiber[i].fv.dimen );
             } /* UNIT_SWITCH == ON */
          } /* i=0, i<fep->no_fiber */
       } /* FIBER_2DS and FIBER_3DS */
    } /* eap->elmt_type == FIBER_** */

     /* DOF of nodes in a elment */

      FLAG = FALSE;
      for (i = 0; i < NO_ELEMENTS_IN_LIBRARY ; i++) {
          if((elmt_library[i].name) != NULL &&
             Streq(frame->eattr[n-1].elmt_type, elmt_library[i].name)) {
             no_dof   = elmt_library[i].no_dof;
             FLAG = TRUE;
             break;
          }
      }

      if(FLAG == FALSE) {
          printf(">>FATAL ERROR:  In assign_properties():\n");
          printf(">>           :  In file elmt.h : elmemt type %s is not defined \n",
                                  frame->eattr[n-1].elmt_type);
          exit(1);
      }

      for(i = 1; i <= no_dof; i++) {
        frame->eattr[n-1].map_ldof_to_gdof[i-1] = eap->map_ldof_to_gdof[i-1];
      }

    return(frame);

#ifdef DEBUG
       printf("*** Leave assign_properties()\n");
#endif
}
