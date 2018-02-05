/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt_set_attr.c : Set Element Attributes
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
#include "fe_functions.h"
#include "elmt.h"

enum { ATTR = 20 };  /* Array of element attributes */

#ifdef __STDC__
EFRAME *Set_Elmt_Attrs(EFRAME *frame)
#else
EFRAME *Set_Elmt_Attrs(frame)
EFRAME *frame;
#endif
{

ELEMENT_ATTR   *eap;
SECTION_ATTR   *sap;
MATERIAL_ATTR  *map;
FIBER_ELMT     *fep;

char      *eap_name;
char    *attr_array[ ATTR ];
char      *new_attr;
int     i, j = 0, k;
int   flag, elmt_attr_no;
int     total_fiber;
int    UNITS_SWITCH;

int     no_integ_pt;
SYMBOL         *slp;

#ifdef DEBUG
       printf("*** Enter Set_Elmt_Attrs() : frame->no_elements = %5d\n", frame->no_elements);
#endif

    UNITS_SWITCH = CheckUnits();

    /* =======================================================*/
    /*  Get Element Attr : section and material properties    */ 
    /*                     for all elements from Hash Nodes   */
    /* =======================================================*/

    for( k=1; k <= ATTR ; k++ )
        attr_array[k-1] = (char *) NULL;

    no_integ_pt = frame->no_integ_pt;
    for(i=1; i<=frame->no_elements; i++) {

       /* [a] : Lookup Element, Section, and Material Attributes */ 

       eap_name = frame->element[i-1].elmt_attr_name;
       eap      = lookup(eap_name)->u.eap;
       if(eap == NULL) {
          FatalError("Elmt_Attribute name not found",(char *)NULL);
       }

#ifdef DEBUG
       printf(" In Set_Elmt_Attrs() : eap->section = %s \n", eap->section);
#endif

       sap = lookup(eap->section)->u.sap;
       if(sap == NULL) {
          FatalError("Section_Attribute name not found",(char *)NULL);
       }

       map = lookup(eap->material)->u.map;
       if(map == NULL) {
          FatalError("Material_Attribute name not found",(char *)NULL);
       }

       fep = (FIBER_ELMT *)NULL;
       if( !(strcmp(eap->elmt_type, "FIBER_2D"))  || !(strcmp(eap->elmt_type, "FIBER_3D"))
       ||  !(strcmp(eap->elmt_type, "FIBER_2DS")) || !(strcmp(eap->elmt_type, "FIBER_3DS")) ) {
          fep = lookup(eap->fiber_attr_name)->u.fep;
          if(fep == NULL) {
              FatalError("Fiber_Elmt_Attr name not found",(char *)NULL);
          }
       }

       /* [b] Check for definition of a new element attribute type */ 

       if (i==1) {
           j++;
           attr_array[0] = SaveString(eap_name);
           elmt_attr_no = j;
           frame->no_element_attr = j;
           frame->element[i-1].elmt_attr_no = elmt_attr_no;
           frame->eattr[elmt_attr_no-1].name = SaveString(eap_name);
           frame->eattr[elmt_attr_no-1].elmt_type = SaveString(eap->elmt_type);
           frame->eattr[elmt_attr_no-1].material  = SaveString(eap->material);
           frame->eattr[elmt_attr_no-1].section   = SaveString(eap->section);
           frame->eattr[elmt_attr_no-1].fiber_attr_name = SaveString(eap->fiber_attr_name);
       }
       else {
           k = 1;
           flag = 0;
           while( k<=ATTR && attr_array[k-1]!=(char *)NULL ) {
               if(strcmp(eap_name,attr_array[k-1])!=0)
                   k++;
               else {
                   flag++;
                   elmt_attr_no = k;
                   frame->element[i-1].elmt_attr_no = k;
                   break;
               }
           }
           if( flag==0 ) {
               j++;
               attr_array[j-1] = SaveString(eap_name);
               elmt_attr_no = j;
               frame->no_element_attr = j;
               frame->element[i-1].elmt_attr_no = elmt_attr_no;
               frame->eattr[elmt_attr_no-1].name = SaveString(eap_name);
               frame->eattr[elmt_attr_no-1].elmt_type = SaveString(eap->elmt_type);
               frame->eattr[elmt_attr_no-1].material  = SaveString(eap->material);
               frame->eattr[elmt_attr_no-1].section   = SaveString(eap->section);
               frame->eattr[elmt_attr_no-1].fiber_attr_name = SaveString(eap->fiber_attr_name);
           }
       }

       /* [c] Assign Section/Material properties to Element Array */ 

       frame = assign_properties(frame, eap, map, sap, fep,i, elmt_attr_no);

       /* allocate space for response and element state in the fiber element */

       if( !(strcmp(eap->elmt_type, "FIBER_2D"))  || !(strcmp(eap->elmt_type, "FIBER_3D")) 
       ||  !(strcmp(eap->elmt_type, "FIBER_2DS")) || !(strcmp(eap->elmt_type, "FIBER_3DS")) ) {
         if( UNITS_SWITCH == ON )
            SetUnitsOff();
         if( !(strcmp(eap->elmt_type, "FIBER_2D")) || !(strcmp(eap->elmt_type, "FIBER_2DS")) ) {
	   frame->element[i-1].rp->Q_saved = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 3, 1);
	   frame->element[i-1].rp->q_saved = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 3, 1);
         }
         else {
	   frame->element[i-1].rp->Q_saved = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 5, 1);
	   frame->element[i-1].rp->q_saved = MatrixAllocIndirect( (char *)NULL, DOUBLE_ARRAY, 5, 1);
         }

         total_fiber = fep->no_fiber + (frame->no_dimen-1)*fep->no_shear;
	 frame->element[i-1].rp->sr_saved = MatrixAllocIndirect( (char *)NULL,
	        DOUBLE_ARRAY, no_integ_pt+2, total_fiber );
	 frame->element[i-1].rp->er_saved = MatrixAllocIndirect( (char *)NULL,
	        DOUBLE_ARRAY, no_integ_pt+2, total_fiber );
	 frame->element[i-1].rp->s0_saved = MatrixAllocIndirect( (char *)NULL,
	        DOUBLE_ARRAY, no_integ_pt+2, total_fiber );
	 frame->element[i-1].rp->e0_saved = MatrixAllocIndirect( (char *)NULL,
	        DOUBLE_ARRAY, no_integ_pt+2, total_fiber );
	 frame->element[i-1].rp->sx_saved = MatrixAllocIndirect( (char *)NULL,
	        DOUBLE_ARRAY, no_integ_pt+2, total_fiber );
	 frame->element[i-1].rp->ex_saved = MatrixAllocIndirect( (char *)NULL,
	        DOUBLE_ARRAY, no_integ_pt+2, total_fiber );

	 frame->element[i-1].esp->yielding_saved =
		(int **)MyCalloc( no_integ_pt+2, sizeof(int *) );
	 frame->element[i-1].esp->pre_range_saved =
		(int **)MyCalloc( no_integ_pt+2, sizeof(int *) );
	 frame->element[i-1].esp->pre_load_saved =
		(int **)MyCalloc( no_integ_pt+2, sizeof(int *) );
	 for( k=1 ; k <= no_integ_pt+2 ; ++k ) {
	   frame->element[i-1].esp->yielding_saved[k-1] =
		  (int *)MyCalloc( total_fiber, sizeof(int) );
	   frame->element[i-1].esp->pre_range_saved[k-1] =
		  (int *)MyCalloc( total_fiber, sizeof(int) );
	   frame->element[i-1].esp->pre_load_saved[k-1] =
		  (int *)MyCalloc( total_fiber, sizeof(int) );
	 }

         if( UNITS_SWITCH == ON )
            SetUnitsOn();
       } /* allocate response and state for fiber element */
    }

    k = 1;
    while( k<=ATTR && attr_array[k-1]!=(char *)NULL ) {
        free((char *) attr_array[k-1] );
        k++;
    }

   /* 
    *  =========================================================
    *  [b] : Get Properties of Rigid Bodies
    *  =========================================================
    */ 

    for(i=1; i<=frame->no_rigid; i++) {

        eap_name = frame->rigid[i-1].rbody_attr_name;
        sap      = lookup(frame->rigid[i-1].rbody_attr_name)->u.sap;
        if(sap == NULL) {
           FatalError("Rigid Body attribute name not found",(char *)NULL);
        }

        if(sap->section_type != 0) 
           frame->rigid[i -1].rb_type = sap->section_type;

        frame->rigid[i -1].prop[1]  = map->density;
        frame->rigid[i -1].prop[2]  = sap->plate_thickness;
        frame->rigid[i -1].prop[3]  = sap->area;
        frame->rigid[i -1].prop[7]  = sap->Ixx;
        frame->rigid[i -1].prop[8]  = sap->Iyy;
        frame->rigid[i -1].prop[9]  = sap->Izz;
        frame->rigid[i -1].prop[10] = sap->Ixy;
        frame->rigid[i -1].prop[11] = sap->Ixz;
        frame->rigid[i -1].prop[12] = sap->Iyz;
    }

#ifdef DEBUG
       printf("*** Leave set_elmt_attrs() \n");
#endif

      return(frame);
}
