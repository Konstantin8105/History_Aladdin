/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  symbol_load.c : Load materials and AISC sections into symbol table
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
 *  Written by: Mark Austin                                          January 1993
 *  ============================================================================= 
 */

#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"

/*
 *  ------------------------------------
 *  Load AISC Sections into Symbol Table
 *  ------------------------------------
 */

void Load_AISC_Sections() { 
char *name;
FILE   *fp;
double weight, area, depth;
double tw, bf, tf, rT, Ix, Iy, Iz;
SYMBOL *sp;
char buffer[MAX_NO_CHAR_PER_LINE];
int i, lenghth1, lenght2;
DIMENSIONS  *dp_weight, *dp_length, *dimen1, *dimen2;
 
    /* Open Input File */
 
    if((fp = fopen("section.h", "r")) == NULL) {
        FatalError("In Load_AISC_Sections() : Can't open 'section.h' data file",
                  (char *) NULL);
    }
 
    /* Read the File Header */
 
    for(i = 1; i <= NO_LINES_IN_HEADER; i++)
        fgets(buffer,MAX_NO_CHAR_PER_LINE,fp);

    name = (char *) calloc(10,sizeof(char));
    if( CheckUnits() == ON ) {
        dp_length = DefaultUnits("in");
        dimen1    = DefaultUnits("lbf");
        dimen2    = DefaultUnits("ft");
        dp_weight = UnitsDiv( dimen1, dimen2, NO );
        free((char *)dimen1->units_name);
        free((char *)dimen1);
        free((char *)dimen2->units_name);
        free((char *)dimen2);
        dimen1 = UnitsPower( dp_length, 2.0, NO );
        dimen2 = UnitsPower( dp_length, 4.0, NO );
    }

    while((fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", name,
           &weight, &area, &depth, &tw, &bf, &tf, &rT, &Iz, &Iy,&Ix)) != EOF) {

        sp = install(name);
	sp->u.sap = (SECTION_ATTR *)Alloc_Section_Attr_Item();

         /* sp->u.sap->tor_const.value = J; */
	if( CheckUnits() == ON ) {
            sp->u.sap->weight.value  = weight*dp_weight->scale_factor;
            sp->u.sap->area.value  = area*dimen1->scale_factor;
            sp->u.sap->depth.value = depth*dp_length->scale_factor;
            sp->u.sap->tw.value  = tw*dp_length->scale_factor;
            sp->u.sap->bf.value  = bf*dp_length->scale_factor;
            sp->u.sap->tf.value  = tf*dp_length->scale_factor;
            sp->u.sap->rT.value  = rT*dp_length->scale_factor;

            sp->u.sap->Ixx.value = Ix*dimen2->scale_factor;
            sp->u.sap->Iyy.value = Iy*dimen2->scale_factor;
            sp->u.sap->Izz.value = Iz*dimen2->scale_factor;

            sp->u.sap->weight.dimen = dp_weight;
            sp->u.sap->area.dimen = dimen1;
            sp->u.sap->depth.dimen = dp_length;
            sp->u.sap->tw.dimen  = dp_length;
            sp->u.sap->bf.dimen  = dp_length;
            sp->u.sap->tf.dimen  = dp_length;
            sp->u.sap->rT.dimen  = dp_length;
            sp->u.sap->Ixx.dimen = dimen2;
            sp->u.sap->Iyy.dimen = dimen2;
            sp->u.sap->Izz.dimen = dimen2;
	}
        else {
            sp->u.sap->weight.value  = weight;
            sp->u.sap->area.value  = area;
            sp->u.sap->depth.value = depth;
            sp->u.sap->tw.value  = tw;
            sp->u.sap->bf.value  = bf;
            sp->u.sap->tf.value  = tf;
            sp->u.sap->rT.value  = rT;

            sp->u.sap->Ixx.value = Ix;
            sp->u.sap->Iyy.value = Iy;
            sp->u.sap->Izz.value = Iz;
        }
   }
 
   fclose(fp);
}


/*
 *  --------------------------------
 *  Load Materials into Symbol Table
 *  --------------------------------
 */
 
void Load_AISC_Material() {
FILE      *fp;
char    *name;
SYMBOL    *sp;
double fy, G, E,nu;
char buffer[MAX_NO_CHAR_PER_LINE];
int i;
DIMENSIONS *dimen;
 
    /* Open Input File */
 
    if((fp = fopen("material.h", "r")) == NULL) {
        FatalError("In Load_AISC_Sections() : Can't open 'material.h'",
                  (char *) NULL);
    }  
 
    /* Read the header */
 
    for(i = 1; i <= NO_LINES_IN_HEADER; i++)
        fgets( buffer, MAX_NO_CHAR_PER_LINE, fp);
 
    name = (char *) malloc(10*sizeof(char));
    if( CheckUnits() == ON )
        dimen = DefaultUnits("ksi");

    while((fscanf(fp, "%s %lf %lf %lf %lf", name, &fy, &G, &E, &nu)) != EOF) {
        sp = install(name);
        sp->u.map = (MATERIAL_ATTR *)Alloc_Material_Attr_Item();

	if( CheckUnits() == ON ) {
            sp->u.map->fy.value   = fy*dimen->scale_factor;
            sp->u.map->G.value    = G*dimen->scale_factor;
            sp->u.map->E.value    = E*dimen->scale_factor;
            sp->u.map->nu         = nu;

	    sp->u.map->fy.dimen = dimen;
	    sp->u.map->E.dimen  = dimen;
	    sp->u.map->G.dimen  = dimen;
	}
        else {
            sp->u.map->fy.value   = fy;
            sp->u.map->G.value    = G;
            sp->u.map->E.value    = E;
            sp->u.map->nu         = nu;
        }
    }  

    fclose(fp);
}
