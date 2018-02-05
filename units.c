/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  units.c : Functions to handle engineering units/dimensions
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
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef  __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"

static enum {UNITS_OFF, UNITS_ON} iFlag = UNITS_ON;
static int UNITS_TYPE;

/* #define DEBUG */

/*
 *  ===========================================================
 *  Utility Functions for Units Operations
 *  
 *  SetUnitsOn()     : Set flag for units "on"
 *  SetUnitsOff()    : Set flag for units "off"
 *  CheckUnits()     : Flag to see if units are "on" or "off"
 *  CheckUnitsType() : Check type of units
 *  ===========================================================
 */

int SetUnitsOn(){
    iFlag = UNITS_ON;  
}

int SetUnitsOff(){
    iFlag = UNITS_OFF;  
}

int CheckUnits(){
    return(iFlag);
}

int CheckUnitsType(){
    return(UNITS_TYPE);
}

#ifdef  __STDC__
int ChangeUnitsType( int index )
#else
int ChangeUnitsType( index )
int index;
#endif
{
    switch( index ) {
       case 1:
          UNITS_TYPE = US;
          break;
       case 2:
          UNITS_TYPE = SI;
          break;
       default:
          UNITS_TYPE = SI_US;
          break;
    }
}


/*
 *  ========================================
 *  Do two quantities have the same units ??
 *  ========================================
 */

#ifdef __STDC__
int SameUnits(DIMENSIONS *d1, DIMENSIONS *d2)
#else
int SameUnits(d1, d2)
DIMENSIONS *d1, *d2;
#endif
{
    if((d1->length_expnt  == d2->length_expnt) &&
       (d1->mass_expnt    == d2->mass_expnt)   &&
       (d1->time_expnt    == d2->time_expnt)   &&
       (d1->temp_expnt    == d2->temp_expnt))  {
         return TRUE;
    }
    else {
         return FALSE;
    }
}


/*
 *  ================================================
 *  Functions for Mutiply/Divide Operations on Units
 *  ================================================
 */

/*
 *  ===============================================================
 *  UnitsMult() : Multiply Dimensional Units.
 * 
 *  Input : DIMENSIONS *d1 -- pointer to dimensions data structure
 *          DIMENSIONS *d2 -- pointer to dimensions data structure
 *  Input : DIMENSIONS  *d -- product of dimensions data structures
 *  ===============================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsMult(DIMENSIONS *d1, DIMENSIONS *d2)
#else
DIMENSIONS *UnitsMult(d1, d2)
DIMENSIONS *d1, *d2;
#endif
{
DIMENSIONS         *d;
DIMENSIONS *dp1, *dp2;
int            length;

#ifdef DEBUG
       printf("Enter UnitsMult() \n");
       UnitsPrint(d1);
       UnitsPrint(d2);
#endif

    if( d1==(DIMENSIONS *)NULL || d2==(DIMENSIONS *)NULL ) {
        printf("Warning: d1 or d2 is NULL,   in UnitsMult()\n");
        return (DIMENSIONS *)NULL;
    }

    d   = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    dp1 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    dp2 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));

    UnitsCopy(dp1,d1);
    UnitsCopy(dp2,d2);

    if(dp1->units_type == dp2->units_type) {
       d->units_type = dp1->units_type;
    }
    else {
       if( dp1->units_type==SI && dp2->units_type==US ) {
           d->units_type = SI;
           if( dp2->units_name != (char *)NULL )
              UnitsTypeConvert( dp2, SI );
       }
       else if( dp1->units_type==US && dp2->units_type==SI ) {
           d->units_type = SI;
           if( dp1->units_name != (char *)NULL )
              UnitsTypeConvert( dp1, SI );
       }
       else if(dp1->units_type == SI_US)
           d->units_type = dp2->units_type;
       else if(dp2->units_type == SI_US)
           d->units_type = dp1->units_type;
    }
	
    if( dp1->units_name!=(char *)NULL ) {
       length = UnitsLength(dp1->units_name,dp2->units_name,"STOP");
       d->units_name = (char *)MyCalloc(length+2,sizeof(char));
       d->units_name = strcpy(d->units_name, dp1->units_name);
       if( dp2->units_name!=(char *)NULL ) {
          d->units_name = strcat(d->units_name, (char *) ".");
          d->units_name = strcat(d->units_name, dp2->units_name);
       }
    }
    else if( dp2->units_name!=(char *)NULL )
       d->units_name = SaveString(dp2->units_name);
    else
       d->units_name = (char *)NULL;

    d->scale_factor = dp1->scale_factor*dp2->scale_factor;
    d->length_expnt = dp1->length_expnt + dp2->length_expnt;
    d->mass_expnt   = dp1->mass_expnt + dp2->mass_expnt;
    d->time_expnt   = dp1->time_expnt + dp2->time_expnt;
    d->temp_expnt   = dp1->temp_expnt + dp2->temp_expnt;
    d->radian_expnt = dp1->radian_expnt + dp2->radian_expnt;

    free((char *)dp1->units_name);
    free((char *)dp1);
    free((char *)dp2->units_name);
    free((char *)dp2);

    return (d);
}

/*
 *  =====================================================================
 *  UnitsMultRep() : Units multiply with replacement
 * 
 *  Input :  DIMENSIONS *d1 -- pointer to dimensions data structure
 *           DIMENSIONS *d2 -- pointer to dimensions data structure
 *           DIMENSIONS  *d -- structure for result of dimensions product
 *  Output : DIMENSIONS  *d -- structure for result of dimensions product
 *  =====================================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsMultRep(DIMENSIONS *d, DIMENSIONS *d1, DIMENSIONS *d2)
#else
DIMENSIONS *UnitsMultRep(d, d1, d2)
DIMENSIONS       *d;
DIMENSIONS *d1, *d2;
#endif
{
DIMENSIONS *dp1, *dp2;
int            length;

#ifdef DEBUG
       printf("ENTER UnitsMultRep() \n");
       UnitsPrint(d);
#endif

    if( d1==(DIMENSIONS *)NULL || d2==(DIMENSIONS *)NULL ) {
        printf("Warning : Dimensions d1 and d2 cannot be NULL\n");
        free((char *)d->units_name);
        free((char *)d);
        d = (DIMENSIONS *) NULL;
        return (DIMENSIONS *) NULL;
    }

    dp1 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    dp2 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));

    UnitsCopy(dp1,d1);
    UnitsCopy(dp2,d2);

    if(dp1->units_type == dp2->units_type) {
       d->units_type = dp1->units_type;
    }
    else {
       if( dp1->units_type==SI && dp2->units_type==US ) {
           d->units_type = SI;
           if( dp2->units_name != (char *)NULL )
              UnitsTypeConvert( dp2, SI );
       }
       else if( dp1->units_type==US && dp2->units_type==SI ) {
           d->units_type = SI;
           if( dp1->units_name != (char *)NULL )
              UnitsTypeConvert( dp1, SI );
       }
       else if(dp1->units_type == SI_US)
           d->units_type = dp2->units_type;
       else if(dp2->units_type == SI_US)
           d->units_type = dp1->units_type;
    }
	
    if( dp1->units_name != (char *)NULL ) {
       free((char *)d->units_name);
       length = UnitsLength(dp1->units_name,dp2->units_name,"STOP");
       d->units_name = (char *) MyCalloc(length+2,sizeof(char));
       d->units_name = strcpy(d->units_name, dp1->units_name);
       if( dp2->units_name != (char *)NULL ) {
          d->units_name = strcat(d->units_name, (char *) ".");
          d->units_name = strcat(d->units_name, dp2->units_name);
       }
    }
    else if( dp2->units_name != (char *)NULL ) {
       free((char *) d->units_name);
       d->units_name = SaveString(dp2->units_name);
    }
    else {
       free((char *) d->units_name);
       d->units_name = (char *)NULL;
    }

    d->scale_factor = dp1->scale_factor*dp2->scale_factor;
    d->length_expnt = dp1->length_expnt + dp2->length_expnt;
    d->mass_expnt   = dp1->mass_expnt + dp2->mass_expnt;
    d->time_expnt   = dp1->time_expnt + dp2->time_expnt;
    d->temp_expnt   = dp1->temp_expnt + dp2->temp_expnt;
    d->radian_expnt = dp1->radian_expnt + dp2->radian_expnt;

    free((char *)dp1->units_name);
    free((char *)dp1);
    free((char *)dp2->units_name);
    free((char *)dp2);

#ifdef DEBUG
       printf("Leave UnitsMultRep() \n");
#endif

    return (d);
}

/*
 *  =====================================================================
 *  UnitsDiv() : Divide units
 * 
 *  Input :  DIMENSIONS *d1 -- pointer to dimensions data structure
 *           DIMENSIONS *d2 -- pointer to dimensions data structure
 *           int FLAG       -- YES or NO flag for "()" brackete.
 *  Output : DIMENSIONS  *d -- structure for result of dimensions product
 *  =====================================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsDiv( DIMENSIONS *d1, DIMENSIONS *d2, int FLAG )
#else
DIMENSIONS *UnitsDiv( d1, d2, FLAG )
DIMENSIONS *d1, *d2;
int            FLAG;           /* YES or NO flag : to determine wheather "()" is needed */
#endif
{
DIMENSIONS         *d;
DIMENSIONS *dp1, *dp2;
int            length;

#ifdef DEBUG
       printf("Enter UnitsDiv() \n");
       UnitsPrint(d1);
       UnitsPrint(d2);
#endif

    if( d2==(DIMENSIONS *)NULL )
        FatalError("In UnitsDiv() : unit divider is a NULL pointer",(char *)NULL);

    if( d1==(DIMENSIONS *)NULL ) {
        printf("Warning: d1 is NULL,   in UnitsDiv()\n");
        return (DIMENSIONS *)NULL;
    }

    d   = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    dp1 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    dp2 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));

    UnitsCopy(dp1,d1);
    UnitsCopy(dp2,d2);

    if(dp1->units_type == dp2->units_type) {
       d->units_type = dp1->units_type;
    }
    else {
       if( dp1->units_type==SI && dp2->units_type==US ) {
           d->units_type = SI;
           if( dp2->units_name != (char *)NULL )
              UnitsTypeConvert( dp2, SI );
       }
       else if( dp1->units_type==US && dp2->units_type==SI ) {
           d->units_type = SI;
           if( dp1->units_name != (char *)NULL )
              UnitsTypeConvert( dp1, SI );
       }
       else if(dp1->units_type == SI_US)
           d->units_type = dp2->units_type;
       else if(dp2->units_type == SI_US)
           d->units_type = dp1->units_type;
    }
	
    if(dp1->units_name != (char *)NULL) {
       if( FLAG == YES) { 
          length = UnitsLength(dp1->units_name,dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+6,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *) "(");
          d->units_name = strcat(d->units_name, dp1->units_name);
          d->units_name = strcat(d->units_name, (char *) ")");
          if( dp2->units_name!=(char *)NULL ) {
              d->units_name = strcat(d->units_name, (char *) "/(");
              d->units_name = strcat(d->units_name, dp2->units_name);
              d->units_name = strcat(d->units_name, (char *) ")");
          }
       }
       else {
          length = UnitsLength(dp1->units_name,dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+2,sizeof(char));
          d->units_name = strcpy(d->units_name, dp1->units_name);
          if( dp2->units_name!=(char *)NULL ) {
              d->units_name = strcat(d->units_name, (char *) "/");
              d->units_name = strcat(d->units_name, dp2->units_name);
          }
       } 
    }
    else if(dp2->units_name != (char *)NULL) {
       if( FLAG == YES) { 
          length = UnitsLength(dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+5,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *)"1/(");
          d->units_name = strcat(d->units_name, dp2->units_name);
          d->units_name = strcat(d->units_name, (char *) ")");
       }
       else {
          length = UnitsLength(dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+3,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *)"1/");
          d->units_name = strcat(d->units_name, dp2->units_name);
       }
    }
    else
       d->units_name = (char *)NULL;

    d->scale_factor = dp1->scale_factor/dp2->scale_factor;
    d->mass_expnt = dp1->mass_expnt - dp2->mass_expnt;
    d->length_expnt = dp1->length_expnt - dp2->length_expnt;
    d->time_expnt   = dp1->time_expnt - dp2->time_expnt;
    d->temp_expnt   = dp1->temp_expnt - dp2->temp_expnt;
    d->radian_expnt = dp1->radian_expnt - dp2->radian_expnt;

    free((char *)dp1->units_name);
    free((char *)dp1);
    free((char *)dp2->units_name);
    free((char *)dp2);

    return (d);
}

/*
 *  =====================================================================
 *  UnitsDivRep() : Units divide with replacement
 * 
 *  Input :  DIMENSIONS *d1 -- pointer to dimensions data structure
 *           DIMENSIONS *d2 -- pointer to dimensions data structure
 *           DIMENSIONS  *d -- structure for result of dimensions product
 *           int FLAG       -- YES or NO flag for "()" brackete.
 *  Output : DIMENSIONS  *d -- structure for result of dimensions product
 *  =====================================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsDivRep( DIMENSIONS *d, DIMENSIONS *d1, DIMENSIONS *d2, int FLAG )
#else
DIMENSIONS *UnitsDivRep( d, d1, d2, FLAG )
DIMENSIONS       *d;
DIMENSIONS *d1, *d2;
int            FLAG;
#endif
{
DIMENSIONS *dp1, *dp2;
int            length;

#ifdef DEBUG
       printf(" Enter UnitsDivRep() : ================\n");
       UnitsPrint(d);
#endif

    if( d2==(DIMENSIONS *)NULL )
        FatalError("In UnitsDivRep() : unit divider is a NULL pointer",(char *)NULL);

    if( d1==(DIMENSIONS *)NULL ) {
        printf("Warning: d1 is NULL,   in UnitsDivRep()\n");
        free((char *)d->units_name);
        free((char *)d);
        d = (DIMENSIONS *)NULL;
        return (DIMENSIONS *)NULL;
    }

    dp1 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    dp2 = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));

    UnitsCopy(dp1,d1);
    UnitsCopy(dp2,d2);

    if(dp1->units_type == dp2->units_type) {
       d->units_type = dp1->units_type;
    }
    else {
       if( dp1->units_type==SI && dp2->units_type==US ) {
           d->units_type = SI;
           if( dp2->units_name != (char *)NULL )
              UnitsTypeConvert( dp2, SI );
       }
       else if( dp1->units_type==US && dp2->units_type==SI ) {
           d->units_type = SI;
           if( dp1->units_name != (char *)NULL )
              UnitsTypeConvert( dp1, SI );
       }
       else if(dp1->units_type == SI_US)
           d->units_type = dp2->units_type;
       else if(dp2->units_type == SI_US)
           d->units_type = dp1->units_type;
    }
	
    if(dp1->units_name != (char *)NULL) {
       if( FLAG == YES) { 
          free((char *)d->units_name);
          length = UnitsLength(dp1->units_name,dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+6,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *)"(");
          d->units_name = strcat(d->units_name, dp1->units_name);
          d->units_name = strcat(d->units_name, (char *) ")");
          if( dp2->units_name!=(char *)NULL ) {
              d->units_name = strcat(d->units_name, (char *) "/(");
              d->units_name = strcat(d->units_name, dp2->units_name);
              d->units_name = strcat(d->units_name, (char *) ")");
          }
       }
       else {
          free((char *)d->units_name);
          length = UnitsLength(dp1->units_name,dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+2,sizeof(char));
          d->units_name = strcpy(d->units_name, dp1->units_name);
          if( dp2->units_name!=(char *)NULL ) {
              d->units_name = strcat(d->units_name, (char *) "/");
              d->units_name = strcat(d->units_name, dp2->units_name);
          }
       } 
    }
    else if(dp2->units_name != (char *)NULL) {
       if( FLAG == YES) { 
          free((char *)d->units_name);
          length = UnitsLength(dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+5,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *)"1/(");
          d->units_name = strcat(d->units_name, dp2->units_name);
          d->units_name = strcat(d->units_name, (char *) ")");
       }
       else {
          free((char *)d->units_name);
          length = UnitsLength(dp2->units_name,"STOP");
          d->units_name = (char *)MyCalloc(length+3,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *)"1/");
          d->units_name = strcat(d->units_name, dp2->units_name);
       }
    }
    else {
       free((char *)d->units_name);
       d->units_name = (char *)NULL;
    }

#ifdef DEBUG
       printf(" === In UnitsDivRep() ============= \n");
       UnitsPrint(dp1);
       UnitsPrint(dp2);
       printf(" === In UnitsDivRep() ============= \n");
#endif

    d->scale_factor = dp1->scale_factor/dp2->scale_factor;
    d->length_expnt = dp1->length_expnt - dp2->length_expnt;
    d->mass_expnt   = dp1->mass_expnt - dp2->mass_expnt;
    d->time_expnt   = dp1->time_expnt - dp2->time_expnt;
    d->temp_expnt   = dp1->temp_expnt - dp2->temp_expnt;
    d->radian_expnt = dp1->radian_expnt - dp2->radian_expnt;

    free((char *)dp1->units_name);
    free((char *)dp1);
    free((char *)dp2->units_name);
    free((char *)dp2);

#ifdef DEBUG
    printf(" At the end of UnitsDivRep() :\n");
    UnitsPrint(d);
    printf(" Leaving UnitsDivRep() :\n");
#endif

    return (d);
}

/*
 *  ==============================================================
 *  UnitsPower() :
 *
 *  Input  : DIMENSIONS *d1 : pointer to "units"
 *         : double value   : value of power "exponent"
 *         : int FLAG       : YES/NO flag for use of "()"
 *  Output : DIMENSIONS *   : pointer to "d1^value" units.
 *  ==============================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsPower( DIMENSIONS *d1, double value, int FLAG )
#else
DIMENSIONS *UnitsPower( d1, value, FLAG )
DIMENSIONS *d1;
double   value;
int       FLAG;
#endif
{
DIMENSIONS *d;
char      *cp;
int    length;

#ifdef DEBUG
       printf(" Enter UnitsPower() :\n");
       UnitsPrint(d1);
#endif

    if( d1==(DIMENSIONS *)NULL ) {
        printf("Warning: d1 is NULL,   in UnitsPower()\n");
        return (DIMENSIONS *)NULL;
    }

    d = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));

    if(d1->units_name != (char *)NULL && value != 0.0) {
       cp = (char *)calloc(20,sizeof(char));
       sprintf(cp,"%g", value);
       length = UnitsLength(d1->units_name,cp,"STOP");
       if( FLAG == YES ) {
          d->units_name = (char *)MyCalloc(length+4,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *) "(");
          d->units_name = strcat(d->units_name, d1->units_name);
          d->units_name = strcat(d->units_name, (char *) ")");
       }
       else {
          d->units_name = (char *)MyCalloc(length+2,sizeof(char));
          d->units_name = strcpy(d->units_name, d1->units_name);
       }
       d->units_name = strcat(d->units_name, (char *) "^");
       d->units_name = strcat(d->units_name, cp);
       free((char *) cp);
    }
    else 
       d->units_name = (char *)NULL;

    if(value > 0) {
       d->scale_factor =  pow(d1->scale_factor, value);
       d->units_type   = d1->units_type;
       d->length_expnt = d1->length_expnt * value;
       d->mass_expnt   = d1->mass_expnt * value;
       d->time_expnt   = d1->time_expnt * value;
       d->temp_expnt   = d1->temp_expnt * value;
       d->radian_expnt = d1->radian_expnt * value;
    }
    else if(value == 0.0) {
       d->scale_factor = 1.0;
       d->units_type   = d1->units_type;
       d->length_expnt = 0.0;
       d->mass_expnt   = 0.0;
       d->time_expnt   = 0.0;
       d->temp_expnt   = 0.0;
       d->radian_expnt = 0.0;
    }
    else {
       d->scale_factor =  1.0/pow(d1->scale_factor, ABS(value));
       d->units_type   = d1->units_type;
       d->length_expnt = d1->length_expnt * value;
       d->mass_expnt   = d1->mass_expnt * value;
       d->time_expnt   = d1->time_expnt * value;
       d->temp_expnt   = d1->temp_expnt * value;
       d->radian_expnt = d1->radian_expnt * value;
    }
	
#ifdef DEBUG
    printf(" Leaving UnitsPower() :\n");
#endif

    return (d);
}

#ifdef __STDC__
DIMENSIONS *UnitsPowerRep( DIMENSIONS *d, DIMENSIONS *d1, double value, int FLAG )
#else
DIMENSIONS *UnitsPowerRep( d, d1, value, FLAG )
DIMENSIONS  *d;
DIMENSIONS *d1;
double   value;
int       FLAG;
#endif
{
char       *cp;
DIMENSIONS *dp;
int     length;

#ifdef DEBUG
       printf(" Enter UnitsPowerRep() :\n");
       UnitsPrint(d);
#endif

    if( d1==(DIMENSIONS *)NULL ) {
        printf("Warning: d1 is NULL,   in UnitsPowerRep()\n");
        free((char *)d->units_name);
        free((char *)d);
        d = (DIMENSIONS *)NULL;
        return (DIMENSIONS *)NULL;
    }

    dp = (DIMENSIONS *) MyCalloc(1, sizeof(DIMENSIONS));
    UnitsCopy(dp,d1);

    free((char *)d->units_name);
    if(dp->units_name != (char *)NULL && value != 0.0) {
       cp = (char *)calloc(20,sizeof(char));
       sprintf(cp,"%g", value);
       length = UnitsLength(dp->units_name,cp,"STOP");
       if( FLAG == YES ) {
          d->units_name = (char *)MyCalloc(length+4,sizeof(char));
          d->units_name = strcpy(d->units_name, (char *) "(");
          d->units_name = strcat(d->units_name, dp->units_name);
          d->units_name = strcpy(d->units_name, (char *) ")");
       }
       else {
          d->units_name = (char *)MyCalloc(length+2,sizeof(char));
          d->units_name = strcpy(d->units_name, dp->units_name);
       }
       d->units_name = strcat(d->units_name, (char *) "^");
       d->units_name = strcat(d->units_name, cp);
       free((char *) cp);
    }
    else
       d->units_name = (char *)NULL;

    if(value > 0) {
       d->scale_factor =  pow(dp->scale_factor, value);
       d->units_type   = dp->units_type;
       d->length_expnt = dp->length_expnt * value;
       d->mass_expnt   = dp->mass_expnt * value;
       d->time_expnt   = dp->time_expnt * value;
       d->temp_expnt   = dp->temp_expnt * value;
       d->radian_expnt = dp->radian_expnt * value;
    }
    else if(value == 0.0) {
       d->scale_factor = 1.0;
       d->units_type   = dp->units_type;
       d->length_expnt = 0.0;
       d->mass_expnt   = 0.0;
       d->time_expnt   = 0.0;
       d->temp_expnt   = 0.0;
       d->radian_expnt = 0.0;
    }
    else {
       d->scale_factor =  1.0/pow(dp->scale_factor, ABS(value));
       d->units_type   = dp->units_type;
       d->length_expnt = dp->length_expnt * value;
       d->mass_expnt   = dp->mass_expnt * value;
       d->time_expnt   = dp->time_expnt * value;
       d->temp_expnt   = dp->temp_expnt * value;
       d->radian_expnt = dp->radian_expnt * value;
    }

    free((char *)dp->units_name);
    free((char *)dp);
	
#ifdef DEBUG
    printf(" Leaving UnitsPowerRep() :\n");
#endif

    return (d);
}

/*
 *  ============================
 *  Copy units from "d2" to "d1"
 *  ============================
 */

#ifdef __STDC__
DIMENSIONS *UnitsCopy(DIMENSIONS *d1, DIMENSIONS *d2)
#else
DIMENSIONS *UnitsCopy(d1, d2)
DIMENSIONS *d1, *d2;
#endif
{
int   length;

#ifdef DEBUG
       printf(" Enter UnitsCopy() :\n");
#endif

   /* [a] : Check for NULL pointer dimensions */

   if( d2==(DIMENSIONS *) NULL ) {
       printf("Warning : d2 is NULL in UnitsCopy(), return NULL\n");
       free((char *) d1->units_name);
       free((char *) d1);
       d1 = (DIMENSIONS *)NULL;
       return d1;
   }

   /* [b] : Copy units data structure */

   if(d2->units_name != (char *)NULL) {
      free((char *)d1->units_name);
      d1->units_name = SaveString(d2->units_name);
   } else {
      free((char *)d1->units_name);
      d1->units_name = (char *)NULL;
   }

   d1->units_type   = d2->units_type;
   d1->scale_factor = d2->scale_factor;
   d1->length_expnt = d2->length_expnt;
   d1->mass_expnt   = d2->mass_expnt;
   d1->time_expnt   = d2->time_expnt;
   d1->temp_expnt   = d2->temp_expnt;
   d1->radian_expnt = d2->radian_expnt;
	
#ifdef DEBUG
       printf(" Leaving UnitsCopy() :\n");
#endif

   return (d1);
}

/*
 *  ======================
 *  Zero out units in "d1"
 *  ======================
 */

#ifdef __STDC__
DIMENSIONS *ZeroUnits(DIMENSIONS *d1)
#else
DIMENSIONS *ZeroUnits(d1)
DIMENSIONS *d1;
#endif
{
#ifdef DEBUG
       printf("Enter ZeroUnits()\n");
#endif

     if( d1==(DIMENSIONS *)NULL ) {
         printf("Warning:  d1 is NULL, in ZeroUnits(), return NULL\n");
         return d1;
     }

     if( d1->units_name != (char *)NULL ) {
         free((char *)d1->units_name);
         d1->units_name   = (char *)NULL;
     }

     d1->scale_factor = 1.0;
     d1->length_expnt = 0.0;
     d1->mass_expnt   = 0.0;
     d1->time_expnt   = 0.0;
     d1->temp_expnt   = 0.0;
     d1->radian_expnt = 0.0;
     if(d1->units_type != SI && d1->units_type != US)
        d1->units_type = SI_US;

#ifdef DEBUG
       printf("Leaving ZeroUnits()\n");
#endif

     return (d1);
}

/*
 *  ===================================
 *  Get Default Units from Symbol Table
 *  ===================================
 */

#ifdef __STDC__
DIMENSIONS *DefaultUnits(char *name)
#else
DIMENSIONS *DefaultUnits(name)
char *name;
#endif
{
DIMENSIONS  *dp, *d;
	
   d = (DIMENSIONS *)NULL;
   d = lookup(name)->u.q->dimen;
   if(d == (DIMENSIONS *)NULL){
      printf("ERROR: \"%s\" is not found in hash_table \n", name);
      FatalError("In DefaultUnits()",(char *)NULL);
   }

   dp    = (DIMENSIONS *) MyMalloc(sizeof(DIMENSIONS));
   dp->units_name   = SaveString(name);
   dp->scale_factor = d->scale_factor;
   dp->units_type   = d->units_type;
   dp->length_expnt = d->length_expnt;
   dp->mass_expnt   = d->mass_expnt;
   dp->time_expnt   = d->time_expnt;
   dp->temp_expnt   = d->temp_expnt;
   dp->radian_expnt = d->radian_expnt;

   return(dp);
}


/*
 *  =================================================================
 *  UnitsSimplify() : Simplify Units Expression
 *
 *  Input  : DIMENSIONS * -- pointer to unsimplified units expression
 *  Output : DIMENSIONS * -- pointer to simplified units expression
 *  =================================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsSimplify(DIMENSIONS *d)
#else
DIMENSIONS *UnitsSimplify(d)
DIMENSIONS *d;
#endif
{
DIMENSIONS *d1, *d2, *d3, *d4, *dp;
DIMENSIONS *dp1, *dp2, *dp3, *dp4;
DIMENSIONS *dd1, *dd2;
SYMBOL     *hp;
double      dlength;

#ifdef DEBUG
       printf(" Enter UnitsSimplify() \n");
       UnitsPrint(d);
#endif

   /* [a] : Check that units are defined */

   if( d == (DIMENSIONS *) NULL ) {
       printf("Warning : d is NULL in UnitsSimplify()\n");
       return (DIMENSIONS *) NULL;
   }

   /* [b] : Don't simplify units if "unit name" already in symbol table */
          
   hp = NULL;
   if(d->units_name != (char *) NULL) {
      hp = lookup(d->units_name);
      if(hp != NULL) {
         return(d);
      }
   } else
      return(d);

   /* [c] : Simplify dimensionless/radians only units */

   if( d->length_expnt == 0 && d->time_expnt == 0 && 
       d->mass_expnt   == 0 && d->temp_expnt == 0 ) {

          if(d->radian_expnt == 0) {
             free((char *)d->units_name);
             d->units_name = (char *) NULL;
             d->scale_factor = 1.0;
          } 

          if(d->radian_expnt == 1) {
             free((char *) d->units_name);
             d->units_name = (char *) SaveString("rad");
          } 

          return(d);
   }

   /* [d] : Simplify LENGTH units alone */

   if(d->length_expnt != 0 && d->time_expnt == 0 && 
      d->mass_expnt   == 0 && d->temp_expnt == 0 ) {

      /* Setup default units for length */

      if(d->units_type == SI)
         dp = DefaultUnits("m");
      else
         dp = DefaultUnits("in");	 

      /* Transfer "radians" units to working data structure */

      if(d->radian_expnt != 0 ) {
         dp->radian_expnt = d->radian_expnt;
      }

      /* Build working data structures */

      if(d->length_expnt == 1) {
         if(strcmp(d->units_name, dp->units_name) != 0) {
            free((char *) d->units_name);
            d->units_name = SaveString(dp->units_name);
            d->scale_factor = dp->scale_factor;
         }
      }
      else if(d->length_expnt < 0) {
         d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(d1);
         d1->units_type = d->units_type;
         if(d->length_expnt == -1) {
            UnitsDivRep(d, d1, dp, NO);
         } else {
            d2 = UnitsPower( dp, ABS(d->length_expnt), NO );
            UnitsDivRep(d, d1, d2, NO);
            free((char *) d2->units_name);
            free((char *) d2);
         }
         free((char *) d1->units_name);
         free((char *) d1);
      } else
         UnitsPowerRep(d, dp, d->length_expnt, NO);

      free((char *) dp->units_name);
      free((char *) dp);

      /* If appropriate, append radian components to units name */

      if(d->radian_expnt != 0) {
         RadUnitsNameExtension( d );
      }

      return(d);
   }

   /* [e] : Simplify MASS units alone */
         
   if(d->length_expnt == 0 && d->time_expnt == 0 && 
      d->mass_expnt   != 0 && d->temp_expnt == 0 ) {

      /* Default units for "m" mass */

      if(d->units_type == SI)
         dp = DefaultUnits("kg");
      else
         dp = DefaultUnits("lb");	 

      /* Transfer "radians" units to working data structure */

      if(d->radian_expnt != 0 ) {
         dp->radian_expnt = d->radian_expnt;
      }

      /* Build working data structures */

      if(d->mass_expnt == 1) {
         if(strcmp(d->units_name, dp->units_name) != 0) {
            free((char *) d->units_name);
            d->units_name = SaveString(dp->units_name);
            d->scale_factor = dp->scale_factor;
         }
      } else if(d->mass_expnt < 0) {
         d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(d1);
         d1->units_type = d->units_type;
         if(d->mass_expnt == -1) {
            UnitsDivRep(d, d1, dp, NO);
         }
         else {
            d2 = UnitsPower( dp, ABS(d->mass_expnt), NO );
            UnitsDivRep(d, d1, d2, NO);
            free((char *) d2->units_name);
            free((char *) d2);
         }
         free((char *) d1->units_name);
         free((char *) d1);
      } else
         UnitsPowerRep(d, dp, d->mass_expnt, NO);

      free((char *) dp->units_name);
      free((char *) dp);

      /* Append radian components to units name */

      if(d->radian_expnt != 0) {
         RadUnitsNameExtension( d );
      }

      return(d);
   }

   /* [f] : Simplify TIME units alone */

   if(d->length_expnt == 0 && d->time_expnt != 0 && 
      d->mass_expnt   == 0 && d->temp_expnt == 0 ) {

#ifdef DEBUG
      printf("In UnitsSimplify() : Start TIME Units\n");
      UnitsPrint(d);
#endif

      /* Special cases of "time" and "radians and time" */

      if(d->radian_expnt == 0 && d->time_expnt == -1) {
         free((char *) d->units_name);
         d->units_name = (char *) SaveString("Hz");
         return(d);
      }

      if(d->radian_expnt == 1 && d->time_expnt == -1) {
         free((char *) d->units_name);
         d->units_name = (char *) SaveString("rad/sec");
         return(d);
      }

      if(d->radian_expnt == 1 && d->time_expnt == -2) {
         free((char *) d->units_name);
         d->units_name = (char *) SaveString("rad/sec^2");
         return(d);
      }

      /* Build default units and transfer "radians" units */
      
      dp = DefaultUnits("sec");
      if(d->radian_expnt != 0 ) {
         dp->radian_expnt = d->radian_expnt;
      }

      /* Build working data structurs */

      if(d->time_expnt == 1) {
         if(strcmp(d->units_name, dp->units_name) != 0) {
            free((char *) d->units_name);
            d->units_name   = SaveString(dp->units_name);
            d->scale_factor = dp->scale_factor;
         }
      } else if( d->time_expnt < 0 ) {
         d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(d1);
         d1->units_type = SI_US;
         if(d->time_expnt == -1) {
            UnitsDivRep(d, d1, dp, NO);
         } else {
            d2 = UnitsPower( dp, ABS(d->time_expnt), NO );
            UnitsDivRep(d, d1, d2, NO);
            free((char *) d2->units_name);
            free((char *) d2);
         }
         free((char *) d1->units_name);
         free((char *) d1);
      } else
         UnitsPowerRep(d, dp, d->time_expnt, NO);

      free((char *) dp->units_name);
      free((char *) dp);

      /* Append radian components to units name */

      if(d->radian_expnt != 0) {
         RadUnitsNameExtension( d );
      }

#ifdef DEBUG
      printf("UnitsSimplify() : At the end of TIME Units\n");
      UnitsPrint(d);
      printf("UnitsSimplify() : Leaving TIME Units\n");
#endif

      return(d);
   }

   /* [g] : Simplify TEMPERATURE units alone */

   if(d->length_expnt == 0 && d->time_expnt == 0 && 
      d->mass_expnt   == 0 && d->temp_expnt != 0 ) {

      /* Build default units */

      if(d->units_type == SI)
         dp = DefaultUnits("deg_C");
      else
         dp = DefaultUnits("deg_F");	 

      /* Transfer "radians" units to working data structure */

      if(d->radian_expnt != 0 ) {
         dp->radian_expnt = d->radian_expnt;
      }

      /* Setup working data structures */

      if(d->temp_expnt == 1) {
         if(strcmp(d->units_name, dp->units_name) != 0) {
            free((char *) d->units_name);
            d->units_name = SaveString(dp->units_name);
            d->scale_factor = dp->scale_factor;
         }
      } else if(d->temp_expnt < 0) {
         d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
         ZeroUnits(d1);
         d1->units_type = d->units_type;
         if(d->temp_expnt == -1) {
            UnitsDivRep(d, d1, dp, NO);
         } else {
            d2 = UnitsPower( dp, ABS(d->temp_expnt), NO );
            UnitsDivRep(d, d1, d2, NO);
            free((char *) d2->units_name);
            free((char *) d2);
         }
         free((char *) d1->units_name);
         free((char *) d1);
       } else
         UnitsPowerRep(d, dp, d->temp_expnt, NO);

       free((char *) dp->units_name);
       free((char *) dp);

       /* Append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
   }

   /* [h] : Simplify groups of related units */

   /* [h.1] : Simplify LENGTH and MASS related units  */

   if( d->length_expnt != 0 && d->time_expnt == 0 && 
       d->mass_expnt   != 0 && d->temp_expnt == 0 ) {

       /* Setup default units names */

       if(d->units_type == SI) {
          dp1 = DefaultUnits("m");
          dp2 = DefaultUnits("kg");
       } else {
          dp1 = DefaultUnits("in");
          dp2 = DefaultUnits("lb");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          dp1->radian_expnt = d->radian_expnt;
          dp2->radian_expnt = d->radian_expnt;
       }

       /* Setup working data structures */

       if(ABS(d->length_expnt) == 1) {
          d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d1,dp1);
       } else
          d1 = UnitsPower( dp1, ABS(d->length_expnt), NO );

       if(ABS(d->mass_expnt) == 1) {
          d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d2,dp2);
       } else
          d2 = UnitsPower( dp2, ABS(d->mass_expnt), NO );

       if(d->length_expnt > 0) {
          if(d->mass_expnt > 0)
             UnitsMultRep(d, d1, d2);	 
          else
             UnitsDivRep(d, d1, d2, NO);	 
       } else if(d->mass_expnt > 0) {
          UnitsDivRep(d, d2, d1, NO);	 
       } else {
          dd1 = UnitsPower( dp1, d->length_expnt, NO );
          dd2 = UnitsPower( dp2, d->mass_expnt, NO );
          UnitsMultRep(d, dd1, dd2);	 
          free((char *) dd1->units_name);
          free((char *) dd1);
          free((char *) dd2->units_name);
          free((char *) dd2);
       }

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) dp1->units_name);
       free((char *) dp1);
       free((char *) dp2->units_name);
       free((char *) dp2);

       /* Append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
   }

   /* [h.2] : Simplify LENGTH and TIME related units */

   if( d->length_expnt != 0 && d->time_expnt != 0 && 
       d->mass_expnt   == 0 && d->temp_expnt == 0 ) {

       /* Setup default units */

       if(d->units_type == SI) {
          dp1 = DefaultUnits("m");
          dp2 = DefaultUnits("sec");
       } else {
          dp1 = DefaultUnits("in");
          dp2 = DefaultUnits("sec");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          dp1->radian_expnt = d->radian_expnt;
          dp2->radian_expnt = d->radian_expnt;
       }

       /* Setup working "units" data structures */

       if(ABS(d->length_expnt) == 1) {
          d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d1,dp1);
       } else
          d1 = UnitsPower( dp1, ABS(d->length_expnt), NO );

       if(ABS(d->time_expnt) == 1) {
          d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d2,dp2);
       } else
          d2 = UnitsPower( dp2, ABS(d->time_expnt), NO );

       if(d->length_expnt > 0) {
          if(d->time_expnt > 0)
             UnitsMultRep(d, d1, d2);	 
          else
             UnitsDivRep(d, d1, d2, NO);	 
       } else if(d->time_expnt > 0) {
           UnitsDivRep(d, d2, d1, NO);	 
       } else {
           dd1 = UnitsPower( dp1, d->length_expnt, NO );
           dd2 = UnitsPower( dp2, d->time_expnt, NO );
           UnitsMultRep(d, dd1, dd2);	 
           free((char *) dd1->units_name);
           free((char *) dd1);
           free((char *) dd2->units_name);
           free((char *) dd2);
       }

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) dp1->units_name);
       free((char *) dp1);
       free((char *) dp2->units_name);
       free((char *) dp2);

       /* If appropriate, append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
   }

   /* [h.3] : Simplify LENGTH and TEMPERATURE related units */

   if( d->length_expnt != 0 && d->time_expnt == 0 && 
       d->mass_expnt   == 0 && d->temp_expnt != 0 ) {

       /* Build default units */

       if(d->units_type == SI) {
          dp1 = DefaultUnits("m");
          dp2 = DefaultUnits("deg_C");
       } else {
          dp1 = DefaultUnits("in");
          dp2 = DefaultUnits("deg_F");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          dp1->radian_expnt = d->radian_expnt;
          dp2->radian_expnt = d->radian_expnt;
       }

       /* Setup working "units" data structures */

       if(ABS(d->length_expnt) == 1) {
          d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d1,dp1);
       } else
          d1 = UnitsPower( dp1, ABS(d->length_expnt), NO );

       if(ABS(d->temp_expnt) == 1) {
          d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d2,dp2);
       } else
          d2 = UnitsPower( dp2, ABS(d->temp_expnt), NO );

       if(d->length_expnt > 0) {
          if(d->temp_expnt > 0)
             UnitsMultRep(d, d1, d2);	 
          else
             UnitsDivRep(d, d1, d2, NO);	 
       } else if(d->temp_expnt > 0) {
          UnitsDivRep(d, d2, d1, NO);	 
       } else {
          dd1 = UnitsPower( dp1, d->length_expnt, NO );
          dd2 = UnitsPower( dp2, d->temp_expnt, NO );
          UnitsMultRep(d, dd1, dd2);	 
          free((char *) dd1->units_name);
          free((char *) dd1);
          free((char *) dd2->units_name);
          free((char *) dd2);
       }

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) dp1->units_name);
       free((char *) dp1);
       free((char *) dp2->units_name);
       free((char *) dp2);

       /* Append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
    }

    /* [h.4] : Simplify TIME and MASS related units */

    if(d->length_expnt == 0 && d->time_expnt == -2 &&
       d->mass_expnt   == 1 && d->temp_expnt ==  0 ) {

       /* Setup default units names */

       if(d->units_type == SI){
          d1 = DefaultUnits("N");
          d2 = DefaultUnits("m");
       } else{
          d1 = DefaultUnits("lbf");
          d2 = DefaultUnits("in");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          d1->radian_expnt = d->radian_expnt;
          d2->radian_expnt = d->radian_expnt;
       }

       UnitsDivRep(d, d1, d2, NO);	 
       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       return(d);
    }

    /* [h.5] : Simplify LENGTH/FORCE related units */

    if(d->length_expnt ==  0 && d->time_expnt == 2 &&
       d->mass_expnt   == -1 && d->temp_expnt == 0 ) {

       /* Setup default units names */

       if(d->units_type == SI){
          d1 = DefaultUnits("m");
          d2 = DefaultUnits("N");
       } else{
          d1 = DefaultUnits("in");
          d2 = DefaultUnits("lbf");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          d1->radian_expnt = d->radian_expnt;
          d2->radian_expnt = d->radian_expnt;
       }

       UnitsDivRep(d, d1, d2, NO);	 
       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);

       /* Append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
    }

    /* [h.6] : Simplify TIME and MASS related units */

    if(d->length_expnt == 0 && d->time_expnt != 0 && 
       d->mass_expnt   != 0 && d->temp_expnt == 0 ) {

       /* Setup default units names */

       if(d->units_type == SI) {
          dp1 = DefaultUnits("sec");
          dp2 = DefaultUnits("kg");
       } else {
          dp1 = DefaultUnits("sec");
          dp2 = DefaultUnits("lb");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          dp1->radian_expnt = d->radian_expnt;
          dp2->radian_expnt = d->radian_expnt;
       }

       /* Setup working "units" data structures */

       if(ABS(d->time_expnt) == 1) {
          d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d1,dp1);
       } else
          d1 = UnitsPower( dp1, ABS(d->time_expnt), NO );

       if(ABS(d->mass_expnt) == 1) {
          d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d2,dp2);
       } else
          d2 = UnitsPower( dp2, ABS(d->mass_expnt), NO );

       if(d->time_expnt > 0) {
          if(d->mass_expnt > 0)
             UnitsMultRep(d, d1, d2);	 
          else
             UnitsDivRep(d, d1, d2, NO);	 
       } else if(d->mass_expnt > 0) {
          UnitsDivRep(d, d2, d1, NO);	 
       } else {
          dd1 = UnitsPower( dp1, d->time_expnt, NO );
          dd2 = UnitsPower( dp2, d->mass_expnt, NO );
          UnitsMultRep(d, dd1, dd2);	 
          free((char *) dd1->units_name);
          free((char *) dd1);
          free((char *) dd2->units_name);
          free((char *) dd2);
       }

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) dp1->units_name);
       free((char *) dp1);
       free((char *) dp2->units_name);
       free((char *) dp2);

       /* Append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
    }

    /* [h.7] : Simplify TIME and TEMPERATURE related units */

    if(d->length_expnt == 0 && d->time_expnt != 0 && 
       d->mass_expnt   == 0 && d->temp_expnt != 0 ) {

       /* Default units */

       if(d->units_type == SI) {
          dp1 = DefaultUnits("sec");
          dp2 = DefaultUnits("deg_C");
       } else {
          dp1 = DefaultUnits("sec");
          dp2 = DefaultUnits("deg_F");
       }

       /* Transfer "radians" units to working data structure */

       if(d->radian_expnt != 0 ) {
          dp1->radian_expnt = d->radian_expnt;
          dp2->radian_expnt = d->radian_expnt;
       }

       /* Setup working "units" data structures */

       if(ABS(d->time_expnt) == 1) {
          d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d1,dp1);
       } else
           d1 = UnitsPower( dp1, ABS(d->time_expnt), NO );

       if(ABS(d->temp_expnt) == 1) {
          d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d2,dp2);
       } else
          d2 = UnitsPower( dp2, ABS(d->temp_expnt), NO );

       if(d->time_expnt > 0) {
          if(d->temp_expnt > 0)
             UnitsMultRep(d, d1, d2);	 
          else
             UnitsDivRep(d, d1, d2, NO);	 
       } else if(d->temp_expnt > 0) {
             UnitsDivRep(d, d2, d1, NO);	 
       } else {
             dd1 = UnitsPower( dp1, d->time_expnt, NO );
             dd2 = UnitsPower( dp2, d->temp_expnt, NO );
             UnitsMultRep(d, dd1, dd2);	 
             free((char *) dd1->units_name);
             free((char *) dd1);
             free((char *) dd2->units_name);
             free((char *) dd2);
       }

       free((char *) d1->units_name);
       free((char *) d1);
       free((char *) d2->units_name);
       free((char *) d2);
       free((char *) dp1->units_name);
       free((char *) dp1);
       free((char *) dp2->units_name);
       free((char *) dp2);

       /* Append radian components to units name */

       if(d->radian_expnt != 0) {
          RadUnitsNameExtension( d );
       }

       return(d);
   }

   /* [h.8] : Simplify MASS and TEMPERATURE related units */

   if(d->length_expnt == 0 && d->time_expnt == 0 && 
      d->mass_expnt   != 0 && d->temp_expnt != 0 ) {

      /* Setup default units */

      if(d->units_type == SI) {
         dp1 = DefaultUnits("kg");
         dp2 = DefaultUnits("deg_C");
      } else {
         dp1 = DefaultUnits("lb");
         dp2 = DefaultUnits("deg_F");
      }

      /* Transfer "radians" units to working data structure */

      if(d->radian_expnt != 0 ) {
         dp1->radian_expnt = d->radian_expnt;
         dp2->radian_expnt = d->radian_expnt;
      }

      /* Setup working "units" data structures */

      if(ABS(d->mass_expnt) == 1) {
         d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
         UnitsCopy(d1,dp1);
      } else
         d1 = UnitsPower( dp1, ABS(d->mass_expnt), NO );

      if(ABS(d->temp_expnt) == 1) {
         d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
         UnitsCopy(d2,dp2);
      } else
         d2 = UnitsPower( dp2, ABS(d->temp_expnt), NO );

      if(d->mass_expnt > 0) {
         if(d->temp_expnt > 0)
            UnitsMultRep(d, d1, d2);	 
         else
            UnitsDivRep(d, d1, d2, NO);	 
      } else if(d->temp_expnt > 0) {
         UnitsDivRep(d, d2, d1, NO);	 
      } else {
         dd1 = UnitsPower( dp1, d->mass_expnt, NO );
         dd2 = UnitsPower( dp2, d->temp_expnt, NO );
         UnitsMultRep(d, dd1, dd2);	 
         free((char *) dd1->units_name);
         free((char *) dd1);
         free((char *) dd2->units_name);
         free((char *) dd2);
      }

      free((char *) d1->units_name);
      free((char *) d1);
      free((char *) d2->units_name);
      free((char *) d2);
      free((char *) dp1->units_name);
      free((char *) dp1);
      free((char *) dp2->units_name);
      free((char *) dp2);

      /* Append radian components to units name */

      if(d->radian_expnt != 0) {
         RadUnitsNameExtension( d );
      }

      return(d);
   }

   /* [i] : Simplify PRESSURE units */

   if(d->length_expnt == -1 && d->time_expnt == -2 && 
      d->mass_expnt   ==  1 && d->temp_expnt ==  0 ) {

      if(d->units_type == SI){
         if(strcmp(d->units_name, "Pa") != 0) {
            dp = DefaultUnits("Pa");
            if(d->radian_expnt != 0 ) {
               dp->radian_expnt = d->radian_expnt;
            }
            UnitsCopy(d,dp);
            free((char *) dp->units_name);
            free((char *) dp);
         }
      } else {
         if(strcmp(d->units_name, "psi") != 0) {
            dp = DefaultUnits("psi"); 
            if(d->radian_expnt != 0 ) {
               dp->radian_expnt = d->radian_expnt;
            }
            UnitsCopy(d,dp);
            free((char *) dp->units_name);
            free((char *) dp);
         }
      }

      return(d);
   }

   /* [j] : Simplify POWER units */

   if( d->length_expnt == 2 && d->time_expnt == -3 && 
       d->mass_expnt   == 1 && d->temp_expnt ==  0 ) {
       if(strcmp(d->units_name, "Watt") != 0) {
          dp = DefaultUnits("Watt");
          UnitsCopy(d,dp);
          free((char *) dp->units_name);
          free((char *) dp);
       }
       return(d);
   }

   /* [k] : Simplify FORCE units */

   if(d->length_expnt != 0 && d->time_expnt != 0 && 
      d->mass_expnt   != 0 && d->temp_expnt == 0 ) {

      dlength = d->length_expnt - d->mass_expnt;
      if( d->time_expnt == -2.0*d->mass_expnt ) {

          /* Setup default units of "N" force */

          if(d->units_type == SI){
             dp1 = DefaultUnits("N");
             dp2 = DefaultUnits("m");
          } else {
             dp1 = DefaultUnits("lbf"); 
             dp2 = DefaultUnits("in");
          }

          /* Transfer "radians" units to working data structure */

          if(d->radian_expnt != 0 ) {
             dp1->radian_expnt = d->radian_expnt;
             dp2->radian_expnt = d->radian_expnt;
          }

          /* Setup dimensions for mass units */

          if(ABS(d->mass_expnt) == 1) {
             d1 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             UnitsCopy(d1,dp1);
          } else
             d1 = UnitsPower( dp1, ABS(d->mass_expnt), NO );

          /* Setup dimensions for mass units */

          if(dlength == 0) {
             d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             ZeroUnits(d2);
             d2->units_type = dp2->units_type;
          } else if(ABS(dlength) == 1) {
             d2 = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
             UnitsCopy(d2,dp2);
          } else
             d2 = UnitsPower( dp2, ABS(dlength), NO );

          if(d->mass_expnt > 0) {
             if(dlength > 0) {
                UnitsMultRep(d, d1, d2);
             }
             else if(dlength == 0) {
                UnitsCopy(d,d1);
             } else
                UnitsDivRep(d, d1, d2, NO);
          }
          else if(dlength >= 0) {
             UnitsDivRep(d, d2, d1, NO);
          }
          else {
             dd1 = UnitsPower( dp1, d->mass_expnt, NO );
             dd2 = UnitsPower( dp2, dlength, NO );
             UnitsMultRep(d, dd1, dd2);	 
             free((char *) dd1->units_name);
             free((char *) dd1);
             free((char *) dd2->units_name);
             free((char *) dd2);
          }

          free((char *)d1->units_name);
          free((char *)d1);
          free((char *)d2->units_name);
          free((char *)d2);
          free((char *)dp1->units_name);
          free((char *)dp1);
          free((char *)dp2->units_name);
          free((char *)dp2);
       } else {

          /* Setup default "force" units */

          if(d->units_type == SI) {
             dp1 = DefaultUnits("m");
             dp2 = DefaultUnits("kg");
             dp3 = DefaultUnits("sec");
          } else {
             dp1 = DefaultUnits("in");
             dp2 = DefaultUnits("lb");
             dp3 = DefaultUnits("sec");
          }

          /* Correction for "radians" units */

          if(d->radian_expnt != 0 ) {
             dp1->radian_expnt = d->radian_expnt;
             dp2->radian_expnt = d->radian_expnt;
          }

          if( d->length_expnt==1 ) {
              d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy(d1,dp1);
          } else
              d1 = UnitsPower(dp1,d->length_expnt, NO);

          if( d->mass_expnt==1 ) {
              d2 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy(d2,dp2);
          } else
              d2 = UnitsPower(dp2,d->mass_expnt, NO);

          if( d->time_expnt==1 ) {
              d3 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
              UnitsCopy(d3,dp3);
          } else
              d3 = UnitsPower(dp3,d->time_expnt, NO);

          dd1 = UnitsMult(d1,d2);
          UnitsMultRep(d,dd1,d3);

          free((char *)d1->units_name);
          free((char *)d1);
          free((char *)d2->units_name);
          free((char *)d2);
          free((char *)d3->units_name);
          free((char *)d3);
          free((char *)dp1->units_name);
          free((char *)dp1);
          free((char *)dp2->units_name);
          free((char *)dp2);
          free((char *)dp3->units_name);
          free((char *)dp3);
          free((char *)dd1->units_name);
          free((char *)dd1);
      }

      /* Append radian components to units name */

      if(d->radian_expnt != 0) {
         RadUnitsNameExtension( d );
      }

      return(d);
   }

   /* [l] : Express units in compact form */

   if(d->units_type == SI) {
      dp1 = DefaultUnits("m");
      dp2 = DefaultUnits("kg");
      dp3 = DefaultUnits("sec");
      dp4 = DefaultUnits("deg_C");
   } else {
      dp1 = DefaultUnits("in");
      dp2 = DefaultUnits("lb");
      dp3 = DefaultUnits("sec");
      dp4 = DefaultUnits("deg_F");
   }

   /* Units of Mass */

   if( d->mass_expnt!=0 ) {
       if( d->mass_expnt==1 ) {
           d2 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
           UnitsCopy(d2,dp2);
       } else
           d2 = UnitsPower(dp2,d->mass_expnt, NO);
   } else {
       d2 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
       ZeroUnits(d2);
       d2->units_type = SI_US;
   }

   /* Units of Length */

   if( d->length_expnt!=0 ) {
      if( d->length_expnt==1 ) {
          d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
          UnitsCopy(d1,dp1);
      } else
          d1 = UnitsPower(dp1,d->length_expnt, NO);
   }
   else {
      d1 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
      ZeroUnits(d1);
      d1->units_type = SI_US;
   }

   /* Units of Time */

   if( d->time_expnt!=0 ) {
       if( d->time_expnt==1 ) {
           d3 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
           UnitsCopy(d3,dp3);
       } else
           d3 = UnitsPower(dp3,d->time_expnt, NO);
   } else {
       d3 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
       ZeroUnits(d3);
       d3->units_type = SI_US;
   }

   if( d->temp_expnt!=0 ) {
       if( d->temp_expnt==1 ) {
           d4 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
           UnitsCopy(d4,dp4);
       } else
           d4 = UnitsPower(dp4,d->temp_expnt, NO);
   } else {
       d4 = (DIMENSIONS *)MyCalloc(1,sizeof(DIMENSIONS));
       ZeroUnits(d4);
       d4->units_type = SI_US;
   }

   dd1 = UnitsMult(d1,d2);
   dd2 = UnitsMult(d3,d4);

   UnitsMultRep(d,dd1,dd2);

   free((char *)d1->units_name);
   free((char *)d1);
   free((char *)d2->units_name);
   free((char *)d2);
   free((char *)d3->units_name);
   free((char *)d3);
   free((char *)d4->units_name);
   free((char *)d4);
   free((char *)dp1->units_name);
   free((char *)dp1);
   free((char *)dp2->units_name);
   free((char *)dp2);
   free((char *)dp3->units_name);
   free((char *)dp3);
   free((char *)dp4->units_name);
   free((char *)dp4);
   free((char *)dd1->units_name);
   free((char *)dd1);
   free((char *)dd2->units_name);
   free((char *)dd2);

#ifdef DEBUG
     printf(" Leaving UnitsSimplify() \n");
#endif

   return(d);
}


/*
 *  ============================================================
 *  RadUnitsNameExtension() : Extend "radians" to units name.
 *  
 *  Input : DIMENSIONS *d : Pointer to dimensions data structure  
 *  Ouput : DIMENSIONS *d : Pointer to "simplified" dimensions
 *  ============================================================
 */

#ifdef __STDC__
RadUnitsNameExtension( DIMENSIONS *d )
#else
RadUnitsNameExtension ( d )
DIMENSIONS *d;
#endif
{

#ifdef DEBUG
       printf("Enter RadUnitsNameExtension() \n");
#endif

   /* [a] : Check for simple extension to name */

   switch ((int) d->radian_expnt) {
       case -1:
            d->units_name = strcat(d->units_name, (char *) "/rad");
            break;
       case 1:
            d->units_name = strcat(d->units_name, (char *) ".rad");
            break;
       default:
            break;
   }

   /* [b] : Energy due to rotation */

   if( d->mass_expnt ==  1 && d->length_expnt == 2 &&
       d->time_expnt == -2 && d->radian_expnt == 1 ) {
       free ( (char *) d->units_name );
       d->units_name = (char *) SaveString ("Joule");
   }

#ifdef DEBUG
       printf("Leave RadUnitsNameExtension() \n");
#endif

}

/*
 *  ============================================================
 *  RadUnitsSimplify() : Simplify Radian Units
 *  
 *  Input : DIMENSIONS *d : Pointer to dimensions data structure  
 *  Ouput : DIMENSIONS *d : Pointer to "simplified" dimensions
 *  ============================================================
 */

#ifdef __STDC__
DIMENSIONS *RadUnitsSimplify( DIMENSIONS *d )
#else
DIMENSIONS *RadUnitsSimplify( d )
DIMENSIONS *d;
#endif
{
SYMBOL     *hp;

#ifdef DEBUG
     printf(" Enter RadUnitsSimplify() \n");
#endif

   /* [a] : Are unit names are in the symbol table */
       
   if( d==(DIMENSIONS *)NULL ) {
       printf("Warning: d is NULL,   in RadUnitsSimplify()\n");
       return (DIMENSIONS *)NULL;
   }

   hp = (SYMBOL *)NULL;
   if(d->units_name != NULL)
      hp = lookup(d->units_name);

   if(hp != NULL)
      return(d);
          
   /* [b] : Check for radian units */
         
   if(d->length_expnt == 0 && d->time_expnt == 0 && 
      d->mass_expnt   == 0 && d->temp_expnt == 0 ) {
      free((char *)d->units_name);
      d->units_name   = SaveString( (char *) "rad" );
      d->scale_factor = 1.0;
      d->units_type   = SI_US;
      return(d);
   }

#ifdef DEBUG
     printf(" Leaving RadUnitsSimplify() \n");
#endif

   return(d);
}

/*
 *  ==========================
 *  UnitsPrint() : Print Units
 *  ==========================
 */

#ifdef __STDC__
void UnitsPrint(DIMENSIONS *d)
#else
void UnitsPrint(d)
DIMENSIONS *d;
#endif
{

#ifdef DEBUG
       printf("\n");
#endif

   if( d==(DIMENSIONS *)NULL ) {
       printf("Unit d is NULL,   in UnitsPrint()\n");
       return;
   }
  
   if( d->units_name != (char *)NULL )
       printf("\n units_name = %s\n", d->units_name);
   else
       printf("\n units_name = %s\n", (char *)"NULL");

   printf(" d->units_type   = %d \n",d->units_type);
   printf(" d->scale_factor = %lg\n",d->scale_factor);
   printf(" d->mass_expnt   = %g \n",d->mass_expnt );
   printf(" d->length_expnt = %g \n",d->length_expnt);
   printf(" d->time_expnt   = %g \n",d->time_expnt  );
   printf(" d->temp_expnt   = %g \n",d->temp_expnt );
   printf(" d->radian_expnt = %g \n",d->radian_expnt );

#ifdef DEBUG
      printf("\n");
#endif

}

/* 
 *  ==================================================================
 *  Usage :
 *  int   : L ; char * name1, * name2, ...                   
 *          L = UnitsLength(name1, "STOP");  
 *  or      L = UnitsLength(name1, name2, ..., "STOP");   
 *  ==================================================================
 */ 

#ifdef __STDC__
int UnitsLength( char *first_name, ... )
#else
int UnitsLength(va_alist)
va_dcl
#endif
{
char       *name;
int       length;

#ifdef __STDC__
va_list  arg_ptr;
#else
va_list  arg_ptr;
char *first_name;
#endif

#ifdef DEBUG
     printf(" Enter UnitsLength() \n");
#endif

#ifdef __STDC__
    va_start(arg_ptr, first_name);
#else
    va_start(arg_ptr);
    first_name = va_arg(arg_ptr, char *);
#endif 

   for (name = first_name, length = (int) 0;    ; name = va_arg(arg_ptr, char * ) ) {
       if(name != (char *)NULL && !strcmp(name, "STOP") )
          break;
        if (name != (char *)NULL) {
              length += strlen(name);
        }
        else 
           length += 0;
   }
   va_end(arg_ptr);

#ifdef DEBUG
     printf(" Leaving UnitsLength() length = %d\n", length);
#endif

   return(length);
}

/* 
 *  ==================================================
 *  BufferPrint() : Print contents of matrix buffer
 *  ==================================================
 */ 

#ifdef __STDC__
void BufferPrint(char *name, DIMENSIONS *array, int no_col)
#else
void BufferPrint(name, array, no_col)
char        *name;
DIMENSIONS *array;
int        no_col;
#endif
{
int i;
#ifdef DEBUG
  printf("\n Enter BufferPrint()\n");
#endif 

   printf("\n In BufferPrint() no_of col = %d \n", no_col);
   printf(" \n buffer_name = %s\n", name);

   if(array != (DIMENSIONS *)NULL) { 
      for (i = 1; i <= no_col; i++) {
           printf("\n col/row[%d] \n", i);
           UnitsPrint(&array[i-1]);
      }
   } else 
      printf("\n Buffer = NULL \n");

#ifdef DEBUG
      printf("\n Leaving BufferPrint()\n");
#endif 
}

/* 
 *  ==================================================
 *  BufferInit() : Allocate and zero Matrix Buffer
 *  
 *  Input  : int no_col : number of elements in buffer  
 *  Output : pointer to array of DIMENSIONS.
 *  ==================================================
 */ 

#ifdef __STDC__
DIMENSIONS *BufferInit(int no_col)
#else
DIMENSIONS *BufferInit(no_col)
int         no_col;
#endif
{
DIMENSIONS *array;
int          i;
#ifdef DEBUG
       printf("\n Enter BufferInit()\n");
#endif 

    array = (DIMENSIONS *) MyCalloc(no_col, sizeof(DIMENSIONS));

    for( i = 1; i <= no_col; i++)
       ZeroUnits(&array[i-1]);

#ifdef DEBUG
       BufferPrint(" no name ", array, no_col);
       printf("\n Leaving BufferInit()\n");
#endif 

  return (array);
}

/*
 *  ========================================================================
 *  ConvertTempUnits() : Convert temperature units
 *
 *  Input :
 *  Output :
 *  ========================================================================
 */

#ifdef __STDC__
double ConvertTempUnits(char *name, double value,int units_type)
#else
double ConvertTempUnits(name, value, units_type)
char      *name;
double    value;
int  units_type;
#endif
{

#ifdef DEBUG
  printf("**** Enter ConvertTempUnits \n");
  printf("  : name = %s \n", name);
  printf("  : value= %g \n", value);
  printf("  : type = %d \n", units_type);
#endif

  switch(units_type) {
     case SI: /* convert to SI */
          if(name != (char *)NULL) {
          if( !strcmp(name, "deg_C") || !strcmp(name, "deg_F") ) 
              value = 5.0*(value - 32.0)/9.0;
          else {
              printf("*** In ConvertTempUnits(): units_name = %s", name);
              FatalError("Try to convert a non temperature units",
                        (char *) NULL);
          }
          }
          break;
     case US: /* convert to US */
          if(name != (char *)NULL) {
          if( !strcmp(name, "deg_C") || !strcmp(name, "deg_F") )
             value = 9.0*value/5.0 + 32.0;
          else {
             printf("*** In ConvertTempUnits(): units_name = %s", name);
             FatalError("Try to convert a non temperature units",
             (char *) NULL);
          }
          }
          break;
  }
  
  return(value);
}

/*
 *  ========================================================================
 *  UnitsScaleConvert() :
 *
 *  Input :
 *  Output :
 *  ========================================================================
 */

#ifdef __STDC__
DIMENSIONS UnitsScaleConvert(DIMENSIONS eng_units, int units_type) 
#else
DIMENSIONS UnitsScaleConvert(eng_units, units_type) 
DIMENSIONS eng_units;
int         units_type;
#endif
{

 switch(units_type) {
    case SI:
         /* convert from US to SI */
         /* [a] length units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("in", eng_units.units_name))
            eng_units.scale_factor = 0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("ft", eng_units.units_name))
            eng_units.scale_factor = 12.0*0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("mil", eng_units.units_name))
            eng_units.scale_factor = 1E-3*0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("micro_in", eng_units.units_name))
            eng_units.scale_factor = 1E-3*0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("yard", eng_units.units_name))
            eng_units.scale_factor = 36.0*0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("mile", eng_units.units_name))
            eng_units.scale_factor = 63360*0.0254;

         /* [b] volume units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("gallon", eng_units.units_name))
            eng_units.scale_factor = 3.785412E-3;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("barrel", eng_units.units_name))
            eng_units.scale_factor = 0.1589873;

         /* [c] mass units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("lb", eng_units.units_name))
            eng_units.scale_factor = 0.4535924;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("grain", eng_units.units_name))
            eng_units.scale_factor = 0.4535924/7E+3;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("ton", eng_units.units_name))
            eng_units.scale_factor = 2E+3*0.4535924;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("klb", eng_units.units_name))
            eng_units.scale_factor = 1E+3*0.4535924;

         /* [d] force units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("lbf", eng_units.units_name))
            eng_units.scale_factor = 0.4535924*9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("kips", eng_units.units_name))
            eng_units.scale_factor = 1E+3*0.4535924*9.80665;

         /* [e] pressure units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("psi", eng_units.units_name))
            eng_units.scale_factor = 0.4535924*9.80665/0.0254/0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("ksi", eng_units.units_name))
            eng_units.scale_factor = 1E+3*0.4535924*9.80665/0.0254/0.0254;

         /* [f] incremental temperature units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("DEG_F", eng_units.units_name))
            eng_units.scale_factor = 5.0/9.0;

         break;
    case US:
         /* convert from SI to US */
         /* [a] length units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("micron", eng_units.units_name))
            eng_units.scale_factor = 1E-6/0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("mm", eng_units.units_name))
            eng_units.scale_factor = 1E-3/0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("cm", eng_units.units_name))
            eng_units.scale_factor = 1E-2/0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("dm", eng_units.units_name))
            eng_units.scale_factor = 1E-1/0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("m", eng_units.units_name))
            eng_units.scale_factor = 1.0/0.0254;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("km", eng_units.units_name))
            eng_units.scale_factor = 1E+3/0.0254;

         /* [b] mass units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("g", eng_units.units_name))
            eng_units.scale_factor = 1E-3/0.4535924;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("kg", eng_units.units_name))
            eng_units.scale_factor = 1.0/0.4535924;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("Mg", eng_units.units_name))
            eng_units.scale_factor = 1E+3/0.4535924;

         /* [c] force units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("N", eng_units.units_name))
            eng_units.scale_factor = 1.0/0.4535924/9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("kN", eng_units.units_name))
            eng_units.scale_factor = 1E+3/0.4535924/9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("kgf", eng_units.units_name))
            eng_units.scale_factor = 1.0/0.4535924;

         /* [d] pressure units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("Pa", eng_units.units_name))
            eng_units.scale_factor = 0.0254*0.0254/0.4535924/9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("kPa", eng_units.units_name))
            eng_units.scale_factor = 1E+3*0.0254*0.0254/0.4535924/9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("MPa", eng_units.units_name))
            eng_units.scale_factor = 1E+6*0.0254*0.0254/0.4535924/9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("GPa", eng_units.units_name))
            eng_units.scale_factor = 1E+9*0.0254*0.0254/0.4535924/9.80665;

         /* [e] energy units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("Jou", eng_units.units_name))
            eng_units.scale_factor = 1.0/0.4535924/0.0254/9.80665;
         else if(eng_units.units_name != (char *)NULL &&
            !strcmp("kJ", eng_units.units_name))
            eng_units.scale_factor = 1E+3/0.4535924/0.0254/9.80665;

         /* [f] incremental temperature units */

         if(eng_units.units_name != (char *)NULL &&
            !strcmp("DEG_C", eng_units.units_name))
            eng_units.scale_factor = 9.0;

         break;
    default:
         printf(" UNITS_TYPE = %d \n", (int) units_type);
         FatalError("Units_type must be SI or US to use this function",
                   (char *) NULL);
         break;
  }

  return (eng_units);
}

/*
 *  ============================================================
 *  Convert Units Type :
 *  ============================================================
 */

#ifdef __STDC__
DIMENSIONS *UnitsTypeConvert( DIMENSIONS *d, int type )
#else
DIMENSIONS *UnitsTypeConvert( d, type )
DIMENSIONS *d;
int      type;
#endif
{
    if( d->units_type!=type || d->units_type!=SI_US ) {
        d->units_type = type;
        free((char *) d->units_name);
        d->units_name = (char *) NULL;
        UnitsSimplify(d);
    }
    return(d);
}
