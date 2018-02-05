/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  units.h : Data structures for engineering quantities
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
 *  Written by: Mark Austin and Xiaoguang Chen                           May 1997
 *  ============================================================================= 
 */

#ifndef UNITS_H
#define UNITS_H

/* Units Type */

#define  SI     100
#define  US     200
#define  SI_US  300

/* Variables Defined for Units Functions */

#define COLUMN       1
#define ROW          2
#define COL_UNITS    3
#define ROW_UNITS    4

/* Engineering Quantities */

typedef struct dimensional_exponents {
        char         *units_name; /* units name              */
        double      scale_factor; /* scale/conversion factor */
        double      length_expnt; /* length                  */
        double        mass_expnt; /* mass                    */
        double        time_expnt; /* time                    */
        double        temp_expnt; /* temperature             */
        double      radian_expnt; /* planar angle            */
        int           units_type; /* US or SI units          */
} DIMENSIONS;

typedef struct engineering_quantity {
        double       value;
        DIMENSIONS  *dimen;
} QUANTITY;

/* Units functions */

#ifdef __STDC__

int            SameUnits( DIMENSIONS *, DIMENSIONS * );
DIMENSIONS    *UnitsMult( DIMENSIONS *, DIMENSIONS * );
DIMENSIONS    *UnitsMultRep( DIMENSIONS *, DIMENSIONS *, DIMENSIONS * );
DIMENSIONS    *UnitsDiv( DIMENSIONS *, DIMENSIONS *, int );
DIMENSIONS    *UnitsDivRep( DIMENSIONS *, DIMENSIONS *, DIMENSIONS *, int );
DIMENSIONS    *UnitsPower( DIMENSIONS *, double , int );
DIMENSIONS    *UnitsPowerRep( DIMENSIONS *, DIMENSIONS *, double , int );
DIMENSIONS    *UnitsCopy( DIMENSIONS *, DIMENSIONS * );
DIMENSIONS    *ZeroUnits( DIMENSIONS * );
DIMENSIONS    *DefaultUnits( char * );
DIMENSIONS    *UnitsSimplify( DIMENSIONS * );
DIMENSIONS    *RadUnitsSimplify( DIMENSIONS * );

void           UnitsPrint( DIMENSIONS * );
int            UnitsLength( char *, ... );

DIMENSIONS    *UnitsTypeConvert( DIMENSIONS * , int );
DIMENSIONS     UnitsScaleConvert( DIMENSIONS , int );
double         ConvertTempUnits( char *, double, int );

DIMENSIONS    *BufferInit( int );
void           BufferPrint( char *, DIMENSIONS *, int );

int       SetUnitsOn();
int       SetUnitsOff();
int       CheckUnits();
int       CheckUnitsType();

QUANTITY          *QuantityUnitsLess( QUANTITY * );

#else  /* start case not STDC */

int            SameUnits();
DIMENSIONS    *UnitsMult();
DIMENSIONS    *UnitsMultRep();
DIMENSIONS    *UnitsDiv();
DIMENSIONS    *UnitsDivRep();
DIMENSIONS    *UnitsPower();
DIMENSIONS    *UnitsPowerRep();
DIMENSIONS    *UnitsCopy();
DIMENSIONS    *ZeroUnits();
DIMENSIONS    *DefaultUnits();
DIMENSIONS    *UnitsSimplify();
DIMENSIONS    *RadUnitsSimplify();

void           UnitsPrint();
int            UnitsLength();

DIMENSIONS    *UnitsTypeConvert();
DIMENSIONS     UnitsScaleConvert();
double         ConvertTempUnits();

DIMENSIONS    *BufferInit();
void           BufferPrint();

int       SetUnitsOn();
int       SetUnitsOff();
int       CheckUnits();
int       CheckUnitsType();

QUANTITY *QuantityUnitsLess();

#endif /* end case not STDC */
#endif /* end case UNITS_H */
