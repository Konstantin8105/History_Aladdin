/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_setflags.c : Finite Element Preprocessor & Base Module
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

#include "defs.h"
#include "units.h"

unsigned int PRINT_PROFILE;
unsigned int PRINT_PLINK;
unsigned int PRINT_MAP_DOF;

/* -------------------------------------------------------- */
/* set_default_values();                                    */
/* function to set default values for problem/analysis      */
/* called in fera_preprocessor()                            */
/*                                                          */
/* These parameters are set as default values for analysis  */
/* and can be overridden by  values given in input file     */
/* -------------------------------------------------------- */

void set_default_values() 
{
	/* set default units type to SI */
	ChangeUnitsType( 2 );
}

/* -------------------- */
/* set_print_output()   */
/* -------------------- */

void set_print_output() {

     PRINT_PROFILE          = ON;
     PRINT_PLINK            = ON;
     PRINT_MAP_DOF          = OFF;
}
