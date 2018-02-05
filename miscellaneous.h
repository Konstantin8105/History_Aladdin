/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  miscellaneous.h : Header file for Miscellaneous Routines
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
 *  Written by: Mark Austin                                         December 1995
 *  ============================================================================= 
 */

#ifndef MISCELLANEOUS_H
#define MISCELLANEOUS_H

/* Miscellaneous Functions */

#ifdef __STDC__

char   *MyMalloc( unsigned int );
char   *MyCalloc( unsigned int, unsigned int );
char   *SaveString( char * );
void    FatalError( char * , ... );

#else  /* case not STDC */

char    *MyMalloc();
char    *MyCalloc();
char    *SaveString();
void    FatalError();

#endif /* end case STDC */

#endif /* end case MISCELLANEOUS_H */
