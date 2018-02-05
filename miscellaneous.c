/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  miscellaneous.c : Miscellaneous Routines
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
 *  Written by: Mark Austin, Xiaoguang Chen, and Wane-Jang Lin      December 1995
 *  ============================================================================= 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#include "miscellaneous.h"

/*
 *  ==========================================================
 *  MyMalloc() : Call malloc() and check for allocation errors
 *  ==========================================================
 */

#ifdef __STDC__
char *MyMalloc( unsigned int uiSize )
#else  /* case not STDC */
char *MyMalloc( uiSize )
unsigned int uiSize;
#endif /* end case not STDC */
{
char *cpTemp;

    cpTemp = (char *) malloc( (unsigned) uiSize );
    if (cpTemp == NULL)
        FatalError("Unable to allocate sufficient memory", (char *) NULL);
    else
        return( cpTemp );
}

/*
 *  ==========================================================
 *  MyCalloc() : Call calloc() and check for allocation errors
 *  ==========================================================
 */

#ifdef __STDC__
char *MyCalloc( unsigned int uiNoItems, unsigned int uiSize)
#else  /* case not STDC */
char *MyCalloc( uiNoItems, uiSize )
unsigned int uiNoItems, uiSize;
#endif /* end case not STDC */
{
char *cpTemp;

    cpTemp = (char *) calloc( (unsigned) uiNoItems, (unsigned) uiSize );
    if (cpTemp == NULL)
        FatalError("Unable to allocate sufficient memory", (char *) NULL);
    else 
        return( cpTemp );
}

/*
 *  =========================================
 *  SaveString() : Allocate Memory for String 
 *  =========================================
 */

#ifdef __STDC__
char *SaveString( char *cpName )
#else  /* case not STDC */
char *SaveString( cpName )
char *cpName;
#endif /* end case not STDC */
{
char *cpTemp;

    if( cpName == (char *)NULL ) {
        cpTemp = (char *)NULL;
        return cpTemp;
    }
    cpTemp = (char *) MyCalloc( strlen(cpName)+1 , sizeof(char) );
    cpTemp = (char *) strcpy(cpTemp, cpName);

    return cpTemp;
}

/*
 *  =============================================================
 *  FatalError : Informs user of a fatal error and aborts program
 *  =============================================================
 */

#ifdef __STDC__

void FatalError( char *cpFirst, ... ) {
va_list  arg_ptr;
char  *cpMessage;

    va_start(arg_ptr, cpFirst );
    for(cpMessage = cpFirst; cpMessage != NULL;
        cpMessage = va_arg(arg_ptr, char *))
        printf("FATAL ERROR >> \"%s\"\n", cpMessage);

    va_end(arg_ptr);

    exit(1);
}

#else  /* Start case not STDC */

void FatalError( va_alist )
va_dcl
{
va_list  arg_ptr;
char  *cpMessage;
char    *cpFirst;

    va_start(arg_ptr);
    cpFirst = va_arg(arg_ptr, char *);

    for(cpMessage = cpFirst; cpMessage != NULL;
        cpMessage = va_arg(arg_ptr, char *))
        printf("FATAL ERROR >> \"%s\"\n", cpMessage);

    va_end(arg_ptr);

    exit(1);
}

#endif /* End case not STDC */
