/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  main.c : Check the input type
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

#include <stdio.h>       
#include <ctype.h>		
#include <signal.h>
#include <setjmp.h>	

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"
#include "code.h"

/* Declarations global frame and working element data structures */

ARRAY     *array;
EFRAME    *frame;

jmp_buf  begin;				/* query this line.	*/

char *progname;

int        lineno = 1;
char   *infile = NULL;
FILE      *fin, *finp;
char          **gargv;
int             gargc;

main(argc, argv)			
char *argv[];
{
int    c;
int    a=0, b=0;
void   fpecatch();
void   set_default_values();
void   set_print_output();
char **p;

   /* [a] : codes for running on FreeBSD platform */

#ifdef __FreeBSD__
   fpsetmask (0);
#endif

   /* [b] : Read flags and options from command input */

   progname = argv[0];
   fin      = stdin;
   finp     = stdin;

   /* [c] : Provide help when command flags are missing */

   if (argc == 1) { 
       ProgramUsage();
       exit(1);
   }

   /* [d] : Read and process command flags */

   while(--argc > 0 && ((*++argv) != NULL)) 
      if((*argv)[0] == '-') 
         while((*argv != NULL) && (c = *++argv[0]) != 0) 
            switch(c) {
	        case 'f':
                     a = c;
                     if(infile == NULL) {

                        if((*(argv-1) != NULL && strcmp(*(argv-1), progname))
                            && ((*(argv-1))[0] != '-'))
                            infile = SaveString(*(argv-1));
                        else {
                            p = argv;
                            while(*++p != NULL)
                               if((*p)[0] != '-')
                                  infile = SaveString(*p);
                         
                           if((*(argv-1) == NULL) && (*p == NULL))
                               FatalError("Specify inputfile name",(char*)NULL);
                        }
                   
		        Open_File(infile);
                     }
		     break;
                case 'k':
                     a = c;
                     finp = fopen("inputfile.std", "w");
                     printf("ALADDIN : ");
                     fflush(stdout);
                     break;
                case 's':
                     b = c;
                     if(a == 0 && infile == NULL) {
                        if((*(argv-1) != NULL && strcmp(*(argv-1), progname))
                            && ((*(argv-1))[0] != '-'))
                            infile = SaveString(*(argv-1));
                        else {
                            p = argv;
                            while(*++p != NULL)
                                 if((*p)[0] != '-')
                                    infile = SaveString(*p);

                        if((*(argv-1) == NULL) && (*p == NULL))
                            FatalError("Specify inputfile name",(char*)NULL);
                        }

                        if(infile != NULL)
		           Open_File(infile);
                     }
		     break;
                case 'h':
                     ProgramUsage();
                     exit(1);
		     break;
	        default:
                     break;
	   }

   /* [e] : Set Output Flags */

      set_default_values();
      set_print_output();

   /* [f] : Load grammar, materials and AISC sections into symbol table */

      Init_Problem();

   /* [g] : Run Problem */
      
      if(b == 's') ScanInput();
      Run(a, b);
}

/*
 *  ================================================
 *  ProgramUsage() : Provide hints on using ALADDIN.
 *  ================================================
 */

ProgramUsage() {
   printf("Usage : ALADDIN -k                 : Input from keyboard            \n");
   printf("        ALADDIN -f <input-file>    : Input from file                \n");
   printf("        ALADDIN -s <input-file>    : Scan input for correct grammar,\n");
   printf("                                   : but do not process.            \n");
}

/*
 *  ================================================
 *  Open_File() : Open Input file for reading
 *  ================================================
 */

#ifdef __STDC_
Open_File( char *filename )
#else
Open_File(filename)
char *filename;
#endif
{
   fin = fopen(filename, "r");
   if(fin == NULL)
      FatalError("Cannot open input file",(char *)NULL);
}

/*
 *  ================================================
 *  ScanInput() : Scan Input for grammatical errors.
 *  ================================================
 */

ScanInput() {
   Init_Code();
   while(yyparse() != 0) {
         Init_Code();
   }
  fclose(fin);
}

/*
 *  ================================================
 *  Run() : Run ALADDIN Program
 *  ================================================
 */


#ifdef __STDC_
Run(int c, int b)
#else
Run(c, b)
int c, b;
#endif
{
void fpecatch();

   setjmp(begin);			
   signal(SIGFPE, fpecatch);

   if(c == 'f')
      Open_File(infile);
   else if (c == 'k' && b == 's') 
      Open_File("inputfile.std");

   Init_Code();
   while(yyparse() != 0) {
         Execute(progbase);
         Init_Code();
   }
}
