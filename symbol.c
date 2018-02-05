/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  symbol.c : Symbol Table : Modified to take (multiple) numbers.
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
#include <string.h>

#include "defs.h"
#include "miscellaneous.h"
#include "units.h"
#include "matrix.h"
#include "fe_database.h"
#include "symbol.h"

#define TABLE_SIZE 101

typedef SYMBOL SymbolNode;

static  SymbolNode *SymbolTable[TABLE_SIZE];
static  SymbolNode *SymbolList;

char	 *malloc(), *calloc();
SymbolNode *AllocSymbolNode();

/* ====================== */
/* (a) Build Symbol table */
/* ====================== */

#ifdef __STDC__
SYMBOL *build_table(char *s, int t, double d)
#else
SYMBOL *build_table(s, t, d)
char  *s;
int    t;
double d;
#endif
{
SYMBOL *sp;

   sp = install(s);
   sp->type    = t;
   sp->u.value = d;

   return (sp);
}

/* =================== */
/* (b) : hash function */
/* =================== */

#ifdef __STDC__
int hash(char *s)
#else
int hash(s)
char	*s;
#endif
{
int mid, len;

    len = strlen(s); 
    mid = (len + 1) / 2;

    return ((s[0] + mid * s[mid - 1] + len * s[len - 1]) % TABLE_SIZE);
}

/* =========================== */
/* (c) : initialize hash table */
/* =========================== */

SymbolInit()
{
SymbolNode *h, *next;
int	   i;

	for (i=0; i < TABLE_SIZE; i++) {
		h = SymbolTable[i];
		while (h != NULL) {
			next = h->next;
			FreeSymbolNode(h);
			h = next;
		}
		SymbolTable[i] = NULL;
	}
}

/* ===================== */
/* (d) : free hash nodes */
/* ===================== */

#ifdef __STDC__
FreeSymbolNode(SymbolNode *h)
#else
FreeSymbolNode(h)
SymbolNode *h;
#endif
{
	h->next = SymbolList;
	SymbolList = h;
}

/* ======================== */
/* (e) : allocate hash node */
/* ======================== */

SYMBOL *AllocSymbolNode()
{
SYMBOL *h;

	while ((h = SymbolList) != NULL) {
		SymbolList = h->next;
	}
	if ((h = (SymbolNode *)calloc(1, sizeof(SymbolNode))) == NULL) {
		fprintf(stderr,"AllocSymbolNode(): malloc failed\n");
		exit(1);
	}
	return h;
}

/* ================================================== */
/* (f) : lookup string name, return NULL if not found */
/* ================================================== */

#ifdef __STDC__
SYMBOL *lookup(char *s)
#else
SYMBOL *lookup(s)
char *s;
#endif
{
SymbolNode  *h;
int hash_val;

#ifdef DEBUG
    printf("****** enter lookup():\n");
#endif

    hash_val = hash(s);

    for( h = SymbolTable[hash_val]; h != NULL; h = h->next) {
     	 if((h->cpSymName != NULL) && strcmp(s, h->cpSymName) == 0) {
	    return h;
         }
    }

#ifdef DEBUG
    printf("Return (SYMBOL *) NULL\n");
    printf("Leaving lookup(): \n");
#endif
    return NULL;
}

/* ================================= */
/* (g) : install node if not present */
/* ================================= */

#ifdef __STDC__
SYMBOL *install(char *s)
#else
SYMBOL *install(s)
char   *s;
#endif
{
SYMBOL *h;
int hash_val;
   
    hash_val = hash(s);

#ifdef DEBUG
    printf("****** In install(): hash_val = %d \n", hash_val);
    printf("****** In install():  s = %s\n", s);
#endif


    /* Allow repetition of numbers in Symbol Table */

    if(strcmp("", s) != 0) {

       if((h = lookup(s)) != NULL)
          return h;
    }


    /* Allocate Symbol Table Node and Name */

    h = AllocSymbolNode();

    if(strcmp("", s) != 0) 

    h->cpSymName = SaveString(s);
    h->next = SymbolTable[hash_val];
    SymbolTable[hash_val] = h;
    
    return h;
}

/* ======================== */
/* (h) : print symbol table */
/* ======================== */

void print_symtable()
{
SYMBOL *h;
int i;

    for(i = 1; i <= TABLE_SIZE; i++) {
        fprintf(stdout," INFO >> i = %3d : ", i);
        for(h = SymbolTable[i-1]; h != NULL; h = h->next)
     	    fprintf(stdout,"%17s ", h->cpSymName);
        fprintf(stdout,"\n");
    }
}
