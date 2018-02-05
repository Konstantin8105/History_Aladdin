/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  grammar.y : YACC grammar
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

%{
#include  <stdio.h>
#include  "defs.h"
#include  "units.h"
#include  "matrix.h"
#include  "fe_database.h"
#include  "symbol.h"
#include  "code.h"
	
#define	  Code2(c1,c2)     Code(c1); Code(c2)
#define   Code3(c1,c2,c3)  Code(c1); Code(c2); Code(c3)
#define   Code5(c1,c2,c3,c4,c5)  Code(c1); Code(c2); Code(c3); Code(c4); Code(c5)


%}

%union {
	SYMBOL     *sym;
	Inst   	  *inst;
	int	   narg;
}

%token    	QUIT END_OF_FILE BREAK 
%token          SET_UNITS_OFF     SET_UNITS_ON	

%token  <sym>   SECTION   	MATERIAL	TYPE	FIBER 

%token  <sym>  	WHILE IF ELSE THEN FOR  NUMBER STRING PRINT UNDEF 
%token  <sym>  	BLTINVar_MATRIX BLTINVar_QUANTITY 
%token  <sym>  	BLTIN_MATRIX	BLTIN1_MATRIX	BLTIN2_MATRIX	BLTIN3_MATRIX
%token  <sym>  	BLTIN_QUANTITY	BLTIN1_QUANTITY BLTIN2_QUANTITY	BLTIN3_QUANTITY	

%token  <sym>  	NODE_QUANT      MESH
%token  <sym>  	LINK_NODE	ADD_ELMT
%token  <sym>  	SECT_ATTR	ELMT_ATTR	MATL_ATTR	FIB_ATTR
%token  <sym>   UNIT
%token  <sym>	MAP  	LDOF 	TO 	GDOF
%token  <sym>	BLTIN_FE_FUNC	BLTIN1_FE_FUNC

%token  <sym>  	VAR		QUAN 	VECT	 	MATX

%token  <sym>  	DIMENSION
%token  <sym>  	LENGTH	MASS	GTIME

%type   <inst>  conditional_quantity conditional_list quantity_dimen
%type   <inst>  quantity  quantity_assign  stmt prlist  stmtlist  exprlist listlist  object 
%type   <inst>  matrix matrix_assign  dimensions  m_to_quantity
%type   <inst>  cond  while  if  end  for 
%type   <inst>  finite_elmt
%type   <inst>  element_property 
%type   <narg>  matrix_seq, matrix_seq_var
%type   <narg>  row_seq dimen_row quantity_row 
%type   <narg>  elmt_object_seq
%type   <narg>  object_seq

%right   '='
%left    OR
%left    AND
%left    GT GE LT LE EQ NE
%left    '+'   '-'	
%left    '*'  '/'  
%left    '%'
%left    UNARYPLUS UNARYMINUS NOT
%right	'^'	
%%

list:	                {;}
    | list '\n'  
    | list QUIT  ';'    { if(fin == stdin && finp != stdin) 
                             { fprintf(finp,"\n"); fclose(finp); }
                          return 0; }
    | list END_OF_FILE  { return 0; }
    | list stmt          { Code(STOP); return 1; }
    | list object        { Code(STOP); return 1; }
    | list error         {  yyerrok; }
    ;

stmt :  matrix       ';' { Code( Pop_Eval ); }
     |  quantity     ';' { Code( Pop_Eval ); }
     |  finite_elmt  ';' { Code( Pop_Eval ); }
     |  PRINT prlist ';' { $$ = $2; }
     |  while '(' cond ')'  stmt end {
                    ($1)[1] = (Inst) $5;
		    ($1)[2] = (Inst) $6;
		    }
     |  if '(' cond  ')' THEN stmt end ELSE stmt end {
		    ($1)[1] = (Inst) $6;
		    ($1)[2] = (Inst) $9;
		    ($1)[3] = (Inst) $10;
		    }
     |  if '(' cond ')' stmt end {
		    ($1)[1] = (Inst) $5;
		    ($1)[3] = (Inst) $6;
		    }
     | for '(' exprlist ';' cond ';' exprlist ')' stmt end {
                    ($1)[1] = (Inst) $3;  /* initiation */
                    ($1)[2] = (Inst) $5;  /* condition  */
                    ($1)[3] = (Inst) $9;  /* body of loop */
                    ($1)[4] = (Inst) $7;  /* increments */
                    ($1)[5] = (Inst) $10; /* end, if cond fails */
                    }
     |  '{' stmtlist '}'  {  $$ = $2; }
     |  BREAK ';'         {  $$ = Code(Bltin_Break); }
     |  SET_UNITS_ON';'   {  $$ = Code(SetUnitsOn); }
     |  SET_UNITS_OFF';'  {  $$ = Code(SetUnitsOff); }
     ;

stmtlist:                 { $$ = progp; } 
        | stmtlist stmt
        ;

cond:      conditional_list { Code(STOP); $$ = $1; }
    ;

conditional_list:                      {$$ = progp; Code(Quantity_Yes);}
                | conditional_quantity {$$ = $1;}
                ;

exprlist:      listlist      {Code(STOP); $$ = $1; } 
        ;

listlist :                       { $$ = progp; } 
         | quantity              {Code(Pop_Eval); $$ = $1;}
         | listlist ',' quantity {Code(Pop_Eval); $$ = $1;} 
         ;

while: WHILE { $$ = Code3(While_Code, STOP, STOP);}
     ;

if:    IF    { $$ = Code(If_Code);  Code3(STOP, STOP, STOP);}
  ;

for:  FOR { $$ = Code5(For_Code, STOP, STOP, STOP, STOP); Code(STOP);}
   ;

end:  { Code(STOP); $$ = progp; } 
   ;

prlist: quantity_dimen 
      | STRING                                { Code2(Print_String, (Inst)$1);}
      | prlist ',' quantity_dimen    
      | prlist ',' STRING                     { Code2(Print_String, (Inst)$3);}
      ;	

quantity_dimen: quantity                    { Code(Print_Expr);  }
              | quantity '(' dimensions ')' { Code(Print_Dimen_Expr);}
      ;

/*   -------------------------------------------------
 *   Finite Element Functions, Objects, and Quantities
 *   ------------------------------------------------- */

finite_elmt : NODE_QUANT '(' quantity ',' matrix ')' {
                    Code2( Bltin_Node_Quant, (Inst) $1->u.voidptr); 
                    }
            | LINK_NODE '(' matrix ',' matrix ')'  {
                    Code2( Bltin_Link_Node, (Inst)$1->u.voidptr);
                    }
            | ADD_ELMT  '(' quantity ',' matrix ',' STRING ')' {
                    Code3( Push_String,    (Inst) $7, String_Eval);
                    Code2( Bltin_Add_Elmt, (Inst) $1->u.voidptr);
                    }
            | MESH '(' ')' {
                    Code2( Bltin_Mesh, (Inst) $1->u.voidptr);
                    }
            | BLTIN_FE_FUNC '(' ')' {
                    Code2( Bltin_Fe_Function, (Inst) $1->u.voidptr);
                    }
            | BLTIN1_FE_FUNC  '(' matrix ')' {
                    Code2( Bltin1_Fe_Function, (Inst) $1->u.voidptr);
                    }
            ;

/* Declaration of Element and Section Objects and Units Type */

object : MATL_ATTR '(' STRING ')' '{' object_seq '}' {
                    Code2( Push_Variable, (Inst) $6);
                    Code3( Push_String,   (Inst) $3, String_Eval);
                    Code( Bltin_Material_Attr);
		    }
       | SECT_ATTR '(' STRING ')' '{' object_seq '}' {
                    Code2( Push_Variable, (Inst) $6);
                    Code3( Push_String,   (Inst) $3, String_Eval);
                    Code( Bltin_Section_Attr);
		    }
       | ELMT_ATTR '(' STRING ')' '{' elmt_object_seq '}' {
                    Code2( Push_Variable, (Inst) $6);
                    Code3( Push_String,   (Inst) $3, String_Eval);
                    Code( Bltin_Element_Attr);
		    }
       | FIB_ATTR '(' VAR ',' STRING ')' '{' object_seq '}' {
                    Code2( Push_Variable, (Inst) $8);
                    Code3( Push_String,   (Inst) $5, String_Eval);
                    Code3( Push_Variable, (Inst)$3, Variable_Eval);
                    Code(  Bltin_Fiber_Attr );
		    }
       | FIB_ATTR '(' NUMBER ',' STRING ')' '{' object_seq '}' {
                    Code2( Push_Variable, (Inst) $8);
                    Code3( Push_String,   (Inst) $5, String_Eval);
                    Code( Push_Dimensionless );
                    Code3( Push_Constant, (Inst)$3, Dimension_Eval);
                    Code(  Bltin_Fiber_Attr );
		    }
       | UNIT '(' STRING ')' ';'  {
                    Code3( Push_String, (Inst) $3, String_Eval );
                    Code( Bltin_Units_Type );
                    }
       ;

elmt_object_seq:                                  { $$ = 0; }
          | element_property ';'                  { $$ = 1; }
          | elmt_object_seq element_property ';'  { $$ = $1 + 1;} 
          ;

object_seq:                                  { $$ = 0; }
          | VAR '=' quantity ';'             { Code2( Push_Variable, (Inst)$1);
                                               $$ = 1; }
          | VAR '=' matrix ';'               { Code2( Push_Variable, (Inst)$1);
                                               $$ = 1; }
          | TYPE '=' STRING  ';'             { Code2(Push_String,   (Inst)$3);
                                               Code2(Push_Variable, (Inst)$1); 
                                               $$ = 1; }
          | object_seq VAR '=' quantity ';'  {Code2( Push_Variable, (Inst)$2);
                                               $$ = $1 + 1;} 
          | object_seq VAR '=' matrix ';'    {Code2( Push_Variable, (Inst)$2);
                                               $$ = $1 + 1;} 
          | object_seq TYPE '=' STRING   ';' { Code2(Push_String,   (Inst)$4); 
                                               Code2(Push_Variable, (Inst)$2); 
                                               $$ = $1 + 1;} 
          ;

element_property : TYPE     '=' STRING   { Code2(Push_String,   (Inst)$3);
                                           Code2(Push_Variable, (Inst)$1); }
                 | SECTION  '=' STRING   { Code2(Push_String,   (Inst)$3);
                                           Code2(Push_Variable, (Inst)$1); }
                 | MATERIAL '=' STRING   { Code2(Push_String,   (Inst)$3);
                                           Code2(Push_Variable, (Inst)$1); }
                 | FIBER    '=' STRING   { Code2(Push_String,   (Inst)$3);
                                           Code2(Push_Variable, (Inst)$1); }
                 | MAP LDOF matrix TO GDOF matrix {;}
                 ;

/* ---------------------------
 * Matrix and Matrix Functions
 * --------------------------- */

matrix : MATX                             { $$ = Code3( Push_Matrix,   (Inst)$1, Matrix_Eval); }
       | BLTIN_MATRIX  '(' ')'            { Code2( Bltin_Matrix,  (Inst)$1->u.matrixptr); }
       | BLTIN1_MATRIX '(' matrix ')'     { Code2( Bltin1_Matrix, (Inst)$1->u.matrixptr); }
       | BLTIN2_MATRIX '(' matrix ',' matrix ')'
					  { Code2( Bltin2_Matrix, (Inst)$1->u.matrixptr);}
       | BLTIN3_MATRIX '(' matrix ',' matrix ',' matrix ')'
					  { Code2( Bltin3_Matrix, (Inst)$1->u.matrixptr);}
       | BLTINVar_MATRIX '(' matrix_seq_var ')' 
                                          {Code2(Push_Variable, (Inst)$3);
                                           $$ = Code2(BltinVar_Matrix, (Inst)$1->u.matrixptr); }
       | '[' matrix_seq ']'               { $$ = Code3( Push_Variable, (Inst)$2, Matrix_Build);}
       | matrix '*' matrix                {  Code( Bltin_Matrix_Mult);   }
       | matrix '+' matrix                {  Code( Bltin_Matrix_Add);    }
       | matrix '-' matrix                {  Code( Bltin_Matrix_Sub);    } 
       | '-' matrix %prec  UNARYMINUS     {  Code( Bltin_Matrix_Negate); }
       | '+' matrix %prec  UNARYPLUS      {  Code( Bltin_Matrix_Affirm); }
       | '(' matrix ')'                   {  $$ = $2; }
       | matrix_assign                    { ; } 
       | quantity '*' matrix              {  Code( Bltin_Quan_Matrix_Mult); $$ = $3;}
       | matrix '*' quantity              {  Code( Bltin_Matrix_Quan_Mult); }  
       | matrix '/' quantity              {  Code( Bltin_Matrix_Quan_Div); }  
       | matrix '^' quantity              {  Code( Bltin_Matrix_Power); }  
       ;

matrix_assign : VAR    '=' matrix { $1->type = MATX;
                                    $1->u.m = (MATRIX *) NULL;
                                    Code3( Push_Matrix, (Inst)$1, Assign_Matrix); }
              | MATX   '=' matrix { Code3( Push_Matrix, (Inst)$1, Assign_Matrix); }
              ;

matrix_seq: row_seq                { Code2( Push_Variable,(Inst)$1); $$ =1; }
          | matrix_seq ';' row_seq { Code2( Push_Variable,(Inst)$3); $$= $1 + 1;}
          ;

row_seq: quantity_row  {$$ = $1;}
       | dimen_row     {$$ = $1;}         
       ;

quantity_row:                                    { $$ = 0; }
            | quantity                           { $$ = 1; }
            | quantity_row ',' quantity          { $$ = $1 + 1;} 
            ;

dimen_row :  dimensions                       { $$ = 1; }
          |  dimen_row ',' dimensions         { $$ = $1 + 1;} 
          ;

matrix_seq_var:                           { $$ = 0; }
              | matrix                    { $$ = 1; }
              | matrix_seq_var ',' matrix { $$ = $1 + 1;} 
              ;


/*  ----------------------
 *  Engineering Quantities 
 *  ---------------------- */

quantity:  NUMBER                           { $$ = Code(Push_Dimensionless); 
                                              Code3(Push_Constant, (Inst)$1, Dimension_Eval); }
        |  NUMBER dimensions                { Code3(Push_Constant, (Inst)$1, Dimension_Eval); $$ = $2;}
        |  VAR                              { $$ =  Code3(Push_Variable, (Inst)$1,  Variable_Eval); }
        |  VAR dimensions                   { Code3(Push_Variable, (Inst)$1, Variable_Eval);
                                              Code(Dimension_Eval); $$ = $2; }
        |  BLTIN_QUANTITY '(' ')'           { Code2( Bltin_Quantity, (Inst)$1->u.quantityptr); }
        |  BLTIN1_QUANTITY '(' quantity ')' { Code2( Bltin1_Quantity, (Inst)$1->u.quantityptr); }
        |  BLTIN2_QUANTITY '(' quantity ','
                               quantity ')' { Code2( Bltin2_Quantity, (Inst)$1->u.quantityptr); }
        |  BLTIN3_QUANTITY '(' matrix   ')' { Code2( Bltin3_Quantity, (Inst)$1->u.quantityptr); }
        |  BLTINVar_QUANTITY '('quantity_row ')'
                                            { Code2(Push_Variable, (Inst)$3);
                                              $$ = Code2(BltinVar_Quantity, (Inst)$1->u.quantityptr); }
	|  '-' quantity  %prec  UNARYMINUS  { Code(Quantity_Negate); }
	|  '+' quantity  %prec  UNARYPLUS   { Code(Quantity_Affirm); }
	|  m_to_quantity	            { $$ = $1; }
        |  quantity_assign                  { ;} 
        |  conditional_quantity             { $$ = $1;} 
        |  '(' quantity ')'                 { $$ = $2; }
        |  quantity '+' quantity            { Code(Quantity_Add); }
        |  quantity '-' quantity            { Code(Quantity_Sub); }
        |  quantity '*' quantity            { Code(Quantity_Mul); }
        |  quantity '/' quantity            { Code(Quantity_Div); }
        |  quantity '^' quantity            { Code(Quantity_Power); }
        |  quantity '%' quantity            { Code(Quantity_Mod); }
        ;

conditional_quantity:  quantity GT  quantity  { Code(Quantity_Gt);   }
	            |  quantity GE  quantity  { Code(Quantity_Ge);   }
	            |  quantity LT  quantity  { Code(Quantity_Lt);   }
	            |  quantity LE  quantity  { Code(Quantity_Le);   }
	            |  quantity EQ  quantity  { Code(Quantity_Eq);   }
	            |  quantity NE  quantity  { Code(Quantity_Ne);   }
	            |  quantity AND quantity  { Code(Quantity_And);  }
	            |  quantity OR  quantity  { Code(Quantity_Or);   }
	            |  NOT quantity           { Code(Quantity_Not); $$ = $2; }
	            ;	

quantity_assign : VAR '=' quantity  {  Code2( Push_Variable, (Inst)$1);
                                       Code( Assign_Quantity); $$ = $3; }
                | MATX '[' quantity ']' '[' quantity ']' '=' quantity { 
                                       Code2( Push_Matrix, (Inst)$1);
                                       Code( Assign_Matrix_Item ); $$ = $9;}
	        ;	

m_to_quantity : MATX '[' quantity ']' '[' quantity ']' {
                                       Code3( Push_Matrix, (Inst)$1, Matrix_Eval);
                                        Code( Quantity_Extract); }
              ;

dimensions: DIMENSION                 { $$ = Code2( Push_Dimension, (Inst)$1); }
	  | dimensions '*' dimensions {  Code( Dimension_Mult ); }
	  | dimensions '/' dimensions {  Code( Dimension_Div); }
	  | dimensions '^' quantity   {  $$ = $1; Code( Dimension_Power); }
	  ;	
%%


/*
 *  ------------------------------
 *  Lexical Analysis
 *  ------------------------------
 */

#include <stdio.h>       
#include <ctype.h>		
#include <signal.h>
#include <setjmp.h>

char   *progname;
extern int      lineno;
extern char    *infile;
extern FILE       *fin;
extern FILE      *finp;
extern jmp_buf   begin;

int indef, c;

int  backslash();
int  follow();
int  yyerror();
void warning();
void EXECUTION_ERROR();
void fpecatch();
void Defn_Only();

yylex()
{
int finished;

  /* [a] : Eat blank space, tab, and newline input -- read comment statements */
       
     while((c = fgetc(fin)) == ' ' || c == '\t' || c == '\n' || c == '/') {

          if(c == '\n') 
             lineno++;

          if(fin == stdin && finp != stdin) {
             if(c == '\n') {
                printf("ALADDIN : ");
                fflush(stdout);
             }
             fputc(c, finp);
          }

          if(c == '/') {
             if(follow('*', TRUE, FALSE) == TRUE) {
                c = fgetc(fin);
                finished = FALSE;
                while(finished == FALSE) {
                    c = fgetc(fin);

                    if(c == '\n') 
                       lineno++;

                    if(c == EOF) 
                       FatalError("ERROR >> ... file ended with unbalanced comment !!\n",
                                  (char *)NULL);
                    if((c == '*') && (follow('/', TRUE, FALSE) == TRUE)) {
                       finished = TRUE;
                    }
               }
             } else {
               return c;
             } 
          }
     }

  /* [b] : return end-of-file */

     if(c == EOF) {
        return END_OF_FILE;
     }

  /* [c] : read and store floating-point numbers */

     if (c == '.' || isdigit(c)) {
 	 double value;
	 ungetc(c, fin);
	 fscanf(fin, "%lf",  &value);

         if(fin == stdin && finp != stdin) fprintf(finp, "%lf", value);

	 yylval.sym = build_table("", NUMBER, value);
	 return NUMBER;
     }

  /* [d] : read and store variable */

     if (isalpha(c) || c == '_') {
 	 SYMBOL *s;
	 char	sbuf[100], *p = sbuf;

	 do {	
	    if(p >= sbuf + sizeof(sbuf) - 1) {
	       p = '\0';
	       EXECUTION_ERROR("ERROR >> name too long");
	    }
	    *p++ = c;
	 } while ((( c = fgetc(fin)) != EOF) && (isalnum(c) || (c == '_' )));

         ungetc(c, fin);
	 *p = '\0';

         if(fin == stdin && finp != stdin) fprintf(finp, "%s", sbuf);

	 if ((s = lookup(sbuf)) == 0) {
	      s = build_table(sbuf, VAR, 0.0);
	 }
	 yylval.sym = s;

         switch((int) s->type) {
 	    case  UNDEF:
            case  VAR:
            case  QUAN:
                  return VAR;
            default:
                  return (s->type);
         }
     }

  /* [e] : Get Quoted String */

     if (c == '"') {
         char sbuf[100], *p;
         if(fin == stdin && finp != stdin) fputc(c, finp);

	 for(p = sbuf; (c = fgetc(fin)) != '"'; p++) {
             if(c == '\n' || c == EOF)
                EXECUTION_ERROR("ERROR >> missing quote");
             if(p >= sbuf + sizeof(sbuf) - 1) {
                p = '\0';
                EXECUTION_ERROR("ERROR >> string name too long");
             }
             if((fin == stdin && finp != stdin))
                 fputc(c, finp);

		*p = backslash(c);
	 }
	 *p = 0;

         yylval.sym = (SYMBOL *) SaveString(sbuf);

         if(fin == stdin && finp != stdin) fputc(c, finp);

	 return STRING;
    }

    if(fin == stdin && finp != stdin) {
       fputc(c, finp);
    }

    switch(c) {
        case '>': 
             return follow('=',  GE,  GT);
        case '<':
             return follow('=',  LE,  LT);
        case '=':
             return follow('=',  EQ, '=');
        case '!':
             return follow('=',  NE, NOT);
        case '|':
             return follow('|',  OR, '|');
        case '&':
             return follow('&', AND, '&');
        case '\n':
             lineno = lineno + 1;
             return '\n';
        default:
	     return c;
    }
}

/* Error Traceback Routines */

int yyerror(s)
char *s;
{
     warning(s, (char *) 0);
}

int backslash(c)
int c;
{
char *strchr();
static char transtab[] = "b\bf\fn\nr\rt\t";

    if(c != '\\')
       return c;

    c = fgetc(fin);

    if(fin == stdin && finp != stdin) fputc(c, finp);

    if(islower(c) && strchr(transtab, c))
       return strchr(transtab, c)[1];

    return c;
}

int follow(expect, ifyes, ifno)
int expect, ifyes, ifno;
{
int c = fgetc(fin);

    if(c == expect) {
       if(fin == stdin && finp != stdin && ( c != '*' && c != '/') )  
          fputc(c, finp);

       return ifyes;
    }

    ungetc(c, fin);
    return ifno;
}

void Defn_Only(s)
char *s;
{
   if(!indef)
      EXECUTION_ERROR("ERROR >> Variable used outside definition", (char *) 0);
}

void fpecatch()
{
     EXECUTION_ERROR("Floating point exception", (char *) 0);
}

void EXECUTION_ERROR(s, t)
char *s, *t;
{
	warning(s, t);
	fseek(fin,0L,2);
	longjmp(begin,0);
}

void warning(s, t)
char *s, *t;
{
	fprintf(stderr, "%s: %s", progname, s);
	if(t)
	  fprintf(stderr,"ERROR >> %s",t);
	if(infile)
	  fprintf(stderr," in file '%s' ",infile);

        fprintf(stderr,"near line %d\n",lineno);
        FatalError("In Warning()",(char *)NULL);

        /* flush rest on file input */
	while(c == '\n' && c != EOF)
	      c = fgetc(fin);

	if(c == '\n')
	   lineno = lineno + 1;
}

