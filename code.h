/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  code.h : Declarations for Stack Machine
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

#ifndef  CODE_H
#define  CODE_H 

typedef union Datum  {
        QUANTITY    *q;
        MATRIX      *m;
	SYMBOL    *sym;
} DATUM;

typedef int (*Inst)();
#define STOP  (Inst) 0

extern Inst *progp, *progbase;
extern Inst  prog[];

#ifdef __STDC__
Inst  *Code( Inst );
int    Execute( Inst * );
int    Push( DATUM );
#else
Inst  *Code();
int    Execute();
int    Push();
#endif

DATUM   Pop();
int     Pop_Eval();
int     Init_Code();

/* Finite State Machine Routines */

int     If_Code();
int     For_Code();
int     While_Code();
int     Call();
int     Arg();
int     Arg_Assign();
int     Proc_Ret();
int     Func_Ret();

/* External Function Declarations */

int     Bltin_Break();
int     Check_Break();
int     After_Break();

int     Variable_Eval();
int     String_Eval();
int     Matrix_Eval();
int     Dimension_Eval();

int     Matrix_Build();
int     Bltin_Matrix();
int     Bltin1_Matrix();
int     Bltin2_Matrix();
int     Bltin3_Matrix();
int     BltinVar_Matrix();

int     Bltin_Quantity();
int     Bltin1_Quantity();
int     Bltin2_Quantity();
int     Bltin3_Quantity();
int     BltinVar_Quantity();

int     BltinVar_String();
int     BltinVar_Void();

int     Assign_Quantity();
int     Assign_Matrix();
int     Assign_Matrix_Item();

int     Push_Variable();
int     Push_Constant();
int     Push_String();
int     Push_Matrix();
int     Push_Argument();

int     Print_Expr();
int     Print_Dimen_Expr();
int     Print_String();

/* Engineering Dimensions/Units */

int     Push_Dimension();
int     Push_Dimensionless();
int     Dimension_Mult();
int     Dimension_Div();
int     Dimension_Div2();
int     Dimension_Power();

/* Engineering Quantities */

int     Quantity_Add();
int     Quantity_Sub();
int     Quantity_Mul();
int     Quantity_Div();
int     Quantity_Affirm();
int     Quantity_Negate();
int     Quantity_Power();
int     Quantity_Mod();
int     Quantity_Gt();
int     Quantity_Lt();
int     Quantity_Eq();
int     Quantity_Ge();
int     Quantity_Le();
int     Quantity_Ne();
int     Quantity_And();
int     Quantity_Or();
int     Quantity_Not();
int     Quantity_Yes();
int     Quantity_Extract();

/* Matrix Functions */

int     Bltin_Matrix_Add();
int     Bltin_Matrix_Sub();
int     Bltin_Matrix_Mult();
int     Bltin_Matrix_Power();
int     Bltin_Matrix_Trans();
int     Bltin_Quan_Matrix_Mult();
int     Bltin_Matrix_Quan_Mult();
int     Bltin_Matrix_Quan_Div();
int     Bltin_Matrix_Negate();
int     Bltin_Matrix_Affirm();

/* Finite Element Functions */

int      Bltin_Node_Quant();
int      Bltin_Mesh();
int      Bltin_Link_Node();
int      Bltin_Add_Elmt();

int      Bltin_Element_Attr();
int      Bltin_Section_Attr();
int      Bltin_Material_Attr();
int      Bltin_Fiber_Attr();

int      Bltin_Fe_Function();
int      Bltin1_Fe_Function();

int      Bltin_Units_Type();

#endif /* end case CODE_H */
