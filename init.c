/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  init.c : Initialize symbol table
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

#include <math.h>
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
#include "code.h"
#include "fe_functions.h"
#include "y.tab.h"

QUANTITY                *Sqrt(), *Fabs(), *Pow();
QUANTITY    *Log(), *Log10(), *Exp(), *Integer();
QUANTITY         *Sin(), *Cos(), *Tan(), *Atan();
QUANTITY                               *Random();

static struct {
	char	*name;
	double	cval;
	}  consts[] = {
		"PI", 		3.14159265358979323846,
		"E",		2.71828182845904523536,
		"GAMMA",	0.57721566490153286060,
		"DEG",	       57.29577951308232087680,
	};

#define NO_CONSTANTS (sizeof(consts)/sizeof(consts[0]))

static struct {
	char	*name;
	int	kval;
	} keywords[] = {
		    "if",       IF,
		  "else",       ELSE,
		  "then",       THEN,
		 "while",       WHILE,
		 "print",       PRINT,
                   "for",       FOR,
		  "quit",       QUIT,
		 "break",       BREAK,
	    "SetUnitsOn",       SET_UNITS_ON,
           "SetUnitsOff",       SET_UNITS_OFF,
	};

#define NO_KEYWORDS (sizeof(keywords)/sizeof(keywords[0]))

/* 
 *  ====================================================
 *  Builtin Quantity Functions : includes math functions
 *  ====================================================
 */

static struct {
	char	       *name;
	QUANTITY  *(*func)();
	} builtin_quantity[] = {
		"random",	Random,
	};

#define NO_BUILTIN_QUANTITY (sizeof(builtin_quantity)/sizeof(builtin_quantity[0]))

/* mathematic functions */

static struct {
	char    	     *name;
	QUANTITY 	*(*func)();
	} builtin1_quantity[] = {
		"sin",		Sin,
		"cos",		Cos,
		"tan",		Tan,
		"atan",		Atan,
		"log",		Log,
		"log10",	Log10,
		"exp",		Exp,
		"int",		Integer,
		"sqrt",		Sqrt,
		"abs",		Fabs,
              "QDimenLess",     QuantityUnitsLess,

	};

#define NO_BUILTIN1_QUANTITY (sizeof(builtin1_quantity)/sizeof(builtin1_quantity[0]))

static struct {
	char	       *name;
	QUANTITY  *(*func)();
	} builtin2_quantity[] = {
		"pow",		Pow,
	};

#define NO_BUILTIN2_QUANTITY (sizeof(builtin2_quantity)/sizeof(builtin2_quantity[0]))

static struct {
	char          *name;
	QUANTITY *(*func)();
	} builtin3_quantity[] = {
              "QuanCast",    QuantityCast,
		   "Det",	MatrixDet,
		   "Max",	MatrixMax,
		   "Min",	MatrixMin,
		"L2Norm",    MatrixL2Norm,
	};

#define NO_BUILTIN3_QUANTITY (sizeof(builtin3_quantity)/sizeof(builtin3_quantity[0]))

/* 
 *  ---------------------------------------------------
 *  Builtin Matrix and FE Functions with 0,1 and 2 args
 *  ---------------------------------------------------
 */

static struct {
	char           *name;
	MATRIX    *(*func)();
	} builtin_matrix[] = {
	              "Stiff",  Form_Stiffness,
	       "ExternalLoad", 	Form_External_Load,
	     "EquivNodalLoad", 	Form_Equiv_Nodal_Load,
	};

#define NO_BUILTIN_MATRIX (sizeof(builtin_matrix)/sizeof(builtin_matrix[0]))

static struct {
	char           *name;
	MATRIX    *(*func)();
	} builtin1_matrix[] = {
		       "Copy", 	MatrixCopy,
		  "Dimension", 	MatrixDimension,
	          "Decompose",  MatrixLU,
		      "Trans", 	MatrixTranspose,
		     "Matrix",  MatrixAllocate,
		       "Diag",  MatrixDiag,
		       "Zero",  MatrixZero,
		        "One",  MatrixOne,
                 "MDimenLess",  MatrixUnitsLess,
                    "Inverse",  MatrixInverse,
	       	       "Mass", 	Form_Mass,
                 "Eigenvalue",  Extract_Eigenvalue,
                "Eigenvector",  Extract_Eigenvector,
                   "GetCoord",  Get_Coord,
                    "GetNode",  Get_Node,
                     "GetDof",  Get_Dof,
               "GetStiffness",  Get_Stiffness,
                 "GetSection",  Get_Section,
                "GetMaterial",  Get_Material,
	};

#define NO_BUILTIN1_MATRIX (sizeof(builtin1_matrix)/sizeof(builtin1_matrix[0]))

static struct {
	char           *name;
	MATRIX    *(*func)();
	} builtin2_matrix[] = {
		      "Solve",    MatrixSolve,
	       "Substitution",    MatrixFB,
                   "GetDispl",    Get_Displ,
                  "GetStress",    Get_Stress,
	};

#define NO_BUILTIN2_MATRIX (sizeof(builtin2_matrix)/sizeof(builtin2_matrix[0]))

static struct {
        char           *name;
        MATRIX   * (*func)();
        } builtin3_matrix[] = {
                "Eigen",  Solve_Eigen,
        };

#define NO_BUILTIN3_MATRIX (sizeof(builtin3_matrix)/sizeof(builtin3_matrix[0]))

static struct {
        char            *name;
        MATRIX     *(*func)();
        } builtinvar_matrix[] = {
	       "InternalLoad",  Form_Internal_Load,
                "PrintStress",  Print_Stress,
                "PrintMatrix",  MatrixPrintVar,
             "PrintMatrixCast", MatrixPrintCast,
	        "ColumnUnits",  MatrixColumnUnits,
	           "RowUnits",  MatrixRowUnits,
                    "Extract",  MatrixExtract,
                        "Put",  MatrixPut,
        };

#define NO_BUILTINVar_MATRIX (sizeof(builtinvar_matrix)/sizeof(builtinvar_matrix[0]))

/* ------------------------------------
 * Finite Element Functions and Objects
 * ------------------------------------ */

static struct {
	char        *name;
	void	(*func)();
	} fe_mesh[] = {
	      "StartMesh",     Start_Mesh,
	        "EndMesh",       End_Mesh,
	      "PrintMesh",     Print_Mesh,
	};

#define NO_FE_MESH (sizeof(fe_mesh)/sizeof(fe_mesh[0]))

static struct {
	char        *name;
	void	(*func)();
	} fe_node[] = {
		"AddNode",	 Add_Node, 
		"FixNode",	 Fix_Node,
               "NodeLoad",      Node_Load,
	};

#define NO_FE_NODE (sizeof(fe_node)/sizeof(fe_node[0]))

static struct {
	char        *name;
	void	(*func)();
	} builtin_fe_function[] = {
         "UpdateResponse",     UpdateResponse,
	};

#define NO_BLTIN_FE_FUNC (sizeof(builtin_fe_function)/sizeof(builtin_fe_function[0]))

static struct {
	char        *name;
	void	(*func)();
	} builtin1_fe_function[] = {
	     "PrintDispl",     Print_Displ,
	   "ElmtStateDet",     Elmt_State_Det,
	};

#define NO_BLTIN1_FE_FUNC (sizeof(builtin1_fe_function)/sizeof(builtin1_fe_function[0]))

static struct {
	char        *name;
	int	     kval;
	} fe_objects[] = {
            "ElementAttr",     ELMT_ATTR,
            "SectionAttr",     SECT_ATTR,
           "MaterialAttr",     MATL_ATTR,
              "FiberAttr",     FIB_ATTR,
           "SetUnitsType",     UNIT,
                   "type",     TYPE,
                "section",     SECTION,
               "material",     MATERIAL,
                  "fiber",     FIBER,
                    "map",     MAP,
                     "to",     TO,
                   "ldof",     LDOF,
                   "gdof",     GDOF,
	       "LinkNode",     LINK_NODE,
		"AddElmt",     ADD_ELMT,
	};

#define NO_OBJECTS (sizeof(fe_objects)/sizeof(fe_objects[0]))

/*
 *  ------------------------------------------------
 *  Adjustable Parameters for Finite Element Problem
 *  ------------------------------------------------
 */

static struct {
	char        *name;
	double 	    value;
	} fe_parameters[] = {
                "NDimension",  (double) UNIT_NDM,
       	       "NDofPerNode",  (double) UNIT_NDF,
        "MaxNodesPerElement",  (double) UNIT_NEN,
           "InPlaneIntegPts",  (double) UNIT_IN_PLANE_INTEG_PTS,
         "ThicknessIntegPts",  (double) UNIT_INTEG_PTS,
             "GaussIntegPts",  (double) UNIT_INTEG_PTS,
        };


#define NO_PARAMETERS (sizeof(fe_parameters)/sizeof(fe_parameters[0]))

/* 
 *  ============================================
 *  Engineering Units : see units.h for details.
 *  ============================================
 */ 

static DIMENSIONS eng_units[] = {

/* LENGTH */         "micron",   1E-6,  1,  0,  0, 0, 0, SI,
                         "mm",   1E-3,  1,  0,  0, 0, 0, SI,
                         "cm",   1E-2,  1,  0,  0, 0, 0, SI,
                         "dm",   1E-1,  1,  0,  0, 0, 0, SI,
                          "m",    1.0,  1,  0,  0, 0, 0, SI,
                         "km",   1E+3,  1,  0,  0, 0, 0, SI,

/* MASS   */              "g",   1E-3,  0,  1,  0, 0, 0, SI,
                         "kg",    1.0,  0,  1,  0, 0, 0, SI,
                         "Mg",   1E+6,  0,  1,  0, 0, 0, SI,

/* TIME   */            "sec",    1.0,  0,  0,  1, 0, 0, SI_US,
                         "ms", 0.0010,  0,  0,  1, 0, 0, SI_US, /* micro sec */
                        "min",   60.0,  0,  0,  1, 0, 0, SI_US,
                         "hr", 3600.0,  0,  0,  1, 0, 0, SI_US, 

/* TEMPERATURE */     "deg_C",    1.0,  0,  0,  0, 1, 0, SI,
                      "DEG_C",    1.0,  0,  0,  0, 1, 0, SI,   /* incr. temp. only */

/* FREQUENCY   */        "Hz",      1.0,  0,  0, -1, 0, 0, SI_US,
/* SPEED       */       "rpm", 1.0/60.0,  0,  0, -1, 0, 0, SI_US,  /* rev per min */
                        "cps",      1.0,  0,  0, -1, 0, 0, SI_US,  /* cycle per sec */

/* FORCE       */         "N",      1.0,  1,  1, -2, 0, 0, SI,
                         "kN",     1E+3,  1,  1, -2, 0, 0, SI,
                        "kgf",  9.80665,  1,  1, -2, 0, 0, SI,

/* PRESSURE    */        "Pa",      1.0, -1,  1, -2, 0, 0, SI,
                        "kPa",     1E+3, -1,  1, -2, 0, 0, SI,
                        "MPa",     1E+6, -1,  1, -2, 0, 0, SI,
                        "GPa",     1E+9, -1,  1, -2, 0, 0, SI,

/* ENERGY      */       "Jou",      1.0,  2,  1, -2, 0, 0, SI,
                         "kJ",     1E+3,  2,  1, -2, 0, 0, SI,

/* POWER       */      "Watt",      1.0,  2,  1, -3, 0, 0, SI,
                         "kW",     1E+3,  2,  1, -3, 0, 0, SI,

/* LENGTH */             "in",      1.0,  1,  0,  0, 0, 0, US,
                        "mil",     1E-3,  1,  0,  0, 0, 0, US,
                   "micro_in",     1E-6,  1,  0,  0, 0, 0, US,
                         "ft",     12.0,  1,  0,  0, 0, 0, US,
                       "yard",     36.0,  1,  0,  0, 0, 0, US,
                       "mile",  63360.0,  1,  0,  0, 0, 0, US,

/* VOLUME */         "gallon",      1.0,  3,  0,  0, 0, 0, US,
                     "barrel",      1.0,  3,  0,  0, 0, 0, US,

/* MASS   */         "lb",          1.0,  0,  1,  0, 0, 0, US,
                  "grain",     1.0/7E+3,  0,  1,  0, 0, 0, US,
                    "ton",         2E+3,  0,  1,  0, 0, 0, US,
                    "klb",         1E+3,  0,  1,  0, 0, 0, US,

/* TEMPERATURE */ "deg_F",          1.0,  0,  0,  0, 1, 0, US,
                  "DEG_F",          1.0,  0,  0,  0, 1, 0, US,/* incre. temp only */

/* FORCE  */        "lbf",          1.0,  1,  1, -2, 0, 0, US,
                   "kips",         1E+3,  1,  1, -2, 0, 0, US,

/* PRESSURE */      "psi",          1.0, -1,  1, -2, 0, 0, US,
                    "ksi",         1E+3, -1,  1, -2, 0, 0, US,

/* PLANE ANGLE */   "deg",   0.01745329,  0,  0,  0, 0, 1, SI_US,
                    "rad",          1.0,  0,  0,  0, 0, 1, SI_US,
};

#define NO_UNITS (sizeof(eng_units)/sizeof(eng_units[0]))

/*
 *  =====================================================
 *  Init_Problem() : Initialize database for symbol table
 *  =====================================================
 */

Init_Problem() {
int            i;
SYMBOL        *s;
int UNITS_SWITCH;

    UNITS_SWITCH = CheckUnits();

    /* Macro Language Keywords and Constants */

    for (i = 0; i < NO_KEYWORDS; i++)
 	build_table(keywords[i].name,  keywords[i].kval, 0.0);

    for (i = 0; i < NO_CONSTANTS; i++) {
	s = build_table(consts[i].name, QUAN, 0.0);

	s->u.q          = (QUANTITY *) MyCalloc(1,sizeof(QUANTITY));
	s->u.q->value   = consts[i].cval;
        if(UNITS_SWITCH==ON) {
           s->u.q->dimen   = (DIMENSIONS *) MyCalloc(1,sizeof(DIMENSIONS));
           ZeroUnits(s->u.q->dimen);
           s->u.q->dimen->units_type = SI_US;
        }
    }

    /* Builtin Quantity Functions */

    for (i = 0; i < NO_BUILTIN_QUANTITY; i++)  {
	s = build_table(builtin_quantity[i].name, BLTIN_QUANTITY, 0.0);
	s->u.quantityptr = builtin_quantity[i].func;
    }

    for (i = 0; i < NO_BUILTIN1_QUANTITY; i++)  {
	s = build_table(builtin1_quantity[i].name, BLTIN1_QUANTITY, 0.0);
	s->u.quantityptr = builtin1_quantity[i].func;
    }

    for (i = 0; i < NO_BUILTIN2_QUANTITY; i++)  {
        s = build_table(builtin2_quantity[i].name, BLTIN2_QUANTITY, 0.0);
        s->u.quantityptr = builtin2_quantity[i].func;
    }

    for (i = 0; i < NO_BUILTIN3_QUANTITY; i++)  {
	s = build_table(builtin3_quantity[i].name, BLTIN3_QUANTITY, 0.0);
	s->u.quantityptr = builtin3_quantity[i].func;
    }

    /* Builtin Matrix Functions */

    for (i = 0; i < NO_BUILTIN_MATRIX; i++)  {
	s = build_table(builtin_matrix[i].name, BLTIN_MATRIX, 0.0);
	s->u.matrixptr = builtin_matrix[i].func;
    }

    for (i = 0; i < NO_BUILTIN1_MATRIX; i++)  {
	s = build_table(builtin1_matrix[i].name, BLTIN1_MATRIX, 0.0);
	s->u.matrixptr = builtin1_matrix[i].func;
    }

    for (i = 0; i < NO_BUILTIN2_MATRIX; i++)  {
       s = build_table(builtin2_matrix[i].name, BLTIN2_MATRIX, 0.0);
       s->u.matrixptr = builtin2_matrix[i].func;
    }

    for (i = 0; i < NO_BUILTIN3_MATRIX; i++)  {
       s = build_table(builtin3_matrix[i].name, BLTIN3_MATRIX, 0.0);
       s->u.matrixptr = builtin3_matrix[i].func;
    }

    for (i = 0; i < NO_BUILTINVar_MATRIX; i++)  {
       s = build_table(builtinvar_matrix[i].name, BLTINVar_MATRIX, 0.0);
       s->u.matrixptr = builtinvar_matrix[i].func;
    }

    /* Engineering Units and Quantities */

    for (i = 0; i < NO_UNITS; i++) {
         if((eng_units[i].units_type != SI_US)        &&
            (eng_units[i].units_type != SI)           &&
            strcmp(eng_units[i].units_name, "deg_C")  &&
            strcmp(eng_units[i].units_name, "deg_F"))
            eng_units[i] = UnitsScaleConvert(eng_units[i], SI);

            s = build_table(eng_units[i].units_name, DIMENSION, 0.0);

	    s->u.q =        (QUANTITY *) MyMalloc(sizeof(QUANTITY));

	    s->u.q->dimen = (DIMENSIONS *)    MyMalloc(sizeof(DIMENSIONS));
	    s->u.q->dimen->units_name   = eng_units[i].units_name;
	    s->u.q->dimen->units_type   = eng_units[i].units_type;
	    s->u.q->dimen->scale_factor = eng_units[i].scale_factor;
	    s->u.q->value               = eng_units[i].scale_factor;
	    s->u.q->dimen->length_expnt = eng_units[i].length_expnt;
	    s->u.q->dimen->mass_expnt   = eng_units[i].mass_expnt;
	    s->u.q->dimen->time_expnt   = eng_units[i].time_expnt;
	    s->u.q->dimen->temp_expnt   = eng_units[i].temp_expnt;
	    s->u.q->dimen->radian_expnt = eng_units[i].radian_expnt;
    }

    /* Finite Element Functions, Objects and Solution Procedures */

    for (i = 0; i < NO_FE_MESH; i++) {
	s = build_table(fe_mesh[i].name, MESH, 0.0);
	s->u.voidptr = fe_mesh[i].func;
    }

    for (i = 0; i < NO_FE_NODE; i++) {
	s = build_table(fe_node[i].name, NODE_QUANT, 0.0);
	s->u.voidptr = fe_node[i].func;
    }

    for (i = 0; i < NO_BLTIN_FE_FUNC; i++) {
	s = build_table(builtin_fe_function[i].name, BLTIN_FE_FUNC, 0.0);
	s->u.voidptr = builtin_fe_function[i].func;
    }

    for (i = 0; i < NO_BLTIN1_FE_FUNC; i++) {
	s = build_table(builtin1_fe_function[i].name, BLTIN1_FE_FUNC, 0.0);
	s->u.voidptr = builtin1_fe_function[i].func;
    }

    for (i = 0; i < NO_OBJECTS; i++) {
        s = build_table(fe_objects[i].name, fe_objects[i].kval, 0.0);
    }

    /* Load Adjustable Parameters for Size of Finite Element Problem */

    for (i = 0; i < NO_PARAMETERS; i++) {
        s = build_table(fe_parameters[i].name, QUAN, 0.0);

        s->u.q =        (QUANTITY *) MyMalloc(sizeof(QUANTITY));
        s->u.q->value               = fe_parameters[i].value;
        if(UNITS_SWITCH==ON) {
           s->u.q->dimen   = (DIMENSIONS *)MyMalloc(sizeof(DIMENSIONS));
           ZeroUnits(s->u.q->dimen);
           s->u.q->dimen->units_type = SI_US;
        }
    }	
}
