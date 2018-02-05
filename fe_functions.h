/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_functions.h : External Function Declarations for Finite Elements
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

#ifndef FE_FUNCTIONS_H
#define FE_FUNCTIONS_H

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

EFRAME    *FrameAlloc();
EFRAME    *profile();
EFRAME    *Set_Elmt_Attrs();
EFRAME    *plink();
EFRAME    *rplink();

ARRAY     *Alloc_p_Array();
ARRAY     *Assign_p_Array();
ARRAY     *Element_Property();
ARRAY     *Mate_Property();
ARRAY     *Eload_Property();

MATRIX    *Element_Matrix();
MATRIX    *Element_Equiv();
QUANTITY  *Element_Vector();

MATRIX    *Assemble_Global();
MATRIX    *Assemble_Global_Load();

int       *Destination_Array();
int        Bound_Disp();
double   **Boundary_Conditions();
double   **MATER_MAT_PLANE();    /* function for elmt_psps.c only */

QUANTITY  *pload();
QUANTITY  *Modify_Load();
QUANTITY  *Addload_Vector();
QUANTITY  *Transform_Force();

QUANTITY  *Assemble_Nodal_Load();
QUANTITY  *Assemble_Gravity_Load();
QUANTITY  *Assemble_Ctrfgl_Load();
double    *Assemble_Equiv_Load(); 

int       *Destination_Array_for_Rigid_Body();
ARRAY     *Assign_p_Array_for_Rigid_Body();
int        Rigidbody_conn();
MATRIX    *Transform_Stiff_Matrix();
MATRIX    *Transform_Rigid_Body_Mass_Matrix();
double   **Transformation_Matrix();
double   **Modify_T_Matrix();

/* functions called in code.c for finite element solution procedures */

void       Start_Mesh();
void       End_Mesh();
void       Print_Mesh();
void       Add_Node();
void       Fix_Node();
void       Node_Load();
void       Link_Node();
void       Add_Elmt();

void       Print_Displ();
#ifdef __STDC__
MATRIX    *Print_Stress( MATRIX *, ... );
#else
MATRIX    *Print_Stress();
#endif

MATRIX    *Form_Stiffness();
MATRIX    *Form_Mass();
MATRIX    *Form_External_Load();
MATRIX    *Form_Equiv_Nodal_Load();
#ifdef __STDC__
MATRIX    *Form_Internal_Load( MATRIX *, ... );
#else
MATRIX    *Form_Internal_Load();
#endif

MATRIX    *Solve_Eigen();
MATRIX    *Velocity_Extract();
MATRIX    *Displacement_Extract();
void       Ldof_to_gdof();

/* Finite Element Allocation Routines */

ELEMENT_ATTR     *Alloc_Element_Attr_Item();
SECTION_ATTR     *Alloc_Section_Attr_Item();
MATERIAL_ATTR    *Alloc_Material_Attr_Item();
FIBER_ELMT       *Alloc_Fiber_Elmt_Attr_Item();

/* functions declarations for Add_Elmt, */

EFRAME    *CheckElementSpace();
EFRAME    *CheckRigidSpace();
EFRAME    *CheckJdiagSpace();   
EFRAME    *CheckNodeSpace();   
EFRAME    *CheckNforcesSpace(); 
EFRAME    *CheckEforcesSpace(); 

/* functions for rule checking/post-processing */

MATRIX    *Get_Coord();
MATRIX    *Get_Node();
MATRIX    *Get_Displ();
MATRIX    *Get_Stress();
MATRIX    *Get_Stiffness();
MATRIX    *Get_Dof();
MATRIX    *Get_Section();
MATRIX    *Get_Material();

/* functions used in elmt_*.c */

MATRIX    *beamst();
MATRIX    *beamms();
MATRIX    *beamms3d();
MATRIX    *beamst3d();
double   **tmat();
double   **rotate();
double   **rotate3d();
int        pstres();
int        gauss();
int        pgauss();
int        shape();
int        shp0();
double   **qushp8();
double   **dktqbm();
int        jacqud();
void       dktb06();
void       dktq06();
void       hshp06();
void       jacq06();
void       jtri06();
void       proj06();
void       rots06();
void       rshp06();
void       stre06();
void       tran06();
double    *pstres06();
void       shp_prt();

/* functions about 4 node shell elmt */

ARRAY     *elmt_shell_4nodes_implicit();
void       Lamina_Sys();
void       Lamina_Sys_Implicit();
void       elmt_shell_shape_4node();
double   **Hourglass_Stress_Rate();
void       Shell_4Node_Mass();
double   **B_MATRIX_4Node();
void       Shell_Stiff_Plane_4node();
double   **Shell_Nodal_Load_Plane();

/* functions about 8 node shell elmt */

ARRAY     *elmt_shell_8nodes_implicit();
void       Lamina_Sys_8node();
void       elmt_shell_shape_8node();
void       Shell_8Node_Mass();
void       Stress_Update_8Node();
double   **B_MATRIX_8Node();
void       Shell_Stiff_Plane_8node();
double   **Shell_Nodal_Load_8Node();

/* functions declarations for 4 Node and 8 Node shell elements */

void       MATER_SHELL_UPDATE();
double   **STRAIN_RATE_SHELL();
double   **Hourglass_Stress_Rate();
double   **Hourglass_stiff();
double   **Strain_Displ_Matrix();
double   **MATER_MAT_SHELL();
void       Load_Curve();
void       Plastic_Deform();

void       DISPL_UPDATE();
void       Stress_Update();
void       BB_Vector();

/* functions declarations for FIBER_2D and FIBER_3D elements */

void       Force_Interpolation_Matrix_2d();
void       Linear_Geometric_Matrix_2d();
void       Section_Tangent_Stiffness_2d();
void       Section_Resisting_Force_2d();
MATRIX    *Rigid_Body_Rotation_2d();
MATRIX    *Element_Transformation_2d();
void       Fiber_Elmt_State_Det_2d();

void       Gauss_Lobatto();
void       Stress_Strain_Relationship();

void       Force_Interpolation_Matrix_3d();
void       Linear_Geometric_Matrix_3d();
void       Section_Tangent_Stiffness_3d();
void       Section_Resisting_Force_3d();
MATRIX    *Rigid_Body_Rotation_3d();
MATRIX    *Element_Transformation_3d();
void       Fiber_Elmt_State_Det_3d();

/* functions for non-linear analysis */

void       SetUpRespondBuffer();
void       UpdateResponse();
void       Elmt_State_Det();

void       SaveRespondBuffer();
void       save_action();

void       SetUpFiberRespondBuffer();
void       SaveFiberRespondBuffer();
HISTORY_DATA    *FiberElmtHistory();

#endif /* end case FE_FUNCTIONS_H */
