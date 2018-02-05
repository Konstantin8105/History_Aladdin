/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  elmt.h : Definitions for finite element library.
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

#ifndef ELMT_H
#define ELMT_H

/* element library functions */

ARRAY     *elmlib();
ARRAY     *elmt_frame_2d();
ARRAY     *elmt_frame_3d();
ARRAY     *elmt_psps();
ARRAY     *elmt_plate();
ARRAY     *elmt_shell_4n_q();   /* Shell Element with Drill DOF */
ARRAY     *elmt_shell_4n();     /* X.G. shell element, 4 node   */
ARRAY     *elmt_shell_8n();     /* X.G. shell element, 8 node   */
ARRAY     *elmt_fiber_2d();     /* 2-D fiber element with one global shear spring */
ARRAY     *elmt_fiber_3d();     /* 3-D fiber element with two global shear springs */

/* element properties handling */

EFRAME    *assign_properties();
void       print_property_frame_2d();
void       print_property_frame_3d();
void       print_property_psps();
void       print_property_plate();
void       print_property_shell_4n();
void       print_property_shell_4n_q();
void       print_property_shell_8n();
void       print_property_fiber_2d();
void       print_property_fiber_3d();

/* Generic Template for Item in Finite Element Library */

static struct {
	char                 *name;         /* name of elment type                      */
        ARRAY  *(*elmt_lib_func)();         /* pointer to elmt library function         */
        void       (*elmt_print)();         /* pointer to elmt library print function   */
        int                 no_dof;         /* No dof per node                          */
        int       no_node_per_elmt;         /* No nodes per element                     */
        int               no_dimen;         /* No dimension of problem                  */
	} elmt_library[] = {
           "FRAME_2D",        elmt_frame_2d,    print_property_frame_2d,   3, 2, 2,
           "FRAME_3D",        elmt_frame_3d,    print_property_frame_3d,   6, 2, 3,
           "SHELL_4N",        elmt_shell_4n,    print_property_shell_4n,   5, 4, 3,
           "SHELL_4NQ",       elmt_shell_4n_q,  print_property_shell_4n_q, 6, 4, 3,
           "SHELL_8N",        elmt_shell_8n,    print_property_shell_8n,   5, 8, 3,
           "PLANE_STRAIN",    elmt_psps,        print_property_psps,       2, 4, 2,
           "PLANE_STRESS",    elmt_psps,        print_property_psps,       2, 4, 2,
           "DKT_PLATE",       elmt_plate,       print_property_plate,      3, 4, 3,
           "FIBER_2D",        elmt_fiber_2d,    print_property_fiber_2d,   3, 2, 2,
           "FIBER_3D",        elmt_fiber_3d,    print_property_fiber_3d,   6, 2, 3,
           "FIBER_2DS",       elmt_fiber_2d,    print_property_fiber_2d,   3, 2, 2,
           "FIBER_3DS",       elmt_fiber_3d,    print_property_fiber_3d,   6, 2, 3,
	};

#define NO_ELEMENTS_IN_LIBRARY (sizeof(elmt_library)/sizeof(elmt_library[0]))

/* ------------------------- */
/* Cases for Element Library */
/* ------------------------- */

#define PROPTY             1 
#define CHERROR            2   
#define STIFF              3 
#define STRESS             4
#define STRESS_MATRIX      5
#define MASS_MATRIX        6 
#define LOAD_MATRIX        7
#define PRESSLD            8
#define STRESS_LOAD        9
#define EQUIV_NODAL_LOAD  10
#define STRESS_UPDATE     11

#endif /* end case ELMT_H */
