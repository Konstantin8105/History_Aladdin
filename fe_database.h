/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_database.h : Database for Finite Element Machine
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

#ifndef FE_DATABASE_H
#define FE_DATABASE_H

typedef struct element_response {
     MATRIX                 *Forces;
     MATRIX                 *stress;
     MATRIX                  *displ; 
     MATRIX              *strain_pl;
     MATRIX         *strain_pl_incr;
     double       *effect_pl_strain;
     double     *eff_pl_strain_incr;
     QUANTITY            min_moment;
     QUANTITY       max_moment, Mzc;
     QUANTITY  min_shear, max_shear;
     MATRIX     *Q_saved,  *q_saved;
     MATRIX  	*sr_saved, *er_saved;    /* matrix[no_integ_pt][ifib] */
     MATRIX     *s0_saved, *e0_saved;
     MATRIX     *sx_saved, *ex_saved;
} RESPONSE;

typedef struct mater_load_curve {
   char                  *name;   /* load curve type                       */
   double                   *R;   /* radius of yield surface               */
   double        **back_stress;   /* back stress vector                    */
   double                   *H;   /* tangent of stress-plastic strain curve*/
   double                alpha;   /* parameter for Ramberg-Osgood relation */
   double                    n;   /* strain hardening exponent             */
   double                 beta;   /* parameter for strain hardening        */
                                  /* beta = 0 kinematic hardening          */
                                  /* beta = 1 isotropic hardening          */
   /*  variables defined for elmt_shell_4n_q  */
   int                   ialph;   /* Rotation constant ; 0=none, 1=hughes  */
   double                  pen;   /* Penalty constant                      */
   double              load[6];   /* Loading parameters, per unit area     */
                                  /* [0]: uniform loading in 1-direction   */
                                  /* [1]: uniform loading in 2-direction   */
                                  /* [2]: pressure loading normal to shell */
                                  /* [3]: gravity loading in x-direction   */
                                  /* [4]: gravity loading in y-direction   */
                                  /* [5]: gravity loading in z-direction   */
} MATER_LOAD_CURVE;

/*
 * ------------------
 * Nodes and Elements  
 * ------------------
 */

typedef struct node {
        QUANTITY       *coord;     /* Coordinates of node            */
        int            *bound_id;  /* Boundary Constratint ID        */
        int             rb_num;    /* component of rigid body number */
        QUANTITY       *disp;
        MATRIX         *TrT;
} NODE;

typedef struct element_state {
        int             state; /* state = 1    elastic-plastic deformation */
                               /* state = 0    elastic deformation         */
	int	**yielding_saved, **pre_range_saved, **pre_load_saved;

} ELEMENT_STATE;

typedef struct element {
        char     *elmt_attr_name;  /* Elmt_Attr name                      */
        int        *node_connect;  /* List of Nodal Connectivities        */
	int         elmt_attr_no;  /* Elmt Attribute No (stored in frame) */
        int             *d_array;  /* Destination array storage           */
        RESPONSE             *rp;  /* Response Pointer                    */
        ELEMENT_STATE       *esp;  /* Deformation State of the Element    */
        MATER_LOAD_CURVE *LC_ptr;  /* properties for describe yield       */
                                   /* surface and stress-strain curve     */

} ELEMENT;

/*
 *  ------------------------------- 
 *  Data Structure For Rigid Bodies 
 *  ------------------------------- 
 */

typedef struct rigid {
        char       *rbody_attr_name;  /* Specified Rbody_Attr name              */
        int              nodes_conn;  /* number of nodes connected              */
        int                 rb_type;  /* assigned aelmt type                    */
        int                     *in;  /* node numbers connected by rigid body   */
	int              *per_nodes;  /* list of perimeter nodes                */
        int               *rest_dof;  /* array of dofs: 1/o-dof un-/ restricted */
        QUANTITY      xcg, ycg, zcg;  /* coordinates of c.g. of rigid body      */
        QUANTITY              *prop;  /* properties input for rigid body        */
} RIGID;


/*
 *  -----------------------------------------------------------------------
 *  Working Arrays for Element, Material and Section Attributes 
 * 
 *  Note : MATERIAL and SECTION Attributes only stored in Symbol Table.
 *       : ELEMENT_ATTR is stored in symbol table and frame data structure.
 *       : Space for "double *d" only used when attached to frame. 
 *  -----------------------------------------------------------------------
 */

typedef struct fiber_attr {
	QUANTITY	y, z;
	QUANTITY	area;
	QUANTITY	Es, Et, fy;
	QUANTITY	Gs, Gt, fv;
} FIBER_ATTR;

typedef struct fiber_elmt {
	int		no_fiber;
	int		no_shear; /* No of uniaxial shear spring in each direction */
	FIBER_ATTR      *fiber;   /* Array of fibers for all fiber properties */
} FIBER_ELMT;

typedef struct elmt_attr {
        char                    *name;
        char               *elmt_type;
        char                *material;
        char                 *section;
        char         *fiber_attr_name;
        int         *map_ldof_to_gdof;
        QUANTITY       *work_material;   /* Working Array for Material Properties    */
        QUANTITY        *work_section;   /* Working Array for Section  Properties    */
	FIBER_ELMT        *work_fiber;
} ELEMENT_ATTR;
                                         /* notes :                                  */
                                         /* For given elmement No, element property  */
                                         /* can be recalled from these two working   */
                                         /* arraries instead of looking up hashing   */
                                         /* table                                    */

typedef struct section_attr {
        char *section_name;
        int   section_type; 
        QUANTITY      Ixx;           /* [0] */
        QUANTITY      Iyy;           /* [1] */
        QUANTITY      Izz;           /* [2] */
        QUANTITY      Ixz;           /* [3] */
        QUANTITY      Ixy;           /* [4] */
        QUANTITY      Iyz;           /* [5] */
        QUANTITY      weight;        /* [6]  Section weight */
        QUANTITY      bf;            /* [7]  Width of flange         */
        QUANTITY      tf;            /* [8]  thickness of flange     */
        QUANTITY      depth;         /* [9]  Section depth           */
        QUANTITY      area;          /* [10] Section area            */
        QUANTITY      plate_thickness; /* [11] */
        QUANTITY      tor_const;     /* [12] Torsional Constant J    */
        QUANTITY      rT;            /* [13] Section radius of gyration */
        QUANTITY      width;         /* [14] Section width           */
        QUANTITY      tw;            /* [15] Thickness of web        */
	double        ks;            /* [16] Shear correction factor */
      } SECTION_ATTR;

 typedef struct materials {
        char                   *name;    /* Material name                    */ 
        QUANTITY                   E;    /* [0] Young's modulus                  */
        QUANTITY                   G;    /* [1] Shear modulus                    */
        QUANTITY                  fy;    /* [2] Yield stress                     */
        QUANTITY                  ET;    /* [3] Tangent Young's Modulus          */
        double                    nu;    /* [4] Poission's ratio                 */
        QUANTITY             density;    /* [5] */
        QUANTITY                  fu;    /* [6] Ultimate stress                  */
        QUANTITY      *alpha_thermal;    /* [7,8,9] thermal expansion coeff(0,1,2) */
	QUANTITY                  Gt;    /* [10] Tangent shear modulus           */
	QUANTITY                  fv;    /* [11] Yielding stress for shear       */
	MATER_LOAD_CURVE     *LC_ptr;    /* parameters for decribing yield   */
                                         /* surface and stress strain curve  */
      } MATERIAL_ATTR;

/* 
 *  ------------------------------------------
 *  Data Structure for Node and Element Arrays
 *  ------------------------------------------
 */ 

typedef struct node_prop {
	int           node_f;   /* node number     */
	QUANTITY         *fn; 
} NODE_LOADS;

typedef struct eload_lib {
        char                 *name;
        short                 type;  /* set to  -1 for dist loads at input stage*/
	QUANTITY       *body_force; 
	QUANTITY       temp_change; 
	QUANTITY      *init_stress; 
	double        *init_strain; 
	QUANTITY         *traction; 
        int                  *nopl;
        MATRIX                 *pr;
        int                face_no;
        int             numnp_face;
        QUANTITY                 P;        /* Load Value in Local y Direction */
        QUANTITY                 a;        /* Start distance of Load (P) from left end */
        QUANTITY                 b;        /* End   distance of Load (P) from left end */
        QUANTITY          px,py,pz;        /* Point Load Value in Local x,y,z Direction */
        QUANTITY          mx,my,mz;        /* Moment Load Value in Local x,y,z Direction */
        QUANTITY          bx,by,bz;        /* Dist Load Value in Local x,y,z Direction */

} ELOAD_LIB;

typedef struct element_loads {
	int            elmt_no;     /* Element number identification          */
	int          elmt_type;     /* Element type   identification          */
        int     no_loads_faces;     /* Number of elmts surface under traction */
        double     *face_direc;     /* Direction cosine of tracted elment     */
	ELOAD_LIB     *elib_ptr;   /* Pointer to element_load struct         */
} ELEMENT_LOADS;


/* ------------------------------------------------ */
/* Data structure for fiber element loading history */
/* ------------------------------------------------ */
typedef struct section_j_1
{
   double   xi;
   double   wi;
   MATRIX  *fx;
   MATRIX  *Dx;
   MATRIX  *rx;
}SECTION_DATA;

typedef struct load_history
{
   int     elmt_no;
   MATRIX  *stress, *strain, *tangent;
   MATRIX  *sr, *er, *s0, *e0;
   int     **yielding, **pre_load, **pre_range, **loading;
}HISTORY_DATA;

typedef struct fiber_respond
{
   int   total_fiber_elmt;
   HISTORY_DATA  *history;
}FIBER_RESPONSE;

/* ------------------------ */
/* Data structure for frame */
/* ------------------------ */

typedef struct frame {
	char                  *name; /* Title */
        int                no_nodes; /* no of nodes                              */
        int                no_rigid; /* no of rigid body elements                */
        int             no_elements; /* no of flexible elements                  */
        int         no_element_attr; /* no of element attributes                 */
        int           no_node_loads; /* no of nodal loads                        */
        int        no_element_loads; /* no of element loads                      */
        int                no_dimen; /* no of dimension                          */ 
        int                  no_dof; /* no of dof per node--in global            */
        int       no_nodes_per_elmt; /* max no of nodes per element              */
	NODE                  *node;
	RIGID                *rigid;
	ELEMENT            *element;
        ELEMENT_ATTR         *eattr;
	NODE_LOADS         *nforces;
	ELEMENT_LOADS      *eforces;
        int                  *jdiag;
        int                   no_eq;
        int                     nre; /* flag for presence of nodal load delet later */
	int             no_integ_pt; /* Gauss_Lobatto integration point */
} EFRAME;


/* 
 * --------------------------
 * Working Array For Elements
 * --------------------------
 */ 
typedef struct Integ_Pts {
    int      surface_pts;  /* No of integration points on surface/x-y plane     */
    int    thickness_pts;  /* No of integration points in thickness/z-direction */ 
    int        integ_pts;  /* User defined inte pts                             */
} INTEG_PTS;

typedef struct array {
        int                   no_dimen;  /* No dimensions in Problem                     */
        int                    elmt_no;  /* Element no                                   */
	int               elmt_attr_no;  /* Elmt Attribute No, as stored in frame        */
        char                *elmt_type;  /* Element Type                                 */
        int             nodes_per_elmt;  /* Max no nodes per element                     */
        int               dof_per_node;  /* Dof per node                                 */
        int              size_of_stiff;  /* Size of stiffness (mass) matrix              */
        int                       type;  /* Flag for type of mass : LUMPED, Consistent   */
        int                   *d_array;  /* Destination Array                            */
        char            *material_name;  /* Material used in the element                 */
        QUANTITY        *work_material;  /* Section/Material working array               */
        QUANTITY         *work_section;  /* Section/Material working array               */
        QUANTITY               **coord;  /* Nodal coordinates                            */
        int              *node_connect;  /* Array of Nodal Connectivities                */
        QUANTITY          *nodal_loads;  /* Array of Point Loads                         */
        QUANTITY           *nodal_temp;  /* Nodal temperature change                     */
        QUANTITY    **nodal_body_force;  /* Nodal body force                             */
        QUANTITY   **nodal_init_stress;  /* Nodal initial stress                         */
        double     **nodal_init_strain;  /* Nodal initial strain                         */
        QUANTITY      **nodal_traction;  /* Nodal surface tracion on load boundary       */
        MATRIX                  *stiff;  /* Stiffness (and Mass) Matrix ???              */
        MATRIX       *equiv_nodal_load;  /* Equivalant nodal loads                       */
        MATRIX           *mater_matrix;  /* Material matrix at given gaussian point      */
        MATRIX                  *displ;  /* Element Level Displacements                  */
        MATRIX             *displ_incr;  /* Element Level Displacements incremental      */
                                         /* for time interval [t,t+dt]                   */
                                         /* they are instored in local coordinate system */
        MATRIX            *strain_rate;  /* Element Level Strain Rate                    */
        MATRIX            *stress_rate;  /* Element Level Stress Rate                    */
        MATRIX                 *stress;  /* Element Level Stress at integration points   */
        MATRIX              *strain_pl;  /* Elemt level plastic strain at integ points   */
        MATRIX         *strain_pl_incr;  /* Elemt level plastic strain incr at integ pts */ 
        double       *effect_pl_strain;  /* effective plastic strain at integration pts  */
        double     *eff_pl_strain_incr;  /* eff. plastic strain incr. at integration pts */
        double              *direc_cos;  /* direction cosine for shell elmt used for     */
                                         /* correction of co-rotational coord            */ 
        int                 elmt_state;  /* Element state : Elastic or Plastic           */
        MATER_LOAD_CURVE       *LC_ptr;  /* material load curve pointer                  */
        double                     eep;  /* Plastic strain ?? (use later)                */
        QUANTITY                eangle;
        QUANTITY                ealpha;
        QUANTITY                length;
        ELEMENT_LOADS   *elmt_load_ptr;
        INTEG_PTS           *integ_ptr;
	FIBER_ELMT          *fiber_ptr;  /* NOTE:  use the same address in frame */
        MATRIX     *Q_saved,  *q_saved;
        MATRIX    *sr_saved, *er_saved, *s0_saved, *e0_saved, *sx_saved, *ex_saved;
	int           **yielding_saved, **pre_range_saved, **pre_load_saved;
} ARRAY;

#endif /* end case FE_DATABASE_H */
