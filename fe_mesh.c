/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_mesh.c : Generate Finite Element Mesh
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
 *  Written by: Mark Austin and Wane-Jang Lin                            May 1997
 *  ============================================================================= 
 */

#include <stdio.h>
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
#include "vector.h"
#include "fe_functions.h"
#include "elmt.h"

/* External Declarations for global frame/working element data structures */

extern ARRAY     *array;
extern EFRAME    *frame;

int    TDOF;
int    TNEQ;
int    MDOF;

/*
#define DEBUG 
*/

/*
 *  -------------------------------------
 *  Add Node to Finite Element Mesh
 *  -------------------------------------
 */
#ifdef __STDC__
void Add_Node(double node_no, MATRIX *m)
#else
void Add_Node(node_no, m)
double node_no;
MATRIX      *m; /* coordinate vector/matrix */
#endif
{
int no = 0, num = 0;
int     i, j, length;
int     UNITS_SWITCH;
NODE             *np;

#ifdef DEBUG
       printf("*** Enter Add_Node() : node_no = %4d :\n ", (int) node_no);
#endif 

    UNITS_SWITCH = CheckUnits();
    frame = CheckNodeSpace(frame, (int) node_no);

    switch( UNITS_SWITCH ) {
      case ON:
         for (i =1; i <= frame->no_dimen; i++) {
            frame->node[(int) (node_no-1)].coord[i-1].value 
            = m->uMatrix.daa[0][i-1];
            UnitsCopy( frame->node[(int)(node_no-1)].coord[i-1].dimen, &(m->spColUnits[i-1]) );
         }
         break;
      case OFF:
         for(i =1; i <= frame->no_dimen; i++) 
             frame->node[(int) (node_no-1)].coord[i-1].value
             = m->uMatrix.daa[0][i-1];
         break;
      default:
         break;
    }
    frame->no_nodes = (int) MAX(node_no, frame->no_nodes); 

#ifdef DEBUG
       printf("*** Leave Add_Node()\n");
#endif 

}


/*
 *  -------------------------------------
 *  Add Element to Finite Element Mesh
 *  -------------------------------------
 */

#ifdef __STDC__
void Add_Elmt(double elmt_no, MATRIX *nodes, char *elmt_attr_name)
#else
void Add_Elmt(elmt_no, nodes, elmt_attr_name)
double        elmt_no;
MATRIX         *nodes;
char  *elmt_attr_name;
#endif
{
int             i, no;
char            *type;

#ifdef DEBUG
       printf("\n*** Enter Add_Elmt(%3d) : nodal connectivity = [ ", (int) elmt_no);
       for(i = 1; i <= (int) nodes->iNoColumns; i++) 
           printf("%2d ", (int) nodes->uMatrix.daa[0][i-1]);
       printf(" ]\n");
       printf("          Add_Elmt(%3d) : elmt attribute  = \"%s\"\n", (int) elmt_no, (char *) elmt_attr_name);
       printf("\n");
#endif 


   frame = CheckElementSpace(frame, (int) elmt_no);

   frame->element[(int) (elmt_no-1)].elmt_attr_name
   = SaveString(elmt_attr_name);

   for(i = 1; i <= nodes->iNoColumns; i++)  {
       frame->element[(int)(elmt_no-1)].node_connect[i-1] 
       = (int) nodes->uMatrix.daa[0][i-1];
   }

   frame->no_elements = (int) MAX(elmt_no, frame->no_elements); 

#ifdef DEBUG
       printf("*** Leave Add_Elmt()\n");
#endif 
}


/*
 *  ----------------------------------------------------
 *  Add Element and Section Attributes to Frame Database
 *  ----------------------------------------------------
 */

#ifdef __STDC__
int  Add_Element_Attr(char *name, ELEMENT_ATTR *ep)
#else
int  Add_Element_Attr(name, ep)
char          *name;
ELEMENT_ATTR    *ep;
#endif
{
SYMBOL          *hp;
int               i;

#ifdef DEBUG
       printf("*** Enter Add_Element_Attr() : name = \"%s\"\n", name);
#endif 

     /* [a] : Install elment attribute name into Symbol Table */

     hp = install(name);
     if(hp == NULL)
        FatalError("In Add_Element_Attr() : hp = NULL !!!\n",(char *)NULL);
     else
        hp->u.eap = ep;

     frame->no_element_attr = frame->no_element_attr + 1;

#ifdef DEBUG
       printf("*** Leave Add_Element_Attr()\n");
#endif 

}

#ifdef __STDC__
int  Add_Fiber_Elmt_Attr(char *name, FIBER_ELMT *fep)
#else
int  Add_Fiber_Elmt_Attr(name, fep)
char          *name;
FIBER_ELMT    *fep;
#endif
{
SYMBOL          *hp;

#ifdef DEBUG
       printf("*** Enter Add_Fiber_Elmt_Attr() : name = \"%s\"\n", name);
#endif 

     /* [a] : Install elment attribute name into Symbol Table */

     hp = install(name);
     if(hp == NULL)
        FatalError("In Add_Fiber_Elmt_Attr() : hp = NULL !!!\n",(char *)NULL);
     else
        hp->u.fep = fep;

#ifdef DEBUG
       printf("*** Leave Add_Fiber_Elmt_Attr()\n");
#endif 
}

#ifdef __STDC__
int Add_Section_Attr(char *name, SECTION_ATTR *sp)
#else
int Add_Section_Attr(name, sp)
char          *name;
SECTION_ATTR    *sp;
#endif
{
SYMBOL          *hp;
int               i;

#ifdef DEBUG
       printf("*** Enter Add_Section_Attr() : name = \"%s\"\n", name);
#endif 

     /* [a] : Install Section into Symbol Table */

     hp = install(name);
     if(hp == NULL)
        FatalError("In Add_Section_Attr() : hp = NULL !!!\n",(char *)NULL);
     if(hp->u.sap == NULL)  hp->u.sap = sp;

#ifdef DEBUG
       printf("*** Leave Add_Section_Attr()\n");
#endif 

}

#ifdef __STDC__
int Add_Material_Attr(char *name, MATERIAL_ATTR *mp)
#else
int Add_Material_Attr(name, mp)
char           *name;
MATERIAL_ATTR    *mp;
#endif
{
SYMBOL              *hp;
int    iInPlaneIntegPts;
int  iThicknessIntegPts;
int       iNO_INTEG_pts;
int                i, j;

#ifdef DEBUG
       printf("*** Enter Add_Material_Attr()\n");
#endif 

      hp = lookup("InPlaneIntegPts");   /* number of integration pts in plane/surface      */
      if(hp == NULL)
         iInPlaneIntegPts = 2*2;        /* 2x2 as default */
      else
         iInPlaneIntegPts = (int) hp->u.q->value;

      hp = lookup("ThicknessIntegPts"); /* number of integration pts in thickness direction*/
      if(hp == NULL)
         iThicknessIntegPts = 2;        /* 2 as default */
      else
         iThicknessIntegPts = (int) hp->u.q->value;

      iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;


     /* [0] calculate non-independent parameters */

     switch(((int) mp->LC_ptr->beta)){
        case 0:  /* Kinematic hardening */
        case 1:  /* Isotropic hardening */
          for(j = 1; j <= iNO_INTEG_pts; j++)
              for(i = 1; i <= 6; i++)
                  mp->LC_ptr->back_stress[i-1][j-1] = 0.0;
        break;
        default:
           printf(" In Add_Material_Attr(): UNDEFINED STRAIN_HARDENING \n");
           FatalError("***** Undefined strain hardening parameter found \n", (char *) NULL);
        break;
     }
     /* For inital yielding */
     for(j = 1; j <= iNO_INTEG_pts; j++)
         mp->LC_ptr->R[j-1] = sqrt(2.0/3.0)*mp->fy.value;
     
     if( mp->LC_ptr->name != (char *)NULL &&
         !strcmp(mp->LC_ptr->name, "Bi-Linear")) {
        for(j = 1; j <= iNO_INTEG_pts; j++) {
            mp->LC_ptr->H[j-1] = mp->ET.value/(1.0 - mp->ET.value/mp->E.value); 
        }
     }
     if(mp->LC_ptr->name != (char *)NULL &&
        !strcmp(mp->LC_ptr->name, "Ramberg-Osgood")) {
         if(mp->LC_ptr->alpha == 0) {
            printf(" **** Incorrect input in MaterialAttr{}\n");
            printf(" **** Zero alpha parameter for Ramberg-Osgood Stress-Strain Curve");
            FatalError(" **** In Add_Material_Attr()", (char *) NULL);
         }

         /* for initilization */
         for(j = 1; j <= iNO_INTEG_pts; j++)
             mp->LC_ptr->H[j-1] = mp->E.value/mp->LC_ptr->alpha/mp->LC_ptr->n;
     }

     /* [a] : Install Material into Symbol Table */

     hp = install(name);
     if(hp == NULL)
        FatalError("In Add_Material_Attr() : hp = NULL !!!\n",(char *)NULL);

     hp->u.map = mp;

#ifdef DEBUG
       printf("*** Leave Add_Material_Attr()\n");
#endif 

}


/*
 *  =============================================
 *  Fix Node in Finite Element Mesh
 *  
 *  Fix Nodes : for boundary conditions.                              
 *            : for master-slave linking of nodes                     
 *            : for rigid bodies                                       
 *  =============================================
 */

#ifdef __STDC__
void Fix_Node(double node_no, MATRIX *m)
#else
void Fix_Node(node_no, m)
double node_no;
MATRIX      *m;
#endif
{
int               i;
int    UNITS_SWITCH;
DIMENSIONS *d1, *d2;

#ifdef DEBUG
       printf("*** Enter Fix_Node()\n");
#endif

    /* [a] : Allocate and initialize space for Displ at Fixed Node */

    UNITS_SWITCH = CheckUnits();
    switch( UNITS_SWITCH ) {
      case ON:
         for(i = 1; i <= frame->no_dof; i++) {
             frame->node[(int)(node_no-1)].bound_id[i-1]   = m->uMatrix.daa[0][i-1];
             frame->node[(int)(node_no-1)].disp[i-1].value = 0.0;
             ZeroUnits( frame->node[(int)(node_no-1)].disp[i-1].dimen );
         }
         break;
      case OFF:
         for(i = 1; i <= frame->no_dof; i++) {
             frame->node[(int)(node_no-1)].bound_id[i-1]   = m->uMatrix.daa[0][i-1];
             frame->node[(int)(node_no-1)].disp[i-1].value = 0.0;
         }
         break;
      default:
         break;
    }

#ifdef DEBUG
       printf("*** Leaving Fix_Node()\n");
#endif

}


/*
 *  ===================================================================
 *  Link Nodal degrees of freedom
 *  
 *  Input :  
 *     spMatrix1 : (1xn) matrix containing list of nodes to be linked.
 *     spMatrix2 : (1x3) or (1x6) matrix containing dofs to be linked.
 *  
 *  Written By; M. Austin                                 October, 1994
 *  ===================================================================
 */

#ifdef __STDC__
void Link_Node( MATRIX *spMatrix1, MATRIX *spMatrix2 )
#else
void Link_Node( spMatrix1, spMatrix2 )
MATRIX *spMatrix1, *spMatrix2;
#endif
{
int ii, ij, ik, iNoNodes, iNoDof;
int iNodeMax, iTemp;

 /* [a] : Check Dimensions of Matrix Input */

    iNoNodes = spMatrix1->iNoColumns;
    iNoDof   = spMatrix2->iNoColumns;

    if(iNoNodes < 2)
       FatalError("Link_Node() : No of nodes < 2 \n",
                 (char *) NULL);

    if(iNoDof != 3 && iNoDof != 6)
       FatalError("Link_Node() : Second arg to LinkNode() must be either (1x3) or (1x6) matrix \n",
                 (char *) NULL);

    if(iNoDof != frame->no_dof)
       FatalError("Link_Node() : No dof's do not match \n",
                 (char *) NULL);

 /* [b] : Find max -ve index in Nodal id[frame->no_dof] array */

    iNodeMax = 0;
    for(ii = 1; ii <= frame->no_nodes; ii++) {
        for(ij = 1; ij <= frame->no_dof; ij++)
	    iNodeMax = MIN( iNodeMax, frame->node[ ii-1 ].bound_id[ ij-1 ]);
    }

 /* [c] : Add linking constraints to profile array */

    for(ii = 1; ii <= iNoNodes; ii = ii + 1) {
        iTemp = iNodeMax - 1;
        ik = spMatrix1->uMatrix.daa[0][ ii-1 ];
        for(ij = 1; ij <= frame->no_dof; ij++) {
            if(spMatrix2->uMatrix.daa[0][ij-1] == 1) {
	       frame->node[ ik-1 ].bound_id[ ij-1 ] = iTemp;
               iTemp = iTemp - 1;
	    }
        }
    }
}


/*
 *  ---------------------------------------
 *  Apply Nodal Load to Finite Element Mesh
 *  ---------------------------------------
 */

#ifdef __STDC__
void Node_Load(double node_no, MATRIX *m)
#else
void Node_Load(node_no, m)
double   node_no;
MATRIX        *m;
#endif
{
double force_val;
int    length, i;
int UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter Node_Load() : node_no = %d\n", (int) node_no);
#endif

    /* [a] : Initialize nodal load parameters */

       frame->nre             = YES ;   /* counter initialized for presence of nodal load */
       frame->no_node_loads   = frame->no_node_loads + 1;
       frame = CheckNforcesSpace(frame, (int) frame->no_node_loads);
       UNITS_SWITCH = CheckUnits();

       frame->nforces[(int)(frame->no_node_loads-1)].node_f 
       = (int) node_no;

    switch( UNITS_SWITCH ) {
      case ON:
         for(i = 1; i <= frame->no_dof; i++) {
           frame->nforces[(int)(frame->no_node_loads-1)].fn[i-1].value
           = m->uMatrix.daa[0][i-1];
           UnitsCopy( frame->nforces[(int)(frame->no_node_loads-1)].fn[i-1].dimen,
                      &(m->spColUnits[i-1]) );
         }
         break;
      case OFF:
         for(i = 1; i <= frame->no_dof; i++)
           frame->nforces[(int)(frame->no_node_loads-1)].fn[i-1].value
           = m->uMatrix.daa[0][i-1];
         break;
      default:
         break;
    }

#ifdef DEBUG
       printf("*** Leave Node_Load()\n");
#endif

}

/*
 *  ===================================================================
 *  Start_Mesh() : Allocate Data Structures before defining Mesh.
 *  End_Mesh()   : Compile Mesh before starting F.E. Solution Procedure
 *  ===================================================================
 */

void Start_Mesh() {

void Load_AISC_Sections();
void Load_AISC_Material();

   /* [a] : Load materials and AISC sections into Symbol Table */

   Load_AISC_Sections();
   Load_AISC_Material();

   /* [b] : Allocate initial memory for a frame */

   frame = FrameAlloc();
   array = Alloc_p_Array();
}

void End_Mesh() {
ELEMENT                *ep;
ELEMENT_ATTR          *eap;
SECTION_ATTR     *sap_temp;
int i, j, k, n,    node_no;
int ii, jj, kk, el_attr_no;
int       no_ldof_per_node;
int              id, index;

   /* [a] : Compute Profile of Equations including Bconds */

   frame = profile(frame, YES, EQNOS);

   /* [b] : Setup Element Attributes */
          
   frame = Set_Elmt_Attrs(frame, array);

   /* [c] : Allocate Memeory for Problem */

   TDOF = frame->no_dof * frame->no_nodes;          /* Total dof for the system              */
   MDOF = frame->no_dof * frame->no_nodes_per_elmt; /* Total dof for an element              */

   TNEQ = frame->no_eq;                             /* Total num of equations for the system */

   /*  [d] : Compute Element Destination Arrays 
    *
    *  In many cases, the no of global dofs and the local dofs are different.
    *  Also for case of mix use of elements, no of nodes per element are also
    *  varies with the element. 
    *
    *  frame->no_nodes_per_elmt = MAX no of nodes per element
    *  elmt_library[n].no_node_per_elmt = actual no of nodes per element
    *  frame->no_dof = MAX no of glocal dofs per node 
    *  elmt_library[n].no_dof = actual local dofs per node
    */

   for(i = 1; i <= frame->no_elements; i++) {  /* loop over elements */

       ep = &frame->element[i-1];
       el_attr_no = ep->elmt_attr_no;
       eap        = &frame->eattr[el_attr_no-1];
       for (ii = 1; ii <= NO_ELEMENTS_IN_LIBRARY; ii++) {
            if(!strcmp(elmt_library[ii-1].name, eap->elmt_type)) {
               n = ii;
               break;
            }
       }

       no_ldof_per_node = elmt_library[n-1].no_dof;           /* not the maximum dofs per elmt */
       NODES_PER_ELMT   = elmt_library[n-1].no_node_per_elmt; /* not the maxmum nodes per elmt */
                                                                     
       LOCAL_DOF        = no_ldof_per_node*NODES_PER_ELMT;

       ep->d_array      = iVectorAlloc((LOCAL_DOF + 1));
       ep->d_array[0]   = 0;

       index = 1;
       for(j = 1; j <= NODES_PER_ELMT; j++) {      /* loop over nodal connectivities */
           node_no = ep->node_connect[j-1];
           if(node_no != 0) {
              for(k = 1; k <= no_ldof_per_node; k++) {  /* loop over dofs */
                  jj = eap->map_ldof_to_gdof[k-1];  
                  if(jj <= 0) {
                     ep->d_array[index] = jj;
                  } else {
                     id = frame->node[node_no - 1].bound_id[jj-1];
                     ep->d_array[index] = id;
                     ep->d_array[0] = ep->d_array[0] + 1;
                  }
                  index = index + 1;
              }
           }
       }
  }

  /* [e] Set up respond temporary buffer for nonlinear analysis responds */

  SetUpRespondBuffer();

}

