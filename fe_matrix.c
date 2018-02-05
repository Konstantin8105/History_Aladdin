/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_matrix.c : Functions to solve (non)linear FE solution procedures
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

#include <stdio.h>
#include <stdlib.h>
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

extern ARRAY     *array;
extern EFRAME    *frame;

/* #define DEBUG */

/* 
 *  ------------------------------ 
 *  Form Stiffness and Mass Matrix
 *  ------------------------------ 
 */ 

MATRIX *Form_Stiffness() {
MATRIX      *stiff;
MATRIX         *Ke; 
int            *Ld;
int        elmt_no;
int            i,j;
int        iMinRow;
int      iColumnNo;
int   *ipColHeight;

   /* [a] : Compute skyline profile for global stiffness matrix */

   ipColHeight = (int *)iVectorAlloc( TNEQ );

   for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {
       Ld    = Destination_Array(frame, elmt_no);
       j = 1;
       while( Ld[j]<0 ) j++;
       iMinRow = Ld[j];
       for( i=j+1 ; i<=Ld[0] ; i++ ) {
            if( Ld[i]>0 )
                iMinRow   = MIN(iMinRow, Ld[i]);
       }

       for( i=j ; i<=Ld[0] ; i++ ) {
            if( Ld[i]>0 )  iColumnNo = Ld[i];
            else           iColumnNo = 1;
            if((1+iColumnNo-iMinRow) > ipColHeight[iColumnNo-1])
                ipColHeight[iColumnNo-1] = 1+iColumnNo-iMinRow;
       }
       free((char *) Ld);
   }

   /* [b] : Allocate memory for global stiffness matrix */

    stiff = MatrixAllocSkyline("stiff", DOUBLE_ARRAY, TNEQ, TNEQ, ipColHeight);
    free((char *)ipColHeight);

   /* [c] : Add element level stiffness matrices into global stiffness matrix */

   for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {

       array = Assign_p_Array(frame, elmt_no, array, STIFF);
       array = Assign_p_Array(frame, elmt_no, array, PROPTY);

       array = Element_Property(array);  /* passing elmt info in             */
                                         /* elmt liberary to working array   */

       Ke = Element_Matrix(array, STIFF);/* set elment stiffness/mass matrix */

       if(Rigidbody_conn(frame,array) == YES)
          Ke = Transform_Stiff_Matrix(frame,array,Ke);

       /* transfer Destination array from frame to Ld */

       Ld    = Destination_Array(frame, elmt_no); 
       stiff = Assemble_Global(stiff, frame, Ld, Ke);

       free((char *) Ld);
    } 

    return(stiff);
}

/*
 *  ========================================================
 *  Form Equivalent Nodal Load due to body force      
 *  inital stress, initial strain and surface loadings
 *  ========================================================
 */

MATRIX *Form_Equiv_Nodal_Load() {
MATRIX  *equiv_nodal_load;
MATRIX        *nodal_load;
MATRIX                *Fe; 
int                   *Ld;
int               elmt_no;
int          dof_per_elmt;
int                   i,j;

   /* [a] : Allocate memory for equivalent nodal load vector */

    dof_per_elmt     = frame->no_nodes_per_elmt*frame->no_dof;
    equiv_nodal_load = MatrixAllocIndirect("equiv_nodal_load", DOUBLE_ARRAY, TNEQ, 1);

   /* [b] : Add element level nodal loads into global load matrix */

    for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {

        array = Assign_p_Array(frame, elmt_no, array, EQUIV_NODAL_LOAD);
        Fe    = Element_Equiv(array, EQUIV_NODAL_LOAD);
        Ld    = Destination_Array(frame, elmt_no); 

        equiv_nodal_load = Assemble_Global_Load(equiv_nodal_load, frame, Ld, Fe);
    }

    return(equiv_nodal_load);
}

/* 
 *  =================================================
 *  Form Global Mass Matrix
 *  
 *  Input  : frame, p, type = CONSISTENT or LUMPED.
 *  Output : MATRIX    *mass
 *  =================================================
 */ 

#ifdef __STDC__
MATRIX *Form_Mass(MATRIX *m)
#else
MATRIX *Form_Mass(m)
MATRIX *m;
#endif
{
MATRIX       *mass;
ELEMENT        *el;
MATRIX         *Me; 
int            *Ld;
int        elmt_no;
int       rigid_no;
int      mass_type;
int       size,i,j;
int        iMinRow;
int      iColumnNo;
int   *ipColHeight;

#ifdef DEBUG
       printf("*** Enter Form_Mass()\n");
#endif

   /* [a] : Setup type of mass[][] to be computed */ 

   mass_type = (int) m->uMatrix.daa[0][0];

   /* [b] : Compute skyline profile for mass matrix */

   ipColHeight = (int *)iVectorAlloc( TNEQ );

   for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {
       Ld    = Destination_Array(frame, elmt_no);
       j = 1;
       while( Ld[j]<0 ) j++;

       iMinRow = Ld[j];
       for( i=j+1 ; i<=Ld[0] ; i++ ) {
            if( Ld[i]>0 )
                iMinRow   = MIN(iMinRow, Ld[i]);
       }

       for( i=j ; i<=Ld[0] ; i++ ) {
            if( Ld[i]>0 )  iColumnNo = Ld[i];
            else           iColumnNo = 1;
            if((1+iColumnNo-iMinRow) > ipColHeight[iColumnNo-1])
                ipColHeight[iColumnNo-1] = 1+iColumnNo-iMinRow;
       }
       free((char *) Ld);
   }

   /* [c] : Allocate memory for skyline mass matrix */

   mass = MatrixAllocSkyline("mass", DOUBLE_ARRAY, TNEQ, TNEQ, ipColHeight);
   free((char *)ipColHeight);

   /* [d] : Form global mass matrix from element mass matrices */

   for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {

      array = Assign_p_Array(frame, elmt_no, array, MASS_MATRIX);
      array = Assign_p_Array(frame, elmt_no, array, PROPTY);
      array->type = mass_type;

      array = Element_Property(array);

      Me = Element_Matrix(array, MASS_MATRIX);

      /* If element connected to rigid body transform  */
      /* mass matrix to working points of Rigid Body   */

      if(Rigidbody_conn(frame,array) == YES)
         Me = Transform_Stiff_Matrix(frame,array,Me);

      /* Me : Local element stiffness/mass matrix, size_of_stiffxsize_of_stiff */
      /* mass : Global stiffness/mass matrix, TNEQxTNEQ    */

      Ld   = Destination_Array(frame, elmt_no);
      mass = Assemble_Global(mass, frame, Ld, Me);

      free((char *) Ld);
   }

   if(frame->no_rigid == 0) {
      return (mass);
   }

   /*
    *  -----------------------------------------------------------------
    *  [e] : Compute Mass Matrix for Rigid Body Components
    * 
    *        Calculate [m] for rigid bodies; add to global mass matrix.
    * 
    *        Transform mass matrix from mass center 'c' to working pt 'p'
    *        by : Mp[][] = Tpc[][].Mc[][].Ttpc[][].
    *  -----------------------------------------------------------------
    */

    /*  ----- fix this section of code later -----------------------

    for(rigid_no = 1; rigid_no <= frame->no_rigid; rigid_no++) {

        array       = Assign_p_Array_for_Rigid_Body(frame, rigid_no, array, MASS_MATRIX);
        array->type = mass_type;                                 * LUMPED or CONSISTENT * 

        Me = Element_Matrix(array,MASS_MATRIX);
        Me = Transform_Rigid_Body_Mass_Matrix(frame, array, Me);

        Ld   = Destination_Array_for_Rigid_Body(frame, array, rigid_no, NO);
        mass = Assemble_Global(mass, frame, Ld, Me);
    }

     ---- fix this block of code later --------------------------*/

    return(mass);
}

/* 
 *  =======================================
 *  Form External and Internal Load Vectors
 *  =======================================
 */ 

MATRIX *Form_External_Load() {
MATRIX                     *load;
NODE_LOADS              *nforces;
DIMENSIONS            *units_buf;
int node_no, i, j, k,jj, counter;
int node_no_1,            length;
int                 UNITS_SWITCH;

#ifdef DEBUG
       printf("*** Enter Form_External_Load() TNEQ = %d\n", TNEQ);
#endif

   /* [a] : Allocate memory for external load vector */

   UNITS_SWITCH = CheckUnits();
   load = MatrixAllocIndirect("load", DOUBLE_ARRAY, TNEQ, 1);    

   if( UNITS_SWITCH == ON ) {
       for (i = 1; i <= TNEQ; i++)
            ZeroUnits(&(load->spRowUnits[i-1]));
   }

   /* [b] : Transfer units into units buffer */

   k = 1;
   nforces = &frame->nforces[k-1];
   node_no_1 = nforces->node_f;

   if( UNITS_SWITCH == ON ) {
       for(i = 1; i <= frame->no_nodes; i++) {
           for(j = 1; j <= frame->no_dof; j++) {
               jj = frame->node[i-1].bound_id[j-1];
               if(jj > 0)
                  UnitsCopy( &(load->spRowUnits[jj-1]), nforces->fn[j-1].dimen );
           }
       }
   }

   for(k = 1; k <= frame->no_node_loads; k++) {
       nforces = &frame->nforces[k-1];
       node_no = nforces->node_f;
       for(j = 1; j <= frame->no_dof; j++) {
           jj = frame->node[node_no-1].bound_id[j-1];
           if(jj > 0)
              load->uMatrix.daa[(int)(jj-1)][0] += nforces->fn[j-1].value;
       }
   }

   if( UNITS_SWITCH == ON ) {
       ZeroUnits(&(load->spColUnits[0]));
       load->spColUnits[0].units_type = nforces->fn[0].dimen->units_type;
   }

#ifdef DEBUG
       MatrixPrintIndirectDouble(load);
       printf("*** Leave Form_External_Load()\n");
#endif

    return(load);
}

/* 
 *  ==============================================================
 *  Form Internal Load Vector
 * 
 *  Usage : load = Form_Internal_Load(displ) 
 *       or load = Form_Internal_Load(displ, displ_incr)         
 *  
 *  Boundary Info : INPUT -- bound_id = 0, no restraint,                
 *                           bound_id > 0, restrained, displ specified  
 *                           bound_id is changed based on INPUT bound_id
 *               : OUTPUT -- bound_id < 0, restrained displ specified    
 *                           bound_id > 0, no restraint                  
 * 
 *  ==============================================================
 */ 
 
#ifdef  __STDC__
MATRIX *Form_Internal_Load( MATRIX *displ, ... )
#else
MATRIX *Form_Internal_Load(va_alist)
va_dcl
#endif
{
#ifndef __STDC__
MATRIX            *displ;
#endif

va_list          arg_ptr;
MATRIX       *displ_incr;
MATRIX             *load;
ELEMENT              *ep;
ELEMENT_ATTR        *eap;
int     node_no, elmt_no;
int     i, j, k, jj, *Ld; 
int         elmt_attr_no;
int     length1, length2;
int         UNITS_SWITCH;
int itemp;

#ifdef DEBUG
       printf("*** Enter Form_Internal_Load()\n");
#endif

   UNITS_SWITCH = CheckUnits();

   /* [a] : read input */

#ifdef __STDC__
   va_start(arg_ptr, displ);
#else
   va_start(arg_ptr);
   displ = va_arg(arg_ptr, MATRIX *);
#endif

   if(array->material_name != NULL && strcmp(array->material_name, "ELASTIC")) {
       displ_incr = va_arg(arg_ptr, MATRIX *);
   } else
       displ_incr = (MATRIX *) NULL;

   va_end(arg_ptr);

   /* [b] : Allocate memory for internal load vector */

   load = MatrixAllocIndirect("Internal Load", DOUBLE_ARRAY, TNEQ, 1);

   /* [c] : Assemble internal load vector from element level loads */

   for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) { 

       ep  = &frame->element[elmt_no-1];
       eap = &(frame->eattr[ep->elmt_attr_no-1]);
       if( !(strcmp(eap->elmt_type, "SHELL_4N")) ||
           !(strcmp(eap->elmt_type, "SHELL_8N")) ) {
             if(displ_incr == (MATRIX *) NULL)
                ep->esp->state = 0; /* Calculate internal load as elastic behaviour */
             else
                ep->esp->state = 1; /* Need check for ELASTIC PLASTIC STATE */
       }

       array = Assign_p_Array(frame, elmt_no, array, STRESS);
       array = elmlib(array, PROPTY);

       /*  Transfer fixed displacements */

       for(i = 1; i <= array->nodes_per_elmt; i++) {

           k = 1;
           node_no = ep->node_connect[i-1];
           for(j = 1; j <= frame->no_dof; j++) {

           switch((int) array->nodes_per_elmt) {
               case 2:
               case 3:
                    jj = frame->node[node_no - 1].bound_id[j-1];  
                    if(jj > 0) {

                       /* displacement  */

                       array->displ->uMatrix.daa[j-1][i-1] = displ->uMatrix.daa[jj-1][0];
                       if( UNITS_SWITCH == ON ) {
                           UnitsCopy(&(array->displ->spRowUnits[j-1]), &(displ->spRowUnits[jj-1]));
                           UnitsCopy(&(array->displ->spColUnits[i-1]), &(displ->spColUnits[0]));
                       }

                       /* displacement increments */

                       if(displ_incr != (MATRIX *) NULL) {
                          array->displ_incr->uMatrix.daa[j-1][i-1] = displ_incr->uMatrix.daa[jj-1][0];
                          if(UNITS_SWITCH == ON) {
                             UnitsCopy(&(array->displ_incr->spRowUnits[j-1]), &(displ->spRowUnits[jj-1]));
                             UnitsCopy(&(array->displ_incr->spColUnits[i-1]), &(displ->spColUnits[0]));
                          }
                       }
                    } else {

                       array->displ->uMatrix.daa[j-1][i-1] = frame->node[node_no -1].disp[j-1].value;

                       if(UNITS_SWITCH == ON ) {

                          /* displacement      */

                          UnitsCopy(&(array->displ->spRowUnits[j-1]), frame->node[node_no -1].disp[j-1].dimen);
                          UnitsCopy(&(array->displ->spColUnits[i-1]), &(displ->spColUnits[0]));

                          /* displcement increments */

                          UnitsCopy(&(array->displ_incr->spRowUnits[j-1]), frame->node[node_no -1].disp[j-1].dimen);
                          UnitsCopy(&(array->displ_incr->spColUnits[i-1]), &(displ->spColUnits[0]));
                       }
                    }
                    break;
               case 4: case 8:
                    jj = frame->node[node_no - 1].bound_id[j-1];
                    if(jj > 0) {

                       /* displacement  */

                       array->displ->uMatrix.daa[k-1][i-1] = displ->uMatrix.daa[jj-1][0];

                       if( UNITS_SWITCH == ON ) {
                           UnitsCopy(&(array->displ->spRowUnits[k-1]), &(displ->spRowUnits[jj-1]));
                           UnitsCopy(&(array->displ->spColUnits[i-1]), &(displ->spColUnits[0]));
                       }

                       /* incremental displacement  */

                       if(displ_incr != (MATRIX *) NULL) {

                          array->displ_incr->uMatrix.daa[k-1][i-1] = displ_incr->uMatrix.daa[jj-1][0];

                          if(UNITS_SWITCH == ON) {
                             UnitsCopy(&(array->displ_incr->spRowUnits[k-1]), &(displ->spRowUnits[jj-1]));
                             UnitsCopy(&(array->displ_incr->spColUnits[i-1]), &(displ->spColUnits[0]));
                          }
                       }
                    } else {
                       array->displ->uMatrix.daa[k-1][i-1] = frame->node[node_no -1].disp[j-1].value;

                       if(UNITS_SWITCH == ON ) {

                       /* displacement    */

                       UnitsCopy(&(array->displ->spRowUnits[k-1]), frame->node[node_no -1].disp[j-1].dimen);
                       UnitsCopy(&(array->displ->spColUnits[i-1]), &(displ->spColUnits[0]));

                       /* displcement increments */

                       UnitsCopy(&(array->displ_incr->spRowUnits[k-1]), frame->node[node_no-1].disp[j-1].dimen);
                       UnitsCopy(&(array->displ_incr->spColUnits[i-1]), &(displ->spColUnits[0]));
                       }
                    }
                    k = k + 1;
                    break;
               default:
                    break;
            }
       }
       }

       array = elmlib(array, LOAD_MATRIX );
       Ld    = Destination_Array(frame, elmt_no);

       if( UNITS_SWITCH == ON ) 
           ZeroUnits(&(load->spColUnits[0]));

       for(j = 1; j <= Ld[0]; j++) {
           k = Ld[j];
           if(k > 0) {
              load->uMatrix.daa[k-1][0] += array->nodal_loads[j-1].value;
              if( UNITS_SWITCH == ON )
                  UnitsCopy(&(load->spRowUnits[k-1]), array->nodal_loads[j-1].dimen);
           }
       }

#ifdef DEBUG
       printf("DIAGNOSTIC : elmt node %4d\n", elmt_no );
       dMatrixPrint("elmt diplacements", array->displ->uMatrix.daa, 3,2);
       printf("destination      :");
       for(itemp = 1; itemp <= Ld[0]; itemp++) {
           printf(" %10d : ", Ld[itemp]);
       }
       printf("\n");
       printf("elmt nodal loads :");
       for(itemp = 1; itemp <= Ld[0]; itemp++) {
           printf(" %10.2e : ", array->nodal_loads[itemp-1].value);
       }
       printf("\n");
       dMatrixPrint("Internal Load", load->uMatrix.daa, TNEQ, 1);
#endif

       free((char *) Ld);
   }

#ifdef DEBUG
       printf("*** Leave Form_Internal_Load()\n");
#endif

    return(load);
}


/* 
 *  =================================================
 *  UTILITY Functions for Stiffness and Mass matrices
 *  =================================================
 */ 

/* 
 *  -----------------------------------------
 *  [a] : Transfer Destination Array to *Ld[]
 *  -----------------------------------------
 */ 

#ifdef __STDC__
int *Destination_Array(EFRAME *frame, int elmt_no)
#else
int *Destination_Array(frame, elmt_no)
EFRAME *frame;
int    elmt_no;
#endif
{
ELEMENT *el;
int     i,no_dof_per_elmt, *Ld;
      
#ifdef DEBUG
    printf(" enter Destination_Array() \n");
    printf(" in Destination_Array() : \n");
    printf("                        : frame->no_dof = %d \n", frame->no_dof);
    printf("                        : frame->no_nodes_per_elmt = %d \n", frame->no_nodes_per_elmt);
#endif

    no_dof_per_elmt  = frame->no_dof*frame->no_nodes_per_elmt;

#ifdef DEBUG
    printf("                        : no_dof_per_elmt = %d \n", no_dof_per_elmt);
#endif

    Ld = iVectorAlloc((no_dof_per_elmt + 1));

    el = &frame->element[elmt_no - 1];
    for(i = 1; i <= no_dof_per_elmt + 1; i++) {
         Ld[i-1] = el->d_array[i-1];
    }

#ifdef DEBUG
    printf(" leaving  Destination_Array()  \n");
#endif
    return(Ld);
}

/* 
 *  ----------------------------------------------------
 *  [c] : Compute Destination Array *Ld[] for Rigid Body
 *  ----------------------------------------------------
 */ 

#ifdef __STDC__
int *Destination_Array_for_Rigid_Body(EFRAME *frp, ARRAY *p, int elmt_no, int dflg)
#else
int *Destination_Array_for_Rigid_Body(frp, p, elmt_no, dflg)
EFRAME       *frp;
ARRAY          *p;
int elmt_no, dflg;
#endif
{
ELEMENT                        *el;
int      no_nodes_per_elmt, no_dof;
int i,j,k,jj,ii,iid,mate_no, rb_ty;
int                           *vld;

       no_nodes_per_elmt = p->nodes_per_elmt;                  /* MIN(p->nodes_per_elmt, frp->no_nodes_per_elmt);*/
       no_dof            = p->dof_per_node;                    /* MIN(p->dof_per_node,frp->no_dof);*/
       vld     = iVectorAlloc(no_dof*no_nodes_per_elmt);

       rb_ty = frp->rigid[elmt_no -1].rb_type ;
       ii    = frp->rigid[elmt_no -1].in[1] ; /* node no of node conn by rbody for dofs */

       for(j=1;j<=no_dof;j++){     /* loop over dofs */
           jj = frp->rigid[elmt_no -1].rest_dof[j];
           if(jj >= 1 && rb_ty != 7) /* dofs not connected */
              vld[j] = 0;           /*k;*/
           else { /* jj=0 i:e; rest dofs & for  rb = 7; storey type rbody */
              k  = frp->node[ii - 1].bound_id[j]; /* equation nos assigned to dofs */
              vld[j] = k;
           }
       } 

       vld[0] = no_dof*no_nodes_per_elmt;

    return(vld);
}


/* 
 *  ------------------------------------------------------------------
 *  Assemble Global  Matrix : Stiffness and Mass  
 *  gmatrix : Global matrix  TNEQxTNEQ
 *  Ge      : Local matrix   array->size_of_stiffxarray->size_of_stiff
 *  ------------------------------------------------------------------
 */ 

#ifdef __STDC__
MATRIX *Assemble_Global(MATRIX *gmatrix, EFRAME *frame, int *Ld, MATRIX *Ge)
#else
MATRIX *Assemble_Global(gmatrix, frame, Ld, Ge)
MATRIX *gmatrix;
EFRAME *frame;
int    *Ld;
MATRIX *Ge;
#endif
{
ELEMENT          *el;
int j,l,k,m, size_ld;
int length1, length2;
int       iMin, iMax;
int     UNITS_SWITCH;
int linked, count, ik;

    size_ld = Ld[0];
    UNITS_SWITCH = CheckUnits();

#ifdef DEBUG
       printf("*** Enter Assemble_Global() : length of Ld = %3d : ", size_ld);
       for(j = 1; j <= size_ld; j++) {
           printf(" %3d ", Ld[j]);
       }
       printf("\n");
#endif

    for(j = 1; j <= size_ld; j++) {
        k = Ld[j];

        count = 0;
        for(ik = 1; ik <= size_ld; ik++) {
            if(k == Ld[ik]) count = count + 1;
        }

        linked = (count > 1);
        if(k > 0 && linked == 0) {
           for(l = 1; l <= size_ld; l++){
               m = Ld[l];
               if(m > 0){
	          switch( gmatrix->eRep ) {
	            case  INDIRECT :
                        gmatrix->uMatrix.daa[k-1][m-1] += Ge->uMatrix.daa[j-1][l-1];
	                break;
	            case  SKYLINE :
	                if( l >= j ) {
	                    iMin = (int) MIN( k, m );
	                    iMax = (int) MAX( k, m );
	                    if( (iMax-iMin+1) <= gmatrix->uMatrix.daa[iMax-1][0] )
                                gmatrix->uMatrix.daa[iMax-1][iMax-iMin+1]
                                += Ge->uMatrix.daa[j-1][l-1];
	                    else
	                        gmatrix->uMatrix.daa[iMax-1][iMax-iMin+1]
                                =  Ge->uMatrix.daa[j-1][l-1];
	                }
	                break;
	            default:
                        FatalError("In Assemble_Global() : Undefined spA->eRep",
                                    (char *) NULL);
	                break;
	          }
                  
                  if( UNITS_SWITCH == ON ) {
                     UnitsCopy( &(gmatrix->spColUnits[m-1]), &(Ge->spColUnits[l-1]) );
                     UnitsCopy( &(gmatrix->spRowUnits[k-1]), &(Ge->spRowUnits[j-1]) );
                  }
               }
           }
#ifdef DEBUG
           printf(" local dof : j = %d \t global dof : k = %d \n", j,k); 
#endif
        }
    }

#ifdef DEBUG
       printf("*** Leave Assemble_Global()\n");
#endif

    return(gmatrix);
}



/* 
 *  -----------------------------------------------
 *  Assemble Global Load: Equiv Nodal Load  
 *  -----------------------------------------------
 */ 

#ifdef __STDC__
MATRIX *Assemble_Global_Load(MATRIX *gmatrix, EFRAME *frame, int *Ld, MATRIX *Ge)
#else
MATRIX *Assemble_Global_Load(gmatrix, frame, Ld, Ge)
MATRIX *gmatrix;
EFRAME *frame;
int    *Ld;
MATRIX *Ge;
#endif
{
ELEMENT          *el;
int j,l,k,m, size_ld;
int           length;
int     UNITS_SWITCH;

    size_ld = Ld[0];
    UNITS_SWITCH = CheckUnits();

#ifdef DEBUG
       printf("*** Enter Assemble_Global_Load() : length of Ld = %3d : ", size_ld);
#endif

    for(j = 1; j <= size_ld; j++) {
        k = Ld[j];
        if(k > 0) {
          gmatrix->uMatrix.daa[k-1][0] += Ge->uMatrix.daa[j-1][0];
          if( UNITS_SWITCH == ON ) {
              UnitsCopy( &(gmatrix->spRowUnits[k-1]),  &(Ge->spRowUnits[j-1]) );
          }
        }
    }
    if( UNITS_SWITCH == ON )
       ZeroUnits( &(gmatrix->spColUnits[0]) );

#ifdef DEBUG
       printf("*** Leave Assemble_Global_Load()\n");
#endif

    return(gmatrix);
}

/* ======================================================= */
/* Set Element Matrix - Stiffness/Mass for Feap problem    */
/* Input  - frame pointer, element number,p array          */
/* Output - Local Element Stiffness/Mass                   */
/* ======================================================= */

#ifdef __STDC__
MATRIX *Element_Matrix(ARRAY *p, int isw)
#else
MATRIX *Element_Matrix(p, isw)
ARRAY    *p;
int     isw;
#endif
{
   p = elmlib(p,isw);
   return(p->stiff);
}

/* ========================================================== */
/* Set Element Equivalent Load -                              */
/* Input  - p working array , task id :isw = EQUIV_NODAL_LOAD */
/* Output - Local Element Equivalent nodal load               */
/* ========================================================== */

#ifdef __STDC__
MATRIX *Element_Equiv(ARRAY *p, int isw)
#else
MATRIX *Element_Equiv(p, isw)
ARRAY    *p;
int     isw;
#endif
{
   p = elmlib(p,isw);
   return(p->equiv_nodal_load);
}

/* ======================================================= */
/* Set Element Vector -  Lumped Mass for Feap problem      */
/* Input  - frame pointer, element number,p array          */
/* Output - Local Lumped  Mass                             */
/* ======================================================= */

#ifdef __STDC__
QUANTITY *Element_Vector(ARRAY *p, int isw)
#else
QUANTITY *Element_Vector(p, isw)
ARRAY     *p;
int      isw;
#endif
{

  /*printf("*** In Element_VECTOR. Flag 1. Elem_# = %d\n",p->elmt_no); */
  p = elmlib(p,isw);
  return(p->nodal_loads);
}

/* ===================================================================== */
/* Rigidbody_conn                                                        */
/* Input - frame pointer frp, array ptr checks if any node of an element */
/*         is connected to a rigid body. It labels it by adding          */
/*         a value of (one thousand * rigid_body no) to such nodes       */
/*                                                                       */
/* NOTE : Nodes attached to rigid body are scaled by a factor 1,000,000  */
/*        in the working-element p->ix[] data structure.                 */
/* ===================================================================== */

#ifdef __STDC__
int Rigidbody_conn(EFRAME *frp, ARRAY *p)
#else
int Rigidbody_conn(frp,p)
EFRAME *frp;
ARRAY  *p;
#endif
{
int i,j,k,nn, no_node,out;
out = NO;

    for(i = 1; i<=frp->no_rigid;i++) {
        no_node = frp->rigid[i-1].nodes_conn;

        for(j = 1; j<=no_node;j++){
            nn = frp->rigid[i-1].in[j];
            for(k = 1; k<=p->nodes_per_elmt;k++) {
                if(p->node_connect[k] == nn){
                   p->node_connect[k] = p->node_connect[k] + i * 1000000;
                   frp->node[nn -1].rb_num = i;
                   out = YES;
                }
            }
        }
    }

    return(out);
}


/* =================================================================== */
/* Transform Stiff Matrix of member connected to Rigid Body            */
/* transforms action at j of member to working pointt f of Rigid Body. */
/*                                                                     */
/* Transforms element S. Matrix Ke to T*Ke*Tt                          */
/* Input - EFRAME *frp;  ARRAY *p;                                     */
/*         Member Stiff Matrix MS                                      */
/* Output- Modified Member Stiff Matrix MS                             */
/* =================================================================== */

#ifdef __STDC__
MATRIX *Transform_Stiff_Matrix(EFRAME *frp, ARRAY *p, MATRIX *MS)
#else
MATRIX *Transform_Stiff_Matrix(frp, p, MS)
EFRAME     *frp;
ARRAY        *p;
MATRIX      *MS;
#endif
{
int    nn,i,j,k,rb_no,node_no,rb_ty;
double **Tfj, **Ti, **TiT, **ktit;
double cgx,cgy,cgz;

    Tfj = (double **)MatrixAllocIndirectDouble(p->dof_per_node,p->dof_per_node);
    Ti  = (double **)MatrixAllocIndirectDouble(p->size_of_stiff,p->size_of_stiff);
    TiT = (double **)MatrixAllocIndirectDouble(p->size_of_stiff,p->size_of_stiff);

    for(i = 1; i<=p->nodes_per_elmt; i++) { 
        nn = p->node_connect[i];
        if(nn > 1000000) {
           node_no = nn - rb_no * 1000000;
           rb_no = nn/1000000; /* find rigid body number */

    /* ---- FIX LATER 
           if(frp->ndm == 3) {
               if(frp->basic_elmt_type == TRUSS_3D)
                 rb_ty = TRUSS_3D;
              else
                 rb_ty = FRAME_3D;
            }
           else 
                 rb_ty = frp->rigid[rb_no -1].rb_type; 
           */ 

           cgx = p->coord[1][i].value - frp->rigid[rb_no -1].xcg.value;
           cgy = p->coord[2][i].value - frp->rigid[rb_no -1].ycg.value;

/******** fix later
           if(frp->ndm >= 3)
              cgz = p->coord[3][i] - frp->rigid[rb_no -1].zcg;
	   else
              cgz = 0.0;
****************************/

           Tfj = (double **) Transformation_Matrix(Tfj,rb_ty, p->dof_per_node,cgx,cgy,cgz);
           Tfj = (double **)       Modify_T_Matrix(Tfj,  frp, p->dof_per_node, rb_no);
       }
       else 
         Tfj = Transformation_Matrix(Tfj,rb_ty, p->dof_per_node,0.0,0.0,0.0);

       /* get other terms by symm          */
       /* printf("node at end = %d\n", i); */

       for(j = 1; j<=p->dof_per_node; j++) {
          for(k = 1; k<=p->dof_per_node; k++) {
               Ti[(i-1)*p->dof_per_node+j][ (i-1)*p->dof_per_node+ k] = Tfj[j][k];
              TiT[(i-1)*p->dof_per_node+k][ (i-1)*p->dof_per_node+ j] = Tfj[j][k];
          }
       }
    }

   /* ------------------ */
   /* Transformed Matrix */
   /* ------------------ */

   ktit = (double **) dMatrixMult(MS->uMatrix.daa,p->size_of_stiff,p->size_of_stiff, TiT, p->size_of_stiff, p->size_of_stiff);
   MS->uMatrix.daa   = (double **) dMatrixMultRep(MS->uMatrix.daa, Ti, p->size_of_stiff, p->size_of_stiff, ktit, p->size_of_stiff, p->size_of_stiff);

   MatrixFreeIndirectDouble(Tfj, p->dof_per_node);
   MatrixFreeIndirectDouble(Ti, p->size_of_stiff);
   MatrixFreeIndirectDouble(TiT, p->size_of_stiff);
   MatrixFreeIndirectDouble(ktit, p->size_of_stiff);

   return(MS);
}


/* ============================================================ */
/* Transform Mass Matrix of RigidBody:                          */
/* Transforms R body Mass Matrix M =p->s to Mp = Tpc.M.Ttpc.    */
/* transforms action from                                       */
/*  mass centres   to  working point of Rigid Body              */
/*                                                              */
/*  here pcx,pcy,pcz are cordinates of mass centre c of RBODY   */
/*                       wrt working point p of RBODY           */
/* Input -       frame pointer frp, Array ptr p.                */
/*                Member Stiff Matrix MS                        */
/* Output-       Modified Member Stiff Matrix MS                */
/* ============================================================ */

#ifdef __STDC__
MATRIX *Transform_Rigid_Body_Mass_Matrix(EFRAME *frp, ARRAY *p, MATRIX *MS)
#else
MATRIX *Transform_Rigid_Body_Mass_Matrix(frp, p, MS)
EFRAME   *frp;
ARRAY    *p;
MATRIX   *MS;
#endif
{
int nn, i,j,k,rb_no,node_no,rb_ty;
double **Tpc, **Ttpc,**MTtpc;
double pcx,pcy,pcz;

     /*----------------------------------------------------------------------*/
     /* Alloc_Matrix                                                         */
     /*----------------------------------------------------------------------*/

        Tpc  = MatrixAllocIndirectDouble(p->dof_per_node,p->dof_per_node);
        Ttpc = MatrixAllocIndirectDouble(p->dof_per_node,p->dof_per_node);

     /*----------------------------------------------------------------------*/
     /* find rigid body number */
     /*----------------------------------------------------------------------*/

     rb_no = p->elmt_no;         
     rb_ty = frp->rigid[rb_no -1].rb_type;

     /*----------------------------------------------------------------------*/
     /* Coordinates of mass centre c  wrt working point p(geom centre)       */
     /*    = global coord of cg c - global coord of wkg pt  p                */ 
     /*----------------------------------------------------------------------*/

     pcx =     p->work_section[4].value - frp->rigid[rb_no -1].xcg.value;
     pcy =     p->work_section[5].value - frp->rigid[rb_no -1].ycg.value;
     pcz = 0.0;
/***********
     if( frp->ndm >= 3)
        pcz =    p->work_section[6] - frp->rigid[rb_no -1].zcg;
***************/
     Tpc = (double **) Transformation_Matrix(Tpc,rb_ty, p->dof_per_node,pcx,pcy,pcz);

     /*----------------------------------------------------------------------*/
     /* Modify depending on the dofs constrained   */
     /*----------------------------------------------------------------------*/

     Tpc = (double **) Modify_T_Matrix(Tpc, frp,p->dof_per_node,rb_no);

    /* Transpose of Tpc matrix */

     for(j = 1; j<=p->dof_per_node;j++) 
       for(k = 1; k<=p->dof_per_node;k++) 
          Ttpc[j][k] = Tpc[k][j];
      
      /* Product : M.Ttpc           */
      /* Product : Mp = Tpc. M.Ttpc */

      MTtpc = (double **) dMatrixMult(MS->uMatrix.daa, p->dof_per_node, p->dof_per_node, Ttpc, p->dof_per_node, p->dof_per_node);
      MS->uMatrix.daa    = (double **) dMatrixMultRep(MS->uMatrix.daa, Tpc, p->dof_per_node, p->dof_per_node, MTtpc, p->dof_per_node, p->dof_per_node);

      MatrixFreeIndirectDouble(Tpc, p->dof_per_node);
      MatrixFreeIndirectDouble(Ttpc, p->dof_per_node);
      MatrixFreeIndirectDouble(MTtpc, p->dof_per_node);

      return(MS);
}


/* ======================================================*/
/* Transformation matrix Tfj (=Tpj) for point j to f(=p) */ 
/* Action:  Ap = Tpj . Aj; Disp:    Dj = Ttpj . Dp       */
/*   -- for framed structures and conn rigid bodies      */
/* function double **Transformation_Matrix               */
/*                    (t,type,no_dof,cx,cy,cz)              */
/*  input;  t         :allocated matrix                  */
/*          type      : type of element                  */
/*          no_dof       : no of dof                        */
/*          cx,cy,cz  : coordinates of  point j wrt      */
/*                                    reference point p  */
/* =======================================================*/

#ifdef __STDC__
double **Transformation_Matrix(double **t, char *type, int no_dof, double cx, double cy, double cz)
#else
double **Transformation_Matrix(t, type,no_dof,cx,cy,cz)
double **t;
char   *type;
int    no_dof;
double cx,cy,cz;
#endif
{
int i,j,n;
     for(n = 0; n < NO_ELEMENTS_IN_LIBRARY; n++) { 

         if(!strcmp(type, elmt_library[n].name))
         break;
       }

       for(i = 1; i <= elmt_library[n].no_dof; i++)
       for(j = 1; j <= elmt_library[n].no_dof; j++)
           t[i][j] = 0.0;

       if(elmt_library[n].no_dof == 3) {
        /*2d beam:2d bar:2d frame: grid*/

             t[1][1] = 1.0;
             t[2][2] = 1.0;
             if(!strcmp(type,"BEAM_2D"))
                t[2][1] = cx;
             else if(!strcmp(type,"TRUSS_2D")) {
                t[3][1] = -cy;
                t[3][2] =  cx;
             }
             else if(!strcmp(type,"FRAME_2D")){
                t[3][1] = -cy;
                t[3][2] =  cx;
                t[3][3] =  1;
             }
             else {                               /*  GRID_XY */
                t[1][3] =  cy;
                t[2][3] = -cx;
                t[3][3] =  1;
             }

          }

       if(elmt_library[n].no_dof == 6) {

             for(i = 1; i <= 6; i++)
                 t[i][i] = 1.0;
     
             if(!strcmp(type,"TABLE_XZ")) {   /* 3d space table o== r storey */
                t[3][1] =  cz;
                t[3][2] = -cx;
             }
             else if(!strcmp(type,"TRUSS_3D")) {   /* space truss */
                t[5][1] =  cz;
                t[6][1] = -cy;
                t[4][2] = -cz;
                t[6][2] =  cx;
                t[4][3] =  cy;
                t[5][3] = -cx;
             }
             else {                               /* space frames */
                for(i = 4; i <= 6; i++)
                    t[i][i] = 1.0;
                t[4][2] = -cz;
                t[4][3] =  cy;
                t[5][1] =  cz;
                t[5][3] = -cx;
                t[6][1] = -cy;
                t[6][2] =  cx;
             }
     }

     return(t);
}


/* 
 *  =====================================================================
 *  Modify Transformation_Matrix(T) : makes rows and columns to zero for
 *  dof not connected to rigid body. Diagonal terms are retained as unity 
 *  
 *  Input - frame pointer frp, prob case  
 *  =====================================================================
 */ 

#ifdef __STDC__
double **Modify_T_Matrix(double **T, EFRAME *frp, int no_dof, int rb_no)
#else
double **Modify_T_Matrix(T,frp, no_dof,rb_no)
double **T;
EFRAME *frp;
int    no_dof, rb_no;
#endif
{
int nn, i,j,k,ndof;
int rb_elmt_type;

    rb_elmt_type  = frp->rigid[rb_no -1].rb_type;

    /* For Storey type rigid body                     */
    /* all dofs(=3) assumed to be connected to R.Body */

    if(rb_elmt_type == 7) { /* Storey type rigid body ;dofs=3 */ 
                            /* all dofs are connected hence returned unaltered     */
                            /* printf("** In Modify T_Mat: Storey Type   \n");  */
       return(T);
    }

    /* general case */

    for(j = 1; j<= no_dof;j++){
        ndof = frp->rigid[rb_no -1].rest_dof[j];
        if(ndof >= 1) { /* dofs not connected by Rigid body are = 1 */
           for(k = 1; k <= no_dof; k++)
               T[j][k] =0.0;
               T[j][j] = 1.0;
        }
    }

    return(T);
}

/* 
 *  ============================================
 *  Set Element Properties for Feap problem        
 *  Input   - p array                               
 *  output  - p array                               
 *  ============================================
 */ 

#ifdef __STDC__
ARRAY *Element_Property(ARRAY *p)
#else
ARRAY *Element_Property(p)
ARRAY *p;
#endif
{
   p = elmlib(p,PROPTY);
   return(p);
}

/* ================================================ */
/* Set Material Properties for Feap problem         */
/* Input  - p array                                 */
/* output  - p array                                */
/* ================================================ */

#ifdef __STDC__
ARRAY *Mate_Property(ARRAY *p)
#else
ARRAY *Mate_Property(p)
ARRAY *p;
#endif
{
int isw = 15;

    p = elmlib(p,isw);
    return(p);
}

/* ================================================== */
/* Assemble Element Property Vector                   */
/* Input  - frame pointer frp, prob case              */
/* ================================================== */

#ifdef __STDC__
ARRAY *Eload_Property(ARRAY *p)
#else
ARRAY *Eload_Property(p)
ARRAY *p;
#endif
{
int  isw = 16;

    p = elmlib(p,isw) ;
    return(p);
}
 
#define PLASTIC 2

/* ======================================================*/
/* function  double   *Transform_Force(frp,F,nn )        */
/*  Calculates  Transform_Force Matrix  Ti   for         */
/*  Relation:        Abi = Ti . Ami                      */
/*  Converts the actions in Ami(at end of the member)    */
/*  to the statically equivalent actions in Abi          */
/*         ( at the working points of the rigid bodies ) */
/* Input - frame pointer frp,                            */
/*      F: Force Vector for the problem which is partially */
/*            modified for node nu 'nn' in transforming  */
/*            node force from this point to wkg pt of RB */
/* ===================================================== */

#ifdef __STDC__
QUANTITY *Transform_Force(EFRAME *frp, QUANTITY *F, int nn )
#else
QUANTITY *Transform_Force(frp,F,nn )
EFRAME   *frp;
QUANTITY *F;
int nn ;
#endif
{
int      i,k,rb_no,rb_elmt_type;
double   **Fj, **Ff, **Tfj;
double   cgx,cgy,cgz ,temp;

   k   = frp->no_dof;       /*for general case*/
   Fj  = MatrixAllocIndirectDouble(k,1);
   Tfj = MatrixAllocIndirectDouble(k,k);

   rb_no = frp->node[nn -1].rb_num;
   rb_elmt_type = 7;
   rb_elmt_type = frp->rigid[rb_no -1].rb_type;

   /* ------------------------------------------------------------------- */
   /* compute coordinates of node wrt ref point- cg(wkg pt) of rigid body */
   /* ------------------------------------------------------------------- */

    cgx = frp->node[nn-1].coord[1].value - frp->rigid[rb_no -1].xcg.value;
    cgy = frp->node[nn-1].coord[2].value - frp->rigid[rb_no -1].ycg.value;

/********
    if(frp->ndm >=3)
       cgz  =   frp->node[nn-1].coord[3] - frp->rigid[rb_no -1].zcg;
    else
       cgz  = 0.0;

**********/
    /* Transformation Matrix */

    Tfj = (double **)Transformation_Matrix(Tfj,(char *)rb_elmt_type,frp->no_dof, cgx,cgy,cgz);
    Tfj = (double **)Modify_T_Matrix(Tfj, frp,frp->no_dof ,rb_no);

    /* Transformation Matrix- store for further use */

    frp->node[nn -1].TrT->uMatrix.daa = (double **) dMatrixTranspose(Tfj,frp->no_dof,frp->no_dof);

    for(i = 1; i<=k;i++)
        Fj[i][1] = F[(nn  - 1) * frp->no_dof +i].value;

    Ff = (double **) dMatrixMult(Tfj, frp->no_dof, frp->no_dof,Fj, frp->no_dof, 1);
    for(i = 1; i<=k;i++)
        F[(nn  - 1) * frp->no_dof+i].value = Ff[i][1];

    MatrixFreeIndirectDouble(Ff, frp->no_dof);
    MatrixFreeIndirectDouble(Fj, frp->no_dof);
    MatrixFreeIndirectDouble(Tfj, frp->no_dof);

    return(F);
}


/* 
 *  ============================================================================ 
 *  Bound_Disp : Checks for displacement at restraint boundary nodes for an elmt             
 *  
 *  Input - frame pointer frp, element no  
 *  ============================================================================ 
 */ 

#ifdef __STDC__
int Bound_Disp(EFRAME *frp, ARRAY *p, int elmt_no)
#else
int Bound_Disp(frp,p, elmt_no)
EFRAME  *frp;
ARRAY   *p;
int     elmt_no;
#endif
{
ELEMENT *el;
double  displ;
int     i,j,k,nen, no_dof,out;

    el  = &frp->element[elmt_no -1];      /* element ptr   */
    p   = Assign_p_Array(frp,elmt_no,p, 10);
    nen = p->nodes_per_elmt;
    no_dof = p->dof_per_node;   /*ed feb 6  MIN(p->dof_per_node,frp->no_dof)*/

    out = NO;         /* default outcome  initialized */

    for(i = 1; i<=nen;i++){
        for(j = 1; j<=no_dof;j++){
            k = el->d_array[(i-1) *no_dof + j];
            if(k<0) {           /* fixed or restrained dofs */
               displ = frp->node[el->node_connect[i]  - 1].disp[j].value; 
               p->nodal_loads[(i-1) * no_dof + j].value = displ;
               if(displ != 0.0) out = YES;
            }
            else
               p->nodal_loads[(i-1) * no_dof + j].value = 0.0;
        }
    }

    return(out);
}

/* =========== */
/* Modify_Load */
/* =========== */

#ifdef __STDC__
QUANTITY *Modify_Load(QUANTITY *b, EFRAME *frp, ARRAY *p, int elmt_no)
#else
QUANTITY *Modify_Load(b, frp, p, elmt_no)
QUANTITY *b;
EFRAME   *frp;
ARRAY    *p;
int      elmt_no;
#endif
{
ELEMENT  *el;
MATRIX   *Ke;
int      *Ld;
int      i,j,k;

    el = &frp->element[elmt_no -1];
    p  = Assign_p_Array(frp,elmt_no,p, LOAD_MATRIX);
    Ke = Element_Matrix(p,STIFF);
    Ld = el->d_array;

    for(i = 1; i<=p->size_of_stiff;i++) {
        k = Ld[i];
        if(k>0) {
           for(j = 1; j<=p->size_of_stiff;j++)
               b[k].value = b[k].value - Ke->uMatrix.daa[i][j] * p->nodal_loads[j].value;
        }
    }

    return(b);
}


/* ============================================== */
/*   Function to Add input loadvector pl          */
/*                      to loadvector fvector     */
/* ============================================== */

#ifdef __STDC__
QUANTITY *Addload_Vector(EFRAME *frp, ARRAY *p, QUANTITY *fvector,  QUANTITY *pl, int elmt_no)
#else
QUANTITY *Addload_Vector(frp,p,fvector, pl, elmt_no)
EFRAME   *frp;
ARRAY    *p;
QUANTITY *fvector, *pl;
int      elmt_no;
#endif
{
int      j,k;
int      i,iid,ii;
ELEMENT  *el;
int      *ldofg;
char     *el_type;
int      elmt_attr_no;

    el           = &frp->element[elmt_no -1];
    elmt_attr_no = el->elmt_attr_no;
    el_type      = frp->eattr[elmt_attr_no-1].elmt_type;
    ldofg        = frp->eattr[elmt_no-1].map_ldof_to_gdof;

    for(i=1;i<=p->nodes_per_elmt;i++){
        ii  = el->node_connect[i];
        iid = (ii-1) * frp->no_dof;
        for(j = 1; j<= p->dof_per_node; j++){
            k = ldofg[j-1];
            if(k > 0)
              fvector[iid + k].value += pl[(i-1)*p->dof_per_node+ j].value;
        }
    }

    return(fvector);
}

/* =============================== */
/* Function to Assemble Nodal Load */
/* =============================== */

#ifdef __STDC__
QUANTITY *Assemble_Nodal_Load(EFRAME *frp, QUANTITY *fv)
#else
QUANTITY *Assemble_Nodal_Load(frp, fv)
EFRAME   *frp;
QUANTITY *fv;
#endif
{
NODE               *np;
NODE_LOADS        *nlp;
int i,node, dof, nload;

 /* Initialize load vector */

    for(dof=1; dof<=frp->no_eq; dof++)
        fv[dof].value = 0.;

    for(nload = 1; nload <= frp->no_node_loads; nload++) {
        nlp  = &frp->nforces[nload-1];
        node = nlp->node_f;
        for(i =1; i <= frp->no_dof; i++)
            fv[(node -1) *  frp->no_dof + i ].value = nlp->fn[i].value;
    }

    return (fv);
}

/* ======================================= */
/*   Function to Assemble  Gravity  Load   */
/* ======================================= */

#ifdef __STDC__
QUANTITY *Assemble_Gravity_Load(EFRAME *frp, QUANTITY *fv)
#else
QUANTITY *Assemble_Gravity_Load(frp, fv)
EFRAME   *frp;
QUANTITY *fv;
#endif
{
NODE               *np;
NODE_LOADS        *nlp;
int i,node, dof, nload;

   /* FILL IN DETAILS LATER */

   return (fv);
}

/* =========================================== */
/*   Function to Apply the Boundary conditions */
/*   to the global stiffness matrix            */
/* =========================================== */

#ifdef __STDC__
double **Boundary_Conditions(EFRAME *frp)
#else
double **Boundary_Conditions(frp)
EFRAME *frp;
#endif
{
NODE                *np;
NODE_LOADS         *nlp;
int  i,node, dof, nload;

   for(node=1;node<=frp->no_nodes;node++) {
       np   = &frp->node[node-1];
       for(i =1; i<= frp->no_dof; i++) {
       if(np->bound_id[i] >0) {
           for(dof=1;dof<=TDOF;dof++) {
               K[dof][(node-1)*frp->no_dof +i] = 0.;
               K[(node-1)*frp->no_dof + 1][dof] = 0.;
           }

           K[(node-1)*frp->no_dof + 1][(node-1)*frp->no_dof+ 1] = 0.;
           F[(node-1)*frp->no_dof + 1]   = 0.;
       }
       }
   }

   return (K);
}

/* =======================================  */
/* Function to Assemble  Centrifugal  Load  */
/* =======================================  */

#ifdef __STDC__
QUANTITY *Assemble_Ctrfgl_Load(EFRAME *frp, QUANTITY *fv)
#else
QUANTITY *Assemble_Ctrfgl_Load(frp,fv)
EFRAME   *frp;
QUANTITY *fv;
#endif
{
NODE               *np;
NODE_LOADS        *nlp;
int i,node, dof, nload;

 /* to be changed */

   for(dof=1; dof<=frp->no_eq; dof++)
      fv[dof].value   = 0.;

   for(nload = 1; nload <= frp->no_node_loads; nload++) {

     nlp  = &frp->nforces[nload-1];
     node = nlp->node_f;
         for(i =1; i <= frp->no_dof; i++)
             fv[(node -1) *  frp->no_dof + i ].value = nlp->fn[i].value;

   }

   return (fv);
}


/* 
 *  =============================================================== 
 *  Functions for design rule checking and analysis post-processing
 *  =============================================================== 
 *
 *  Get_Coord() : Retrieve coordinates for a particular node.
 *
 *  Input  : Matrix *m     -- node number.
 *  Output : Matrix *coord -- matrix of nodal coordinates.
 *  =============================================================== 
 */ 

#ifdef __STDC__
MATRIX *Get_Coord(MATRIX *m)
#else
MATRIX *Get_Coord(m)
MATRIX *m;
#endif
{
MATRIX *coord;
int    nodeno;
int        ii;

       nodeno = (int) m->uMatrix.daa[0][0];
       coord  = MatrixAllocIndirect("Node Coord", DOUBLE_ARRAY, 1, frame->no_dimen);

       if( CheckUnits() == ON ) {
          for( ii=0 ; ii<frame->no_dimen ; ii++ ) {
             coord->uMatrix.daa[0][ii] = frame->node[nodeno-1].coord[ii].value;
             UnitsCopy( &(coord->spColUnits[ii]), frame->node[nodeno-1].coord[ii].dimen );
          }
          ZeroUnits( &(coord->spRowUnits[0]) );
       }
       else {
          for( ii=0 ; ii<frame->no_dimen ; ii++ )
             coord->uMatrix.daa[0][ii] = frame->node[nodeno-1].coord[ii].value;
       }

    return(coord);
}

/* 
 *  ================================================================= 
 *  Get_Node() : Retrieve list of nodes attached to an element.
 *
 *  Input  : Matrix *m       -- element no.
 *  Output : Matrix *connect -- matrix of nodes connected to element.
 *  ================================================================= 
 */ 

#ifdef __STDC__
MATRIX *Get_Node(MATRIX *m)
#else
MATRIX *Get_Node(m)
MATRIX *m;
#endif
{
MATRIX *connect;
int      elmtno;
int          ii;

       elmtno   = (int) m->uMatrix.daa[0][0];
       connect  = MatrixAllocIndirect("Node Connect", DOUBLE_ARRAY, 1, frame->no_nodes_per_elmt);

       if( CheckUnits() == ON ) {
          for( ii=0 ; ii<frame->no_nodes_per_elmt ; ii++ ) {
             connect->uMatrix.daa[0][ii] = (double) frame->element[elmtno-1].node_connect[ii];
             ZeroUnits( &(connect->spColUnits[ii]) );
          }
          ZeroUnits( &(connect->spRowUnits[0]) );
       }
       else {
          for( ii=0 ; ii<frame->no_nodes_per_elmt ; ii++ )
             connect->uMatrix.daa[0][ii] = (double) frame->element[elmtno-1].node_connect[ii];
       }

    return(connect);
}

/* 
 *  ================================================================= 
 *  Get_Displ() : Return displacements at a particular node.
 *
 *  Input  : Matrix *m     -- node no.
 *  Output : Matrix *displ -- matrix of nodes connected to element.
 *  ================================================================= 
 */ 

#ifdef __STDC__
MATRIX *Get_Displ(MATRIX *m1, MATRIX *m2)
#else
MATRIX *Get_Displ(m1,m2)
MATRIX *m1, *m2;
#endif
{
MATRIX *displ;
int    nodeno;
int     ii,jj;
int UNITS_SWITCH, UnitsType;

#ifdef DEBUG
       printf("*** Enter Get_Displ()\n");
#endif

       nodeno = (int) m1->uMatrix.daa[0][0];
       displ  = MatrixAllocIndirect("Node Displ", DOUBLE_ARRAY, 1, frame->no_dof);
       UNITS_SWITCH = CheckUnits();
       UnitsType    = CheckUnitsType();

       switch( UNITS_SWITCH ) {
          case ON:
               for( ii=0 ; ii < frame->no_dof ; ii++ ) {
                    jj = frame->node[nodeno-1].bound_id[ii];
                    if(jj > 0) {
                       if(m2->spRowUnits[jj-1].units_name != NULL ) {
                          UnitsTypeConvert(&(m2->spRowUnits[jj-1]), UnitsType);
                       }
                       RadUnitsSimplify( &(m2->spRowUnits[jj-1]) );
                       frame->node[nodeno-1].disp[ii].value = m2->uMatrix.daa[jj-1][0];
                       UnitsCopy( frame->node[nodeno-1].disp[ii].dimen, &(m2->spRowUnits[jj-1]) );
                    }
                    else {
                       if ( ii < frame->no_dimen ) {
                            UnitsCopy( frame->node[nodeno-1].disp[ii].dimen, 
                                       frame->node[nodeno-1].coord[ii].dimen );
                       }
                    }
               }
               break;
          case OFF:
               for( ii=0 ; ii<frame->no_dof ; ii++){
                    jj = frame->node[nodeno-1].bound_id[ii];
                    if(jj > 0)
                       frame->node[nodeno-1].disp[ii].value = m2->uMatrix.daa[jj-1][0];
               }
               break;
          default:
               break;
       } 

       if( UNITS_SWITCH == ON ) {
          for( ii=0 ; ii<frame->no_dof ; ii++ ) {
             displ->uMatrix.daa[0][ii] = frame->node[nodeno-1].disp[ii].value;
             UnitsCopy( &(displ->spColUnits[ii]), frame->node[nodeno-1].disp[ii].dimen );
          }
          ZeroUnits( &(displ->spRowUnits[0]) );
       }
       else {
          for( ii=0 ; ii<frame->no_dof ; ii++ )
             displ->uMatrix.daa[0][ii] = frame->node[nodeno-1].disp[ii].value;
       }

#ifdef DEBUG
       printf("*** Leave Get_Displ()\n");
#endif

    return(displ);
}

/* 
 *  ================================================================= 
 *  Get_Stress() : Return stresses within an element.
 *
 *  Input  : Matrix *m      -- node no.
 *  Output : Matrix *stress -- matrix of stresses within element.
 *  ================================================================= 
 */ 

#ifdef __STDC__
MATRIX *Get_Stress(MATRIX *m1, MATRIX *m2)
#else
MATRIX *Get_Stress(m1,m2)
MATRIX *m1, *m2;
#endif
{
MATRIX        *stress;
ELEMENT           *ep;
ELEMENT_ATTR     *eap;
DIMENSIONS *dp_length;
int UNITS_SWITCH;
int      elmt_no;
int node_no, elmt_attr_no;
int  i,j,k,ii,jj, iNoCols;

    elmt_no = (int) m1->uMatrix.daa[0][0];

#ifdef DEBUG
       printf("*** Enter Get_Stress() : elmt_no = %d\n", elmt_no );
#endif

    /* Allocate working array for elmt_no and assign element properties */

    UNITS_SWITCH = CheckUnits();
    array = Assign_p_Array(frame, elmt_no, array, STRESS);
    array = elmlib(array, PROPTY);

    /* Transfer Fixed Displacements */

    ep           = &frame->element[elmt_no-1];
    elmt_attr_no = ep->elmt_attr_no;  
    eap          = &frame->eattr[elmt_attr_no-1];

    for(i=1; i <= array->nodes_per_elmt; i++) {
        k = 1; 
        node_no = ep->node_connect[i-1];
        for(j = 1; j <= array->dof_per_node; j++) {
            switch( (int) array->nodes_per_elmt ) {
                case 2:
                case 3:
                     ii = eap->map_ldof_to_gdof[j-1];
                     jj = frame->node[node_no - 1].bound_id[ii-1];
                     if(jj > 0) {
                        array->displ->uMatrix.daa[j-1][i-1] = m2->uMatrix.daa[jj-1][0];
                        if( UNITS_SWITCH == ON ) {
                            UnitsCopy(&(array->displ->spRowUnits[j-1]), &(m2->spRowUnits[jj-1]));
                            ZeroUnits(&(array->displ->spColUnits[i-1]));
                        }
                     } else {
                        array->displ->uMatrix.daa[j-1][i-1]
                           = frame->node[node_no -1].disp[ii-1].value;
                        if( UNITS_SWITCH == ON ) {
                            UnitsCopy(&(array->displ->spRowUnits[j-1]),
                                      frame->node[node_no -1].disp[ii-1].dimen);
                            ZeroUnits(&(array->displ->spColUnits[i-1]));
                        }
                     }
                     break;
                case 4:
                case 8:
                     ii = eap->map_ldof_to_gdof[k-1];
                     jj = frame->node[node_no - 1].bound_id[ii-1];
                     if(jj > 0) {
                        array->displ->uMatrix.daa[k-1][i-1] = m2->uMatrix.daa[jj-1][0];
                        if( UNITS_SWITCH == ON ) {
                            UnitsCopy(&(array->displ->spRowUnits[k-1]), &(m2->spRowUnits[jj-1]));
                            ZeroUnits( &(array->displ->spColUnits[i-1]) );
                        }
                     } else {
                        array->displ->uMatrix.daa[k-1][i-1]
                           = frame->node[node_no -1].disp[ii-1].value;
                        if( UNITS_SWITCH == ON ) {
                            UnitsCopy( &(array->displ->spRowUnits[k-1]),
                                       frame->node[node_no -1].disp[ii-1].dimen);
                            ZeroUnits( &(array->displ->spColUnits[i-1]) );
                        }
                     }
                     k = k + 1;
                     break;
                 default:
                     break;
           }
       }
   }

   /* Allocate memory for the array of element stresses/internal forces. */

   switch ( frame->no_dimen ) {
       case 1:
            iNoCols = 2;
            break;
       case 2:
            iNoCols = 5;
            break;
       case 3:
            iNoCols = 9;
            break;
       default:
            break;
   }

   stress = MatrixAllocIndirect("Element Stress", DOUBLE_ARRAY,
                                array->nodes_per_elmt, iNoCols );

   /* Retrieve element level stresses/forces */

#ifdef DEBUG
       printf("*** Go to elmlib( array, STRESS_MATRIX )\n");
#endif

   array = elmlib(array, STRESS_MATRIX );

#ifdef DEBUG
       printf("*** Back from elmlib( array, STRESS_MATRIX )\n");
#endif

   /* Transfer "units" to "stress" matrix */

   if( UNITS_SWITCH == ON ) {
       for( i = 1; i <= 4*frame->no_dimen - 3 ; i++ ) {
            UnitsCopy( &(stress->spColUnits[ i - 1 ]),
                       &(array->stress->spColUnits[i-1]) );
       }

       for( j = 1; j <= array->nodes_per_elmt ; j++ ) {
            ZeroUnits( &(stress->spRowUnits[j-1]) );
       }
   }

   /* Transfer "coordinate" and "force/stress values" to "stress" matrix */

   for( i = 1; i <= array->nodes_per_elmt; i++ )
   for( j = 1; j <= 4*frame->no_dimen - 3 ; j++ )
        stress->uMatrix.daa[i-1][j-1] = array->stress->uMatrix.daa[i-1][j-1];

#ifdef DEBUG
       printf("*** Leave Get_Stress()\n");
#endif

    return(stress);
}

/* 
 *  ================================================================= 
 *  Get_Dof() : Return dof associated with a particular node.
 *
 *  Input  : Matrix *m      -- node no.
 *  Output : Matrix *gdof   -- matrix of degrees of freedom.
 *  ================================================================= 
 */ 

#ifdef __STDC__
MATRIX *Get_Dof(MATRIX *m)
#else
MATRIX *Get_Dof(m)
MATRIX *m;
#endif
{
MATRIX  *gdof;
int    nodeno;
int        ii;

   nodeno = (int) m->uMatrix.daa[0][0];
   gdof   = MatrixAllocIndirect("Node Global Dof", DOUBLE_ARRAY, 1, frame->no_dof);

   if( CheckUnits() == ON ) {
       for( ii=0 ; ii<frame->no_dof ; ii++ ) {
            gdof->uMatrix.daa[0][ii] = (double) frame->node[nodeno-1].bound_id[ii];
            ZeroUnits( &(gdof->spColUnits[ii]) );
       }
       ZeroUnits( &(gdof->spRowUnits[0]) );
   } else {
       for( ii=0 ; ii<frame->no_dof ; ii++ )
             gdof->uMatrix.daa[0][ii] = (double) frame->node[nodeno-1].bound_id[ii];
   }

   return(gdof);
}


/* 
 *  ================================================================= 
 *  Get_Stiffness() : Return element stiffness matrix 
 *
 *  Input  : Matrix *m      -- element no.
 *  Output : Matrix *stiff  -- element stiffness matrix.
 *
 *  Note : This routine is a real hack (need to tidy up later).
 *  ================================================================= 
 */ 

#ifdef __STDC__
MATRIX *Get_Stiffness(MATRIX *m)
#else
MATRIX *Get_Stiffness(m)
MATRIX *m;
#endif
{
MATRIX      *stiff;
MATRIX         *Ke; 
MATRIX      *Ksave; 
int        elmt_no;
int        iElmtNo;

   iElmtNo = (int) m->uMatrix.daa[0][0];

   /* [a] : Compute element level stiffness matrices */

   array = Assign_p_Array(frame, iElmtNo, array, STIFF);
   array = Assign_p_Array(frame, iElmtNo, array, PROPTY);
   array = Element_Property(array); 

   /* [b] : Compute and save element level stiffness matrices */

   Ke    = Element_Matrix(array, STIFF);
   Ksave = MatrixCopy( Ke );

   return( Ksave );
}


/* 
 *  ==================================================================== 
 *  Get_Section() : Return section properties of a particular element.
 *
 *  Input  : Matrix *m       -- elementno.
 *  Output : Matrix *section -- matrix of section properties.
 *  ==================================================================== 
 *  The details are:
 *
 *      QUANTITY  Ixx;              [0]
 *      QUANTITY  Iyy;              [1]
 *      QUANTITY  Izz;              [2]
 *      QUANTITY  Ixz;              [3]
 *      QUANTITY  Ixy;              [4]
 *      QUANTITY  Iyz;              [5]
 *      QUANTITY  weight;           [6]  Section weight
 *      QUANTITY  bf;               [7]  Width of flange
 *      QUANTITY  tf;               [8]  thickness of flange
 *      QUANTITY  depth;            [9]  Section depth
 *      QUANTITY  area;             [10] Section area
 *      QUANTITY  plate_thickness;  [11]
 *      QUANTITY  tor_const;        [12] Torsional Constant J
 *      QUANTITY  rT;               [13] Section radius of gyration
 *      QUANTITY  width;            [14] Section width
 *      QUANTITY  tw;               [15] Thickness of web
 *      double    ks;               [16] Shear Section Correction Factor
 *  ==================================================================== 
 */ 

#ifdef __STDC__
MATRIX *Get_Section(MATRIX *m)
#else
MATRIX *Get_Section(m)
MATRIX *m;
#endif
{
MATRIX  *section;
int       elmtno;
int elmt_attr_no;
char       *name;
SECTION_ATTR *sap;

#ifdef DEBUG
       printf("*** Enter Get_Section()\n");
#endif

   elmtno       = (int) m->uMatrix.daa[0][0];
   elmt_attr_no = frame->element[elmtno-1].elmt_attr_no;
   name         = frame->eattr[elmt_attr_no-1].section;
   sap          = lookup(name)->u.sap;
   section      = MatrixAllocIndirect("Elmt Section", DOUBLE_ARRAY, 16, 1);

   section->uMatrix.daa[0][0]  = sap->Ixx.value;
   section->uMatrix.daa[1][0]  = sap->Iyy.value;
   section->uMatrix.daa[2][0]  = sap->Izz.value;
   section->uMatrix.daa[3][0]  = sap->Ixz.value;
   section->uMatrix.daa[4][0]  = sap->Ixy.value;
   section->uMatrix.daa[5][0]  = sap->Iyz.value;
   section->uMatrix.daa[6][0]  = sap->weight.value;
   section->uMatrix.daa[7][0]  = sap->bf.value;
   section->uMatrix.daa[8][0]  = sap->tf.value;
   section->uMatrix.daa[9][0]  = sap->depth.value;
   section->uMatrix.daa[10][0] = sap->area.value;
   section->uMatrix.daa[11][0] = sap->plate_thickness.value;
   section->uMatrix.daa[12][0] = sap->tor_const.value;
   section->uMatrix.daa[13][0] = sap->rT.value;
   section->uMatrix.daa[14][0] = sap->width.value;
   section->uMatrix.daa[15][0] = sap->tw.value;

   if( CheckUnits() == ON ) {
       ZeroUnits( &(section->spColUnits[0]) );
       if( sap->Ixx.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[0]), sap->Ixx.dimen );
       else
           ZeroUnits( &(section->spRowUnits[0]) );
       if( sap->Iyy.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[1]), sap->Iyy.dimen );
       else
           ZeroUnits( &(section->spRowUnits[1]) );
       if( sap->Izz.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[2]), sap->Izz.dimen );
       else
           ZeroUnits( &(section->spRowUnits[2]) );
       if( sap->Ixz.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[3]), sap->Ixz.dimen );
       else
           ZeroUnits( &(section->spRowUnits[3]) );
       if( sap->Ixy.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[4]), sap->Ixy.dimen );
       else
           ZeroUnits( &(section->spRowUnits[4]) );
       if( sap->Iyz.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[5]), sap->Iyz.dimen );
       else
           ZeroUnits( &(section->spRowUnits[5]) );
       if( sap->weight.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[6]), sap->weight.dimen );
       else
           ZeroUnits( &(section->spRowUnits[6]) );
       if( sap->bf.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[7]), sap->bf.dimen );
       else
           ZeroUnits( &(section->spRowUnits[7]) );
       if( sap->tf.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[8]), sap->tf.dimen );
       else
           ZeroUnits( &(section->spRowUnits[8]) );
       if( sap->depth.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[9]), sap->depth.dimen );
       else
           ZeroUnits( &(section->spRowUnits[9]) );
       if( sap->area.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[10]), sap->area.dimen );
       else
           ZeroUnits( &(section->spRowUnits[10]) );
       if( sap->plate_thickness.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[11]), sap->plate_thickness.dimen );
       else
           ZeroUnits( &(section->spRowUnits[11]) );
       if( sap->tor_const.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[12]), sap->tor_const.dimen );
       else
           ZeroUnits( &(section->spRowUnits[12]) );
       if( sap->rT.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[13]), sap->rT.dimen );
       else
           ZeroUnits( &(section->spRowUnits[13]) );
       if( sap->width.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[14]), sap->width.dimen );
       else
           ZeroUnits( &(section->spRowUnits[14]) );
       if( sap->tw.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(section->spRowUnits[15]), sap->tw.dimen );
       else
           ZeroUnits( &(section->spRowUnits[15]) );
   }

#ifdef DEBUG
       printf("*** Leave Get_Section()\n");
#endif

    return(section);
}

/* 
 *  ===================================================================== 
 *  Get_Material() : Return material properties for a particular element.
 *
 *  Input  : Matrix *m        -- elementno.
 *  Output : Matrix *material -- matrix of section properties.
 *  ===================================================================== 
 *  The details are:
 * 
 *  QUANTITY               E;    [0] Young's modulus
 *  QUANTITY               G;    [1] Shear modulus
 *  QUANTITY              fy;    [2] Yield stress
 *  QUANTITY              ET;    [3] Tangent Young's Modulus
 *  double                nu;    [4] Poission's ratio
 *  QUANTITY         density;    [5] 
 *  QUANTITY              fu;    [6] Ultimate stress
 *  QUANTITY  *alpha_thermal;    [7] thermal expansion coefficient
 *  ===================================================================== 
 */ 

#ifdef __STDC__
MATRIX *Get_Material(MATRIX *m)
#else
MATRIX *Get_Material(m)
MATRIX *m;
#endif
{
MATRIX *material;
int       elmtno;
int elmt_attr_no;
int           ii;
char       *name;
MATERIAL_ATTR *map;

#ifdef DEBUG
       printf("*** Enter Get_Material()\n");
#endif

   elmtno = (int) m->uMatrix.daa[0][0];
   elmt_attr_no = frame->element[elmtno-1].elmt_attr_no;
   name = frame->eattr[elmt_attr_no-1].material;
   map = lookup(name)->u.map;
   material  = MatrixAllocIndirect("Elmt Material", DOUBLE_ARRAY, 10, 1);

   material->uMatrix.daa[0][0] = map->E.value;
   material->uMatrix.daa[1][0] = map->G.value;
   material->uMatrix.daa[2][0] = map->fy.value;
   material->uMatrix.daa[3][0] = map->ET.value;
   material->uMatrix.daa[4][0] = map->nu;
   material->uMatrix.daa[5][0] = map->density.value;
   material->uMatrix.daa[6][0] = map->fu.value;
   material->uMatrix.daa[7][0] = map->alpha_thermal[0].value;
   material->uMatrix.daa[8][0] = map->alpha_thermal[1].value;
   material->uMatrix.daa[9][0] = map->alpha_thermal[2].value;

   if( CheckUnits() == ON ) {
       ZeroUnits( &(material->spColUnits[0]) );
       if( map->E.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[0]), map->E.dimen );
       else
           ZeroUnits( &(material->spRowUnits[0]) );
       if( map->G.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[1]), map->G.dimen );
       else
           ZeroUnits( &(material->spRowUnits[1]) );
       if( map->fy.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[2]), map->fy.dimen );
       else
           ZeroUnits( &(material->spRowUnits[2]) );
       if( map->ET.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[3]), map->ET.dimen );
       else
           ZeroUnits( &(material->spRowUnits[3]) );
       ZeroUnits( &(material->spRowUnits[4]) );
       if( map->density.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[5]), map->density.dimen );
       else
           ZeroUnits( &(material->spRowUnits[5]) );
       if( map->fu.dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[6]), map->fu.dimen );
       else
           ZeroUnits( &(material->spRowUnits[6]) );
       if( map->alpha_thermal[0].dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[7]), map->alpha_thermal[0].dimen );
       else
           ZeroUnits( &(material->spRowUnits[7]) );
       if( map->alpha_thermal[1].dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[8]), map->alpha_thermal[1].dimen );
       else
           ZeroUnits( &(material->spRowUnits[8]) );
       if( map->alpha_thermal[2].dimen != (DIMENSIONS *)NULL )
           UnitsCopy( &(material->spRowUnits[9]), map->alpha_thermal[2].dimen );
       else
           ZeroUnits( &(material->spRowUnits[9]) );
   }

#ifdef DEBUG
       printf("*** Leave Get_Material()\n");
#endif

   return(material);
}
