/*
 *  ============================================================================= 
 *  ALADDIN Version 2.1.
 *                                                                     
 *  fe_print.c : Print FE Mesh and Solution
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

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#include "defs.h"
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

/* #define DEBUG */


/*
 *  =========================
 *  Print Finite Element Mesh
 *  =========================
 */

void Print_Mesh() {
MATRIX  *m = (MATRIX *)NULL;
NODE                    *np;
ELEMENT               *elmt;
ELEMENT_ATTR           *eap;
RIGID                   *rb;
NODE_LOADS         *nforces;
ELEMENT_LOADS         *elsp;
ELOAD_LIB              *elp;
void           (*voidptr)();
int          iel_no,i,j,k,l;
int            UNITS_SWITCH;

#ifdef DEBUG
       printf("*** In Print_Mesh()\n");
#endif

    UNITS_SWITCH = CheckUnits();

    printf("\n");
    printf("===========================================================================\n");
    printf("Title : DESCRIPTION OF FINITE ELEMENT MESH                                 \n");
    printf("===========================================================================\n");
    printf("\n");

    printf("=======================\n");
    printf("Profile of Problem Size\n");
    printf("=======================\n\n");

    printf("   Dimension of Problem        = %6d\n", frame->no_dimen);
    printf("\n");
    printf("   Number Nodes                = %6d\n", frame->no_nodes);
    printf("   Degrees of Freedom per node = %6d\n", frame->no_dof);
    printf("   Max No Nodes Per Element    = %6d\n", frame->no_nodes_per_elmt);
    printf("\n");
    printf("   Number Elements             = %6d\n", frame->no_elements);
    printf("   Number Element Attributes   = %6d\n", frame->no_element_attr);
    printf("   Number Loaded Nodes         = %6d\n", frame->no_node_loads);       
    printf("   Number Loaded Elements      = %6d\n", frame->no_element_loads);       
    printf("\n");

    /* [a] : Print Nodes */
     
    switch((int) frame->no_dimen ) {
        case 3:
             printf("\n");
             printf("============================================================================================\n");
             printf("Node#      X_coord          Y_coord          Z_coord        Dx    Dy    Dz    Rx    Ry    Rz\n");
             printf("============================================================================================\n");
             printf("\n");

             switch( UNITS_SWITCH ) {
                 case ON:
                      for(i=1; i <= frame->no_nodes; i++) {
                         np = &frame->node[i-1];
                         printf("%5d %13.4e %.3s %13.4e %.3s %13.4e %.3s ", i,
                                 np->coord[0].value/np->coord[0].dimen->scale_factor,
                                 np->coord[0].dimen->units_name,
                                 np->coord[1].value/np->coord[1].dimen->scale_factor,
                                 np->coord[1].dimen->units_name,
                                 np->coord[2].value/np->coord[2].dimen->scale_factor,
                                 np->coord[2].dimen->units_name);

                         for(j=1; j <= frame->no_dof; j++)
                             printf("%5d ", np->bound_id[j-1]);

                         printf("\n");
                      }
                      break;
                 case OFF:
                      for(i=1; i <= frame->no_nodes; i++) {
                         np = &frame->node[i-1];
                         printf("%5d %13.4e     %13.4e     %13.4e     ", i,
                                 np->coord[0].value,
                                 np->coord[1].value,
                                 np->coord[2].value);

                         for(j=1; j <= frame->no_dof; j++)
                             printf("%5d ", np->bound_id[j-1]);

                         printf("\n");
                      }
                      break;
                  default:
                      break;
             }
             break;
        case 2:
             printf("\n");
             printf("===================================================\n");
             printf("Node#      X_coord         Y_coord        Tx     Ty\n");
             printf("===================================================\n");
             printf("\n");

             if( UNITS_SWITCH == ON ) {
                 for(i=1;i<=frame->no_nodes; i++) {
                     np = &frame->node[i-1];
                     printf("%5d %13.4e %.3s %13.4e %.3s ", i,
                             np->coord[0].value/np->coord[0].dimen->scale_factor,
                             np->coord[0].dimen->units_name,
                             np->coord[1].value/np->coord[1].dimen->scale_factor,
                             np->coord[1].dimen->units_name);
                     for(j=1; j<= frame->no_dof;j++)
                         printf(" %5d ", np->bound_id[j-1]);
                     printf("\n");
                 }
             } else {
                 for(i=1;i<=frame->no_nodes; i++) {
                     np = &frame->node[i-1];
                     printf("%5d %13.4e     %13.4e     ", i,
                            np->coord[0].value, np->coord[1].value);
                     for(j=1; j<= frame->no_dof;j++)
                         printf("%5d ", np->bound_id[j-1]);
                     printf("\n");
                 }
             }
             break;
         default:
             printf(">>ERROR : Ndm = %d; should be 2 or 3   \n", frame->no_dimen); 
             break;
   }

   /* [b] : Print Element Data : frame->node_per_elmt = max nodes per element */

   printf("\n\n");
   printf("=======================");
   for(i = 1; i <= frame->no_nodes_per_elmt; i++)
       printf("===========");
   printf("=======================\n");
   printf("Element#     Type      ");
   for(i=1;i <= frame->no_nodes_per_elmt; i++)
       printf("   node[%1d]", i);
   printf("      Element_Attr_Name\n");
   printf("======================================");
   for(i=1;i <= frame->no_nodes_per_elmt;i++)
       printf("==========");
   printf("============\n\n");

   for(i=1; i <= frame->no_elements; i++) {

           printf("%8d", i);
           elmt = &frame->element[i-1];
           j =  elmt->elmt_attr_no;
           printf("   %-12s", frame->eattr[j-1].elmt_type);

           for(j = 1; j <= frame->no_nodes_per_elmt; j++) {
	       if(elmt->node_connect[j-1] != 0)
                  printf("%10d", elmt->node_connect[j-1]);
	       else
                  printf("          ");
	   }

           printf("      %s\n", elmt->elmt_attr_name);
       }
       printf("\n");

    /* [c] : Print rigid body data */

       if(frame->no_rigid >= 1) {
          printf("\n\n");
          printf("------------------------ \n");
          printf("RBody# : Rigid Body Data \n");
          printf("------------------------ \n");

          for(i=1; i <=frame->no_rigid; i++) {
              rb = &frame->rigid[i-1];
          
              printf("\n");
              printf("%3d : name     = \"%s\" \n", i, rb->rbody_attr_name);
              printf("      type     = %2d\n", rb->rb_type);
              printf("      no nodes = %4d\n", rb->nodes_conn);
              printf("      connectivity : ");

              for(j=1;j <= rb->nodes_conn;j++) 
                  printf("%2d ", rb->in[j]);  

              printf("\n");  

              printf("          rest dof : ");

              for(j=1; j <= 6; j++) 
                  printf("%2d ", rb->rest_dof[j]);  

              printf("\n");  
              printf("      Centre Mass : (x,y,z) = (%12.4f, %12.4f, %12.3f)\n",rb->xcg.value,rb->ycg.value,rb->zcg.value);  
              printf("      Density   = %14.5f\n",  rb->prop[1].value);
              printf("      Thickness = %14.5f\n",  rb->prop[2].value);
              printf("      Area      = %14.5f\n",  rb->prop[3].value);  
              printf("      Ixx       = %14.5f\n",  rb->prop[7].value);
              printf("      Iyy       = %14.5f\n",  rb->prop[8].value);
              printf("      Izz       = %14.5f\n",  rb->prop[9].value);
              printf("      Ixy       = %14.5f\n",  rb->prop[10].value);
              printf("      Iyz       = %14.5f\n",  rb->prop[11].value);
              printf("      Izx       = %14.5f\n",  rb->prop[12].value);
              printf("\n");  
          }
     }

    /* [c] : Print material data */

    if(frame->no_element_attr >= 1 ) {
       printf("====================== \n");
       printf("Element Attribute Data \n");
       printf("====================== \n");

       for(i = 1; i <= frame->no_element_attr; i++) {
           eap = &frame->eattr[i-1];
           printf("\n\n"); 
           printf("Element Attribute No %3d : name     = \"%s\" \n ", i, eap->name); 
           printf("                        : section  = \"%s\" \n", eap->section); 
           printf("                         : material = \"%s\" \n ", eap->material); 
           printf("                        : type     = %s\n", eap->elmt_type); 

           for(j = 0; j <= NO_ELEMENTS_IN_LIBRARY; j++){
               if(!strcmp(eap->elmt_type, elmt_library[j].name)) {
                  (*(void (*)()) *(elmt_library[j].elmt_print) )(frame, i);
                  break;
               }
           }
       }
    }

    /* [d] : Print nodal forces */

    if(frame->no_node_loads >= 1 && (frame->no_dimen == 2)) {

       printf("                        \n");
       printf("EXTERNAL NODAL LOADINGS \n");

       nforces = &frame->nforces[0];
       if(frame->no_dof == 3) {
          if(UNITS_SWITCH == ON ) {
             for( i=1 ; i<=frame->no_dof ; i++ )
                  UnitsSimplify( nforces->fn[i-1].dimen );

             printf("================================================\n");
             printf("Node#        Fx (%s)        Fy (%s)      Mz (%s)\n",
                     nforces->fn[0].dimen->units_name,
                     nforces->fn[1].dimen->units_name,
                     nforces->fn[2].dimen->units_name);
             printf("================================================\n");

             for(i = 1;i <= frame->no_node_loads; i++) {
                 nforces = &frame->nforces[i-1];
                 for( j=1 ; j<=frame->no_dof ; j++ )
                      UnitsSimplify( nforces->fn[j-1].dimen );
                 printf("%5d ",(int) nforces->node_f); 
                 printf(" %12.2f ",   nforces->fn[0].value/nforces->fn[0].dimen->scale_factor);
                 printf(" %12.2f ",   nforces->fn[1].value/nforces->fn[1].dimen->scale_factor);
                 printf(" %12.2f\n",   nforces->fn[2].value/nforces->fn[2].dimen->scale_factor);
             }
          } else {
             printf("=================================== \n");
             printf("Node#        Fx         Fy       Mz \n");
             printf("=================================== \n");

             for(i = 1;i <= frame->no_node_loads; i++) {
                 nforces = &frame->nforces[i-1];
                 printf("%5d ",(int) nforces->node_f); 
                 printf(" %12.2f ",   nforces->fn[0].value);
                 printf(" %12.2f ",   nforces->fn[1].value);
                 printf(" %12.2f\n ",   nforces->fn[2].value);
             }
          }
      }

      if(frame->no_dof == 2) {
          if(UNITS_SWITCH == ON ) {
             for( i=1 ; i<=frame->no_dof ; i++ )
                 UnitsSimplify( nforces->fn[i-1].dimen );

             printf("=====================================\n");
             printf("Node#        Fx (%s)        Fy (%s)\n",
                     nforces->fn[0].dimen->units_name,
                     nforces->fn[1].dimen->units_name);
             printf("=====================================\n");

             for(i = 1;i <= frame->no_node_loads; i++) {
                 nforces = &frame->nforces[i-1];
                 for( j=1 ; j<=frame->no_dof ; j++ )
                      UnitsSimplify( nforces->fn[j-1].dimen );
                 printf("%5d ",(int) nforces->node_f);
                 printf(" %12.2f ",   nforces->fn[0].value/nforces->fn[0].dimen->scale_factor);
                 printf(" %12.2f\n",   nforces->fn[1].value/nforces->fn[1].dimen->scale_factor);
             }
          } else {
             printf("=============================\n");
             printf("Node#        Fx            Fy  \n");
             printf("=============================\n");
             for(i = 1;i <= frame->no_node_loads; i++) {
                 nforces = &frame->nforces[i-1];
                 printf("%5d ",(int) nforces->node_f); 
                 printf(" %12.2f ",   nforces->fn[0].value);
                 printf(" %12.2f\n",   nforces->fn[1].value);
             }
         }
      }
   }

   if(frame->no_node_loads >= 1 && (frame->no_dimen == 3)) {

      printf("                        \n");
      printf("EXTERNAL NODAL LOADINGS \n");

      nforces = &frame->nforces[0];
      if(frame->no_dof == 6) {
         if( UNITS_SWITCH == ON ) {
            for( i=1 ; i<=frame->no_dof ; i++ )
                 UnitsSimplify( nforces->fn[i-1].dimen );

            printf("============================================================================================\n");
            printf("Node#       Fx (%s)      Fy (%s)     Fz (%s)     Mx (%s)    My (%s)     Mz (%s)\n",
                   nforces->fn[0].dimen->units_name,
                   nforces->fn[1].dimen->units_name,
                   nforces->fn[2].dimen->units_name,
                   nforces->fn[3].dimen->units_name,
                   nforces->fn[4].dimen->units_name,
                   nforces->fn[5].dimen->units_name);
            printf("============================================================================================\n");
         } else {
            printf("========================================== \n");
            printf("Node#       Fx      Fy     Fz     Mx    My \n");
            printf("========================================== \n");
         }
      }

      if(frame->no_dof == 5) {
         if( UNITS_SWITCH == ON ) {
            for( i=1 ; i<=frame->no_dof ; i++ )
                 UnitsSimplify( nforces->fn[i-1].dimen );

            printf("===================================================================================\n");
            printf("Node#       Fx (%s)     Fy (%s)    Fz (%s)    Mx (%s)   My (%s)\n",
                    nforces->fn[0].dimen->units_name,
                    nforces->fn[1].dimen->units_name,
                    nforces->fn[2].dimen->units_name,
                    nforces->fn[3].dimen->units_name,
                    nforces->fn[4].dimen->units_name);
            printf("===================================================================================\n");
         } else {
            printf("========================================== \n");
            printf("Node#       Fx      Fy     Fz     Mx    My \n");
            printf("========================================== \n");
         }
      }

      if( UNITS_SWITCH == ON ) {
         for(i = 1;i <= frame->no_node_loads; i++) {
             nforces = &frame->nforces[i-1];

             for( j=1 ; j<=frame->no_dof ; j++ )
                  UnitsSimplify( nforces->fn[j-1].dimen );

             printf("%5d ",(int) nforces->node_f); 
             printf("%12.3f ",   nforces->fn[0].value/nforces->fn[0].dimen->scale_factor);
             printf("%12.3f ",   nforces->fn[1].value/nforces->fn[1].dimen->scale_factor);
             printf("%12.3f ",   nforces->fn[2].value/nforces->fn[2].dimen->scale_factor);
             printf("%12.3f ",   nforces->fn[3].value/nforces->fn[3].dimen->scale_factor);
             printf("%12.3f ",   nforces->fn[4].value/nforces->fn[4].dimen->scale_factor);
             if(frame->no_dof == 6)
                printf("%12.3f",  nforces->fn[5].value/nforces->fn[5].dimen->scale_factor); 
             printf("\n"); 
         }
      }

      if( UNITS_SWITCH == OFF ) {
         for(i = 1;i <= frame->no_node_loads; i++) {
             nforces = &frame->nforces[i-1];
             printf("%5d ",(int) nforces->node_f); 
             printf("%12.3f ",   nforces->fn[0].value);
             printf("%12.3f ",   nforces->fn[1].value);
             printf("%12.3f ",   nforces->fn[2].value);
             printf("%12.3f ",   nforces->fn[3].value);
             printf("%12.3f ",   nforces->fn[4].value);

             if(frame->no_dof == 6)
                printf("%12.3f",  nforces->fn[5].value);
             printf("\n"); 
         }
      }
   }

   /* [e] : Print Element Forces */

   if(frame->no_element_loads >= 1) {
       printf("\n");
       printf("================================================\n");
       printf("Element Loads : a,b = (in), Px,Py,Pz = (kips/in)\n");
       printf("================================================\n");

       printf("\n");
       printf("Elmt#\n");
       printf("=====\n\n");

       for(i=1; i <= frame->no_element_loads; i++) {
           elsp = &frame->eforces[i-1];
           elmt = &frame->element[elsp->elmt_no -1];
           elp    = elsp->elib_ptr;

           if(iel_no == 1) {
              for(j=1;j<=elsp->no_loads_faces;j++) {
                  printf("%5d %4d %5d %12.5f %12.5f %5d %13.5f %13.5f\n", elsp->elmt_no,
                         elp->face_no,elp->nopl[1],elp->pr->uMatrix.daa[1][1],   
                         elp->pr->uMatrix.daa[1][2],elp->nopl[2],
                         elp->pr->uMatrix.daa[2][1],elp->pr->uMatrix.daa[2][2]);
              }
           }

           if(iel_no == 5) {
           if(elp->type == -1) { /* dist loading */
               printf("%5d :  a = %13.5f  b = %13.5f\n",
                      elsp->elmt_no, elp->a.value, elp->b.value);
               printf("      : Px = %13.5f Py = %13.5f Pz = %13.5f\n",
                      elp->bx.value,elp->by.value,elp->bz.value);
           } else { /* point load  & moments*/
               printf("%5d :  a = %13.5f  b = %13.5f\n",elsp->elmt_no,
                      elp->a.value, elp->b.value);
               printf("      : Px = %13.5f Py = %13.5f Pz = %13.5f\n",
                      elp->bx.value,elp->by.value,elp->bz.value);
               printf("      : Mx = %13.5f My = %13.5f Mz = %13.5f\n",
                      elp->mx.value,elp->my.value,elp->mz.value);
           }
           }

           if(iel_no == 7 || iel_no == 2) {
           if(elp->type == -1) { /* distributed loading */
              printf("%5d : a = %13.5f b = %13.5f Py = %13.5f\n",
                      elsp->elmt_no, elp->a.value, elp->b.value, elp->by.value);
           } else
              printf("%5d : a = %13.5f b = %13.5f P  = %13.5f\n",
                      elsp->elmt_no, elp->a.value, elp->b.value, elp->P.value);
           }
       } 
    }

    /* [f] : End Message */

    printf("\n\n");
    printf("============= End of Finite Element Mesh Description ==============\n");
    printf("\n\n");
}


/* 
 *  =========================================================
 *  Print_Stress() : Print Stresses.
 *
 *  Input  : Matrix * : Matrix of Stresses.....
 *  Output : void
 *  =========================================================
 */ 

#ifdef  __STDC__
MATRIX *Print_Stress( MATRIX *m1, ... )
#else
MATRIX *Print_Stress(va_alist)
va_dcl
#endif
{
#ifdef __STDC__
va_list        arg_ptr;
#else
va_list        arg_ptr;
MATRIX             *m1;
#endif

ELEMENT            *ep;
ELEMENT_ATTR      *eap;
double       *elvector;
double       mzc,   mz; 
int            node_no;
int            elmt_no;
int       elmt_attr_no;
int    i, j, k, jj, ii;
int             length;
int       UNITS_SWITCH;

DIMENSIONS      *dp_length_temp;
DIMENSIONS       *dp_force_temp;
double    xi[17],eta[17],wg[17];
double              *gamma, *wt;
int                     l, lint;

SYMBOL                     *slp;
int            iInPlaneIntegPts;
int          iThicknessIntegPts;
int               iNO_INTEG_pts;

#ifdef __STDC__
       va_start(arg_ptr, m1);
#else
       va_start(arg_ptr);
       m1 = va_arg(arg_ptr, MATRIX *)
#endif
       va_end(arg_ptr);

    UNITS_SWITCH = CheckUnits();

    /* [a] : Print Header */

    printf("\n\n");
    printf("MEMBER FORCES \n");
    printf("------------------------------------------------------------------------------\n");

    /* [b] :Calculate and Print Element Stresses (linear and nonlinear cases) */

    /* Stresses corresponding to linear response */

    if(m1 != (MATRIX *) NULL)
    for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {

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
                switch((int) array->nodes_per_elmt) {
                  case 2:
                  case 3:
                       ii = eap->map_ldof_to_gdof[j-1];
                       jj = frame->node[node_no - 1].bound_id[ii-1];
                       if(jj > 0) {
                          array->displ->uMatrix.daa[j-1][i-1] = m1->uMatrix.daa[jj-1][0];
                          if( UNITS_SWITCH == ON ) {
                              UnitsCopy(&(array->displ->spRowUnits[j-1]), &(m1->spRowUnits[jj-1]));
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
                          array->displ->uMatrix.daa[k-1][i-1] = m1->uMatrix.daa[jj-1][0];
                          if( UNITS_SWITCH == ON ) {
                              UnitsCopy(&(array->displ->spRowUnits[k-1]), &(m1->spRowUnits[jj-1]));
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

        array = elmlib(array, STRESS);

        for( i=1 ; i <= frame->no_dof ; i++ )
        for( j=1 ; j <= frame->no_nodes_per_elmt ; j++ )
            frame->element[elmt_no-1].rp->Forces->uMatrix.daa[i-1][j-1]
            = array->nodal_loads[frame->no_dof*(j-1)+i-1].value;

        if( UNITS_SWITCH == ON ) {
            for( i=1 ; i <= frame->no_dof ; i++ )
                 UnitsCopy( &(frame->element[elmt_no-1].rp->Forces->spRowUnits[i-1]),
                              array->nodal_loads[i-1].dimen );
            for( j=1 ; j <= frame->no_nodes_per_elmt ; j++ )
                 ZeroUnits( &(frame->element[elmt_no-1].rp->Forces->spColUnits[j-1]) );
        }
    }

    /* 
     *  ====================================================================
     *  Print stresses corresponding to "nonlinear" response of structure
     *  These stresses are stored internally in Aladdin database       
     *                                                                 
     *  Aladdin command is : PrintStress() with no matrix argument.
     *                                                                 
     *  Note : Code in this section still needs to be checked          
     *  ====================================================================
     */ 

    if(m1 == (MATRIX *) NULL) {

       slp = lookup("InPlaneIntegPts");  /* no of integration pts in plane/surface */
       if(slp == NULL)
          iInPlaneIntegPts = 2*2;        /* 2x2 as default */
       else
          iInPlaneIntegPts = (int) slp->u.q->value;

       slp = lookup("ThicknessIntegPts"); /* no of integration pts in thickness direction */
       if(slp == NULL)
          iThicknessIntegPts = 2;        /* 2 as default */
       else
          iThicknessIntegPts = (int) slp->u.q->value;
  
       iNO_INTEG_pts = iInPlaneIntegPts*iThicknessIntegPts;

       lint = (int) 0;
       l = (int) iInPlaneIntegPts;
       for(i = 1; i<= l*l; i++) {
           xi[i-1]  = 0.0;
           eta[i-1] = 0.0;
           wg[i-1]  = 0.0;
       }

       if(l*l != lint)
          pgauss(l,&lint,xi,eta,wg);

       gamma = dVectorAlloc(iThicknessIntegPts+1);
       wt    = dVectorAlloc(iThicknessIntegPts+1);
       gauss(gamma,wt,iThicknessIntegPts);

       for(elmt_no = 1; elmt_no <= frame->no_elements; elmt_no++) {
           ep           = &frame->element[elmt_no-1];
           elmt_attr_no = ep->elmt_attr_no;
           eap          = &frame->eattr[elmt_attr_no-1];
           if(elmt_no == 1) 
              printf(" Element : %s \n Material : %s \n\n", eap->elmt_type, eap->material);

           printf("\n STRESS in  Element No  %d \n",elmt_no);
           printf(" =============================================================================================================== \n");
           printf(" Gaussion    xi       eta        gamma    Stre-xx         Stre-yy         Stre-xy         Stre-yz        Stre-zx \n");
           if(UNITS_SWITCH == OFF)
              printf("  Points \n");
           if(UNITS_SWITCH == ON)   {
              printf("  Points                                    %s             %s             %s              %s             %s \n",
                          ep->rp->stress->spRowUnits[0].units_name,
                          ep->rp->stress->spRowUnits[1].units_name,
                          ep->rp->stress->spRowUnits[2].units_name,
                          ep->rp->stress->spRowUnits[3].units_name,
                          ep->rp->stress->spRowUnits[4].units_name);
              printf(" ---------------------------------------------------------------------------------------------------------------\n \n");

              for(i = 1; i <= iThicknessIntegPts; i++)
              for(j = 1; j <= lint; j++) {
                  k = lint*(i-1)+j;
                  printf("   %d  %10.4f %10.4f %10.4f", k, xi[j-1], eta[j-1], gamma[i]);
                  printf("\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n", 
                             ep->rp->stress->uMatrix.daa[0][k-1]/ep->rp->stress->spRowUnits[0].scale_factor,
                             ep->rp->stress->uMatrix.daa[1][k-1]/ep->rp->stress->spRowUnits[1].scale_factor,
                             ep->rp->stress->uMatrix.daa[2][k-1]/ep->rp->stress->spRowUnits[2].scale_factor,
                             ep->rp->stress->uMatrix.daa[3][k-1]/ep->rp->stress->spRowUnits[3].scale_factor,
                             ep->rp->stress->uMatrix.daa[4][k-1]/ep->rp->stress->spRowUnits[4].scale_factor);
              }
              }
              if(UNITS_SWITCH == OFF)
                for(i = 1; i <= iThicknessIntegPts; i++)
                    for(j = 1; j <= lint; j++) {
                        k = lint*(i-1)+j;
                        printf("   %d  %10.4f %10.4f %10.4f", k, xi[j-1], eta[j-1], gamma[i]);
                        printf("\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",
                             ep->rp->stress->uMatrix.daa[0][k-1],
                             ep->rp->stress->uMatrix.daa[1][k-1],
                             ep->rp->stress->uMatrix.daa[2][k-1],
                             ep->rp->stress->uMatrix.daa[3][k-1],
                             ep->rp->stress->uMatrix.daa[4][k-1]);
                }
        }
    }
}

/* 
 *  =========================================================
 *  Print_Displ() : Print Displacements
 *
 *  Input  : Matrix * : Matrix of Displacements to be printed
 *  Output : void
 *  =========================================================
 */ 

#ifdef __STDC__
void Print_Displ(MATRIX *m1)
#else
void Print_Displ(m1)
MATRIX *m1;
#endif
{
int node_no, j, k, ii, jj;
int          UNITS_SWITCH;
double                 da;

    UNITS_SWITCH = CheckUnits();

    if(frame->no_dof == 6) {
    printf("\n");
    printf("==================================================================================================================\n");
    printf(" Node                                             Displacement\n");
    printf("   No            displ-x           displ-y           displ-z             rot-x             rot-y             rot-z\n");
    printf("==================================================================================================================\n");

    if( UNITS_SWITCH == ON ) {
        printf(" units");
        for( ii=1 ; ii<=6 ; ii++ ) {
           if( ii <= 3 ) {
               if((CheckUnitsType()==SI) || (CheckUnitsType()==SI_US) )
                   printf("                m ");
               else 
                   printf("               in ");
           } else
                   printf("              rad ");
        }
    }
    printf("\n");
    }

    if(frame->no_dof == 5) {
    printf("\n");
    printf("================================================================================================\n");
    printf(" Node                                   Displacement\n");
    printf("   No            displ-x           displ-y           displ-z             rot-x             rot-y\n");
    printf("================================================================================================\n");

    if( UNITS_SWITCH == ON ) {
          printf(" units");
          for( ii=1 ; ii<=5 ; ii++ ) {
             if( ii <= 3 ) {
                if( (CheckUnitsType()==SI) || (CheckUnitsType()==SI_US) )
                   printf("                m ");
                else 
                   printf("               in ");
             }
             else
                   printf("              rad ");
          }
    }
    printf("\n");
    }

    if(frame->no_dof == 3) {
    printf("\n");
    printf("============================================================\n");
    printf(" Node                           Displacement                \n");
    printf("   No            displ-x           displ-y             rot-z\n");
    printf("============================================================\n");

    if( UNITS_SWITCH == ON ) {
        printf(" units");
        for( ii=1 ; ii<=3 ; ii++ ) {
             if( ii <= 2 ) {
                if( (CheckUnitsType()==SI) || (CheckUnitsType()==SI_US) )
                   printf("                m ");
                else 
                   printf("               in ");
             }
             else
                   printf("              rad ");
          }
    }
    printf("\n");
    }

    if(frame->no_dof == 2) {
    printf("\n");
    printf("==========================================\n");
    printf(" Node                  Displacement       \n");
    printf("   No            displ-x           displ-y\n");
    printf("==========================================\n");

    if( UNITS_SWITCH == ON ) {
        printf(" units");
        for( ii=1; ii <= 2; ii++ ) {
             if( (CheckUnitsType()==SI) || (CheckUnitsType()==SI_US) )
                 printf("                m ");
             else 
                 printf("               in ");
        }
    }
    printf("\n");
    } 

 /*
  *  ==================================================
  *  Print displacements for UNITS_SWITCH == ON and OFF
  *  ==================================================
  */

  if( UNITS_SWITCH == ON ) {

      for(node_no = 1; node_no <= frame->no_nodes; node_no++) {

      printf("%4d  ", node_no);
      for(j  = 1; j <= frame->no_dof; j++){

          jj = frame->node[node_no - 1].bound_id[j-1];
          if(jj > 0) {
             if(m1->spRowUnits[jj-1].units_name != NULL ) {
                UnitsTypeConvert(&(m1->spRowUnits[jj-1]), CheckUnitsType());
             }
             RadUnitsSimplify( &(m1->spRowUnits[jj-1]) );
             frame->node[node_no-1].disp[j-1].value = m1->uMatrix.daa[jj-1][0];
             UnitsCopy( frame->node[node_no-1].disp[j-1].dimen, &(m1->spRowUnits[jj-1]) );
             da =  MatrixContentScaleIndirectDouble(m1,jj, 1);
             printf("     %13.5e", da);
          } else {

             switch(CheckUnitsType()) {
                case SI:
                     if( (frame->node[node_no-1].disp[j-1].dimen->units_name != NULL) &&
                        !strcmp( frame->node[node_no-1].disp[j-1].dimen->units_name, "deg_F"))

                     frame->node[node_no -1].disp[j-1].value
                     = ConvertTempUnits(frame->node[node_no-1].disp[j-1].dimen->units_name,
                                        frame->node[node_no -1].disp[j-1].value, US);
                     break;
                case US:
                     if((frame->node[node_no-1].disp[j-1].dimen->units_name != NULL) &&
                         !strcmp(frame->node[node_no-1].disp[j-1].dimen->units_name, "deg_C"))

                         frame->node[node_no -1].disp[j-1].value
                         = ConvertTempUnits(frame->node[node_no-1].disp[j-1].dimen->units_name,
                                            frame->node[node_no -1].disp[j-1].value, SI);
                     break;
             }

             RadUnitsSimplify( frame->node[node_no-1].disp[j-1].dimen );
             printf("     %13.5e", frame->node[node_no -1].disp[j-1].value/
                          frame->node[node_no-1].disp[j-1].dimen->scale_factor);
         }
      }
      printf("\n");
      }
  }

  /* Print displacements when units are OFF */

  if( UNITS_SWITCH == OFF ) {
   
      for(node_no = 1; node_no <= frame->no_nodes; node_no++) {
          printf("%4d  ", node_no);
          for(j  = 1; j <= frame->no_dof; j++){
              jj = frame->node[node_no - 1].bound_id[j-1];
              if(jj > 0) {
                 frame->node[node_no-1].disp[j-1].value = m1->uMatrix.daa[jj-1][0];
                 printf("     %13.5e", m1->uMatrix.daa[jj-1][0]);
              }
              else
                 printf("     %13.5e", frame->node[node_no -1].disp[j-1].value);
          }
          printf("\n");
      }
  }

}
