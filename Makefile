# ==========================================================
# Makefile for ALADDIN
#
# To compile program, type: make aladdin
#
# Written By: M. Austin                       September 1995
# ==========================================================

CFLAGS = -g
YFLAGS = -d
CC = gcc
# CC = /usr/local/SUNWspro/SC2.0.1/cc

.cc.o:	        
		c++ -c $(CFLAGS) -O  $<
.c.o:
		$(CC) $(CFLAGS) -c $<

# Macros

MISC =		miscellaneous.o

MATRIX =	matrix.o 		matrix_double.o	\
		matrix_indirect.o	matrix_skyline.o \
		vector.o		vector_double.o	\
		vector_integer.o

KERNEL =  	grammar.o 	code.o 		init.o		main.o \
		symbol.o	symbol_load.o	units.o		math.o		

FINITE =	fe_allocate.o 	fe_mesh.o       fe_matrix.o \
        	fe_nonlinear.o	fe_profile.o	fe_setflags.o \
        	fe_print.o	fe_checkmemory.o

ELEMENTS =  	elmt_set_attr.o	elmt_library.o	elmt_lamina_sys.o \
        	elmt_frame2d.o	elmt_frame3d.o	elmt_psps.o 	elmt_plate.o \
		elmt_shell_4n.o	elmt_shell_8n.o elmt_shell_4n_q.o \
        	elmt_fiber2d.o	elmt_fiber3d.o

# Make ALADDIN main program

ALADDIN: $(MATRIX) $(MISC) $(KERNEL) $(FINITE) $(ELEMENTS)
	$(CC) $(MATRIX) $(MISC) $(KERNEL) $(FINITE) $(ELEMENTS) -o ALADDIN -lm

# ALADDIN header file dependencies

$(MISC)   $(KERNEL) $(FINITE) $(MATRIX) $(ELEMENTS) : miscellaneous.h
$(KERNEL) $(FINITE) $(MATRIX) $(ELEMENTS) : defs.h
$(KERNEL) $(FINITE) $(MATRIX) $(ELEMENTS) : matrix.h
$(KERNEL) $(FINITE) $(MATRIX) $(ELEMENTS) : units.h
$(KERNEL) $(FINITE) $(ELEMENTS)           : symbol.h
code.o fe_allocate.o fe_checkmemory.o init.o : x.tab.h
grammar.o  code.o  fe_allocate.o  fe_checkmemory.o  init.o  main.o:  code.h
$(ELEMENTS)  fe_matrix.o  fe_mesh.o  fe_print.o  fe_profile.o  init.o:  elmt.h
grammar.o code.o init.o main.o symbol.o symbol_load.o units.o:  fe_database.h
$(ELEMENTS) $(FINITE) : fe_database.h
$(ELEMENTS) fe_matrix.o fe_mesh.o fe_print.o fe_setflag.o code.o init.o main.o: fe_functions.h
$(MATRIX) elmt_frame3d.o elmt_psps.o elmt_shell_4n.o elmt_shell_4n_q.o elmt_shell_8n.o \
fe_checkmemory.o fe_matrix.o fe_mesh.o fe_print.o fe_profile.o : vector.h

x.tab.h:   y.tab.h
	cmp -s x.tab.h y.tab.h || cp y.tab.h x.tab.h

# Remove all object files.

clean:
	/bin/rm -f *.o [xy].tab.[ch]
