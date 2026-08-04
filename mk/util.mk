SFILE = mod_catch.f90
DEP = NONE
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_stack_unlim.f90
DEP = NONE
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_swap.f90
DEP = NONE 
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_rand_perm.f90
DEP = mod_swap.o 
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_string.f90
DEP = mod_catch.o 
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_sort.f90
DEP = mod_catch.o 
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_graph.f90
DEP = NONE 
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_pertab.f90
DEP = NONE
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_xinfo.f90
DEP = NONE
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_euclid.f90
DEP = mod_pertab.o mod_catch.o mod_string.o \
		read_xyz_inc.f90 lapack_iface_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_point_group.f90
DEP = mod_pertab.o mod_euclid.o mod_string.o
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_rbalign.f90
DEP = mod_euclid.o mod_sort.o mod_graph.o \
		matchtool_iface_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_block_input.f90
DEP = mod_catch.o \
		mod_string.o \
		generic_append_inc.f90 \
		append_block_input_inc.f90 \
		gather_block_input_inc.f90 \
		popul_block_input_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_quasi_newton.f90
DEP = mod_catch.o  func_iface_inc.f90 \
		lapack_iface_inc.f90 \
		grad_iface_inc.f90 \
		diag_hess_iface_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

