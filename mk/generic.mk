SFILE = mod_global.f90
DEP = mod_catch.o 
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_lattice.f90
DEP = mod_global.o \
		mod_lpo_inc.f90 \
		prior_queue_inc.f90 \
		ffhash_inc.f90 \
		mod_lpo_lock_set_inc.f90 \
		mod_lpo_cost_set_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_proj_plot.f90
DEP = mod_lattice.o \
	   mod_xinfo.o \
		generic_append_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_crystal.f90
DEP = mod_block_input.o \
		mod_proj_plot.o \
		mod_sys.o
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = main.f90
DEP = mod_crystal.o mod_stack_unlim.o
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

