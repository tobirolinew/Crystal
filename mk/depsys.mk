SFILE = mod_sys.f90
DEP = mod_lattice.o \
		mod_user_rigidMC.o
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

include mk/sys/rigidMC.mk
