SFILE = mod_f2b.f90
DEP = mod_euclid.o mod_string.o \
		type_f2b_inc.f90 read_f2b_pars_inc.f90 \
		add_OA_sites_inc.f90 f2b_func_inc.f90 \
		pot_f2b_inc.f90 pot_f2b_homo_inc.f90 \
		grad_f2b_homo_inc.f90 \
		cpot_f2b_inc.f90 cpot_f2b_homo_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_engine_rigidMC.f90
DEP = mod_lattice.o mod_euclid.o mod_point_group.o mod_catch.o \
		mod_swap.o mod_rand_perm.o mod_string.o mod_sort.o mod_graph.o \
		mod_quasi_newton.o mod_f2b.o \
		generic_append_inc.f90 \
      type_control_rigidMC_inc.f90 type_layout_rigidMC_inc.f90 \
      type_refframe_rigidMC_inc.f90 type_sentinel_rigidMC_inc.f90 \
      type_landscape_rigidMC_inc.f90 type_geo_rigidMC_inc.f90 \
      type_monbox_rigidMC_inc.f90 type_coorbox_rigidMC_inc.f90 \
      init_control_rigidMC_inc.f90 init_layout_rigidMC_inc.f90 \
      init_refframe_rigidMC_inc.f90 init_sentinel_rigidMC_inc.f90 \
      init_landscape_rigidMC_inc.f90 init_monbox_rigidMC_inc.f90 \
      init_coorbox_rigidMC_inc.f90 \
		init_geo_rigidMC_inc.f90 manip_geo_rigidMC_inc.f90 \
		reparam_geo_rigidMC_inc.f90 sync_geo_rigidMC_inc.f90 \
      calc_geo_rigidMC_inc.f90 assess_geo_rigidMC_inc.f90 \
      screen_geo_rigidMC_inc.f90 stream_geo_rigidMC_inc.f90 \
      pot_geo_rigidMC_inc.f90 grad_geo_rigidMC_inc.f90 \
      hess_geo_rigidMC_inc.f90 disp_geo_rigidMC_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

SFILE = mod_user_rigidMC.f90
DEP = mod_global.o mod_rbalign.o update_wcoors_inc.f90 \
		update_lcoors_inc.f90 exchange_wcoors_inc.f90 \
		resolve_geo_uptr_inc.f90 mod_engine_rigidMC.o \
		generic_append_inc.f90 globalize_rigidMC_inc.f90 \
		load_rigidMC_inc.f90 customize_rigidMC_inc.f90 \
		cost_rigidMC_inc.f90 task_rigidMC_inc.f90
CMD = $(FC) $(FFLAGS) -c $$< -o $$@
$(eval $(call obgen, $(SFILE), $(DEP), $(CMD)))

