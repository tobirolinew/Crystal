define obgen
	$(if $(strip $(2)),,$(error obgen: empty DEP for $(1)))
   $(eval SFILE0 = $(1))
   $(eval SFILE_BASE = $(basename $(SFILE0)))
   $(eval OBJ = $(BUILD)/$(SFILE_BASE).o)
   $(eval FIND_SFILE = $(shell find $(SRC) -type f -name $(SFILE0)))
   $(eval SFILE = $(subst ./,,$(FIND_SFILE)))
   $(eval ODEP0 = $(filter %.o,$(2)))
   $(eval ODEP = $(foreach f,$(ODEP0),$(BUILD)/$(f)))
   $(eval SDEP0 = $(filter-out %.o,$(2)))
   $(foreach f,$(SDEP0), \
      $(eval CURR_SRC = $(f)) \
      $(eval CURR_FOUND = \
         $(shell find $(SRC) -type f -name $(CURR_SRC))) \
      $(eval $(CURR_SRC)_FOUND = $(CURR_FOUND)) \
   )
   $(eval FIND_SDEP = $(foreach f,$(SDEP0),$($(f)_FOUND)))
   $(eval SDEP = $(subst ./,,$(FIND_SDEP)))
   $(eval DEP = $(ODEP) $(SDEP))
   $(OBJ): $(SFILE) $(DEP)
		$(3)
   $(eval SFILES := $(SFILES) $(SFILE))
endef
