#FC = ifort
FC = ifx
#FC = gfortran
#MODE = debug
MODE = release
SRC = src
BUILD = build

ifneq ($(filter gfortran%,$(FC)),)
  TC = gnu
  OFLAGS = -O2
  DFLAGS = -O0 -g -fcheck=all
  GFLAGS = -cpp -fopenmp -J$(BUILD) -Wunused-variable \
			  -fbacktrace -ffpe-summary=none -fno-stack-arrays \
			  -no-pie
else ifneq ($(filter ifort ifx,$(FC)),)
  TC = intel
  ifeq ($(FC),ifort)
    CHFLAG = -check all
    DIAGFLAG = -diag-disable=10448
  else ifeq ($(FC),ifx)
    CHFLAG = -check all,nouninit
  endif
  OFLAGS = -O3 -fp-model strict -fPIC
  DFLAGS = -O0 $(CHFLAG) -debug -warn all
  GFLAGS = -traceback -cpp -qopenmp -module $(BUILD) \
			  $(DIAGFLAG)
else
  $(error Unsupported Fortran compiler: FC=$(FC))
endif
LDFLAGS = -Wl,-z,execstack,--no-warn-execstack

LAPACK_PATH = ~/.local/lib/lapack-3.12.0/lib/$(TC)
LAPACK = $(LAPACK_PATH)/liblapack.a
BLAS = $(LAPACK_PATH)/librefblas.a

IFLAGS = -I$(SRC)/util

ifeq ($(MODE),release)
  FFLAGS = $(strip $(OFLAGS) $(GFLAGS) $(IFLAGS))
else ifeq ($(MODE), debug)
  FFLAGS = $(strip $(DFLAGS) $(GFLAGS) $(IFLAGS))
else
  $(error Unsupported MODE: MODE=$(MODE))
endif

all: crystal.x

include mk/obgen.mk
include mk/util.mk
include mk/depsys.mk
include mk/generic.mk

OBJS = $(patsubst %.f90, $(BUILD)/%.o, \
		 $(notdir $(SFILES)))

$(OBJS): | $(BUILD)

$(BUILD):
	@mkdir -p $(BUILD)

crystal.x: $(OBJS)
	@echo Creating crystal.x ...
	@$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LAPACK) $(BLAS) -o $@

.PHONY: all clean

clean:
	rm -f $(BUILD)/*.* crystal.x
