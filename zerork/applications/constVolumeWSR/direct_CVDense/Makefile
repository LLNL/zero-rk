ZERORK_ROOT:=$(realpath ../../../)
include ../../../Makefile.defs

TARGETS  := cv_wsr_direct_CVDense.x
OBJS     := cvIDT_sweep.o ode_funcs.o \
            sweep_util.o utility_funcs.o 

CXXINCL:= -I$(SUPERLU_ROOT)/SRC -I$(SUNDIALS_ROOT)/include \
          -I$(ZERORK_ROOT)/include

CXXLIBS:= -L$(SUPERLU_ROOT)/lib -L$(SUNDIALS_ROOT)/lib \
          -L$(ZERORK_ROOT)/lib


CXXLIBS += -lderivative -lzerork -lckconverter \
           -l$(SUPERLU_LIBNAME) -lsundials_nvecserial -lsundials_cvode


ifeq ($(USE_MKL),Y)
  CXXLIBS+= -L$(MKLROOT)/lib -Wl,-Bstatic -Wl,--start-group \
            -lmkl_core -lmkl_sequential \
            -lmkl_intel_lp64 -lmkl_blas95_lp64 \
            -lmkl_lapack95_lp64 -Wl,--end-group -Wl,-Bdynamic
endif

ifeq ($(USE_AMD_LIBM),Y)
  CXXLIBS+= -L$(AMDLIBMROOT)/lib/static -lamdlibm
endif


ifneq ($(BIGRED2_CRAY),Y)
  CXXLIBS += -llapack -lblas
endif

CLINK    := $(CXX)
END_LIBS := -lm -ldl

all: $(TARGETS)

cv_wsr_direct_CVDense.x: $(OBJS)
	$(CLINK) $(CXXFLAGS) $(CXXINCL) $(OBJS) $(CXXLIBS) $(END_LIBS) \
          -o $@

clean:
	rm -f *.o *~ $(TARGETS)
