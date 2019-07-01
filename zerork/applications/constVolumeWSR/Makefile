ZERORK_ROOT:=$(realpath ../../)
include ../../Makefile.defs
include ../../../build.config

TARGETS  := constVolumeWSR.x direct_CVDense
OBJS     := cvIDT_sweep.o matrix_funcs.o ode_funcs.o \
            cv_param_sparse.o sweep_util.o utility_funcs.o

CXXINCL:= -I$(SUPERLU_ROOT)/SRC -I$(SUNDIALS_ROOT)/include \
          -I../../include

CXXLIBS:= -L$(SUPERLU_ROOT)/lib -L$(SUNDIALS_ROOT)/lib \
          -L../../lib \
          -L$(LAPACK_LIB_DIR) $(LAPACK_LIBS) \
          -L$(BLAS_LIB_DIR) $(BLAS_LIBS)

ZERORKLIBS:=$(ZERORK_ROOT)/lib/libzerork.a $(ZERORK_ROOT)/lib/libckconverter.a

ifeq ($(strip $(SUPERLU_VER)),4_3_static)
  CXXLIBS += -Wl,-Bstatic -lsundials_nvecserial -lsundials_cvode \
             -l$(SUPERLU_LIBNAME) -lzerork -lckconverter \
             -Wl,-Bdynamic
else
  CXXLIBS += -Wl,-Bstatic -lsundials_nvecserial -lsundials_cvode \
             -l$(SUPERLU_LIBNAME) -lzerork -lckconverter \
             -Wl,-Bdynamic
endif

ifeq ($(USE_MKL),Y)
  CXXLIBS+= -L$(MKLROOT)/lib -Wl,-Bstatic -Wl,--start-group \
            -lmkl_core -lmkl_sequential \
            -lmkl_intel_lp64 -lmkl_blas95_lp64 \
            -lmkl_lapack95_lp64 -Wl,--end-group -Wl,-Bdynamic
else
  CXXLIBS += -L$(LAPACK_LIB_DIR) $(LAPACK_LIBS) \
             -L$(BLAS_LIB_DIR) $(BLAS_LIBS)
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

constVolumeWSR.x: $(OBJS) $(ZERORKLIBS)
	$(CLINK) $(CXXFLAGS) $(CXXINCL) $(OBJS) $(CXXLIBS) $(END_LIBS) \
          -o $@
.PHONY: direct_CVDense
direct_CVDense:
	cd ./direct_CVDense ; make all


clean:
	rm -f *.o *~ *.x
	cd ./direct_CVDense ; make clean
