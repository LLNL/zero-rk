ZERORK_ROOT:=$(realpath ../../)
include ../../Makefile.defs

TARGETS = functionTester.x

CXXINCL  := -I../../include
CXXLIBS  := -L../../lib
CXXLIBS  +=  -lzerork -lckconverter

ifeq ($(USE_MKL),Y)
  CXXLIBS+= -L$(MKLROOT)/lib -Wl,-Bstatic -Wl,--start-group \
            -lmkl_core -lmkl_sequential \
            -lmkl_intel_lp64 -lmkl_blas95_lp64 \
            -lmkl_lapack95_lp64 -Wl,--end-group -Wl,-Bdynamic
endif

ifeq ($(USE_AMD_LIBM),Y)
  CXXLIBS+= -L$(AMDLIBMROOT)/lib/static -lamdlibm
endif

CLINK    := $(CXX)
END_LIBS := -lm -ldl

all: $(TARGETS)

functionTester.x: functionTester.o functionTester.h
	$(CLINK) $(CXXFLAGS) $(CXXINCL) $< $(CXXLIBS) $(END_LIBS) \
          -o $@

clean:
	rm -f *.o *~ *.x
