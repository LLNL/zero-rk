#Makefile for zerork_cvIDTsweep

#USER Specifiable Options and External Package Directories

include ./Makefile.defs
include ../build.config
##
##


EXE_EXT :=.x
LIB_NAME:=zerork_reactor

#TODO: Move to defs
LIBTYPE:=shared

PROG_LIB_NAME:=zerork_reactor

# the object files to be linked together. 
COMMON_OBJS:= matrix_funcs.o ode_funcs.o \
       cv_param_sparse.o utility_funcs.o \
       multizone.o


PROG_LIB_DRIVER:=zerork_reactor.o
PROG_DRIVER:=zerork_cvIDT_multireactor.o
PROG_DTEST_DRIVER:=zerork_cvIDT_multireactor_dtest.o
IFP_OBJ:=zerork_cvIDT_multireactorIFP.o
IFP_LIB_OBJ:=zerork_reactorIFP.o
LIB_DRIVER:=zerork_reactor_lib.o
LIB_DRIVER_SRC:=zerork_reactor_lib.cpp

PROG_LIB_OBJS:=$(IFP_LIB_OBJ) $(PROG_LIB_DRIVER)
PROG_OBJS:=$(COMMON_OBJS) $(IFP_OBJ) $(PROG_DRIVER)
PROG_DTEST_OBJS:=$(COMMON_OBJS) $(IFP_OBJ) $(PROG_DTEST_DRIVER)
LIB_OBJS:=$(COMMON_OBJS) $(LIB_DRIVER)

ALL_OBJS:=$(COMMON_OBJS) $(PROG_DRIVER) $(PROG_DRIVER) $(PROG_LIB_DRIVER) \
          $(LIB_DRIVER) $(IFP_OBJ) $(PROG_DTEST_DRIVER) $(IFP_LIB_OBJ)

COMMON_OBJS:= $(patsubst %.o,src/%.o,$(COMMON_OBJS))
PROG_LIB_OBJS:= $(patsubst %.o,src/%.o,$(PROG_LIB_OBJS))
PROG_OBJS:= $(patsubst %.o,src/%.o,$(PROG_OBJS))
PROG_DTEST_OBJS:= $(patsubst %.o,src/%.o,$(PROG_DTEST_OBJS))
LIB_OBJS:= $(patsubst %.o,src/%.o,$(LIB_OBJS))
ALL_OBJS:= $(patsubst %.o,src/%.o,$(ALL_OBJS))



CXXINCL:= -I$(SUNDIALSROOT)/include \
          -I$(ZERORKROOT)/include \
          -I$(SUPERLUROOT)/SRC \
          -I$(SPIFYROOT)/include

CXXLIBS:= -L$(SUPERLUROOT)/lib -L$(SUNDIALSROOT)/lib \
          -Wl,-Bstatic -lsundials_nvecserial  -lsundials_cvode \
          -lsuperlu_4.3


ifeq ($(USE_ARKODE),Y)
  DEFS+= -DARKODE
  CXXLIBS+= -lsundials_arkode  #Sundials-2.6.0
endif




ifeq ($(LIBTYPE), shared)
ZERORKLIB:=-L$(ZERORKROOT)/lib -lcvode_user -lzerork -lckconverter
else
ZERORKLIB:=$(ZERORKROOT)/lib/libzerork.a
ZERORKLIB+= $(ZERORKROOT)/lib/libckconverter.a
ZERORKLIB+= $(ZERORKROOT)/lib/libcvode_user.a
ZERORKLIB+= $(SUPERLUROOT)/lib/libsuperlu_4.3.a
ZERORKLIB+= $(SUNDIALSROOT)/lib/libsundials_nvecserial.a
ZERORKLIB+= $(SUNDIALSROOT)/lib/libsundials_cvode.a

endif

PARSER_LIBS:= -L$(SPIFYROOT)/lib -lspify -lyaml-cpp


ifeq ($(USE_METIS),Y)
  CXXINCL+= -I$(METISROOT)/include
  CXXLIBS+= -L$(METISROOT)/lib -lmetis
  DEFS+= -DCVIDT_USE_METIS
endif

ifeq ($(USE_MKL),Y)
  CXXLIBS+= -L$(MKLROOT)/lib -Wl,--start-group \
            -lmkl_core -lmkl_sequential \
            -lmkl_intel_lp64 -lmkl_blas95_lp64 \
            -lmkl_lapack95_lp64 -Wl,--end-group
endif

ifeq ($(USE_AMD_LIBM),Y)
  CXXLIBS+= -L$(AMDLIBMROOT)/lib/static -lamdlibm
endif

ifeq ($(USE_MPI),Y)
  DEFS+= -DUSE_MPI
endif

CXXLIBS+= -Wl,-Bdynamic -lpthread -ldl

ifneq ($(USE_MKL),Y)
  CXXLIBS += -L$(LAPACK_LIB_DIR) $(LAPACK_LIBS) \
             -L$(BLAS_LIB_DIR) $(BLAS_LIBS)
endif

ifeq ($(LIBTYPE), shared)
LIBRARY:= lib$(LIB_NAME).so
else
LIBRARY:= lib$(LIB_NAME).a
endif

PROGRAM_LIB:= $(PROG_LIB_NAME)$(EXE_EXT)

DEPENDS:= $(OBJS:.o=.d)


IFP_LIB_SRC:= src/zerork_reactorIFP.cpp src/zerork_reactorIFP.h



# flags passed to the C++ compiler/linker for the linking step
LCXXFLAGS:= $(CXX_FLAGS) $(CXXLIBS) -Wl,-Bsymbolic -Wl,-Bsymbolic-functions


# RULES BELOW


all: $(PROGRAM_LIB) $(LIBRARY)

$(IFP_LIB_SRC) : src/generate_zerork_reactorIFP.py
	cd src; ./generate_zerork_reactorIFP.py;

# how to compile C++ source files to object files
src/%.o : src/%.cpp
	$(CXX) -c $< $(CXXINCL) $(CXX_FLAGS) $(DEFS) -o $@

src/%.o : src/%.c
	$(CC) -c $< $(CXXINCL) $(CXX_FLAGS) $(DEFS) -o $@

$(PROGRAM_LIB): $(PROG_LIB_OBJS) $(LIBRARY)
	$(CXX) -o $(PROGRAM_LIB) $^ $(ZERORKLIB) $(LCXXFLAGS) $(PARSER_LIBS)

$(PROGRAM): $(PROG_OBJS)
	$(CXX) -o $(PROGRAM) $^ $(ZERORKLIB) $(LCXXFLAGS) $(PARSER_LIBS)


ifeq ($(LIBTYPE),shared)
$(LIBRARY): $(LIB_OBJS)
	$(CXX) -o $(LIBRARY) $^ $(ZERORKLIB) $(LCXXFLAGS) -$(LIBTYPE)

else
$(LIBRARY): $(LIB_OBJS)
	@echo Making static lib $(LIBRARY)
	@objdirs="";  objs=""; orig_pwd=$$PWD; \
	for lib in $(ZERORKLIB); do \
	  tmpdir=`/bin/mktemp -d`;  \
	  echo $${tmpdir}; \
	  cd $${tmpdir}; \
	  ar x $${lib}; \
	  cd $$orig_pwd; \
	  objs="$${objs} $${tmpdir}/*.o"; \
	  objdirs="$${objdirs} $${tmpdir}"; \
 	done; \
	ar rvs $(LIBRARY) $^ $${objs}; \
	rm -rf $${objdirs};

endif

%.d:
	g++ -MM -I$(SUPERLU_INCDIR) $*.cpp > $*.d

clean:
	$(RM) $(ALL_OBJS) $(PROGRAM) $(LIBRARY) $(PROGRAM_LIB) $(IFP_LIB_SRC)

depends: $(DEPENDS)
	cat *.d > .depends
	$(RM) $(DEPENDS) 

TAGS: 
	etags *.h *.cpp

ifeq ($(wildcard .depends), .depends)
include .depends
endif

