# makefile for the zerork/ext/cvode_user directory

ZERORK_ROOT:=$(realpath ../../)
include ../../Makefile.defs

TARG_LIB = libcvode_user.a

# list of object files to be created
OBJS = cvode_user_sparselu.o cvode_user_timed_dense.o

# test executable
TEST = test_user_sparselu.x


# destinations for headerS and library
INC_DIR:=$(ZERORK_ROOT)/include
LIB_DIR:=$(ZERORK_ROOT)/lib
LOC_INC_DIR = $(INC_DIR)/ext

CINCL += -I$(SUNDIALS_ROOT)/include -I$(SUPERLU_ROOT)/SRC
CLIBS += -L$(SUNDIALS_ROOT)/lib -L$(SUPERLU_ROOT)/lib -L$(LIB_DIR)

CLIBS += -l$(SUPERLU_LIBNAME) -lsundials_cvode -lsundials_nvecserial

ifneq ($(BIGRED2_CRAY),Y)
  CLIBS += -llapack -lblas
endif

CLIBS += -lm

CLINK    := $(CC)

all: $(TEST)

$(TARG_LIB): $(OBJS) | $(LIB_DIR) $(LOC_INC_DIR)
	ar ruv $(TARG_LIB) $^
	ranlib $(TARG_LIB)
	cp *.h $(LOC_INC_DIR)/.
	cp $(TARG_LIB) $(LIB_DIR)/.

$(TEST): test_user_sparselu.c $(TARG_LIB)
	$(CLINK) $^ $(CFLAGS) $(CINCL) $(CLIBS) -o $@

clean:
	rm -f *.o *~ $(TARG_LIB) $(TEST)


$(INC_DIR):
	mkdir $(INC_DIR)

$(LOC_INC_DIR): | $(INC_DIR)
	mkdir $(LOC_INC_DIR)

$(LIB_DIR):
	mkdir $(LIB_DIR)

