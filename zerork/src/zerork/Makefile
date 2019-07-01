# makefile for the zerork directory

ZERORK_ROOT:=$(realpath ../../)
include ../../Makefile.defs

TARG_LIB:=libzerork.a

LOC_INC_DIR=$(INC_DIR)/zerork

# list of object files to be created
OBJS:= element.o species.o mechanism.o utilities.o \
       nasa_poly.o info_net.o rate_const.o \
       perf_net.o fast_exps.o external_funcs_blank.o \
       plog_reaction.o 

include ../Makefile.src.defs

fast_exps.o: fast_exps.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCL) -c fast_exps.cpp

#
