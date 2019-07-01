ZERORK_ROOT:=$(realpath ../../)
include ../../Makefile.defs

TARGETS = randomStateGen.x

CLINK    := $(CXX)
END_LIBS := -lm

all: $(TARGETS)

randomStateGen.x: randomStateGen.o
	$(CLINK) $(CXXFLAGS) $(CXXINCL) $< $(CXXLIBS) $(END_LIBS) \
          -o $@

clean:
	rm -f *.o *~ *.x
