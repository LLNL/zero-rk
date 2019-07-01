# main makefile for Zero-RK

TARGETS = CKconverter zerork derivative matrix utilities transport

all: $(TARGETS)

.PHONY: CKconverter
CKconverter:
	cd ./CKconverter ; make all

.PHONY:zerork
zerork:
	cd ./zerork ; make all

.PHONY:derivative
derivative:
	cd ./derivative ; make all

.PHONY:matrix
matrix:
	cd ./matrix ; make all

.PHONY:utilities
utilities:
	cd ./utilities ; make all


.PHONY:transport
transport:
	cd ./transport ; make all

.PHONY:clean
clean:
	cd ./CKconverter ; make clean
	cd ./zerork ; make clean
	cd ./derivative ; make clean
	cd ./matrix ; make clean
	cd ./utilities ; make clean
	cd ./transport ; make clean
# .PHONY forces make to forces to make to follow a rule unconditionally, for
# some reason the TARGETS weren't being processed
