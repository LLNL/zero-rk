# ext directory makefile for Zero-RK

TARGETS = cvode_user

all: $(TARGETS)

.PHONY: cvode_user
cvode_user:
	cd ./cvode_user ; make all

.PHONY:clean
clean:
	cd ./cvode_user ; make clean

# .PHONY forces make to forces to make to follow a rule unconditionally, for
# some reason the TARGETS weren't being processed
