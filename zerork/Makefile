# main makefile for the zerork project

TARGETS = src tools applications ext

all: $(TARGETS)

.PHONY: src
src:
	cd ./src ; make all

.PHONY: tools
tools:
	cd ./tools ; make all

.PHONY: applications
applications:
	cd ./applications ; make all

.PHONY: ext
ext:
	cd ./ext ; make all

clean:
	cd ./src ; make clean
	cd ./tools ; make clean
	cd ./applications ; make clean
	cd ./ext ; make clean
	rm -f lib/* include/zerork/*.h include/CKconverter/*.h
	rm -rf lib include
