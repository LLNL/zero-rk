#!/bin/bash


wget http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_4.3.tar.gz

tar -xzf superlu_4.3.tar.gz 

SUPERLUINST=${PWD}/../opt/SuperLU_4.3
mkdir -p ${SUPERLUINST}/{lib,SRC}

cd SuperLU_4.3


cat <<EOF > make.inc
PLAT = _linux
SUPERLULIB   	= ${SUPERLUINST}/lib/libsuperlu_4.3.a
TMGLIB       	= libtmglib.a

BLASDEF 	= -DUSE_VENDOR_BLAS
BLASLIB 	= -lblas

LIBS		= \$(SUPERLULIB) \$(BLASLIB)

ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc
CFLAGS       = -DPRNTlevel=0 -O3 -fPIC
NOOPTS       = -fPIC
FORTRAN	     = gfortran
FFLAGS       = -O2 -fPIC
LOADER       = \$(CC)
LOADOPTS     = -fPIC

CDEFS        = -DAdd_
EOF

make

cp SRC/*.h ${SUPERLUINST}/SRC/

cd ../
rm -rf SuperLU_4.3
