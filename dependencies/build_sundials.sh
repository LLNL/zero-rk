#!/bin/bash


tar -xzvf sundials-2.5.0.tar.gz

cd sundials-2.5.0
source ../../build.config

./configure CC=gcc CFLAGS="-m64 -O3 -march=native -mfpmath=sse -mmmx -msse -msse2 -m3dnow -fPIC" \
            F77=gfortran FFLAGS="-m64 -O3 -mfpmath=sse -mmmx -msse -msse2 -m3dnow" \
            --with-ldflags="-L${LAPACK_LIB_DIR} -L${BLAS_LIB_DIR}" \
            --with-precision=double \
            --enable-static --enable-shared \
            --disable-examples \
            --disable-mpi \
            --prefix=${PWD}/../../opt/sundials-2.5.0

make

make install

cd ../
rm -rf sundials-2.5.0
