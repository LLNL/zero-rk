#!/bin/bash


export PATH=${PWD}/opt/bin:${PATH}
export LD_LIBRARY_PATH=${PWD}/opt/lib:${LD_LIBRARY_PATH}

cd dependencies
#./build_openmpi.sh
./build_yaml-cpp.sh
./build_spify.sh
./build_superlu.sh
./build_sundials.sh
cd ../

cd zerork
make all
cd ..

cd zerork_reactor
make all
cd ..

