#!/bin/bash

mkdir -p ../opt/{include,lib}
ln -s ./lib ../opt/lib64

wget https://github.com/LLNL/spify/archive/v1.0.2.tar.gz
tar -xzf v1.0.2.tar.gz
cd spify-1.0.2

#git clone ~/src/spify
#cd spify/

mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=`pwd`/../../../opt
make install

cd ../../
rm -rf spify-1.0.2
rm v1.0.2.tar.gz
