#!/bin/bash

wget https://api.github.com/repos/jbeder/yaml-cpp/tarball/release-0.5.1 -O yaml-cpp-0.5.1.tar.gz
tar -xzvf yaml-cpp-0.5.1.tar.gz

mv jbeder-yaml-cpp-fa6a71e yaml-cpp-0.5.1


mkdir -p ../opt/{include,lib}

cd yaml-cpp-0.5.1

for f in src/*.cpp
do
  gcc -c -O3 -fPIC $f -o $f.o -I./include
done

ar ruv libyaml-cpp.a src/*.o
ranlib libyaml-cpp.a
cp -r ./include/* ../../opt/include
cp libyaml-cpp.a ../../opt/lib

cd ../
rm -rf yaml-cpp-0.5.1

