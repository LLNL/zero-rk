#!/bin/bash

mkdir -p ../opt/{include,lib}

tar -xzvf spify.tar.gz


cd spify

sed -i 's|^YAMLCPPDIR.*$|YAMLCPPDIR:=../../opt|' Makefile
make

cp -rL ./include/spify ../../opt/include
cp lib/* ../../opt/lib
cp src/SpifyParserGenerator.py ../../opt/include/spify

cd ../
rm -rf spify

