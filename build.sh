#!/bin/bash

CXXFLAGS="-std=c++20 -isystem lib/sdsl-lite/include"
#CXXFLAGS="$CXXFLAGS -O3 -march=native -DNDEBUG"
CXXFLAGS="$CXXFLAGS -O3 -mpopcnt -mbmi -DNDEBUG"
#CXXFLAGS="$CXXFLAGS -O0 -ggdb"


mkdir -p build
g++ -c $CXXFLAGS src/main.cpp -o build/main.o
g++ -c $CXXFLAGS src/utils/utils.cpp -o build/utils.o
g++ -c $CXXFLAGS src/oss/expand.cpp -o build/expand.o
g++ -c $CXXFLAGS src/oss/isValid.cpp -o build/isValid.o

g++ -c $CXXFLAGS src/oss/generator/h2.cpp -o build/h2.o
g++ -c $CXXFLAGS src/oss/generator/pigeon.cpp -o build/pigeon.o
g++ -c $CXXFLAGS src/oss/generator/suffixFilter.cpp -o build/suffixFilter.o
g++ -c $CXXFLAGS src/oss/generator/zeroOnesZero.cpp -o build/zeroOnesZero.o


g++ build/main.o build/utils.o build/expand.o build/isValid.o build/h2.o build/pigeon.o build/suffixFilter.o build/zeroOnesZero.o -ldivsufsort64 -lfmt
