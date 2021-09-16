#!/bin/bash

CXXFLAGS="-std=c++20 -isystem lib/sdsl-lite/include"
CXXFLAGS="$CXXFLAGS -O3 -march=native -DNDEBUG"
#CXXFLAGS="$CXXFLAGS -O3 -mpopcnt -mbmi -DNDEBUG"
#CXXFLAGS="$CXXFLAGS -O0 -ggdb"


CXX="ccache g++"

mkdir -p build
$CXX -c $CXXFLAGS src/main.cpp -o build/main.o
$CXX -c $CXXFLAGS src/utils/utils.cpp -o build/utils.o
$CXX -c $CXXFLAGS src/oss/expand.cpp -o build/expand.o
$CXX -c $CXXFLAGS src/oss/isValid.cpp -o build/isValid.o

$CXX -c $CXXFLAGS src/oss/generator/h2.cpp -o build/h2.o
$CXX -c $CXXFLAGS src/oss/generator/pigeon.cpp -o build/pigeon.o
$CXX -c $CXXFLAGS src/oss/generator/suffixFilter.cpp -o build/suffixFilter.o
$CXX -c $CXXFLAGS src/oss/generator/zeroOnesZero.cpp -o build/zeroOnesZero.o


$CXX build/main.o build/utils.o build/expand.o build/isValid.o build/h2.o build/pigeon.o build/suffixFilter.o build/zeroOnesZero.o -ldivsufsort64 -lfmt
