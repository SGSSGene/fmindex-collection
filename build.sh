#!/bin/bash

g++ -std=c++20 src/divsufsort.cpp -ldivsufsort64 -isystem lib/sdsl-lite/include -lfmt -O3 -mpopcnt -mbmi
