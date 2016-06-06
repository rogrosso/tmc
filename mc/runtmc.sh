#!/bin/bash

clang++ -o tmc tmc.cpp MarchingCubes.cpp -std=c++14
./tmc skull128.bin skull128.obj skull128.off
