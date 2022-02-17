#!/bin/sh

cd ..

mkdir build
cd build
cmake ..
make -j8

mkdir xcode
cd xcode
cmake -G Xcode ../..

#cd ..
#cmake --build . --config Release
