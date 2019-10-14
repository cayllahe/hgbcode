#!/bin/bash
echo "Compiling Area Simplification"
cd AreaSimplification
make clean
make
cd ..
echo "Compiling graph Image Conversion"
cd graphImageConversion
make clean
./makelin
cd ..
echo "Compiling HGB"
cd HGB
make clean
./makelin
cd ..
echo "Done..."


