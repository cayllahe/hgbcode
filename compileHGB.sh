#!/bin/bash
echo "Compiling area filtering"
cd UpperBoundBenjamin
make clean
make
cd ..
echo "Compiling saliency map"
cd sm
make clean
./makelin
cd ..
echo "Compiling image to graph"
cd image2graph
make clean
./makelin
cd ..
echo "Compiling HGB"
cd HGB
make clean
./makelin
cd ..
echo "Done..."


