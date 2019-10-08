#!/bin/bash
IN=$1
OUT=$2
MOD=2
TYPE=$3
TH=$4
PRANK=$5

UNAME=$(uname)


filenameIN=$(basename "$IN")
fnameIN="${filename%.*}" 

filenameOUT=$(basename "$OUT")
fnameOUT="${filename%.*}" 


#tmpFile1=`mktemp`

graphIn=`basename graphIn_`
tmpGraph=`mktemp /tmp/${graphIn}XXXXXX`

tempfoo=`basename hgb_`
tmpFile1=`mktemp /tmp/${tempfoo}XXXXXX`

#echo $fname

#	echo "Linux"
./image2graph/linux/bin/ppm2graph $IN 4 $tmpGraph.graph
./HGB/linux/bin/HgbSegInterval $tmpGraph.graph $tmpFile1.graph $TYPE $TH $PRANK
./UpperBoundBenjamin/linux/bin/areaFilter $tmpFile1.graph 0.001 $tmpFile1_filteredC.graph
./sm/linux/bin/graph2pgm $tmpFile1_filteredC.graph 0 $tmpFile1_filteredC.pgm
./sm/linux/bin/float2byte $tmpFile1_filteredC.pgm $MOD $OUT
rm $tmpGraph.graph $tmpFile1.graph $tmpFile1_filteredC.graph $tmpFile1_filteredC.pgm $tmpFile1 $tmpGraph
display $OUT & 

