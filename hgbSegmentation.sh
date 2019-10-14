#!/bin/bash
IN=$1
OUT=$2
MOD=2
TYPE=$3

PARAM=$4
AREASIMP=$5

if [ $# -ne 5 ]
then
	AREASIMP=0.004	
fi 	

UNAME=$(uname)


filenameIN=$(basename "$IN")
fnameIN="${filename%.*}" 

filenameOUT=$(basename "$OUT")
fnameOUT="${filename%.*}" 


graphIn=`basename graphIn_`
tmpGraph=`mktemp /tmp/${graphIn}XXXXXX`

tempfoo=`basename hgb_`
tmpFile1=`mktemp /tmp/${tempfoo}XXXXXX`

./graphImageConversion/linux/bin/ppm2graph $IN 4 $tmpGraph.graph
./HGB/linux/bin/HgbSegInterval $tmpGraph.graph $tmpFile1.graph $TYPE $PARAM
./AreaSimplification/linux/bin/areaFilter $tmpFile1.graph $AREASIMP $tmpFile1_filteredC.graph
./graphImageConversion/linux/bin/graph2pgm $tmpFile1_filteredC.graph 0 $tmpFile1_filteredC.pgm
./graphImageConversion/linux/bin/float2byte $tmpFile1_filteredC.pgm $MOD $OUT
rm $tmpGraph.graph $tmpFile1.graph $tmpFile1_filteredC.graph $tmpFile1_filteredC.pgm $tmpFile1 $tmpGraph
display $OUT & 

