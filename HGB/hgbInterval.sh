
IN=$1
OUT=$2
TYPE=$3
AREA=0.005
MOD=2

filenameIN=$(basename "$IN")
fnameIN="${filename%.*}" 

#HgbSegInterval image.graph outfile.graph TYPE
# prueba con los tipos:
#TYPE 1 = original HGB
#TYPE 2 = max HGB
#TYPE 12 = rank HGB

linux/bin/HgbSegInterval $1 /tmp/$fnameIN_tmpSegC.graph $TYPE
../AreaSimplification/linux/bin/areaFilter /tmp/$fnameIN_tmpSegC.graph $AREA /tmp/$fnameIN_filteredC.graph
/Users/edus/Dropbox/MatMorphology/codigo/sm/linux/bin/graph2pgm /tmp/$fnameIN_filteredC.graph 0 /tmp/$fnameIN_filteredC.pgm
/Users/edus/Dropbox/MatMorphology/codigo/sm/linux/bin/float2byte /tmp/$fnameIN_filteredC.pgm $MOD $OUT
