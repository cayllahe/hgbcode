
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mcweightgraph.h>
#include <lhierarchie.h>
#include <time.h>
#include <lcomptree.h>
#include <string.h>


int main(int argc,char ** argv){
    graphe *g;
    double *Av;
    
    g = ReadGraphe(argv[1],&Av);
    
    if( argc==3 && ( strcmp(argv[2],"MST") == 0) ) {
        IncrementalQFZfromMST(g, argv[1]);
    }
    else{
        IncrementalQFZ(g, argv[1]);
    }
    
    
    
    //SaveGraphe(g, argv[2],Av);
    terminegraphe(g);
    free(Av);
}
