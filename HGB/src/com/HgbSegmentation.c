//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
//#include <mcweightgraph.h>
//#include <lhierarchie.h>
//#include <dealpng.h>

#include <graphSegmentation.h>
#include <time.h>

void usage(char *arg){
    printf("#################################################################\n\n");
    printf("USAGE: %s IN.graph  OUT.graph  \n",arg);
    printf("#################################################################\n\n");
}

void loadGraph(char* argv,  graphe **g, double **Av ){

    *g = ReadGraphe(argv , Av);
    
}
/**
 Efficient scale computation with incremental QFZ
*/
int main(int argc,char ** argv){
    
    
    
    graphe *g;
    double *Av;
    
    if(argc!=3){
        usage(argv[0]);
        exit(1);
    }
    loadGraph(argv[1], &g , &Av);
    
    incrementalSegmentationBranch(g);
    
    SaveGraphe(g, argv[2],Av) ;
    
    terminegraphe(g);
    free(Av);
    return 1;
}

