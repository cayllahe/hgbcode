//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#include <graphSegmentation.h>
#include <time.h>

void usage(char *arg){
    printf("#################################################################\n\n");
    printf("USAGE: %s IN.graph OUT.graph  Option\n",arg);
    //printf("Types \n");
    //printf("1  Naive dual HGB  \n");
    printf("#################################################################\n\n");
}

void loadGraph(char* argv,  graphe **g, double **Av ){
    *g = ReadGraphe(argv , Av);
}

int main(int argc,char ** argv){
    
    graphe *g;
    double *Av;
    
    int op = -1;
    if(argc!=4){
        usage( argv[0] );
        exit(1);
    }
    op = atoi(argv[3]);
    
    loadGraph(argv[1], &g , &Av);

    switch(op){
        case 1:
            dualhgb(g);
            break;
        case 2:
            //dualhgbInverse(g);
            dualhgbMin(g);
            break;
        case 3:
            dualhgbParentMin(g);
            break;
            
        case 4:
            dualhgbMean(g);
            break;
        case 5:
            dualhgbMeanNormal(g);
            break;
            
        case 6:
            dualhgbTruncated(g);
            break;
            
    }
    
    
    SaveGraphe(g, argv[2],Av) ;
    
    terminegraphe(g);
    free(Av);
    return 1;
}



