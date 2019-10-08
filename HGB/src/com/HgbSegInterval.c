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
    printf("USAGE: %s IN.graph  OUT.graph TYPE Thresholding PRANK \n",arg);
    printf("Types \n");
    printf("1  Min  \n");
    printf("2  Max  \n");
    printf("3  Mean  \n");
    printf("4  Median \n");
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
    int op = -1;
    int th = -1;
    double param = -1;
    
    if(argc<4 || argc > 6){
        usage(argv[0]);
        exit(1);
    }
    op = atoi(argv[3]);
    
    if(argc==4){
        th = 500;
        param = 0.01;
    }
    else if(argc==5){
        th = atoi(argv[4]);
        param = 0.01;
    }
    else if(argc==6){
        th = atoi(argv[4]);
        param = atof(argv[5]);
    }
    else
        exit(1);
    /*
     
    if(!(op>=1&&op<=25)) {
        printf("Wrong type!, valid types\n");
        printf("1 Min \n");
        printf("2 Max  \n");
        printf("3   \n");
        printf("4  \n");
        printf("5   \n");
        printf("6    \n");
        printf("7  \n");
        printf("8   \n");
        printf("...   \n");
        exit(1);
    }*/
    
    loadGraph(argv[1], &g , &Av);
   
    //incrementalSegmentationInterval(g, op);
    intervalSegmentation(g, op, th, param  );
    
    SaveGraphe(g, argv[2],Av) ;
    
    terminegraphe(g);
    free(Av);
    return 1;
}



