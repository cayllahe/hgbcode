/*
Copyright ESIEE (2019) 

edward.cayllahua@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/* authors : Edward Cayllahua - J. Cousty and Silvio Guimarães */

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#include <graphSegmentation.h>
#include <time.h>

void usage(char *arg){
    printf("#################################################################\n\n");
    printf("USAGE: %s IN.graph  OUT.graph TYPE Parameter \n",arg);
    printf("Types \n");
    printf("1: Use Min: select minimum value on positive observation intervals\n"); 
    printf("2: Use Max: select the last upper bound on negative intervals\n");     
    printf("3: Lower-length: On positive intervals, apply length threshold and min-rule\n");
    printf("4: Upper-length: apply length threshold and max-rule\n");
    printf("5: Lower-area: On positive intervals, apply area and  min-rule\n");
    printf("6: Upper-Narea: On negative intervals, apply area and  max-rule\n");            
    printf("7: Lower-depth: On positive intervals, apply depth filter and min-rule\n");
    printf("8: Upper-Ndepth: On negative intervals, apply depth filter and max-rule\n");            
    printf("9: Lower p-rank: On positive intervals, apply rank filter and  min-rule\n");
    printf("10: Upper p-rank: On negative intervals, apply rank filter and  max-rule\n");
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
        param = atof(argv[4]);
    }
    else
        exit(1);
     
    if(!(op>=1&&op<=10)) {
        usage(argv[0]);
        exit(1);
    }    
    loadGraph(argv[1], &g , &Av);
    
    intervalSegmentation(g, op, param, param  );
    
    SaveGraphe(g, argv[2],Av) ;
    
    terminegraphe(g);
    free(Av);
    return 1;
}



