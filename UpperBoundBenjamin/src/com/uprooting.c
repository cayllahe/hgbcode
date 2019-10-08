/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

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

/* authors : J. Cousty - L. Najman and M. Couprie */


#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mcweightgraph.h>
#include <lwatershedcut.h>
//#include <mcmesh.h>
//#include <mciomesh.h>
//#include <mccomptreemesh.h>
#include <lhierarchie.h>

#define MAX_LPE 255

#define MINJ(x,y) ((x)>(y)?(y):(x))
#define MAXJ(x,y) ((x)<(y)?(y):(x))

void usage(){
  printf("#################################################################\n\n");
  printf("USAGE: uprooting graphIN  graphOUT \n\n");
  printf("graphIN: graph filename   \n");
  printf("graphOUT: graph out filename  \n");
 	
  printf("#################################################################\n\n");
}


/* =============================================================== */
int32_t main(argc, argv) 
     /* =============================================================== */
     int32_t argc; char **argv; 

{
  /* Pour la LPE */
  graphe *g;
  double *Av,*Fv/*,*F*/;
  int32_t M,N,nblabels,i, nbUprootings;
  u_int32_t *LABELv;
  int32_t * seq;


  if(argc != 3){
    usage();
    fprintf(stderr, "\n", argv[0]);
    exit(1);
  }

  g = ReadGraphe(argv[1],&Av);
			
  M = g->ind;
  N = g -> nbsom;
  
  nblabels = mBorderWshed2drapide(g, &LABELv);
  bKernelToRootedMSF(g, LABELv);
  /* The zero level of g is then an MSF rooted in the minima of the input */

#ifdef DEBUG
  SaveGraphe(g, "/tmp/RootedMSF",Av);
#endif

  /* From the loaded input graph, extract the sequence S such that we
     afterward compute an uprooting by S.  

     In the loaded graph, the value Av[x] of a vertex x gives the
     orderded extinction value of x.
 */
  
  seq = (int32_t *)malloc(sizeof(int32_t)*(nblabels+1));
  if(seq == NULL)
    fprintf(stderr,"Erreur de malloc seq\n");
  
  nbUprootings = 0;
  for(i=0; i <N;i++){
    if(Av[i] > 0){
      seq[(int32_t)(Av[i])] = i;
      nbUprootings = MAXJ(nbUprootings,(int32_t)Av[i]);
    }
  }

  free(LABELv);
  
  uprooting(g, seq, MAXJ(nblabels, nbUprootings));
  SaveGraphe(g, argv[2],Av); 

  terminegraphe(g);
  free(Av);   
  return 0;
} /* main */

