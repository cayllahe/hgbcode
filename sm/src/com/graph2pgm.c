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
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mccodimage.h>
#include <mcimage.h>
#include <mcweightgraph.h>
#include <mcutil.h>

void usage(char *arg){
  printf("#################################################################\n\n");
  printf("USAGE: %s OUT.graph mode IN.pgm \n",arg); 	
  printf("#################################################################\n\n");
}

/***********************************************************/
/* Input:                                                  */

/* Output:                                                 */

/*                                                         */

/***********************************************************/


int32_t main(argc, argv) 
     /* =============================================================== */
     int32_t argc; char **argv; 
{
  struct xvimage *Im;
  graphe *g;
  int32_t rs, cs, rsIm, csIm, ds, N, mode,i,j,x,y;
  pcell p;
  double *Fv;
  float *F, tmp;
  
  if(argc != 4){
    usage(argv[0]);
    exit(1);
  }
  
  mode = atoi(argv[2]);
  g = ReadGraphe(argv[1], &Fv);

  rs = g->rs;
  cs = g->cs;
  N = rs * cs;
  if( (rs == -1) || (cs == -1) ) {
    fprintf(stderr, "%s: %s is not the graph of an image\n", argv[0], argv[1]);
    exit(1);
  }
  
  switch(mode){
  case 0:
    /* Case of a 4-connected graph where each edge is weighted by the
       absolute difference of intensity between its extremity
       pixels */
    rsIm = 2*rs - 1;
    csIm = 2*cs - 1;
    Im = allocimage(NULL, rsIm, csIm, 1, VFF_TYP_FLOAT); 
    F = FLOATDATA(Im);

    for(x = 0; x < N; x++){
      for(p = g->listar[x]; p != NULL; p = p->next){
	y = g->queue[p->index];
	i = x%rs;
	j = x/rs;
	if(y == x) y = g->tete[p->index];
	if( (y - x) == 1)
	  F[2*i+1 + (2*j)*rsIm] =  (float)getweight(g,p->index);
	if( (y - x) == rs)
	  F[2*i + (2*j +1) * rsIm] = (float)getweight(g,p->index);
      }
    }

    for(i = 0; i < rs -1; i++)
      for(j = 0; j < cs - 1; j++)
      {
	tmp = max(F[2*i + (2*j+1)*rsIm], F[2*i + 1 + (2*j)*rsIm] );
	tmp = max(tmp,F[2*i+2 + (2*j+1)*rsIm]);
	F[2*i +1 + (2*j+1)*rsIm] = max(tmp, F[2*i+1 + (2*j + 2)*rsIm]);
      }
    writeimage(Im,argv[3]);
    freeimage(Im);
    break;
  default :
    fprintf(stderr, "%s: mode %d not available\n",argv[0], mode);
    exit(0);
  }
  terminegraphe(g);
  free(Fv);
  return 0;
}

