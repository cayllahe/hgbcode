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

void usage(char *arg){
  printf("#################################################################\n\n");
  printf("USAGE: %s IN.pgm mode OUT.graph \n",arg); 	
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
  int32_t rs, cs, ds, N, mode,i,j,x,y;
  uint8_t *F;
  double *Fv; 
  
  if(argc != 4){
    usage(argv[0]);
    exit(1);
  }
  
  mode = atoi(argv[2]);
  Im = readimage(argv[1]);
  if(depth(Im) != 1) {
    fprintf(stderr, "%s: 3D not yet implemented\n", argv[0]);
    exit(1);
  } 
  rs = rowsize(Im);
  cs = colsize(Im);
  N = rs*cs;

  switch(datatype(Im)){
  case 	VFF_TYP_1_BYTE: 
    F = UCHARDATA(Im);
    break;
  default : 
    fprintf(stderr,"%s: not yet implemented for this kind of data - only byte images implemented \n");
    exit(0);
  }
  
  switch(mode){
  case 0:
    /* Case of a 4-connected graph where each edge is weighted by the
       absolute difference of intensity between its extremity
       pixels */
    g = initgraphe(N, 2*(2*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++)
      for (i = 0; i < rs; i++)
      {
	if (i < rs - 1)
	{
	  x = j * rs + i;
	  y = j * rs + i + 1;
	  addarete(g, x, y, (double)(abs( (int)F[x] - (int)F[y])));
	}
	if (j < cs - 1)
	{
	  x = j * rs + i;
	  y = (j+1) * rs + i;
	  addarete(g, x, y, (double)(abs( (int)F[x] - (int)F[y])));
	}
      }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
      fprintf(stderr,"%s: cannot allocate (%d bytes) memory\n", argv[0], N*sizeof(double));
      terminegraphe(g);
      exit(1);
    }
    /* The area of any pixel is 1 */

    for(i = 0; i < N; i++){
      Fv[i] = 1.0;
    }

    SaveGraphe(g, argv[3], Fv);
    terminegraphe(g);
    free(Fv);
    break;
  default :
    fprintf(stderr, "%s: mode %d not available\n",argv[0], mode);
    exit(0);
  }
  freeimage(Im);
  return 0;
}
