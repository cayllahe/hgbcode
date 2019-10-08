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

/* authors : B. Perret */

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
  printf("USAGE: %s IN.pgm OUT.graph \n",arg); 	
  printf("Transform a pgm image ucm in byte to a 4 connected graph.\n");
  printf("The pgm is assumed to have a 1 pixel border padding.\n");  	
  printf("Hence with a pgm image of size (2*h+1, 2*w+1)\n"); 
  printf("we obtain a graph with h lines and w colonnes.\n"); 
  printf("#################################################################\n\n");
}




int32_t main(argc, argv) 
     /* =============================================================== */
     int32_t argc; char **argv; 
{
  struct xvimage *Im;
  graphe *g;
  int32_t  ds, i,j,x,y, pos;
  uint8_t *F;
  double *Fv; 
  
  if(argc != 3){
    usage(argv[0]);
    exit(1);
  }
  

  Im = readimage(argv[1]);
  if(depth(Im) != 1) {
    fprintf(stderr, "%s: 3D not yet implemented\n", argv[0]);
    exit(1);
  } 
  int32_t rs = rowsize(Im);
  int32_t cs = colsize(Im);
  int32_t N = rs*cs;
  
  int32_t rsD = rs/2;
  int32_t csD = cs/2;
  int32_t ND = rsD*csD;

  switch(datatype(Im)){
  case 	VFF_TYP_1_BYTE: 
    F = UCHARDATA(Im);
    break;
  default : 
    fprintf(stderr,"%s: not yet implemented for this kind of data - only byte images implemented \n");
    exit(0);
  }
  


    g = initgraphe(ND, 2*(2*ND-rsD-csD));
    setSize(g,rsD,csD);
    for (j = 0; j < csD; j++)
      for (i = 0; i < rsD; i++)
      {
	if (i < rsD - 1)
	{
	  pos = (2*j+1) * rs + 2*(i+1);
	  x = j * rsD + i;
	  y = j * rsD + i+1;
	  addarete(g, x, y, F[pos]);
	}
	if (j < csD - 1)
	{
	  pos = 2*(j+1) * rs + 2*i+1;
	  x = j * rsD + i;
	  y = (j+1) * rsD + i;
	  addarete(g, x, y, F[pos]);
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

    SaveGraphe(g, argv[2], Fv);
    terminegraphe(g);
    free(Fv);

  
  freeimage(Im);
  return 0;
}
