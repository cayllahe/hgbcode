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
#include <lhierarchie.h>


void usage(char *arg){
  printf("#################################################################\n\n");
  printf("USAGE: %s IN.graph  OUT.graph [Bitmap.txt] \n",arg);
  printf("#################################################################\n\n");
}


int32_t main(argc, argv) 
     /* =============================================================== */
     int32_t argc; char **argv; 
{
  graphe *g;
  double *Av;
  FILE *BitmapFile;
  double *Bitmap;
  int32_t nbmin,u,i;

  if( (argc != 3) && (argc != 4) ){
    usage(argv[0]);
    exit(1);
  }
  
  g = ReadGraphe(argv[1],&Av);

  ultramopen(g);
  
  if(argc == 4){
    /* Then the bitmap is loaded and applied to the edge weights */
    if( (BitmapFile = fopen(argv[3],"r")) != NULL){
      fscanf(BitmapFile, "#NbMin: %d\n", &nbmin);  
      Bitmap = calloc(nbmin+1, sizeof(double));
      if(Bitmap == NULL){
	fprintf(stderr, "%d cannod allocate (%d byte) memory\n", argv[0], nbmin);
	exit(1);
      }
      Bitmap[0] = 0;
      for(i = 0; i < nbmin; i++){
	fscanf(BitmapFile, "%lf\n", &(Bitmap[i+1]));
      }
      fclose(BitmapFile);
      
      for(u = 0; u < g->ind; u++){
	i = (int32_t)getweight(g, u);
	  setweight(g,u, Bitmap[i]);
      }
      free(Bitmap);
    } else 
      fprintf(stderr,"%s cannot open %s\n", argv[0], argv[3]);
  }

  SaveGraphe(g, argv[2],Av);
  terminegraphe(g);
  free(Av);
}
