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

void usage(char *arg){
  printf("#################################################################\n\n");
  printf("USAGE: %s IN.graph mode OUT.graph [bitmap.txt]\n",arg); 
  
  printf("If mode = 0 or 1: areas are computed. If mode = 2 or 3: dynamics are computed. If mode = 4 or 5: volumes are computed.\n");

  printf("If mode is odd then a bitmap (bitmap.txt) is produced to provide the mapping between order of extinctions and values of attributes.\n");
/*   printf("#################################################################\n\n"); */
}

/***********************************************************/
/* Input:                                                  */
/*   IN.graph an edge (and possibly vertex weighted graph) */
/*   mode : an integer                                     */
/* Output:                                                 */
/*   OUT.graph an edge and vertex weighted graph           */
/*   where the edge weights are the same as IN.graph and   */
/*   where the vertex weights represent the extinction     */
/*   values of the vertices in the minima of the edge      */
/*   graph IN.graph                                        */
/*   The optional file Bitmap.txt stores a mapping from    */ 
/*   {1, ..., N} into the measure of extinction            */  
/*   where the vertices of OUT.graph are weighted 1,..., N */
/* The measure used for the extinction depends of the      */
/* mandatory parameter mode                                */
/*                                                         */
/* For the moment only area is implemented                 */
/* Any even value of mode leads to absolute values in the  */
/* output graph, whereas any odd value leads to relative   */
/* values (vertices weighted with 1,...,N)                 */
/***********************************************************/

int32_t main(argc, argv) 
     /* =============================================================== */
     int32_t argc; char **argv; 
{
  graphe *g;
  double *Av;
  double *Bitmap;
  FILE *BitmapFile;
  int32_t mode, nbmin, i;

  if( (argc != 4) && (argc != 5)){
    usage(argv[0]);
    exit(1);
  }
  
  g = ReadGraphe(argv[1],&Av);
  mode = atoi(argv[2]);
  nbmin = lextinctionvalues(g,mode,Av, &Bitmap);
  SaveGraphe(g, argv[3],Av);
  if( (mode%2 != 0) && (argc == 5) ){
    if( (BitmapFile = fopen(argv[4],"w")) != NULL){    
      fprintf(BitmapFile, "#NbMin: %d\n", nbmin);
      for(i = 0; i < nbmin; i++)
	fprintf(BitmapFile, "%lg\n", Bitmap[i]);
      fclose(BitmapFile);
    } else
      fprintf(stderr, "%s, Bitmap saving failed: could not open %s\n", argv[0], argv[4]);
  } 

  terminegraphe(g);
  free(Av);
  free(Bitmap);
}

