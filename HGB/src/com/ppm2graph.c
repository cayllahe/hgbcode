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
  printf("USAGE: %s IN.ppm mode OUT.graph \n",arg);
  printf("When mode = 1 Euclidean distance in Lab space is used as gradient\n");
  printf("#################################################################\n\n");
}
/*
int32_t max(int32_t a, int32_t b){
  if (a<b) return b;
  else return a;
}*/

double decompand(double V){
  if (V <= 0.04045) return V/12.92;
  else return pow( (V+0.055) / 1.055, 2.4);		     
}

#define KAPPA 903.3
#define EPSILON  0.008856

double simplified(double v, double ref){
  double vv;
  vv = v / ref;

  if(vv > EPSILON) vv = pow(vv, 1.0/3.0);
  else vv = (KAPPA * vv + 16) / 116;

  return vv;
}

void lab(uint8_t R, uint8_t G, uint8_t B, double *l, double* a, double *b){
  double dr, dg, db, x, y, z;

  dr = decompand( ((double) R)/255 );
  dg = decompand( ((double) G)/255 );
  db = decompand( ((double) B)/255 );
  
  x = 0.4124564*dr +  0.3575761*dg  + 0.1804375*db;
  y = 0.2126729*dr +  0.7151522*dg  + 0.0721750*db;
  z = 0.0193339*dr +  0.1191920*dg  + 0.9503041*db;

  x = simplified(x, 0.95047);
  y = simplified(y, 1.0);
  z = simplified(z, 1.08883);
  


  (*l) = 116 * y - 16;
  (*a) = 500 * (x - y);
  (*b) = 200 * (y - z);
}

double labgrad(uint8_t Rx, uint8_t Gx, uint8_t Bx, uint8_t Ry, uint8_t Gy, uint8_t By){
  double lx, ax, bx, ly, ay, by;
  lab(Rx, Gx, Bx, &lx, &ax, &bx);
  lab(Ry, Gy, By, &ly, &ay, &by);
  
  return sqrt( (lx-ly) * (lx-ly) + (ax-ay) * (ax-ay) + (bx-by)*(bx-by) );   
}



int32_t main(argc, argv) 
     /* =============================================================== */
     int32_t argc; char **argv; 
{
  struct xvimage *R, *G, *B;
  graphe *g;
  int32_t rs, cs, ds, N, mode,i,j,x,y;
  uint8_t *fR, *fG, *fB;
  double *Fv; 
  double l, a, b;

  if(argc != 4){
    usage(argv[0]);
    exit(1);
  }
  
  mode = atoi(argv[2]);
  readrgbimage(argv[1], &R, &G, &B);
  if(depth(R) != 1) {
    fprintf(stderr, "%s: 3D not yet implemented\n", argv[0]);
    exit(1);
  } 
  rs = rowsize(R);
  cs = colsize(R);
  N = rs*cs;


  fR = UCHARDATA(R);
  fG = UCHARDATA(G);
  fB = UCHARDATA(B);


  switch(mode){
  case 1:
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
  	  addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
  	}
  	if (j < cs - 1)
  	{
  	  x = j * rs + i;
  	  y = (j+1) * rs + i;
  	  addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
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
  freeimage(R);
  freeimage(G);
  freeimage(B);
  return 0;
}
