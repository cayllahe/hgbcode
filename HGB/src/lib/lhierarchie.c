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
#include <mcweightgraph.h>
#include <mcunionfind.h>
#include <lcomptree.h>
#include <llca.h>


/*********************************/
/*                               */
/* the output weighted graph g   */
/* is the ultrametric opening of */
/* the input weighted graph g    */
/*                               */
/*********************************/
//int32_t ultramopen(graphe *g)
//#define F_NAME "ultramopen"
int32_t ultramopen(graphe *g)
{  
  //int32_t logn, nbRepresent, u, x, y, c1;
 // int32_t nbarcs = g->ind;
  double *STAltitudes;
  JCctree *ST;

  /* Data structure for fast computation of LCA (least common ancestor) */
  //int32_t *Euler, *Depth, *Represent, *Number, **Minim;

  STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
  if (STAltitudes == NULL){
    fprintf(stderr,"ultramopen cannot allocate memory for STAltitudes\n");
    exit(0);
  }
  
    
    
  MST* mst;
  mst = (MST*)malloc(sizeof(MST));
  // To obtain the binary partition tree by altitude ordering
    
    
    
  BPTAO(g, &ST, STAltitudes, mst);
    
  // Computing the Quasi Flat Zones
  JCctree *QCT;
  CanonizeQBT(ST, STAltitudes, &QCT);
  
    
    
 
  //componentTreePrint(ST, STAltitudes, ST->root, NULL);
    
 // componentTreePrint(QCT, STAltitudes, QCT->root, NULL);
    
    
  double sm ;
  sm = computeSM(QCT, g, STAltitudes);
    
    
    
    printf("sum out: %f\n", sm);
    
  componentTreeFree(ST);
  free(STAltitudes);
  return 0;
}
#undef F_NAME

/*
JCctree* getQBT(graphe *g)
#define F_NAME "getQBT"
{
    
    double *STAltitudes;
    JCctree *ST;
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr,"getQBT cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    
    // To obtain the binary partition tree by altitude ordering
    BPTAO(g, &ST, STAltitudes);
    
    
    //componentTreePrint(ST, STAltitudes, ST->root, NULL);
    
    componentTreeFree(ST);
    free(STAltitudes);
    return 0;
}
#undef F_NAME

*/


