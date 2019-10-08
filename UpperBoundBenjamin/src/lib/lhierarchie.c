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
#include <fib.h>
#include <lcomptree.h>
#include <llca.h>

int32_t uprooting(
		  graphe *g,    /* weighted graph (G,F) st MSF(M_F) is the graph induced by the edges of altitude -1 */
		  int32_t *seq, /*  sequence S of the minima such that an uprooting by S is returned, each minimum is identified by one of its verices */
		  int32_t l     /* length of seq */
)
{

  Tarjan *T;
  struct fibheap **heap;
  void **clefs; // le numero des aretes (ici 2 aretes differente encodent (x,y) et (y,x))
  int32_t *H; // Indicator of the uprooting
  int32_t xp,yp,v,y,zp,i,M,N;
  void* vv;

  M = g->ind;
  N = g -> nbsom;

  T = CreeTarjan(N);
  if(T == NULL)
    fprintf(stderr,"Ne peut allouer Tarjan\n");

  heap = (struct fibheap **) malloc(sizeof(struct fibheap*)*N);
  if(heap == NULL)
    fprintf(stderr,"Ne peut allouer Heap\n");

  clefs = (void**) malloc(sizeof(void*)*M*2);
  if(clefs == NULL)
    fprintf(stderr,"Ne peut allouer clefs\n");

  H = (int32_t *)calloc(1,M * sizeof(int32_t));	
  if(H == NULL)
    fprintf(stderr,"Ne peut allouer H\n");

  // Uprooting Init

  for(i=0; i <N; i++){
    TarjanMakeSet(T, i);
    heap[i] = fh_makekeyheap();
  }
  
  for(i = 0; i <M; i++) 
    if(getweight(g,i) == -1)
      H[i] =0;
    else H[i] = l +1;

  for(i = 0; i < M; i++)
    if(getweight(g,i) == -1){ /* Si l'arete i appartient a la foret de poids min */
      xp = TarjanFind(T,g->tete[i]);
      yp = TarjanFind(T,g->queue[i]);
     
      if(xp != yp)
	TarjanLinkSafe(T,xp,yp);
    }

  for(i = 0; i < M; i++)
    if(getweight(g,i) > 0){
      xp = TarjanFind(T,g->tete[i]);
      yp = TarjanFind(T,g->queue[i]);
      if(xp != yp){
	/* The edge i is direct for xp  since it starts in xp. Thus i is inserted in the heap of xp. */
	clefs[i] = fh_insertkey(heap[xp], (int)getweight(g,i), (void*)i);
	/* The edge i is indirect for yp since it ends in yp. Thus M+i is inserted in the heap of xp. */
	clefs[M+i] = fh_insertkey(heap[yp], (int)getweight(g,i), (void*)(M+i));
      }
    }

  // Incremental Uprooting
  for(i=1; i < l; i++){
    xp = TarjanFind(T, seq[i]);
    do{
      vv = fh_extractmin(heap[xp]);
      v = (int)vv;
      if( (v/M) == 0){
	/* The edge v is direct : it starts in a vertex of the component labelled seq[i] */
	y = g->queue[v];
      } else {
	/* The edge v is indirect : its symetric v%M ends in a vertex of the component labelled seq[i] */
	y = g->tete[v%M];
      }
      yp = TarjanFind(T,y);

#ifdef DEBUG
      fprintf(stderr,"Arete min extraite (x,y) = (%d,%d) et representante (xp,yp) = (%d,%d) a la valeur %d\n",  g->queue[v%M], g->tete[v%M] ,xp,yp, (int)getweight(g,v%M));
#endif

    }while(yp == xp);

#ifdef DEBUG 
    fprintf(stderr,"Edge v is outgoing: xp actually disctinct from yp\n");
#endif 

    /* The weight of edge v (and of its symmetric) in the indicator of
       the uprooting is i. We only store the weight of the direct
       edge */

    H[v%M] = i;

#ifdef DEBUG
    fprintf(stderr,"H[ %d = (%d,%d) ] = %d et v =%d \n",v%M,g->tete[v%M], g->queue[v%M],i,v);
#endif     

    zp = TarjanLinkSafe(T,xp,yp);

    if(zp == xp){
    } else
      if (zp == yp){
	yp = xp;
	xp = zp;
      } else fprintf(stderr,"Gros probleme 2\n");
    
    heap[zp] = fh_union(heap[xp],heap[yp]);
    heap[yp] = NULL;
  }
   
  for (i = 0; i < M; i++){
      setweight(g,i,H[i]);
  }

  
  for(i=0; i <N; i++){
    if(heap[i] != NULL) fh_deleteheap(heap[i]);
  }


  TarjanTermine(T);
  
  free(seq);
  free(clefs);
  free(H);
  free(heap);
  return 1;
}




/*********************************/
/*                               */
/* the output weighted graph g   */
/* is the ultrametric opening of */
/* the input weighted graph g    */
/*                               */
/*********************************/
int32_t ultramopen(graphe *g)
#define F_NAME "ultramopen"
{  
  int32_t logn, nbRepresent, u, x, y, c1;
  int32_t nbarcs = g->ind;
  double *STAltitudes;
  JCctree *ST;

  /* Data structure for fast computation of LCA (least common ancestor) */
  int32_t *Euler, *Depth, *Represent, *Number, **Minim;

  STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
  if (STAltitudes == NULL){
    fprintf(stderr,"ultramopen cannot allocate memory for STAltitudes\n");
    exit(0);
  }
 
  saliencyTree(g, &ST, STAltitudes);

  Euler = (int32_t *)calloc(2*ST->nbnodes-1, sizeof(int32_t));
  Represent = (int32_t *)calloc(2*ST->nbnodes-1, sizeof(int32_t));
  Depth = (int32_t *)calloc(ST->nbnodes, sizeof(int32_t));
  Number = (int32_t *)calloc(ST->nbnodes, sizeof(int32_t));
  
  if ((Euler == NULL) || (Represent == NULL)
      || (Depth == NULL) || (Number == NULL)) {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return 0;
  }
  
  Minim = LCApreprocess(ST, Euler, Depth, Represent, Number, &nbRepresent, &logn);

  // For any edge of g
  for(u = 0; u < nbarcs; u++){
    x = g->tete[u]; y = g->queue[u];
    c1 = Represent[LowComAncFast(x, y, Euler, Number, Depth, Minim)];
    setweight(g,u,STAltitudes[c1]);
  }  
  free(Euler);
  free(Represent);
  free(Depth);
  free(Number);
  free(Minim[0]);
  free(Minim);
  componentTreeFree(ST);
  free(STAltitudes);
  return 0;
}
#undef F_NAME
