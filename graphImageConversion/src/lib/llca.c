#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mcweightgraph.h>
#include <lcomptree.h>
#include <assert.h>
 
/////////////////////////////////////////
// lca : Nearest (Lowest) Common Ancestor
// 
// From: The LCA Problem Revisited 
// M.A. Bender - M. Farach-Colton
//
// from lwshedtopo.c

// Depth-first preprocessing
int32_t LCApreprocessDepthFirst(JCctree *CT, int32_t node, int32_t depth, int32_t *nbr, int32_t *rep, int32_t *Euler, int32_t *Represent, int32_t *Depth, int32_t *Number)
{
  int32_t son;
  JCsoncell *sc;
  if (CT->tabnodes[node].nbsons > -1) {
    (*nbr)++;
    Euler[*nbr] = node;
    Number[node] = *nbr;
    Depth[node] = depth;
    Represent[*nbr] = node;
    (*rep)++;
    for (sc = CT->tabnodes[node].sonlist; sc != NULL; sc = sc->next)    {
      son = sc->son;
      LCApreprocessDepthFirst(CT, son, depth+1, nbr, rep, Euler, Represent, Depth, Number);
      Euler[++(*nbr)] = node;
    }
  }
  return *nbr;
}

int32_t ** LCApreprocess(JCctree *CT, int32_t *Euler, int32_t *Depth, int32_t *Represent, int32_t *Number, int32_t *nbR, int32_t *lognR)
{
  //O(n.log(n)) preprocessing
  int32_t nbr, rep, nbNodes;

  nbr = -1; // Initialization number of euler nodes
  rep = 0;
  nbr = LCApreprocessDepthFirst(CT, CT->root, 0, &nbr, &rep, Euler, Represent, Depth, Number);
  nbNodes = rep;

  // Check that the number of nodes in the tree was correct
  assert((nbr+1) == (2*nbNodes-1));

  int32_t nbRepresent = 2*nbNodes-1;
  int32_t logn = (int32_t)(ceil(log((double)(nbRepresent))/log(2.0)));
  *nbR = nbRepresent;
  *lognR = logn; 

  int32_t i,j,k1,k2;
  int32_t *minim = (int32_t *)calloc(logn*nbRepresent, sizeof(int32_t));
  int32_t **Minim = (int32_t **)calloc(logn, sizeof(int32_t*));
  Minim[0] = minim;

  for (i=0; i<nbRepresent-1; i++) {
    if (Depth[Euler[i]] < Depth[Euler[i+1]]) {
      Minim[0][i] = i;
    } else {
      Minim[0][i] = i+1;
    }
  }
  Minim[0][nbRepresent-1] = nbRepresent-1;

  for (j=1; j<logn; j++) {
    k1 = 1<<(j-1);
    k2 = k1<<1;
    Minim[j] = &minim[j*nbRepresent];
    for (i=0; i<nbRepresent; i++) {
      if ((i+ k2) >= nbRepresent) {
	Minim[j][i] = nbRepresent-1;
      } else {
	if (Depth[Euler[Minim[j-1][i]]] <= Depth[Euler[Minim[j-1][i+k1]]]) {
	  Minim[j][i] = Minim[j-1][i];
	} else {
	  Minim[j][i] = Minim[j-1][i+k1];
	}
      }
    }
  }
  return Minim;
}

int32_t LowComAncFast(int32_t n1, int32_t n2, int32_t *Euler, int32_t *Number, int32_t *Depth, int32_t **Minim)
#undef F_NAME
#define F_NAME "LowComAncFast"
{
  int32_t ii, jj, kk, k;

  ii = Number[n1];
  jj = Number[n2];
  if (ii == jj)
    return ii;

  if (ii > jj) {
    kk = jj;
    jj = ii;
    ii = kk;
  }

  k = (int32_t)(log((double)(jj - ii))/log(2.));

  if (Depth[Euler[Minim[k][ii]]] < Depth[Euler[Minim[k][jj-(1<<(k))]]]) {
    return Number[Euler[Minim[k][ii]]];
  } else {
    return Number[Euler[Minim[k][jj-(1<<k)]]];
  }
}
