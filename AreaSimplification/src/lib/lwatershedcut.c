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
#include <string.h>
#include <float.h>
#include <mclifo.h>
#include <mcweightgraph.h>
#include <lwatershedcut.h> 


#define MINJ(x,y) ((x)>(y)?(y):(x))
#define MAXJ(x,y) ((x)<(y)?(y):(x))

/*  
    ga en sortie est une M-border watershed de ga en entree.
    De plus retourne une carte de labels des sommets ds un minimum de F */
/*  mBorderWshed2d was previously called flowLPE2d*/
/* Codé par JC */
u_int32_t mBorderWshed2drapide(graphe *g, u_int32_t ** Res)
#undef F_NAME
#define F_NAME mBorderWshed2drapide
{
  u_int32_t x,y,z,u,w, nlabels, label;
  pcell p;
  u_int32_t *Vminima; 

  int32_t N = g->nbsom;                    /* taille image */
  int32_t M = g->ind;
  double *VF;                           /* fonction indicatrice sur les
					   sommets */
  double *F = g->weight;
  int32_t nbminimas;
  Lifo *L;
  
  /******************INITIALISATION*********************/ 
  
  Vminima = malloc(sizeof(u_int32_t) * N);
  if(Vminima == NULL){
    fprintf(stderr,"Cannot allocate memory for Vminima\n");
    exit(0);
  }

  if( (VF = malloc(sizeof(double) * N)) == NULL) {
    fprintf(stderr,"%s ne peut allouer VF \n", F_NAME);
    exit(1);
  }  

  /* Valuation des sommets */
  for(x = 0; x < N; x++) {
    VF[x] = FLT_MAX;
    Vminima[x] = -1;
    for(p = g->listar[x]; p != NULL; p= p->next)
      if (getweight(g,p->index) < VF[x]) VF[x] = getweight(g,p->index);
  }
  
  /* Initialisation de la LIFO */
  L = CreeLifoVide(M);
  /* Calcul des sommets qui sont dans des minima de F */
 
  nlabels = 0;

  for(x = 0; x < N; x++){
    if(Vminima[x] == -1){ /* On trouve un sommet non encore etiquete */ 
      nlabels ++;
      Vminima[x] = nlabels;
      LifoPush(L,x);
   
      while(!LifoVide(L)) {
	w = LifoPop(L);
	label = Vminima[w];
	for(p = g->listar[w]; p != NULL; p = p->next){
	  if(getweight(g,p->index) == VF[w]){
	    y = voisin(g,w,p->index);
	    if( (label > 0) && (VF[y] < VF[w]) ){
	      label = 0;
	      nlabels --;
	      Vminima[w] = label;
	      LifoPush(L,w);
	    }
	    else if(VF[y] == VF[w]){
	      if( ( (label > 0) && (Vminima[y] == -1) ) ||
		  ( (label == 0) && (Vminima[y] != 0)) ){
		Vminima[y] = label;
		LifoPush(L,y);
	      }
	    }
	  }
	}
      }
    }
  }
  //  printf("nlabels %d \n",nlabels);




  /* Les aretes adjacentes a un minimum sont inserees dans L */
  for(u = 0; u < M; u++){
    x = g->tete[u];
    y = g->queue[u];
    if( ( MINJ(Vminima[x],Vminima[y]) == 0) && /* un des deux sommets non ds un minima */
	(MAXJ(Vminima[x],Vminima[y]) > 0) ) /* et l'autre dans un minima */ 
      LifoPush(L,u); 
  }
   
  /*************BOUCLE PRINCIPALE*****************/
  while(!LifoVide(L)) {
    u = LifoPop(L);
    x = g->tete[u];
    y = g->queue[u];
    if (VF[x] > VF[y]) {z = y; y = x; x = z;}
    if((VF[x] < F[u]) && (VF[y] == F[u])){      /* u est une arete de bord */
      F[u] = VF[x];
      VF[y] = F[u];
      Vminima[y] = Vminima[x];
      for(p = g->listar[y]; p != NULL; p = p->next){
	z = voisin(g,y,p->index); 
	if(Vminima[z] == 0)
	  LifoPush(L, p->index);
      } /* for(p = g-> .. */
    } /* if( (VF[x] < ... */
  }/* while(!LifoVide... */

  free(VF);
  LifoTermine(L);
  (*Res) = Vminima;
  return nlabels;
}// mBorderWshed2drapide(...)

/* Input: a graph g which is a B-Kernel and a map Lab which labels the
   vertex sets of the minima of this B-kernel 

   Output: a graph g whose 0 level set is a forest such that each
   tree of this forest spans exactly one minimum of the input graph
   g */
void bKernelToRootedMSF(graphe *g, u_int32_t *Lab){
  int32_t M,N,x,y,z,curLab;
  double curVal;
  Lifo *L;
  pcell p;

  M =  g->ind;
  N = g -> nbsom;
  L = CreeLifoVide(M);
  
  for(x = 0; x < N; x++){
    if(Lab[x] != 0){
      curLab = Lab[x];
      curVal = FLT_MAX;
      for(p = g->listar[x]; p != NULL; p= p->next)
	if (getweight(g,p->index) < curVal) curVal = getweight(g,p->index);
      /* The CC that contains x is then depth-first explored */
      Lab[x] = 0;
      LifoPush(L,x);
      while(!LifoVide(L)){
	y = LifoPop(L);
	for(p = g->listar[y]; p != NULL; p= p->next){
	  z = voisin(g,y,p->index);
	  if( (Lab[z] == curLab) && (getweight(g,p->index) == curVal))
	  {
	    setweight(g, p->index ,-1);
	    LifoPush(L,z);
	    Lab[z] = 0;
	  }
	}
      }
    }
  }
  LifoTermine(L);
}
