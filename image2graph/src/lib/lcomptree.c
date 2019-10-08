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
#include <string.h>
#include <mcweightgraph.h>
#include <mcunionfind.h>
#include <mcsort.h>
#include <lcomptree.h>


/* ==================================== */
JCctree * componentTreeAlloc(int32_t N)
/* ==================================== */
#undef F_NAME
#define F_NAME "componentTreeAlloc"
{
  JCctree *CT;
  CT = (JCctree *)malloc(sizeof(JCctree)); 
  CT->tabnodes = (JCctreenode *)malloc(N * sizeof(JCctreenode));
  CT->tabsoncells = (JCsoncell *)malloc(2*N * sizeof(JCsoncell));
  CT->flags = NULL;
  // CT->flags = (uint8_t *)calloc(N, sizeof(char));
  // If needed, must be allocated separatly
  //  memset(CT->flags, 0, N);
  if ((CT == NULL) || (CT->tabnodes == NULL) || (CT->tabsoncells == NULL))
  { 
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  CT->nbnodes = N;
  CT->nbsoncells = 0;
  return CT;
} // componentTreeAlloc()

/* ==================================== */
void componentTreeFree(JCctree * CT)
/* ==================================== */
{
  free(CT->tabnodes);
  free(CT->tabsoncells);
  if(CT->flags != NULL)free(CT->flags);
  free(CT);
} // ComponentTreeFree()


/* ==================================== */
void componentTreePrint(JCctree *CT, double *Alt, int32_t root, int32_t *Att1)
/* ==================================== */
{
  JCsoncell *s; 
  char* c;
  if(root <= CT->nbnodes/2) c = "(leaf)"; else c = "(branch)";
  printf("---------------------\nnode %d %s\nAltitude: %lf\n", root, c, Alt[root]);
  if( Att1 != NULL )
    printf("Attribut 1: %d\n",Att1[root]);
  
  printf("children: ");

  for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
    if(s->son <= CT->nbnodes/2) c = 'S'; else c = 'C';
    printf("%c%d (altitude %lf) ",c,s->son, Alt[s->son]);
  }
  printf("\n---------------------\n");
  for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next)
    componentTreePrint(CT, Alt, s->son, Att1);
}

/* ==================================== */
void calculReversePointer(JCctree *CT, int32_t root)  
/* ==================================== */
{
  JCsoncell *s; 
  for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next) 
  {
    calculReversePointer(CT, s->son);
    CT->tabnodes[s->son].father = root;
  }
}

/* The algorithm is an original adaptation/improvment to edge-weighted
   graph of the algorithm of Najman and Couprie published in IEEE IP
   in 2006 / Its description is provided in Cousty, Najman, Perret
   ISMM2013 and in Najman, Cousty, Perret ISMM 2013. */
int32_t BPTAO(/* INPUTS */
		     graphe *g, 
		     /* OUTPUTS */
		     JCctree ** SaliencyTree, /* Component tree - BPTAO  */
		     double *STaltitude  /* weights associated to the
					    nodes of the BPTAO: Must
					    be allocated elsewhere */
		     )
{
  int32_t i,x1,x2,n1,n2,z,k, STx, STy, nbsoncellsloc;
  JCctree *ST;
  int32_t *clefs; 
  int32_t *STmap;
  JCsoncell * newsoncell1;
  JCsoncell * newsoncell2;
  Tarjan *T;
  int32_t taille = g->nbsom;
  int32_t narcs = g->ind;
  
  if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de malloc\n"); 
  }
   
  if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
    fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
    exit(0);
  }
  
  for(i = 0; i < narcs; i++)
    clefs[i] = i; 
  d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
  if( (ST = componentTreeAlloc((2*taille))) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
    exit(0);
  }
  ST->nbnodes = taille;
  ST->nbsoncells = 0;
  T = CreeTarjan(taille);
  for(i = 0; i < taille; i++)
  {
    STmap[i] = i;
    ST->tabnodes[i].nbsons = 0;
    ST->tabnodes[i].sonlist = NULL;
    TarjanMakeSet(T,i);
  }

  for(i = 0; i < narcs; i++){ 
    // for each edge of the graph or MST taken in increasing order of altitude
    n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
    if(n1 != n2) {
      /* If the two vertices do not belong to a same connected
	 component at a lowest level */
      // Which component of ST n1 and n2 belongs to ?
      x1 = STmap[n1]; x2 = STmap[n2];    
      // Create a new component
      z = ST->nbnodes; ST->nbnodes++; 
      nbsoncellsloc = ST->nbsoncells;
      ST->tabnodes[z].nbsons = 2;
      // the altitude of the new component is the altitude of the edge
      // under consideration
      STaltitude[z] = getweight(g, clefs[i]);
      // add x1 and x2 to the lists of sons of the new component
      newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);  
      newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
      ST->tabnodes[z].sonlist = newsoncell1;
      newsoncell1->son = x1;
      newsoncell1->next = newsoncell2;
      newsoncell2->son = x2;
      newsoncell2->next = NULL;
      ST->tabnodes[z].lastson = newsoncell2;
      ST->nbsoncells += 2
;
    
      // then link n1 and n2
      k = TarjanLink(T, n1, n2);
      STmap[k] = z;
    }//if(...)
  }//for(...)

  
  for(i = 0; i < taille; i++)
    STaltitude[i] = 0; /* altitudes of the leaves is 0 */

  /* Construct father relationship */
  ST->tabnodes[ST->nbnodes-1].father = -1;
  ST->root = ST->nbnodes - 1;
  calculReversePointer(ST, ST->root); 


  (*SaliencyTree) = ST;

  /* liberation de la memoire */
  free(STmap); 
  free(clefs);
  TarjanTermine(T);
  return ST->nbnodes - 1;
}

void compactRec(JCctree *CT, double *Val, int32_t root){
  JCsoncell *tmp;
  JCsoncell *f1 = CT->tabnodes[root].sonlist;
  JCsoncell *f2 = CT->tabnodes[root].lastson;
  
  if (f1 == NULL){ /* A leaf of the tree, that is a vertex of the initial graph */
    Val[root] = Val[CT->tabnodes[root].father];
  } else {
    compactRec(CT, Val, f1->son);
    compactRec(CT, Val, f2->son);

    if( (Val[root] == Val[f1->son]) && (f1->son > CT->nbnodes/2) ){
      /* Then f1 must be compacted */
      CT->tabnodes[root].sonlist = CT->tabnodes[f1->son].sonlist;
      tmp = CT->tabnodes[f1->son].lastson;
    } else tmp = f1;

    if( (Val[root] == Val[f2->son]) && (f2->son > CT->nbnodes/2) ){
      /* Then f1 must be compacted */
      tmp->next = CT->tabnodes[f2->son].sonlist;
      CT->tabnodes[root].lastson = CT->tabnodes[f2->son].lastson;
    } else  tmp->next = f2;
  }
}


int32_t QFZ(/* INPUTS */
		     graphe *g, 
		     /* OUTPUTS */
		     JCctree **CompTree, /* Component tree */
		     double *CTaltitude  /* weights associated to the
					    nodes of the component
					    tree : Must be allocated
					    elsewhere */
		      )
{
  QFZ(g, CompTree, CTaltitude);
  compactRec((*CompTree), CTaltitude, (*CompTree)->root);
  calculReversePointer((*CompTree), (*CompTree)->root);
  return 1;
}

