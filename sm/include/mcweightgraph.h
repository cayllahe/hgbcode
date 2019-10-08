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

/* $Id: mcgraphe.h,v 1.3 2006/02/28 07:49:12 michel Exp $ */

/* ================================================ */
/* types publics */
/* ================================================ */
#define MAXEDGES 2000000
#define MAXDOUBLE FLOAT_MAX

#define MINDOUBLE (1e-999)

typedef struct cell {
  //  u_int32_t val;
  u_int32_t index;
  struct cell * next;
} cell;
typedef cell * pcell;

typedef struct {
  u_int32_t v_sommets;
  u_int32_t nbsom;
  u_int32_t nbmaxar;
  u_int32_t nbar; 
  u_int32_t ind;  /*number of edges*/
  double * weight; /*edge weight indexed by cell index*/
  u_int32_t *tete;
  u_int32_t *queue;
  u_int32_t cs;
  u_int32_t rs;
  cell * tasar;    /* tableau des cellules (tas) */
  cell * librear;  /* liste des cellules libres geree en pile lifo */
  pcell * listar;  /* tableau des listes d'aretes indexe par les sommets */
} graphe;

/* ================================================ */
/* prototypes */
/* ================================================ */

extern graphe * initgraphe(u_int32_t nbsom, u_int32_t nbmaxar);
extern void terminegraphe(graphe * g);
extern pcell allouecell(graphe * g);
extern void liberecell(graphe * g, pcell p);
extern void retiretete(graphe * g, pcell * pliste);
extern void retirearete(graphe * g, u_int32_t som, u_int32_t a);
extern int32_t estarete(graphe * g, u_int32_t som, u_int32_t a);
extern int32_t indexarete(graphe * g, u_int32_t som, u_int32_t a, u_int32_t ind);
extern void ajoutearete(graphe * g, u_int32_t som, u_int32_t a);
extern void addarete(graphe * g, u_int32_t som, u_int32_t a, double weight );
extern void setweight(graphe * g, u_int32_t index, double value);
extern void maille4graphe(graphe * g, u_int32_t rs, u_int32_t cs);
//extern void nettoiegraphe(graphe * g);
extern void aretesommets(u_int32_t a, u_int32_t N, u_int32_t rs, u_int32_t * s1, u_int32_t * s2);
extern int32_t estsuccesseur(graphe * g, u_int32_t som, u_int32_t a);
extern int32_t estsymetrique(graphe * g);
extern double getweight(graphe * g, u_int32_t index);
extern u_int32_t voisin(graphe *g, u_int32_t x, u_int32_t u);
extern graphe * ReadGraphe(char * filename,double **Fv);
extern void SaveGraphe(graphe * g, char *filename,double *Fv );
extern void setSize(graphe *g, int32_t rs , int32_t cs);
