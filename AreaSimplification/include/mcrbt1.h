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


/* $Id: mcrbt1.h,v 1.4 2006/02/28 07:49:12 michel Exp $ */
#define RBT_Black 0
#define RBT_Red   1

typedef struct {
  double x;
  double y;
  double z;
} TypRbtKey;
typedef int32_t TypRbtAuxData;

#define EQUALKEY(j,k) ((j.x==k.x)&&(j.y==k.y)&&(j.z==k.z))
#define LESSKEY(j,k) ((j.x<k.x) || ((j.x==k.x)&&(j.y<k.y)) || ((j.x==k.x)&&(j.y==k.y)&&(j.z<k.z)))
#define COPYKEY(j,k) {j.x=k.x;j.y=k.y;j.z=k.z;}
#define PRINTKEY(j) printf("%g %g %g",j.x,j.y,j.z)

typedef struct RBTELT {
  TypRbtAuxData auxdata;
  TypRbtKey key;
  char color;
  struct RBTELT * left;
  struct RBTELT * right;
  struct RBTELT * parent;
} RbtElt;

typedef struct {
  int32_t max;             /* taille max du rbt (en nombre de points) */
  int32_t util;            /* nombre de points courant dans le rbt */
  int32_t maxutil;         /* nombre de points utilises max (au cours du temps) */
  RbtElt *root;        /* racine de l'arbre */
  RbtElt *nil;         /* sentinelle et element dont l'adresse joue le role de NIL */
  RbtElt *libre;       /* pile des cellules libres */
  RbtElt elts[1];      /* tableau des elements physiques */
} Rbt;

/* ============== */
/* prototypes     */
/* ============== */

extern Rbt * CreeRbtVide(
  int32_t taillemax);

extern void RbtFlush(
  Rbt * T);

extern int32_t RbtVide(
  Rbt * T);

extern void RbtTermine(
  Rbt * T);

extern void RbtPrint(
  Rbt * T);

extern RbtElt * RbtSearch(
  Rbt * T, TypRbtKey k);

extern RbtElt * RbtMinimum(
  Rbt * T, RbtElt * x);

extern RbtElt * RbtMaximum(
  Rbt * T, RbtElt * x);

extern RbtElt * RbtSuccessor(
  Rbt * T, RbtElt * x);

extern RbtElt * RbtInsert(
  Rbt ** T, TypRbtKey k, TypRbtAuxData d);

extern void RbtDelete(
  Rbt * T, RbtElt * z);

extern TypRbtAuxData RbtPopMin(
  Rbt * T);

extern TypRbtAuxData RbtPopMax(
  Rbt * T);

extern TypRbtKey RbtMinLevel(
  Rbt * T);
