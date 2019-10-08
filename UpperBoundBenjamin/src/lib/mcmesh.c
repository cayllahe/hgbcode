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

/* $Id: mcmesh.c,v 1.10 2006/10/18 11:51:04 michel Exp $ */
/* 
  Gestion d'une triangulation
  Michel Couprie  -  Mai 2001
  Update Fevrier 2002 : Edges et lissage courbures
  Update Fevrier 2004 : Addition de bruit gaussien
  Update Fevrier 2004 : RegulMeshHamam
  Update F�rier 2006 : mesures
*/

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mcmesh.h>
#include <mcrbt1.h>
#include <mcgeo.h>
#include <mcprobas.h>
#include <mcutil.h>

#ifndef EPSILON
#define EPSILON 1E-20
#endif

/*
#define DEBUG
*/
#define PARANO
#define VERBOSE
#define MESURE

#define NITERMAX 20000

meshtabvertices *Vertices = NULL;
meshtabfaces *Faces = NULL;
meshtabedges *Edges = NULL;
meshtablinks *Links = NULL;
Rbt * RBT;

double normevect(double x, double y, double z)
{
  return sqrt(x*x + y*y + z*z);
}

/* ==================================== */
meshtabvertices * AllocVertices(int32_t taillemax)
/* ==================================== */
#undef F_NAME
#define F_NAME "AllocVertices"
{
  meshtabvertices * T = (meshtabvertices *)calloc(1,sizeof(meshtabvertices) + taillemax*sizeof(meshvertex));
  if (T == NULL) 
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  T->lab = (u_int8_t *)calloc(taillemax, sizeof(char));
  if (T->lab == NULL) 
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  T->tmp = (u_int8_t *)calloc(taillemax, sizeof(char));
  if (T->tmp == NULL) 
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  T->max = taillemax;
  T->cur = 0;
  return T;
} /* AllocVertices() */

/* ==================================== */
void ReAllocVertices(meshtabvertices **A)
/* ==================================== */
{
  int32_t i, taillemax;
  meshtabvertices * T, *Tmp;

  //printf("ReAllocVertices: ancienne taille %d nouvelle taille %d\n", (*A)->max, 2 * (*A)->max);

  taillemax = 2 * (*A)->max;  /* alloue le double de l'ancienne taille */ 
  T = AllocVertices(taillemax);
  T->cur = (*A)->cur;
  memcpy(T->v, (*A)->v, T->cur * sizeof(meshvertex));
  memcpy(T->lab, (*A)->lab, T->cur * sizeof(char));
  memcpy(T->tmp, (*A)->tmp, T->cur * sizeof(char));
  Tmp = *A;
  *A = T;
  free(Tmp->lab);
  free(Tmp->tmp);
  free(Tmp);
} /* ReAllocVertices() */

/* ==================================== */
meshtabfaces * AllocFaces(int32_t taillemax)
/* ==================================== */
#undef F_NAME
#define F_NAME "AllocFaces"
{
  meshtabfaces * T = (meshtabfaces *)calloc(1,sizeof(meshtabfaces) + taillemax*sizeof(meshface));
  if (T == NULL) 
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  T->max = taillemax;
  T->cur = 0;
  return T;
} /* AllocFaces() */

/* ==================================== */
void ReAllocFaces(meshtabfaces **A)
/* ==================================== */
{
  int32_t i, taillemax;
  meshtabfaces * T, *Tmp;

  //printf("ReAllocFaces: ancienne taille %d nouvelle taille %d\n", (*A)->max, 2 * (*A)->max);

  taillemax = 2 * (*A)->max;  /* alloue le double de l'ancienne taille */ 
  T = AllocFaces(taillemax);
  T->cur = (*A)->cur;
  memcpy(T->f, (*A)->f, T->cur * sizeof(meshface));
  Tmp = *A;
  *A = T;
  free(Tmp);
} /* ReAllocFaces() */

/* ==================================== */
meshtabedges * AllocEdges(int32_t taillemax)
/* ==================================== */
#undef F_NAME
#define F_NAME "AllocEdges"
{
  meshtabedges * T = (meshtabedges *)calloc(1,sizeof(meshtabedges) + taillemax*sizeof(meshedge));
  if (T == NULL) 
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  T->max = taillemax;
  T->cur = 0;
  return T;
} /* AllocEdges() */

/* ==================================== */
void ReAllocEdges(meshtabedges **A)
/* ==================================== */
{
  int32_t i, taillemax;
  meshtabedges * T, *Tmp;

  //printf("ReAllocEdges: ancienne taille %d nouvelle taille %d\n", (*A)->max, 2 * (*A)->max);

  taillemax = 2 * (*A)->max;  /* alloue le double de l'ancienne taille */ 
  T = AllocEdges(taillemax);
  T->cur = (*A)->cur;
  memcpy(T->e, (*A)->e, T->cur * sizeof(meshedge));
  Tmp = *A;
  *A = T;
  free(Tmp);
} /* ReAllocEdges() */

/* ==================================== */
meshtablinks * AllocLinks(int32_t nvert, int32_t nedge)
/* ==================================== */
#undef F_NAME
#define F_NAME "AllocLinks"
{
  meshtablinks * T = (meshtablinks *)calloc(1,sizeof(meshtablinks));
  if (T == NULL) 
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  T->nvert = nvert;
  T->nedge = nedge;
  T->lastneigh = (int32_t *)calloc(1,nvert * sizeof(int32_t));
  T->neigh = (int32_t *)calloc(1,nedge * sizeof(int32_t));
  if ((T->lastneigh == NULL) || (T->neigh == NULL))
  {
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  return T;
} /* AllocLinks() */

/* ==================================== */
void InitMesh(int32_t taillemax)
/* ==================================== */
{
  Vertices = AllocVertices(taillemax);
  Faces = AllocFaces(taillemax);
  RBT = CreeRbtVide(taillemax);
} /* InitMesh() */

/* ==================================== */
void TermineMesh()
/* ==================================== */
{
  free(Vertices->lab);
  free(Vertices->tmp);
  free(Vertices);
  free(Faces);
  if (Edges) free(Edges);
  if (Links) { free(Links->lastneigh); free(Links->neigh); free(Links); }
  RbtTermine(RBT);
} /* TermineMesh() */

/* ==================================== */
int32_t NotIn(
  int32_t e,
  int32_t *list,                   
  int32_t n)                       
/* ==================================== */
{
/* renvoie 1 si e n'est pas dans list, 0 sinon */
/* e : l'element a rechercher */
/* list : la liste (tableau d'entiers) */
/* n : le nombre d'elements dans la liste */
  while (n > 0)
    if (list[--n] == e) return 0;
  return 1;
} /* NotIn() */

/* ==================================== */
int32_t AddVertex(double x, double y, double z, int32_t indface)
/* ==================================== */
/* modifie les var. globales Faces, Vertices */
#undef F_NAME
#define F_NAME "AddVertex"
{
  int32_t i;
  TypRbtKey point;
  RbtElt * re;

  /* cherche si le sommet est deja la */
  point.x = x;
  point.y = y;
  point.z = z;
  re = RbtSearch(RBT, point); 
  if (re != RBT->nil)
  {
    i = re->auxdata; /* index du vertex */
    /* il est la : on lui ajoute la face si elle n'y est pas deja */
    if (NotIn(indface, Vertices->v[i].face, Vertices->v[i].nfaces)) 
    {
      if (Vertices->v[i].nfaces >= MAXADJFACES)
      {
        fprintf(stderr, "%s : WARNING: more than %d faces\n", F_NAME, MAXADJFACES);
        fprintf(stderr, "x=%g, y=%g, z=%g\n", x, y, z);
        goto skipadd;
      }
      Vertices->v[i].face[ Vertices->v[i].nfaces++ ] = indface;
    }
skipadd:
    return i;
  } /* if (re != RBT->nil) */
  if (Vertices->cur >= Vertices->max) ReAllocVertices(&Vertices);
  i = Vertices->cur;
  Vertices->cur += 1;
  Vertices->v[i].x = x;
  Vertices->v[i].y = y;
  Vertices->v[i].z = z;
  Vertices->v[i].face[ 0 ] = indface;
  Vertices->v[i].nfaces = 1;
  (void)RbtInsert(&RBT, point, i);
  return i;
} /* AddVertex() */


/* ==================================== */
int32_t AddVertexFixe(double x, double y, double z, int32_t indface)
/* ==================================== */
/* modifie les var. globales Faces, Vertices */
{
  int32_t i;
  i = AddVertex(x, y, z, indface);
  Vertices->lab[i] = 1;
  return i;
} /* AddVertexFixe() */

/* ==================================== */
void AddFace(double x1, double y1, double z1, 
             double x2, double y2, double z2, 
             double x3, double y3, double z3
            )
/* ==================================== */
/* modifie les var. globales Faces, Vertices */
{
  int32_t iv1, iv2, iv3, i, indface;
  if (Faces->cur >= Faces->max) ReAllocFaces(&Faces);
  indface = Faces->cur;
  Faces->cur += 1;
  iv1 = AddVertex(x1, y1, z1, indface);
  iv2 = AddVertex(x2, y2, z2, indface);
  iv3 = AddVertex(x3, y3, z3, indface);
  Faces->f[indface].vert[0] = iv1;
  Faces->f[indface].vert[1] = iv2;
  Faces->f[indface].vert[2] = iv3;
  Faces->f[indface].xn = Faces->f[indface].yn = Faces->f[indface].zn = 0.0;
} /* AddFace() */

/* ==================================== */
void AddFaceFixe(double x1, double y1, double z1, 
                 double x2, double y2, double z2, 
                 double x3, double y3, double z3,
                 int32_t fix1, int32_t fix2, int32_t fix3
                )
/* ==================================== */
/* modifie les var. globales Faces, Vertices */
{
  int32_t iv1, iv2, iv3, i, indface;
  if (Faces->cur >= Faces->max) ReAllocFaces(&Faces);
  indface = Faces->cur;
  Faces->cur += 1;
  if (fix1)
    iv1 = AddVertexFixe(x1, y1, z1, indface);
  else 
    iv1 = AddVertex(x1, y1, z1, indface);
  if (fix2)
    iv2 = AddVertexFixe(x2, y2, z2, indface);
  else 
    iv2 = AddVertex(x2, y2, z2, indface);
  if (fix3)
    iv3 = AddVertexFixe(x3, y3, z3, indface);
  else 
    iv3 = AddVertex(x3, y3, z3, indface);
  Faces->f[indface].vert[0] = iv1;
  Faces->f[indface].vert[1] = iv2;
  Faces->f[indface].vert[2] = iv3;
  Faces->f[indface].xn = Faces->f[indface].yn = Faces->f[indface].zn = 0.0;
} /* AddFaceFixe() */

/* ==================================== */
int32_t AddEdge(int32_t v1, int32_t v2, int32_t f1, int32_t f2)
/* ==================================== */
/* modifie la var. globales Edges */
{
  int32_t indedge;
  if (Edges->cur >= Edges->max) ReAllocEdges(&Edges);
  indedge = Edges->cur;
  Edges->cur += 1;
  Edges->e[indedge].v1 = v1;
  Edges->e[indedge].v2 = v2;
  Edges->e[indedge].f1 = f1;
  Edges->e[indedge].f2 = f2;
  return indedge;
} /* AddEdge() */

/* ==================================== */
void SaveCoords()
/* ==================================== */
{
  int32_t i;
  for (i = 0; i < Vertices->cur; i++)
  {
    Vertices->v[i].xp = Vertices->v[i].x;
    Vertices->v[i].yp = Vertices->v[i].y;
    Vertices->v[i].zp = Vertices->v[i].z;
  }
} /* SaveCoords() */

/* ==================================== */
void SaveOriginalCoords()
/* ==================================== */
{
  int32_t i;
  for (i = 0; i < Vertices->cur; i++)
  {
    Vertices->v[i].xo = Vertices->v[i].x;
    Vertices->v[i].yo = Vertices->v[i].y;
    Vertices->v[i].zo = Vertices->v[i].z;
  }
} /* SaveOriginalCoords() */

/* ==================================== */
void RestoreCoords()
/* ==================================== */
{
  int32_t i;
  for (i = 0; i < Vertices->cur; i++)
  {
    Vertices->v[i].x = Vertices->v[i].xp;
    Vertices->v[i].y = Vertices->v[i].yp;
    Vertices->v[i].z = Vertices->v[i].zp;
  }
} /* RestoreCoords() */

/* ==================================== */
void ComputeEdges()
/* ==================================== */
/*
  Contruit le tableau des cotes (edges),
  et met a jour le champ 'edge' des sommets.
*/
#undef F_NAME
#define F_NAME "ComputeEdges"
{
  int32_t i, j, k, n, e, nvertices;
  meshvertex V;
  meshface F;
  int32_t link[MAXADJFACES];
  int32_t f1, f2;

  if (Edges == NULL)
  {
    fprintf(stderr, "%s : Edges array must be allocated\n", F_NAME);
    exit(0);
  }
  nvertices = Vertices->cur;
  for (i = 0; i < nvertices; i++) Vertices->v[i].nedges = 0;
/*end initialisation*/

  for (i = 0; i < nvertices; i++)
  {
    V = Vertices->v[i];
    n = 0;
	if (V.nfaces > 25 )printf("Too much adjacent faces in vertex = %d  nfaces %d\n", i,  V.nfaces);
    for (j = 0; j < V.nfaces; j++) /* parcourt les faces adjacentes */
    {
	//printf("vertex = %d cares %d nfaces %d\n", i, j, V.nfaces);                                     /* et calcule le link */
      F = Faces->f[V.face[j]];
      k = F.vert[0]; if ((k != i) && NotIn(k, link, n)) link[n++] = k;
      k = F.vert[1]; if ((k != i) && NotIn(k, link, n)) link[n++] = k;
      k = F.vert[2]; if ((k != i) && NotIn(k, link, n)) link[n++] = k;
    } /* for j */

    for (k = 0; k < n; k++)   /* parcourt le link et cree les cotes */
    {
      if (link[k] > i) /* pour ne compter un cote qu'une seule fois */
      {
        f2 = f1 = -1;
        for (j = 0; j < V.nfaces; j++) /* parcourt les faces adjacentes */
        {                    /* et trouve les 2 qui contiennent link[k] */
          F = Faces->f[V.face[j]];
          if ((F.vert[0] == link[k]) || (F.vert[1] == link[k]) || (F.vert[2] == link[k]))
	  {
            if (f1 == -1) f1 = V.face[j]; else { f2 = V.face[j]; break; }
	  }
	
        } /* for j */
        e = AddEdge(i, link[k], f1, f2);
        Vertices->v[i].edge[Vertices->v[i].nedges++] = e;
        Vertices->v[link[k]].edge[Vertices->v[link[k]].nedges++] = e;
      } /* if */
    } /* for k */
  } /* for i */
} /* ComputeEdges() */
