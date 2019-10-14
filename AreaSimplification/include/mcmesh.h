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


/* $Id: mcmesh.h,v 1.6 2006/10/09 12:06:37 michel Exp $ */
#define MAXADJFACES 25

typedef struct {
  double x, y, z;    /* coordonnees */
  double xp, yp, zp; /* coordonnees bis (utilisees aussi pour stocker la normale) */
  double xo, yo, zo; /* coordonnees ter (pour memoriser la position originale) */
                     // utilisees seulement par RegulMeshHC
  int32_t nfaces;
  int32_t face[MAXADJFACES]; /* indices des faces adjacentes (pas plus de MAXADJFACES) */
  int32_t nedges;
  int32_t edge[MAXADJFACES]; /* indices des cotes adjacents (pas plus de MAXADJFACES) */
  double value; 
  float curv1, curv2; // Pour les courbures
  float dcurv[4]; //pour les derivatives des courbures
  int32_t red, green, blue; // pour la couleur
} meshvertex;

typedef struct {
  int32_t vert[3];       /* indices des sommets adjacents */
  double xn, yn, zn; /* normale a la face */
} meshface;

typedef struct {
  int32_t v1, v2;        /* indices des sommets adjacents */
  int32_t f1, f2;        /* indices des faces adjacentes */
  double curv;       /* pour stocker la courbure (angle entre -pi et pi) */
} meshedge;

typedef struct {
  int32_t max;            /* taille max du tableau de sommets */
  int32_t cur;            /* taille courante du tableau de sommets */
  u_int8_t *lab; /* tableau de labels associes aux sommets */
  u_int8_t *tmp; /* tableau de valeurs associes aux sommets */
  meshvertex v[1];    /* tableau des elements physiques */
} meshtabvertices;

typedef struct {
  int32_t max;         /* taille max du tableau de faces */
  int32_t cur;         /* taille courante du tableau de faces */
  meshface f[1];   /* tableau des elements physiques */
} meshtabfaces;

typedef struct {
  int32_t max;         /* taille max du tableau de cotes */
  int32_t cur;         /* taille courante du tableau de cotes */
  meshedge e[1];   /* tableau des elements physiques */
} meshtabedges;

typedef struct {
  int32_t nvert;       /* nombre de sommets */
  int32_t nedge;       /* nombre de cotes */
  int32_t *lastneigh;  /* tableau des index des derniers successeurs */
  int32_t *neigh;      /* tableau des successeurs */
} meshtablinks;

/* boite englobante */
typedef struct {
  double bxmin;
  double bxmax;
  double bymin;
  double bymax;
  double bzmin;
  double bzmax;
} meshbox;

extern meshtabvertices *Vertices;
extern meshtabfaces *Faces;
extern meshtabedges *Edges;
extern meshtablinks *Links;

/* ==================================== */
/* prototypes */
/* ==================================== */

extern meshtabvertices * AllocVertices(int32_t taillemax);
extern meshtabfaces * AllocFaces(int32_t taillemax);
extern meshtabedges * AllocEdges(int32_t taillemax);
/* extern meshtablinks * AllocLinks(int32_t nvert, int32_t nedge); */
extern void ReAllocVertices(meshtabvertices **A);
extern void ReAllocFaces(meshtabfaces **A);
extern void ReAllocEdges(meshtabedges **A);
extern void InitMesh(int32_t taillemax);
extern void TermineMesh();
extern void AddFace(double x1, double y1, double z1, 
             double x2, double y2, double z2, 
             double x3, double y3, double z3
	     );
extern void AddFaceFixe(double x1, double y1, double z1, 
                 double x2, double y2, double z2, 
                 double x3, double y3, double z3,
                 int32_t fix1, int32_t fix2, int32_t fix3
	         );
extern void SaveCoords();
extern void RestoreCoords();
extern void ComputeEdges();



