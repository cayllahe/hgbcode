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


/* $Id: mcgeo.h,v 1.7 2006/02/28 07:49:12 michel Exp $ */
#ifndef EPSILON
#define EPSILON 1E-20
#endif

typedef double vec2[2];
typedef double vec3[3];
typedef double mat22[2][2];
typedef double mat33[3][3];

typedef struct {
  double xa, ya, xb, yb, xc, yc;  /* coordonnees des sommets */
  double a1, a2, a3, b1, b2, b3, c1, c2, c3;  
                                  /* coef. des equations des cotes (ai*x+bi*y+ci = 0) */
  double s1, s2, s3;              /* si = ai*xj+bi*yj+ci avec j sommet oppose au cote i */
  double xmin, ymin, xmax, ymax;  /* coord. du rectangle englobant */
} triangle;

typedef struct {
  double x, y;
} point;

typedef struct vecteur2Dint
	{ int32_t x,y; } vect2Dint;   /* structure representant un vecteur entier */

typedef struct {
  double x, y, z;
} point3;

typedef struct {
  double a, b, c; /* coef. de l'equation (a*x+b*y+c = 0) */
} droite;

typedef struct {
  double x0, y0, r; /* coef. de l'equation (x-x0)2 + (y-y0)2 = r2 */
} cercle;

typedef struct {
  double x1, y1;
  double x2, y2;  
  double a, b, c; /* coef. de l'equation (a*x+b*y+c = 0) de la droite support */
} segment;

typedef struct {
  segment s1, s2, s3, s4;
  point centre;
  double h, w, angle;
} rectangle;

typedef struct {
  segment s1, s2;
  cercle c1, c2;
  point centre;
  double h, w, angle;
} rectanglearrondi;

typedef struct {
  int32_t rs;                     /* taille d'une ligne */
  int32_t cs;                     /* taille d'une colonne */
  int32_t *tab;
} mesh;

#define OBJTYPE_LINE 0
#define OBJTYPE_CLOSEDLINE 1
#define OBJTYPE_SPLINE 2
#define OBJTYPE_CLOSEDSPLINE 3

typedef struct {
  int32_t objtype;
  int32_t npoints;
  point3 *points;
} object;

typedef struct {
  int32_t nobj;
  object **tabobj;
} scene;

/* ============================================= */
/* prototypes pour mcgeo.c */
/* ============================================= */

extern double detpq_qr(point p, point q, point r);
extern int32_t estsitue(point p, point q, point r);
extern int32_t estagauche(point p, point q, point r);
extern int32_t estadroite(point p, point q, point r);
extern double carrenorme(point p);
extern double carredistance(point p, point q);
extern double cosangle(point p, point q, point r);
extern int32_t sontsecants(point p1, point p2, point p3, point p4);
extern double distpointdroite(point p, droite d);
extern double distpointsegment(point p, segment s);
extern double distpointcercle(point p, cercle c);
extern double distpointrectangle(point p, rectangle *r);
extern int32_t initrectangle(rectangle *r);
extern double distpointdroite3(point3 p, point3 p1, point3 p2);

extern void printmat(double * mat, int32_t n, int32_t m);
extern void printvec(double * v, int32_t n);
extern double det2(mat22 m);
extern double det3(mat33 m);
extern int32_t solsyst2(mat22 m, vec2 b, vec2 sol);
extern int32_t solsyst3(mat33 m, vec3 b, vec3 sol);
extern int32_t invmat2(mat22 ma, mat22 mr);
extern int32_t invmat3(mat33 ma, mat33 mr);
extern void multmat2vec2(mat22 m, vec2 v, vec2 sol);
extern void multmat3vec3(mat33 m, vec3 v, vec3 sol);
extern void multmat2mat2(mat22 m1, mat22 m2, mat22 mr);
extern void multmat3mat3(mat33 m1, mat33 m2, mat33 mr);
extern int32_t equdroite(double x1, double y1, double x2, double y2, double * a, double * b, double * c);
extern int32_t inittriangle(triangle *t);
extern int32_t danstriangle(triangle *t, double x, double y);
extern void printtriangle(triangle *t);
extern int32_t dansdisque(double x, double y, double xc, double yc, double r);
extern mesh * readmesh(char *filename);
extern int32_t writemesh(mesh *msh, char *filename);
extern mesh * allocmesh(int32_t rs, int32_t cs);
extern void printmesh(mesh *msh);
extern int32_t checkmesh(mesh *msh, int32_t rs, int32_t cs);
extern void freemesh(mesh *msh);
extern double dist2(double x1, double y1, double x2, double y2);
extern double dist3(double x1, double y1, double z1, double x2, double y2, double z2);
extern int32_t identcercle(int32_t n, double *x, double *y, double *x0, double *y0, double *r);
extern scene * readscene(char *filename);
extern int32_t writescene(scene *scn, char *filename);
extern scene * copyscene(scene * s);

#define _mcgeo_h
