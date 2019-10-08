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

/* $Id: mciomesh.c,v 1.9 2006/02/28 07:49:16 michel Exp $ */
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <mcmesh.h>
#include <mcrbt1.h>
#include <mciomesh.h>

//#define VERBOSE

extern Rbt * RBT;

u_int16_t LE_ReadUnsignedShort(FILE *fd)
{
  u_int16_t tmp, tmp1;
  fread(&tmp, 2 * sizeof(char), 1, fd);
  tmp1 = tmp & 0x00ff;
  tmp = tmp >> 8;
  tmp = tmp | (tmp1 << 8);      
  return tmp;
}

u_int32_t LE_ReadUnsignedLong(FILE *fd)
{
  u_int32_t tmp;
  char tmp1[4];
  fread(tmp1, 4*sizeof(char), 1, fd);
  tmp = (u_int32_t)tmp1[0];
  tmp |= ((u_int32_t)tmp1[1]) << 8;
  tmp |= ((u_int32_t)tmp1[2]) << 16;
  tmp |= ((u_int32_t)tmp1[3]) << 24;
  return tmp;
}

u_int32_t ReadUnsignedLong(FILE *fd)
{
  u_int32_t tmp;
  fread((char *)(&tmp), sizeof(int32_t), 1, fd);
  return tmp;
}

float ReadFloat32(FILE *fd)
{
  float tmp;
  fread((char *)(&tmp), sizeof(float), 1, fd);
  return tmp;
}

/* =============================================================== */
void genheaderPOV(FILE *fileout, int32_t obj_id, meshbox MB)
/* =============================================================== */
{
  fprintf(fileout, "#include \"colors.inc\"\n");
  fprintf(fileout, "#include \"shapes.inc\"\n");
  fprintf(fileout, "#include \"textures.inc\"\n");
  fprintf(fileout, "#include \"stones.inc\"\n");

  fprintf(fileout, "\n");
  fprintf(fileout, "camera { location <0,0,%g> up <0,1,0> right <1,0,0> look_at <0,0,0> }\n",
                   2.7*MB.bzmax);
  fprintf(fileout, "\n");
  fprintf(fileout, "light_source { <%g,%g,%g> color White }\n", 
                   2*MB.bxmin, 2*MB.bymax, 2*MB.bzmax);
  fprintf(fileout, "light_source { <%g,%g,%g> color Gray80 }\n", 
                   2*MB.bxmax, 2*MB.bymax, 2*MB.bzmax);
  fprintf(fileout, "light_source { <%g,%g,%g> color Gray80 }\n", 
                   2*MB.bxmin, 2*MB.bymin, 2*MB.bzmax);
  fprintf(fileout, "light_source { <%g,%g,%g> color Gray80 }\n", 
                   2*MB.bxmax, 2*MB.bymin, 2*MB.bzmax);
  fprintf(fileout, "light_source { <%g,%g,%g> color White }\n", 
                   2*MB.bxmax, 2*MB.bymax, 2*MB.bzmin);
  fprintf(fileout, "light_source { <%g,%g,%g> color White }\n", 
                   0, 0, 3*MB.bzmax);
  fprintf(fileout, "\n");
  fprintf(fileout, "#declare mytexture = Bright_Bronze\n");
  fprintf(fileout, "\n");
  fprintf(fileout, "mesh { // obj_id %d\n", obj_id);
} /* genheaderPOV() */

/* =============================================================== */
void genheaderAC(FILE *fileout, meshbox MB, double red, double green, double blue, char *name)
/* =============================================================== */
{
    fprintf(fileout, "AC3Db\n");
    fprintf(fileout, "MATERIAL \"ac3dmat0\" rgb %g %g %g  amb 0.2 0.2 0.2  emis 0 0 0  spec 0.5 0.5 0.5  shi 10  trans 0\n", red, green, blue);
    fprintf(fileout, "OBJECT world\n");
    fprintf(fileout, "kids 1\n");
    fprintf(fileout, "OBJECT poly\n");
    fprintf(fileout, "name \"%s\"\n", name);
    fprintf(fileout, "loc %g %g %g\n", (MB.bxmax - MB.bxmin) / 2, (MB.bymax - MB.bymin) / 2, 
                     (MB.bzmax - MB.bzmin) / 2);
} /* genheaderAC() */

/* =============================================================== */
void genheaderVTK(FILE *fileout, char *name)
/* =============================================================== */
{
    fprintf(fileout, "# vtk DataFile Version 3.0\n");
    fprintf(fileout, "%s\n", name);
    fprintf(fileout, "ASCII\n");
    fprintf(fileout, "DATASET POLYDATA\n");
} /* genheaderVTK() */

/* =============================================================== */
void genheaderCOL(FILE *fileout, int32_t nbfaces)
/* =============================================================== */
{
  fprintf(fileout, "%d\n", nbfaces);
} /* genheaderCOL() */

/* =============================================================== */
void genheaderDXF(FILE *fileout)
/* =============================================================== */
{
    fprintf(fileout, "0\nSECTION\n2\nENTITIES\n");
} /* genheaderDXF() */

/* =============================================================== */
void genfooterPOV(FILE *fileout)
/* =============================================================== */
{
  fprintf(fileout, " texture{ mytexture }\n");
  fprintf(fileout, " rotate <0,0,0>\n");
  fprintf(fileout, "}\n");
} /* genfooterPOV() */

/* =============================================================== */
void genfooterAC(FILE *fileout)
/* =============================================================== */
{
} /* genfooterAC() */

/* =============================================================== */
void genfooterDXF(FILE *fileout)
/* =============================================================== */
{
  fprintf(fileout, "0\nENDSEC\n0\nEOF\n");
} /* genfooterDXF() */

/* =============================================================== */
void genfacePOV(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
             double x3, double y3, double z3
            )
/* =============================================================== */
{
    fprintf(fileout, "triangle{");
    fprintf(fileout, "<%g,%g,%g>,", x1, y1, z1);
    fprintf(fileout, "<%g,%g,%g>,", x2, y2, z2);
    fprintf(fileout, "<%g,%g,%g>", x3, y3, z3);
    fprintf(fileout, "}\n");
} /* genfacePOV() */

/* =============================================================== */
void genfaceSPOV(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
	     double x3, double y3, double z3,
             double nx1, double ny1, double nz1, 
             double nx2, double ny2, double nz2, 
             double nx3, double ny3, double nz3
            )
/* =============================================================== */
{
    fprintf(fileout, "smooth_triangle{");
    fprintf(fileout, "<%g,%g,%g>,", x1, y1, z1);
    fprintf(fileout, "<%g,%g,%g>,", nx1, ny1, nz1);
    fprintf(fileout, "<%g,%g,%g>,", x2, y2, z2);
    fprintf(fileout, "<%g,%g,%g>,", nx2, ny2, nz2);
    fprintf(fileout, "<%g,%g,%g>,", x3, y3, z3);
    fprintf(fileout, "<%g,%g,%g>" , nx3, ny3, nz3);
    fprintf(fileout, "}\n");
} /* genfaceSPOV() */

/* =============================================================== */
void gencarrePOV(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
             double x3, double y3, double z3, 
             double x4, double y4, double z4
            )
/* =============================================================== */
{
    fprintf(fileout, "polygon { 5, ");
    fprintf(fileout, "<%g,%g,%g>, ", x1, y1, z1);
    fprintf(fileout, "<%g,%g,%g>, ", x2, y2, z2);
    fprintf(fileout, "<%g,%g,%g>, ", x3, y3, z3);
    fprintf(fileout, "<%g,%g,%g>, ", x4, y4, z4);
    fprintf(fileout, "<%g,%g,%g> ", x1, y1, z1);
    fprintf(fileout, "}\n");
} /* gencarrePOV() */

/* =============================================================== */
void genfaceCOL(FILE *fileout, int32_t face_id, int32_t object_id,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
	     double x3, double y3, double z3,
             double xp1, double yp1, double zp1, 
             double xp2, double yp2, double zp2, 
             double xp3, double yp3, double zp3
            )
/* =============================================================== */
{
    fprintf(fileout, "3DFACE ");
    fprintf(fileout, "%d %d\n", face_id, object_id);
    fprintf(fileout, "%g, %g, %g, %g, %g, %g, \n", x1, y1, z1, xp1, yp1, zp1);
    fprintf(fileout, "%g, %g, %g, %g, %g, %g, \n", x2, y2, z2, xp2, yp2, zp2);
    fprintf(fileout, "%g, %g, %g, %g, %g, %g\n", x3, y3, z3, xp3, yp3, zp3);
} /* genfaceCOL() */

/* =============================================================== */
void genfaceDXF(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
	     double x3, double y3, double z3,
             double xp1, double yp1, double zp1, 
             double xp2, double yp2, double zp2, 
             double xp3, double yp3, double zp3
            )
/* =============================================================== */
{
#define ALPHA 0.05
    fprintf(fileout, "0\n3DFACE\n"); // en DXF une 3DFACE est definie par quatre points
    fprintf(fileout, "10\n%g\n20\n%g\n30\n%g\n", x1, y1, z1);
    fprintf(fileout, "11\n%g\n21\n%g\n31\n%g\n", x2, y2, z2);
    fprintf(fileout, "12\n%g\n22\n%g\n32\n%g\n", x3, y3, z3);
    fprintf(fileout, "13\n%g\n23\n%g\n33\n%g\n", x3, y3, z3); // le 3eme est duplique
} /* genfaceDXF() */

/* ==================================== */
void SaveMeshPOV(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
{
  int32_t i;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3 ;
  meshvertex V;

  for (i = 0; i < Faces->cur; i++)
  {
    V = Vertices->v[ Faces->f[i].vert[0] ];
    x1 = V.x; y1 = V.y; z1 = V.z;
    V = Vertices->v[ Faces->f[i].vert[1] ];
    x2 = V.x; y2 = V.y; z2 = V.z;
    V = Vertices->v[ Faces->f[i].vert[2] ];
    x3 = V.x; y3 = V.y; z3 = V.z;
    genfacePOV(fileout, x1, y1, z1, x2, y2, z2, x3, y3, z3);
  }
} /* SaveMeshPOV() */

/* ==================================== */
void SaveMeshSPOV(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
{
  int32_t i;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3 ;
  double nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3 ;
  meshvertex V;

  for (i = 0; i < Faces->cur; i++)
  {
    V = Vertices->v[ Faces->f[i].vert[0] ];
    x1 = V.x; y1 = V.y; z1 = V.z;
    nx1 = V.xp; ny1 = V.yp; nz1 = V.zp;
    V = Vertices->v[ Faces->f[i].vert[1] ];
    x2 = V.x; y2 = V.y; z2 = V.z;
    nx2 = V.xp; ny2 = V.yp; nz2 = V.zp;
    V = Vertices->v[ Faces->f[i].vert[2] ];
    x3 = V.x; y3 = V.y; z3 = V.z;
    nx3 = V.xp; ny3 = V.yp; nz3 = V.zp;
    genfaceSPOV(fileout, x1, y1, z1, x2, y2, z2, x3, y3, z3,
                nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3);
  }
} /* SaveMeshSPOV() */

/* ==================================== */
void SaveMeshCOL(FILE *fileout, int32_t obj_id)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
{
  int32_t i;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3 ;
  double xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3 ;
  meshvertex V;

  for (i = 0; i < Faces->cur; i++)
  {
    V = Vertices->v[ Faces->f[i].vert[0] ];
    x1 = V.x; y1 = V.y; z1 = V.z;
    xp1 = V.xp; yp1 = V.yp; zp1 = V.zp;
    V = Vertices->v[ Faces->f[i].vert[1] ];
    x2 = V.x; y2 = V.y; z2 = V.z;
    xp2 = V.xp; yp2 = V.yp; zp2 = V.zp;
    V = Vertices->v[ Faces->f[i].vert[2] ];
    x3 = V.x; y3 = V.y; z3 = V.z;
    xp3 = V.xp; yp3 = V.yp; zp3 = V.zp;
    genfaceCOL(fileout, i+1, obj_id, x1, y1, z1, x2, y2, z2, x3, y3, z3, 
                        xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3);
  }

} /* SaveMeshCOL() */

/* ==================================== */
void SaveMeshMCM(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
   V %d     (Sommets - champ obligatoire)
   %g %g %g (coord. vertex)
   ...

   v        (normales au sommets - champ optionnel)
   %g %g %g (vect. normal au vertex)
   ...

   f        (faces adjacentes au sommets - champ optionnel)
   %d (nb. faces adj) %d %d ... %d (ind. faces)
   ...

   x %d     (indices des sommets fixes - champ optionnel)
   %d (ind. vertex)
   ...

   F %d     (Faces - champ obligatoire)
   %d %d %d (face: ind. vertices)
   ...

   n        (normales au faces - champ optionnel)
   %g %g %g (vect. normal a la face)
   ...

   e %d     (edges - champ optionnel)
   %d %d (indices sommets) %d %d (indices faces, -1 indique l'absence)

*/
{
  int32_t i, j, nfixes = 0;

  // SOMMETS
  fprintf(fileout, "V %d\n", Vertices->cur);
  for (i = 0; i < Vertices->cur; i++)
  {
    if (Vertices->lab[i]) nfixes++;
    fprintf(fileout, "%g %g %g", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // NORMALES AUX SOMMETS
  fprintf(fileout, "v\n");
  for (i = 0; i < Vertices->cur; i++)
  {
    fprintf(fileout, "%g %g %g", Vertices->v[i].xp, Vertices->v[i].yp, Vertices->v[i].zp);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // FACES ADJACENTES AUX SOMMETS
  fprintf(fileout, "f\n");
  for (i = 0; i < Vertices->cur; i++)
  {
    fprintf(fileout, "%d ", Vertices->v[i].nfaces);
    for (j = 0; j < Vertices->v[i].nfaces; j++)
      fprintf(fileout, "%d ", Vertices->v[i].face[j]);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // SOMMETS FIXES
  fprintf(fileout, "x %d\n", nfixes);
  for (i = 0; i < Vertices->cur; i++)
    if (Vertices->lab[i]) fprintf(fileout, "%d\n", i);
  fprintf(fileout, "\n");

  // FACES
  fprintf(fileout, "F %d\n", Faces->cur);
  for (i = 0; i < Faces->cur; i++)
    fprintf(fileout, "%d %d %d\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  fprintf(fileout, "\n");

  // NORMALES AUX FACES
  fprintf(fileout, "n\n");
  for (i = 0; i < Faces->cur; i++)
  {
    fprintf(fileout, "%g %g %g", Faces->f[i].xn, Faces->f[i].yn, Faces->f[i].zn);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // EDGES
  if (Edges != NULL) 
  {
    fprintf(fileout, "e %d\n", Edges->cur);
    for (i = 0; i < Edges->cur; i++)
      fprintf(fileout, "%d %d %d %d\n", Edges->e[i].v1, Edges->e[i].v2, Edges->e[i].f1, Edges->e[i].f2);

  }


} /* SaveMeshMCM() */

/* ==================================== */
void SaveMeshPLY(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
ply
format ascii 1.0
element vertex 8
property float x
property float y
property float z
property red uchar
property green uchar
property blue uchar
element face 7
property list uchar int vertex_index
end_header
0 0 0 255 0 0
0 0 1 255 0 0
0 1 1 255 0 0
0 1 0 255 0 0
1 0 0 0 0 255
1 0 1 0 0 255
1 1 1 0 0 255
1 1 0 0 0 255
3 0 1 2
3 0 2 3
3 7 6 5
3 0 4 5
3 1 5 6
3 2 6 7
3 3 7 4

*/
{
  int32_t i = 0;


   fprintf(fileout, "ply\n");
   fprintf(fileout, "format ascii 1.0\n");
   fprintf(fileout, "element vertex %d\n", Vertices->cur);
   fprintf(fileout, "property float x\n");
   fprintf(fileout, "property float y\n");
   fprintf(fileout, "property float z\n");
   fprintf(fileout, "property float curv1\n");
   fprintf(fileout, "property float curv2\n");
   fprintf(fileout, "property float dcurv1\n");
   fprintf(fileout, "property float dcurv2\n");
   fprintf(fileout, "property float dcurv3\n");
   fprintf(fileout, "property float dcurv4\n");
   fprintf(fileout, "property uchar red\n");
   fprintf(fileout, "property uchar green\n");
   fprintf(fileout, "property uchar blue\n");
   fprintf(fileout, "element face %d\n", Faces->cur);
   fprintf(fileout, "property list uchar int vertex_indices\n");
   fprintf(fileout, "end_header\n");

   for (i = 0; i < Vertices->cur; i++)
  {
    
    fprintf(fileout, "%g %g %g %g %g %g %g %g %g %d %d %d", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z, Vertices->v[i].curv1, Vertices->v[i].curv2, Vertices->v[i].dcurv[0], Vertices->v[i].dcurv[1], Vertices->v[i].dcurv[2], Vertices->v[i].dcurv[3], Vertices->v[i].red, Vertices->v[i].green, Vertices->v[i].blue );
    fprintf(fileout, "\n");
  }
  

  for (i = 0; i < Faces->cur; i++)
    fprintf(fileout, "3 %d %d %d\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  fprintf(fileout, "\n");


  

} /* SaveMeshPLY() */

/* ==================================== */
void SaveMeshPLYsimple(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
ply
format ascii 1.0
element vertex 8
property float x
property float y
property float z
property red uchar
property green uchar
property blue uchar
element face 7
property list uchar int vertex_index
end_header
0 0 0 255 0 0
0 0 1 255 0 0
0 1 1 255 0 0
0 1 0 255 0 0
1 0 0 0 0 255
1 0 1 0 0 255
1 1 1 0 0 255
1 1 0 0 0 255
3 0 1 2
3 0 2 3
3 7 6 5
3 0 4 5
3 1 5 6
3 2 6 7
3 3 7 4

*/
{
  int32_t i = 0;


   fprintf(fileout, "ply\n");
   fprintf(fileout, "format ascii 1.0\n");
   fprintf(fileout, "element vertex %d\n", Vertices->cur);
   fprintf(fileout, "property float x\n");
   fprintf(fileout, "property float y\n");
   fprintf(fileout, "property float z\n");
   fprintf(fileout, "property uchar red\n");
   fprintf(fileout, "property uchar green\n");
   fprintf(fileout, "property uchar blue\n");
   fprintf(fileout, "element face %d\n", Faces->cur);
   fprintf(fileout, "property list uchar int vertex_indices\n");
   fprintf(fileout, "end_header\n");

   for (i = 0; i < Vertices->cur; i++)
  {
    
    fprintf(fileout, "%g %g %g %d %d %d", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z, Vertices->v[i].red, Vertices->v[i].green, Vertices->v[i].blue );
    fprintf(fileout, "\n");
  }
  

  for (i = 0; i < Faces->cur; i++)
    fprintf(fileout, "3 %d %d %d\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  fprintf(fileout, "\n");


  

} /* SaveMeshPLYsimple() */

/* ==================================== */
void SaveMeshPLYvalue(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
ply
format ascii 1.0
element vertex 8
property float x
property float y
property float z
property red uchar
property green uchar
property blue uchar
element face 7
property list uchar int vertex_index
end_header
0 0 0 255 0 0
0 0 1 255 0 0
0 1 1 255 0 0
0 1 0 255 0 0
1 0 0 0 0 255
1 0 1 0 0 255
1 1 1 0 0 255
1 1 0 0 0 255
3 0 1 2
3 0 2 3
3 7 6 5
3 0 4 5
3 1 5 6
3 2 6 7
3 3 7 4

*/
{
  int32_t i = 0;


   fprintf(fileout, "ply\n");
   fprintf(fileout, "format ascii 1.0\n");
   fprintf(fileout, "element vertex %d\n", Vertices->cur);
   fprintf(fileout, "property float x\n");
   fprintf(fileout, "property float y\n");
   fprintf(fileout, "property float z\n");
   fprintf(fileout, "property float value\n");	
   fprintf(fileout, "property uchar red\n");
   fprintf(fileout, "property uchar green\n");
   fprintf(fileout, "property uchar blue\n");
   fprintf(fileout, "element face %d\n", Faces->cur);
   fprintf(fileout, "property list uchar int vertex_indices\n");
   fprintf(fileout, "end_header\n");

   for (i = 0; i < Vertices->cur; i++)
  {
    
    fprintf(fileout, "%g %g %g %g %d %d %d", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z,Vertices->v[i].value, Vertices->v[i].red, Vertices->v[i].green, Vertices->v[i].blue );
    fprintf(fileout, "\n");
  }
  

  for (i = 0; i < Faces->cur; i++)
    fprintf(fileout, "3 %d %d %d\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  fprintf(fileout, "\n");


  

} /* SaveMeshPLYvalue() */

/* ==================================== */
void SaveMeshPLYlabels(FILE *fileout, int32_t *FaceLabels)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
ply
format ascii 1.0
element vertex 8
property float x
property float y
property float z
property red uchar
property green uchar
property blue uchar
element face 7
property list uchar int vertex_index
end_header
0 0 0 255 0 0
0 0 1 255 0 0
0 1 1 255 0 0
0 1 0 255 0 0
1 0 0 0 0 255
1 0 1 0 0 255
1 1 1 0 0 255
1 1 0 0 0 255
3 0 1 2
3 0 2 3
3 7 6 5
3 0 4 5
3 1 5 6
3 2 6 7
3 3 7 4

*/
{
  int32_t i = 0;


   fprintf(fileout, "ply\n");
   fprintf(fileout, "format ascii 1.0\n");
   fprintf(fileout, "element vertex %d\n", Vertices->cur);
   fprintf(fileout, "property float x\n");
   fprintf(fileout, "property float y\n");
   fprintf(fileout, "property float z\n");
   fprintf(fileout, "property float curv1\n");
   fprintf(fileout, "property float curv2\n");
   fprintf(fileout, "property float dcurv1\n");
   fprintf(fileout, "property float dcurv2\n");
   fprintf(fileout, "property float dcurv3\n");
   fprintf(fileout, "property float dcurv4\n");
   fprintf(fileout, "property uchar red\n");
   fprintf(fileout, "property uchar green\n");
   fprintf(fileout, "property uchar blue\n");
   fprintf(fileout, "element face %d\n", Faces->cur);
   fprintf(fileout, "property list uchar int vertex_indices\n");
   fprintf(fileout, "end_header\n");

   for (i = 0; i < Vertices->cur; i++)
  {
    
    fprintf(fileout, "%g %g %g %g %g %g %g %g %g %d %d %d", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z, Vertices->v[i].curv1, Vertices->v[i].curv2, Vertices->v[i].dcurv[0], Vertices->v[i].dcurv[1], Vertices->v[i].dcurv[2], Vertices->v[i].dcurv[3], Vertices->v[i].red, Vertices->v[i].green, Vertices->v[i].blue );
    fprintf(fileout, "\n");
  }
  

  for (i = 0; i < Faces->cur; i++)
    fprintf(fileout, "3 %d %d %d %d\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2], FaceLabels[i]);
  fprintf(fileout, "\n");


  


} /* SaveMeshPLYsimple() */



/* ==================================== */
void SaveMeshYAMS(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format:

MeshVersionFormatted	1
Dimension	3

#set of vertices
Vertices
num of vertices
-24.98689 -66.48193 -5.470881 0
-0.554481 -66.45267 -10.20065 0
-7.827967 -64.71899 -25.157 0
......
.....
..
....
Triangles
num of triangles
71545 71546 71547 0
12979 71547 71499 0
13013 71548 71547 0
  
*/

{
  int32_t i = 0;


   fprintf(fileout, "MeshVersionFormatted  1\n");
   fprintf(fileout, "Dimension  3\n");
   fprintf(fileout, "\n");
   fprintf(fileout, "Vertices\n");
   fprintf(fileout, "%d\n", Vertices->cur);
   for (i = 0; i < Vertices->cur; i++)
  {
    
    fprintf(fileout, "%g %g %g 0", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
    fprintf(fileout, "\n");
  }
   fprintf(fileout, "\n");
   fprintf(fileout, "Triangles\n");
   fprintf(fileout, "%d\n", Faces->cur);
   for (i = 0; i < Faces->cur; i++){
    fprintf(fileout, "%d %d %d 0\n", Faces->f[i].vert[0] + 1, Faces->f[i].vert[1] + 1, Faces->f[i].vert[2] + 1);
  fprintf(fileout, "\n");
   }
fprintf(fileout, "\n");
fprintf(fileout, "End\n");

 
} /* SaveMeshYAMS() */
 
/* ==================================== */
void SaveMeshVTK(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
POINTS %d float   // Sommets - champ obligatoire
%g %g %g          // coord. vertex
   ...
POLYGONS %d %d    // Faces - champ obligatoire 
                  // arg1: nb polygones; arg2: nb valeurs (=4*arg1 pour des triangles)
3 %d %d %d        // face: ind. vertices
   ...
*/
{
  int32_t i, j;

  // SOMMETS
  fprintf(fileout, "POINTS %d float\n", Vertices->cur);
  for (i = 0; i < Vertices->cur; i++)
  {
    fprintf(fileout, "%g %g %g", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // FACES
  fprintf(fileout, "POLYGONS %d %d\n", Faces->cur, 4*Faces->cur);
  for (i = 0; i < Faces->cur; i++)
    fprintf(fileout, "3 %d %d %d\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  fprintf(fileout, "\n");

} /* SaveMeshVTK() */

 
/* ==================================== */
void SaveCurveVTK(FILE *fileout, u_int8_t* markEdge, int32_t nbMarked) 
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
POINTS %d float   // Sommets - champ obligatoire
%g %g %g          // coord. vertex
   ...
POLYGONS %d %d    // Faces - champ obligatoire 
                  // arg1: nb polygones; arg2: nb valeurs (=4*arg1 pour des triangles)
3 %d %d %d        // face: ind. vertices
   ...
*/
{
  int32_t i, j;

  // SOMMETS
  fprintf(fileout, "POINTS %d float\n", Vertices->cur);
  for (i = 0; i < Vertices->cur; i++)
  {
    fprintf(fileout, "%g %g %g", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // LINES
  fprintf(fileout, "LINES %d %d\n", nbMarked, 3*nbMarked);
  //fprintf(fileout, "LINES %d %d\n", Edges->cur, 3*Edges->cur);
  for (i = 0; i < Edges->cur; i++)
    if(markEdge[i]){
      fprintf(fileout, "2 %d %d\n", Edges->e[i].v1, Edges->e[i].v2);
    }
  fprintf(fileout, "\n");

} /* SaveLineVTK() */


/* ==================================== */
void SaveWeightedCurveVTK(FILE *fileout, u_int8_t* markEdge, int32_t nbMarked, double *Weight) 
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
POINTS %d float   // Sommets - champ obligatoire
%g %g %g          // coord. vertex
   ...
POLYGONS %d %d    // Faces - champ obligatoire 
                  // arg1: nb polygones; arg2: nb valeurs (=4*arg1 pour des triangles)
3 %d %d %d        // face: ind. vertices
   ...
*/
{
  int32_t i, j;

  // SOMMETS
  fprintf(fileout, "POINTS %d float\n", Vertices->cur);
  for (i = 0; i < Vertices->cur; i++)
  {
    fprintf(fileout, "%g %g %g", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");

  // LINES
  fprintf(fileout, "LINES %d %d\n", nbMarked, 3*nbMarked);
  //fprintf(fileout, "LINES %d %d\n", Edges->cur, 3*Edges->cur);
  for (i = 0; i < Edges->cur; i++)
    if(markEdge[i]){
      fprintf(fileout, "2 %d %d\n", Edges->e[i].v1, Edges->e[i].v2);
    }
  fprintf(fileout, "\n");

  // We now save the weights
  
  fprintf(fileout,"CELL_DATA %d\n",nbMarked);
  fprintf(fileout,"SCALARS cell_scalars float 1\n");
  fprintf(fileout,"LOOKUP_TABLE default\n");

  for (i = 0; i < Edges->cur; i++)
    if(markEdge[i]){
      fprintf(fileout, "%lf\n", Weight[i]);
    }
  fprintf(fileout, "\n");

} /* SaveLineVTK() */

/* ==================================== */
void SaveMeshAC(FILE *fileout)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
{
  int32_t i, j;

  fprintf(fileout, "numvert %d\n", Vertices->cur);
  for (i = 0; i < Vertices->cur; i++)
  {
    fprintf(fileout, "%g %g %g ", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
    fprintf(fileout, "\n");
  }
  fprintf(fileout, "\n");
  fprintf(fileout, "numsurf %d\n", Faces->cur);
  for (i = 0; i < Faces->cur; i++)
  {
    fprintf(fileout, "SURF 0x30\n");
    fprintf(fileout, "mat 0\n");
    fprintf(fileout, "refs 3\n");
    fprintf(fileout, "%d 0 0\n%d 0 1\n%d 1 0\n", Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  }
} /* SaveMeshAC() */

/* ==================================== */
void SaveMeshGL(FILE *fileout, int32_t obj_id)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
/* format: 
static GLfloat normals<obj_id>[<nvertices>][3] = {
  {%g, %g, %g},
  {%g, %g, %g},
   ...
};
static GLfloat vertices<obj_id>[<nvertices>][3] = {
  {%g, %g, %g},
  {%g, %g, %g},
   ...
};
static GLfloat face_indiciess<obj_id>[<nfaces>][6] = {
  {%d, %d, %d, %d, %d, %d}, // indices points 1,2,3 normales 1,2,3
  {%d, %d, %d, %d, %d, %d},
   ...
};
*/
{
  int32_t i, j;

  fprintf(fileout, "static GLfloat normals%d[%d][3] = {\n", obj_id, Vertices->cur);
  for (i = 0; i < Vertices->cur-1; i++)
    fprintf(fileout, "  {%g, %g, %g},\n", Vertices->v[i].xp, Vertices->v[i].yp, Vertices->v[i].zp);
  fprintf(fileout, "  {%g, %g, %g}\n};\n\n", Vertices->v[i].xp, Vertices->v[i].yp, Vertices->v[i].zp);

  fprintf(fileout, "static GLfloat vertices%d[%d][3] = {\n", obj_id, Vertices->cur);
  for (i = 0; i < Vertices->cur-1; i++)
    fprintf(fileout, "  {%g, %g, %g},\n", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);
  fprintf(fileout, "  {%g, %g, %g}\n};\n\n", Vertices->v[i].x, Vertices->v[i].y, Vertices->v[i].z);

  fprintf(fileout, "static GLint face_indicies%d[%d][6] = {\n", obj_id, Faces->cur);
  for (i = 0; i < Faces->cur-1; i++)
    fprintf(fileout, "  {%d, %d, %d, %d, %d, %d},\n", 
                     Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2], 
                     Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
  fprintf(fileout, "  {%d, %d, %d, %d, %d, %d}\n};\n\n", 
                   Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2], 
                   Faces->f[i].vert[0], Faces->f[i].vert[1], Faces->f[i].vert[2]);
} /* SaveMeshGL() */

/* ==================================== */
void SaveMeshDXF(FILE *fileout, int32_t obj_id)
/* ==================================== */
/* fileout doit avoir ete ouvert en ecriture */
{
  int32_t i;
  double x1, y1, z1, x2, y2, z2, x3, y3, z3 ;
  double xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3 ;
  meshvertex V;

  for (i = 0; i < Faces->cur; i++)
  {
    V = Vertices->v[ Faces->f[i].vert[0] ];
    x1 = V.x; y1 = V.y; z1 = V.z;
    xp1 = V.xp; yp1 = V.yp; zp1 = V.zp;
    V = Vertices->v[ Faces->f[i].vert[1] ];
    x2 = V.x; y2 = V.y; z2 = V.z;
    xp2 = V.xp; yp2 = V.yp; zp2 = V.zp;
    V = Vertices->v[ Faces->f[i].vert[2] ];
    x3 = V.x; y3 = V.y; z3 = V.z;
    xp3 = V.xp; yp3 = V.yp; zp3 = V.zp;
    genfaceDXF(fileout, x1, y1, z1, x2, y2, z2, x3, y3, z3, 
                        xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3);
  }

} /* SaveMeshDXF() */

/* ==================================== */
void LoadMeshMCM(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/* format: 
   V %d     (Sommets - champ obligatoire)
   %g %g %g (coord. vertex)
   ...

   v        (normales au sommets - champ optionnel)
   %g %g %g (vect. normal au vertex)
   ...

   f        (faces adjacentes au sommets - champ optionnel)
   %d (nb. faces adj) %d %d ... %d (ind. faces)
   ...

   x %d     (indices des sommets fixes - champ optionnel)
   %d (ind. vertex)
   ...

   F %d     (Faces - champ obligatoire)
   %d %d %d (face: ind. vertices)
   ...

   n        (normales au faces - champ optionnel)
   %g %g %g (vect. normal a la face)
   ...
  
*/
{
  int32_t i, j, n, nvert=-1, nfaces=-1, nvertfix;
  double x, y, z;
#define sbuf 1024
  char buf[sbuf];
  char *ret;

  while (1)
  {
    do { ret=fgets(buf, sbuf, filein); }
    while (ret && (buf[0]!='V') && (buf[0]!='v') && (buf[0]!='f') && 
                  (buf[0]!='x') && (buf[0]!='F') && (buf[0]!='n'));
    if (ret == NULL) goto end;;
    if (buf[0]=='V')
    {
      sscanf(buf+1, "%d", &nvert);
      Vertices = AllocVertices(nvert);
      Vertices->cur = nvert;
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf", &x, &y, &z);
        Vertices->v[i].x = x; Vertices->v[i].y = y; Vertices->v[i].z = z;
      } // for i
    }
    else if (buf[0]=='F')
    {
      sscanf(buf+1, "%d", &nfaces);
      Faces = AllocFaces(nfaces);
      Faces->cur = nfaces;
      for (i = 0; i < nfaces; i++)
        fscanf(filein, "%d%d%d", &(Faces->f[i].vert[0]), 
               &(Faces->f[i].vert[1]), &(Faces->f[i].vert[2]));
    }
    else if (buf[0]=='x')
    {
      sscanf(buf+1, "%d", &nvertfix); // les vertex fixes
      for (i = 0; i < nvertfix; i++)
      {
        fscanf(filein, "%d", &j);
        Vertices->lab[j] = 1;
      }
    }
    else if (buf[0]=='v')
    {
      if (nvert == -1)
      {
        fprintf(stderr, "LoadMeshMCM: bad file format\n");
        exit(0);
      }
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf", &x, &y, &z);
        Vertices->v[i].xp = x; Vertices->v[i].yp = y; Vertices->v[i].zp = z;
      } // for i
    }
    else if (buf[0]=='f')
    {
      if (nvert == -1)
      {
        fprintf(stderr, "LoadMeshMCM: bad file format\n");
        exit(0);
      }
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%d", &n);
        Vertices->v[i].nfaces = n;
        for (j = 0; j < n; j++)
          fscanf(filein, "%d", &(Vertices->v[i].face[j]));
      } // for i
    }
    else if (buf[0]=='n')
    {
      if (nfaces == -1)
      {
        fprintf(stderr, "LoadMeshMCM: bad file format\n");
        exit(0);
      }
      for (i = 0; i < nfaces; i++)
      {
        fscanf(filein, "%lf%lf%lf", &x, &y, &z);
        Faces->f[i].xn = x; Faces->f[i].yn = y; Faces->f[i].zn = z;
      } // for i
    }
  } // while (1)

 end:
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshMCM: bad file format\n");
    exit(0);
  }
} /* LoadMeshMCM() */

/* ==================================== */
void LoadMeshPLY(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/* format: 
ply
format ascii 1.0
element vertex 141362 ---------------------->number of vertices
property float x
property float y
property float z
property curv1
property curv2
property dcurv1
property dcurv2
property dcurv3
property dcurv4
property uchar diffuse_red
property uchar diffuse_green
property uchar diffuse_blue
element face 282712 -------------------->number of faces
property list uchar int vertex_indices
end_header
-24.98689 -66.48193 -5.470881 255 197 3
-0.554481 -66.45267 -10.20065 255 253 208
-7.827967 -64.71899 -25.157
......
.....
..
....
3 71545 71546 71547
3 12979 71547 71499
3 13013 71548 71547
3 12980 71499 71548
3 71547 71548 71499
3 12980 71548 71500
3 13013 71549 71548
3 12995 71500 71549
3 71548 71549 71500
3 12996 71550 71506
3 12997 71551 71550  ...
  
*/

{

 int32_t i, j, t, n, nf, nvert=-1, nfaces=-1, f1, f2, f3;
  int32_t r, g, b;
  double x, y, z;
  float c1,c2, d1, d2, d3, d4;
#define sbuf 1024
  char buf[sbuf];
  char *ret;

  while (1)
  {
    //printf("LoadPLY: loading \n");
    do { ret=fgets(buf, sbuf, filein); }
    while (ret && strncmp(buf,"ply",3) && strncmp(buf,"element vertex",14)&& strncmp(buf,"property float",14) && strncmp(buf,"element face",12) && strncmp(buf,"end_header",10));
    if (ret == NULL) goto end;

if (strncmp(buf,"ply",3) == 0) printf("LoadPLY");

if (strncmp(buf,"element vertex",14) == 0)
    {
      sscanf(buf+15, "%d", &nvert);

      printf("LoadPLY: loading %d points\n", nvert);

      Vertices = AllocVertices(nvert);
      
      Vertices->cur = nvert;
	}

if (strncmp(buf,"property float",14) == 0)
    {
      if(!((!strncmp(buf+15,"x",1)) ||(!strncmp(buf+15,"y",1))|| (!strncmp(buf+15,"z",1)) || (!strncmp(buf+15,"curv1",5)) || (!strncmp(buf+15,"curv2",5)) || (!strncmp(buf+15,"dcurv1",6))|| (!strncmp(buf+15,"dcurv2",6))|| (!strncmp(buf+15,"dcurv3",6))|| (!strncmp(buf+15,"dcurv4",6)))) goto end;
	}


if (strncmp(buf,"element face",12) == 0)
    {
      sscanf(buf+13, "%d", &nfaces);
      printf("LoadPLY: loading %d faces\n", nfaces);
      Faces = AllocFaces(nfaces);
      Faces->cur = nfaces;
	}


if (strncmp(buf,"end_header",10) == 0)
    {
	sscanf(buf+10, "%d", &t);
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf%f%f%f%f%f%f%d%d%d", &x, &y, &z, &c1, &c2, &d1, &d2, &d3, &d4, &r, &g, &b );
        Vertices->v[i].x = x; Vertices->v[i].y = y; Vertices->v[i].z = z;
	Vertices->v[i].curv1 = c1; Vertices->v[i].curv2 = c2;
	Vertices->v[i].dcurv[0] = d1; Vertices->v[i].dcurv[1] = d2 ; Vertices->v[i].dcurv[2] = d3 ; Vertices->v[i].dcurv[3] = d4;
	Vertices->v[i].red = r; Vertices->v[i].green = g; Vertices->v[i].blue = b;
	Vertices->v[i].nfaces = 0;
	
	}
	t = nfaces + nvert;
	j=i;
	i=0;
	
	
      for (; j < t; j++){
        fscanf(filein, "%d%d%d%d",&nf, &f1, &f2, &f3);
	
    /*
	We store vertices of each face, increase the counter of the faces that
	the vertex belongs to, and we add at the array of adjacent faces of the vertex the face index

 */

	Faces->f[i].vert[0] = f1;
	Vertices->v[f1].nfaces ++;
	Vertices->v[f1].face[Vertices->v[f1].nfaces - 1] = i;
	Faces->f[i].vert[1] = f2;
        Vertices->v[f2].nfaces ++;
	Vertices->v[f2].face[Vertices->v[f2].nfaces - 1] = i;
	Faces->f[i].vert[2] = f3;
	Vertices->v[f3].nfaces ++;
	Vertices->v[f3].face[Vertices->v[f3].nfaces - 1] = i;
	i ++;
	}
	}
  } // while (1)

 end:
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshPLY: bad file format\n");
    exit(0);
  }
} /* LoadMeshPLY() */

/* ==================================== */
void LoadMeshPLYvalue(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/* format: 
ply
format ascii 1.0
element vertex 141362 ---------------------->number of vertices
property float x
property float y
property float z
property float value
property uchar diffuse_red
property uchar diffuse_green
property uchar diffuse_blue
element face 282712 -------------------->number of faces
property list uchar int vertex_indices
end_header
-24.98689 -66.48193 -5.470881 255 197 3
-0.554481 -66.45267 -10.20065 255 253 208
-7.827967 -64.71899 -25.157
......
.....
..
....
3 71545 71546 71547
3 12979 71547 71499
3 13013 71548 71547
3 12980 71499 71548
3 71547 71548 71499
3 12980 71548 71500
3 13013 71549 71548
3 12995 71500 71549
3 71548 71549 71500
3 12996 71550 71506
3 12997 71551 71550  ...
  
*/

{

 int32_t i, j, t, n, nf, nvert=-1, nfaces=-1, f1, f2, f3;
  int32_t r, g, b;
  double x, y, z;
  float c1,c2, d1, d2, d3, d4;
#define sbuf 1024
  char buf[sbuf];
  char *ret;

  while (1)
  {
    //printf("LoadPLY: loading \n");
    do { ret=fgets(buf, sbuf, filein); }
    while (ret && strncmp(buf,"ply",3) && strncmp(buf,"element vertex",14)&& strncmp(buf,"property float",14) && strncmp(buf,"element face",12) && strncmp(buf,"end_header",10));
    if (ret == NULL) goto end;

if (strncmp(buf,"ply",3) == 0) printf("LoadPLY");

if (strncmp(buf,"element vertex",14) == 0)
    {
      sscanf(buf+15, "%d", &nvert);

      printf("LoadPLY: loading %d points\n", nvert);

      Vertices = AllocVertices(nvert);
      
      Vertices->cur = nvert;
	}

if (strncmp(buf,"property float",14) == 0)
    {
      if(!((!strncmp(buf+15,"x",1)) ||(!strncmp(buf+15,"y",1))|| (!strncmp(buf+15,"z",1)) || (!strncmp(buf+15,"value",5)) )) goto end;
	}


if (strncmp(buf,"element face",12) == 0)
    {
      sscanf(buf+13, "%d", &nfaces);
      printf("LoadPLY: loading %d faces\n", nfaces);
      Faces = AllocFaces(nfaces);
      Faces->cur = nfaces;
	}


if (strncmp(buf,"end_header",10) == 0)
    {
	sscanf(buf+10, "%d", &t);
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf%f%d%d%d", &x, &y, &z, &c1, &r, &g, &b );
        Vertices->v[i].x = x; Vertices->v[i].y = y; Vertices->v[i].z = z;
	Vertices->v[i].value = c1; 
	Vertices->v[i].red = r; Vertices->v[i].green = g; Vertices->v[i].blue = b;
	Vertices->v[i].nfaces = 0;
	}
	t = nfaces + nvert;
	j=i;
	i=0;
	
	
      for (; j < t; j++){
        fscanf(filein, "%d%d%d%d",&nf, &f1, &f2, &f3);
	
    /*
	We store vertices of each face, increase the counter of the faces that
	the vertex belongs to, and we add at the array of adjacent faces of the vertex the face index

 */

	Faces->f[i].vert[0] = f1;
	Vertices->v[f1].nfaces ++;
	Vertices->v[f1].face[Vertices->v[f1].nfaces - 1] = i;
	Faces->f[i].vert[1] = f2;
        Vertices->v[f2].nfaces ++;
	Vertices->v[f2].face[Vertices->v[f2].nfaces - 1] = i;
	Faces->f[i].vert[2] = f3;
	Vertices->v[f3].nfaces ++;
	Vertices->v[f3].face[Vertices->v[f3].nfaces - 1] = i;
	i ++;
	}
	}
  } // while (1)

 end:
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshPLY: bad file format\n");
    exit(0);
  }
} /* LoadMeshPLY() */


/* ==================================== */
void LoadMeshYAMS(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/* format:

MeshVersionFormatted	1
Dimension	3

#set of vertices
Vertices
num of vertices
-24.98689 -66.48193 -5.470881 0
-0.554481 -66.45267 -10.20065 0
-7.827967 -64.71899 -25.157 0
......
.....
..
....
Triangles
num of triangles
71545 71546 71547 0
12979 71547 71499 0
13013 71548 71547 0
  
*/

{

 int32_t i, j, t, n, nf, nvert=-1, nfaces=-1, f1, f2, f3;
  int32_t r, g, b, res;
  double x, y, z;
  float c1,c2, d1, d2, d3, d4;
#define sbuf 1024
  char buf[sbuf];
  char *ret;

  printf("loading YAMS \n" );
     while (1)
  {
    do { ret=fgets(buf, sbuf, filein);
	 }
    while (ret && (buf[0]!='V') && (buf[0]!='T') );
    if (ret == NULL) goto end;;
    if (strncmp(buf,"Vertices",8) == 0)
    {
	printf("vertices \n" );
	//sscanf(buf+10, "%d", &nvert);
      fscanf(filein, "%d",&nvert );
      printf("num vertices: %d\n",nvert );
      Vertices = AllocVertices(nvert);
      Vertices->cur = nvert;

      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf%d", &x, &y, &z, &res);
        Vertices->v[i].x = x; Vertices->v[i].y = y; Vertices->v[i].z = z;
	 Vertices->v[i].red = 200; Vertices->v[i].green = 200; Vertices->v[i].blue = 200;
      } // for i
    }
    else if (strncmp(buf,"Triangles",9) == 0)
    {
      //sscanf(buf+11, "%d", &nfaces);
      fscanf(filein, "%d", &nfaces);
      printf("num faces: %d\n",nfaces );
      Faces = AllocFaces(nfaces);
      Faces->cur = nfaces;

      for (i = 0; i < nfaces; i++){
        fscanf(filein, "%d%d%d%d", &f1, &f2, &f3 , &res);
	Faces->f[i].vert[0] = f1 -1;
 	Faces->f[i].vert[1] = f2 - 1;
	Faces->f[i].vert[2] = f3 -1;
	}
    }
    
    
  } // while (1)

 end:
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshYAMS: bad file format\n");
    exit(0);
  }//printf("LoadPLY: loading \n");
   
} /* LoadMeshPLY() */




/* ==================================== */
void LoadMeshIFS(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/*
The little-endian IFS (indexed face set) file format is:

file :=
  fileheader +
  vertexheader +
  vertex* +
  triheader +
  tri*

fileheader   := (string32) "IFS" + (float32)1.0 + (string32)modelname
vertexheader := (string32) "VERTICES" + (uint32)numVertices
vertex       := (float32)x + (float32)y + (float32)z
triheader    := (string32)"TRIANGLES" + (uint32)numFaces
tri          := (uint32)v0 + (uint32)v1 + (uint32)v2

Where string32 is an unsigned 32-bit integer length (which includes null termina
tor) followed by a string of that many u_int8_t's.  The strings above all i
nclude null termination.

Vertex indices are zero based.

Triangles have counter-clockwise winding direction in a right handed coordinate
system.  Triangles are "one sided".

Models "face" along the negative z-direction.  For example, the P51 Mustang
should fly along its object space -z axis.
*/
{
  int32_t i, j, n, nvert=-1, nfaces=-1, nvertfix;
  float version;
#define sbuf 1024
  char buf[sbuf];
  char *ret;
  u_int32_t s;

  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
  if (strcmp(buf, "IFS") != 0)
  {
    fprintf(stderr, "LoadMeshIFS: bad file format\n");
    exit(0);
  }
  version = ReadFloat32(filein);
#ifdef VERBOSE
  printf("IFS version=%f\n", version);
#endif
  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
#ifdef VERBOSE
  printf("model name: %s\n", buf);
#endif
  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
  nvert = (int32_t)ReadUnsignedLong(filein);
#ifdef VERBOSE
  printf("%s: %ld\n", buf, nvert);
#endif
  Vertices = AllocVertices(nvert);
  Vertices->cur = nvert;
  for (i = 0; i < nvert; i++)
  {
    Vertices->v[i].x = (double)ReadFloat32(filein);
    Vertices->v[i].y = (double)ReadFloat32(filein);
    Vertices->v[i].z = (double)ReadFloat32(filein);
  } // for i
  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
  nfaces = (int32_t)ReadUnsignedLong(filein);
#ifdef VERBOSE
  printf("%s: %ld\n", buf, nfaces);
#endif
  Faces = AllocFaces(nfaces);
  Faces->cur = nfaces;
  for (i = 0; i < nfaces; i++)
  {
    Faces->f[i].vert[0] = (int32_t)ReadUnsignedLong(filein);
    Faces->f[i].vert[1] = (int32_t)ReadUnsignedLong(filein);
    Faces->f[i].vert[2] = (int32_t)ReadUnsignedLong(filein);
  }
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshIFS: bad file format\n");
    exit(0);
  }
} /* LoadMeshIFS() */

/* ==================================== */
void LoadBuildMCM(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/* format: 
   V %d     (Sommets - champ obligatoire)
   %g %g %g (coord. vertex)
   ...

   F %d     (Faces - champ obligatoire)
   %d %d %d (face: ind. vertices)
   ...
  
*/
{
  int32_t i, j, n, nvert=-1, nfaces=-1, nvertfix;
  double x, y, z;
#define sbuf 1024
  char buf[sbuf];
  char *ret;
  double *vx, *vy, *vz;
  int32_t f1, f2, f3;
  while (1)
  {
    do { ret=fgets(buf, sbuf, filein); }
    while (ret && (buf[0]!='V') && (buf[0]!='F'));
    if (ret == NULL) goto end;
    if (buf[0]=='V')
    {
      sscanf(buf+1, "%d", &nvert);
      Vertices = AllocVertices(nvert);
      RBT = CreeRbtVide(nvert);
      vx = (double *)calloc(1,nvert * sizeof(double));
      vy = (double *)calloc(1,nvert * sizeof(double));
      vz = (double *)calloc(1,nvert * sizeof(double));
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf", &x, &y, &z);
        vx[i] = x; vy[i] = y; vz[i] = z;
      } // for i
    }
    else if (buf[0]=='F')
    {
      sscanf(buf+1, "%d", &nfaces);
      Faces = AllocFaces(nfaces);
      for (i = 0; i < nfaces; i++)
      {
        fscanf(filein, "%d%d%d", &f1,  &f2,  &f3);
        AddFace(vx[f1], vy[f1], vz[f1], vx[f2], vy[f2], vz[f2], vx[f3], vy[f3], vz[f3]);
      }
    }
  } // while (1)

 end:
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadBuildMCM: bad file format\n");
    exit(0);
  }
  RbtTermine(RBT);

} /* LoadBuildMCM() */

/* ==================================== */
void LoadBuildVTK(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/* format: 
# vtk DataFile Version 3.0
commentaire
ASCII
DATASET POLYDATA
POINTS %d float    (Sommets - champ obligatoire)
%g %g %g           (coord. vertex)
   ...

POLYGONS %d %d     (Faces - champ obligatoire)
3 %d %d %d         (face: ind. vertices)
   ...
*/
{
  int32_t i, j, n, nvert=-1, nfaces=-1, nvertfix;
  double x, y, z;
#define sbuf 1024
  char buf[sbuf];
  char *ret;
  double *vx, *vy, *vz;
  int32_t nf, f1, f2, f3;
  while (1)
  {
    do { ret=fgets(buf, sbuf, filein); }
    while (ret && strncmp(buf,"POINTS",6) && strncmp(buf,"POLYGONS",8));
    if (ret == NULL) goto end;
    if (strncmp(buf,"POINTS",6) == 0)
    {
      sscanf(buf+7, "%d", &nvert);
#ifdef VERBOSE
      printf("LoadBuildVTK: loading %d points\n", nvert);
#endif
      Vertices = AllocVertices(nvert);
      RBT = CreeRbtVide(nvert);
      vx = (double *)calloc(1,nvert * sizeof(double));
      vy = (double *)calloc(1,nvert * sizeof(double));
      vz = (double *)calloc(1,nvert * sizeof(double));
      for (i = 0; i < nvert; i++)
      {
        fscanf(filein, "%lf%lf%lf", &x, &y, &z);
        vx[i] = x; vy[i] = y; vz[i] = z;
      } // for i
    }
    else if (strncmp(buf,"POLYGONS",8) == 0)
    {
      sscanf(buf+9, "%d", &nfaces);
#ifdef VERBOSE
      printf("LoadBuildVTK: loading %d faces\n", nfaces);
#endif
      Faces = AllocFaces(nfaces);
      for (i = 0; i < nfaces; i++)
      {
        fscanf(filein, "%d%d%d%d", &nf, &f1,  &f2,  &f3);
        if (nf != 3)
        {
          fprintf(stderr, "LoadMeshVTK: faces must be triangles\n");
          exit(0);
        }        
        AddFace(vx[f1], vy[f1], vz[f1], vx[f2], vy[f2], vz[f2], vx[f3], vy[f3], vz[f3]);
      }
    }
  } // while (1)

 end:
  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshVTK: bad file format\n");
    exit(0);
  }
  RbtTermine(RBT);
  free(vx); 
  free(vy); 
  free(vz); 
} /* LoadBuildVTK() */

/* ==================================== */
void LoadBuildIFS(FILE *filein)
/* ==================================== */
/* filein doit avoir ete ouvert en lecture */
/*
The little-endian IFS (indexed face set) file format is:

file :=
  fileheader +
  vertexheader +
  vertex* +
  triheader +
  tri*

fileheader   := (string32) "IFS" + (float32)1.0 + (string32)modelname
vertexheader := (string32) "VERTICES" + (uint32)numVertices
vertex       := (float32)x + (float32)y + (float32)z
triheader    := (string32)"TRIANGLES" + (uint32)numFaces
tri          := (uint32)v0 + (uint32)v1 + (uint32)v2

Where string32 is an unsigned 32-bit integer length (which includes null termina
tor) followed by a string of that many u_int8_t's.  The strings above all i
nclude null termination.

Vertex indices are zero based.

Triangles have counter-clockwise winding direction in a right handed coordinate
system.  Triangles are "one sided".

Models "face" along the negative z-direction.  For example, the P51 Mustang
should fly along its object space -z axis.
*/
{
  int32_t i, j, n, nvert=-1, nfaces=-1, nvertfix;
  float version;
#define sbuf 1024
  char buf[sbuf];
  char *ret;
  u_int32_t s;
  double *vx, *vy, *vz;
  int32_t f1, f2, f3;

  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
  if (strcmp(buf, "IFS") != 0)
  {
    fprintf(stderr, "LoadMeshIFS: bad file format\n");
    exit(0);
  }
  version = ReadFloat32(filein);
#ifdef VERBOSE
  printf("IFS version=%f\n", version);
#endif
  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
#ifdef VERBOSE
  printf("model name: %s\n", buf);
#endif
  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
  nvert = (int32_t)ReadUnsignedLong(filein);
#ifdef VERBOSE
  printf("%s: %ld\n", buf, nvert);
#endif
  Vertices = AllocVertices(nvert);
  RBT = CreeRbtVide(nvert);
  vx = (double *)calloc(1,nvert * sizeof(double));
  vy = (double *)calloc(1,nvert * sizeof(double));
  vz = (double *)calloc(1,nvert * sizeof(double));
  for (i = 0; i < nvert; i++)
  {
    vx[i] = (double)ReadFloat32(filein);
    vy[i] = (double)ReadFloat32(filein);
    vz[i] = (double)ReadFloat32(filein);
  } // for i
  s = LE_ReadUnsignedLong(filein);
  fread(buf, sizeof(char), s, filein);
  nfaces = (int32_t)ReadUnsignedLong(filein);
#ifdef VERBOSE
  printf("%s: %ld\n", buf, nfaces);
#endif
  Faces = AllocFaces(nfaces);
  Faces->cur = nfaces;
  for (i = 0; i < nfaces; i++)
  {
    f1 = (int32_t)ReadUnsignedLong(filein);
    f2 = (int32_t)ReadUnsignedLong(filein);
    f3 = (int32_t)ReadUnsignedLong(filein);
    AddFace(vx[f1], vy[f1], vz[f1], vx[f2], vy[f2], vz[f2], vx[f3], vy[f3], vz[f3]);
  }

  if ((nvert == -1) || (nfaces == -1))
  {
    fprintf(stderr, "LoadMeshIFS: bad file format\n");
    exit(0);
  }
  free(vx); 
  free(vy); 
  free(vz); 
} /* LoadBuildIFS() */




