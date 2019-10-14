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


/* $Id: mciomesh.h,v 1.5 2006/02/28 07:49:12 michel Exp $ */
#define UNKNOWN -1
#define POV  0  /* PovRay 3.1 */ 
#define COL  1  /* private format (ESIEE - IFC students) */
#define MCM  2  /* private format (MC) */
#define AC   3  /* AC3D format */
#define POVB 4  /* PovRay 3.1 without header and footer */ 
#define GL   5  /* Structures de donnees C pour OpenGL */
#define DXF  6  /* DXF */
#define VTK  7  /* VTK PolyData */
#define IFS  8  /* Indexed Face Set (Brown University) */
#define PGM  9  /* 3D pgm extension (point clouds only) */
#define PLY  10  /* Princeton Format */

extern void genheaderPOV(FILE *fileout, int32_t obj_id, meshbox MB);
extern void genheaderAC(FILE *fileout, meshbox MB, double red, double green, double blue, char *name);
extern void genheaderCOL(FILE *fileout, int32_t nbfaces);
extern void genheaderDXF(FILE *fileout);
extern void genheaderVTK(FILE *fileout, char *name);
extern void genfooterPOV(FILE *fileout);
extern void genfooterAC(FILE *fileout);
extern void genfooterDXF(FILE *fileout);
extern void genfacePOV(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
             double x3, double y3, double z3
		       );
extern void gencarrePOV(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
             double x3, double y3, double z3, 
             double x4, double y4, double z4
		 );
extern void genfaceCOL(FILE *fileout, int32_t face_id, int32_t object_id,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
	     double x3, double y3, double z3,
             double xp1, double yp1, double zp1, 
             double xp2, double yp2, double zp2, 
             double xp3, double yp3, double zp3);
extern void genfaceDXF(FILE *fileout,
             double x1, double y1, double z1, 
             double x2, double y2, double z2, 
	     double x3, double y3, double z3,
             double xp1, double yp1, double zp1, 
             double xp2, double yp2, double zp2, 
             double xp3, double yp3, double zp3);
extern void SaveMeshPOV(FILE *fileout);
extern void SaveMeshCOL(FILE *fileout, int32_t obj_id);
extern void SaveMeshMCM(FILE *fileout);
extern void SaveMeshPLY(FILE *fileout);
extern void SaveMeshPLYlabels(FILE *fileout, int32_t *FaceLabels);
extern void SaveMeshPLYsimple(FILE *fileout);
extern void SaveMeshYAMS(FILE *fileout);
extern void SaveMeshVTK(FILE *fileout);
extern void SaveMeshAC(FILE *fileout);
extern void SaveMeshGL(FILE *fileout, int32_t obj_id);
extern void SaveMeshDXF(FILE *fileout, int32_t obj_id);
extern void LoadMeshMCM(FILE *filein);
extern void LoadMeshPLY(FILE *filein);
extern void LoadMeshYAMS(FILE *filein);
extern void LoadMeshIFS(FILE *filein);
extern void LoadMeshAC(FILE *filein);
extern void LoadBuildMCM(FILE *filein);
extern void LoadBuildVTK(FILE *filein);
extern void LoadBuildIFS(FILE *filein);
extern void SaveCurveVTK(FILE *fileout, u_int8_t* markEdge, int32_t nbMarked);
