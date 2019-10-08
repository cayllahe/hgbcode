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
#include <math.h>
#include <float.h>
#include <string.h>
#include <mcweightgraph.h>
#include <mcmesh.h>
#include <mciomesh.h>
#include <lgraphmesh.h>

#define EPS 0.000000001 



//****************************************************************************80

double triangle_area_3d_3 ( double *t )

  //****************************************************************************80
  //
  //  Purpose:
  ////    TRIANGLE_AREA_3D_3 computes the area of a triangle in 3D.
  ////  Discussion:
  ////    This routine uses Heron's formula
  ////  Modified:
  ////   31 July 2005
  ////  Author:
  ////    John Burkardt
  ////  Reference:
  ////    Adrian Bowyer, John Woodwark,
  //    A Programmer's Geometry,
  //    Butterworths, 1983.
  ////  Parameters:
  //
  //    Input, double T[3*3], the triangle vertices.
  //
  //    Output, double AREA, the area of the triangle.
  //
{
# define DIM_NUM 3

  double area;
  int i;
  int j;
  int jp1;
  double s[DIM_NUM];

  for ( j = 0; j < DIM_NUM; j++ )
  {
    jp1 = ( j + 1 ) % DIM_NUM;
    s[j] = 0.0;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      s[j] = s[j] + pow ( t[i+j*DIM_NUM] - t[i+jp1*DIM_NUM], 2 );
    }
    s[j] = sqrt ( s[j] );
  }

  area = (   s[0] + s[1] + s[2] ) 
    * ( - s[0] + s[1] + s[2] ) 
    * (   s[0] - s[1] + s[2] ) 
    * (   s[0] + s[1] - s[2] );

  if ( area < 0.0 )
  {
    area = -1.0;
    return area;
  }

  area = 0.25 * sqrt ( area );
  return area;
# undef DIM_NUM
}

double signe(double a)
{

  if( a < 0 ){ 
    return -1.0;
  }else {
    return 1.0;
  }
}
double maxdoub(double a, double b)
{

  if( a >= b ){ 
    return a;
  }else {
    return b;
  }
}
double mindoub(double a, double b)
{

  if( a <= b ){ 
    return a;
  }else {
    return b;
  }
}

double minarray(double *a, int32_t size)

{
  int32_t i;
  double tmp=1000;
  for(i=0;i<size;i++){
    if( a[i] <= tmp ){ 
      tmp=a[i];
    }
  }
  if(size!=0)return tmp;
}

double maxarray(double *a, int32_t size)

{
  int32_t i;
  double tmp=-FLT_MAX;
  for(i=0;i<size;i++){
    if( a[i] >= tmp ){ 
      tmp=a[i];
    }
  }
  if(size!=0)return tmp;
}

double meanarray(double *a, int32_t size)

{
  int32_t i;
  double tmp=0;
	
  for(i=0;i<size;i++){
		 
    tmp=tmp+a[i];
  }
  tmp=tmp/size;
  if(size!=0)return tmp;
}







//****************************************************************************80
computeTrianglesArea(double *Atriangles)
  //****************************************************************************80
  /*COMPUTE TRIANGLES AREA*/

{
  int32_t i;

  int32_t numfaces = Faces->cur;



  double *Tvert= (double *)calloc(1,9 * sizeof(double));
  for(i = 0; i < numfaces; i++){
	
    Tvert[0]=Vertices->v[Faces->f[i].vert[0]].x;
    Tvert[1]=Vertices->v[Faces->f[i].vert[0]].y;
    Tvert[2]=Vertices->v[Faces->f[i].vert[0]].z;

    Tvert[3]=Vertices->v[Faces->f[i].vert[1]].x;
    Tvert[4]=Vertices->v[Faces->f[i].vert[1]].y;
    Tvert[5]=Vertices->v[Faces->f[i].vert[1]].z;

    Tvert[6]=Vertices->v[Faces->f[i].vert[2]].x;
    Tvert[7]=Vertices->v[Faces->f[i].vert[2]].y;
    Tvert[8]=Vertices->v[Faces->f[i].vert[2]].z;
		
    Atriangles[i]= triangle_area_3d_3 ( Tvert );
    if(Atriangles[i]==-1)printf("BAD AREA :  face :%d\n", i);
  }

  /*TRIANGLES AREA COMPUTED*/
  free(Tvert);

}


//****************************************************************************80
void computecurv(char *curvOption)
  //****************************************************************************80
{
  int32_t i,h;
  double value, k1,k2,H,K; 
  double DerCurv[4];
  int32_t numvert = Vertices->cur;
	
  for (i = 0; i < numvert ; i++){

    if(!strcmp(curvOption, "dcurv")){

      /*DERIVATIVES*/
      for (h = 0; h < 4 ; h++){
	DerCurv[h] = Vertices->v[i].dcurv[h];
		
      }
				
      value = DerCurv[0]*DerCurv[0] + DerCurv[1]*DerCurv[1] + DerCurv[2]*DerCurv[2] + DerCurv[3]*DerCurv[3];
		
      /*fin DERIVATIVES*/

    }else{

      k1 = Vertices->v[i].curv1;
      k2 = Vertices->v[i].curv2;
		
		

      if(!strcmp(curvOption, "mean")){

	H = (k1 + k2)/2;
	value = 256 / M_PI * (atanf(-H*100) + M_PI_2) ;			

      }else if(!strcmp(curvOption, "max")){

	value =  maxdoub (k1*k1, k2*k2) ;

      }else if(!strcmp(curvOption, "meanConvex")){
	H = (k1 + k2)/2;
	value = 1 / M_PI * (atanf(H) + M_PI_2) ;

      }else if(!strcmp(curvOption, "gauss")){
	K = k1*k2;
	value = 1 / M_PI * (atanf(K) + M_PI_2) ;
      }else if(!strcmp(curvOption, "deviation")){
	K = k1*k2;
	H = (k1 + k2)/2;
	value = (4*H*H - 2*K*K) ;
	value = 1 / M_PI * (atanf(value) + M_PI_2) ;
      }
    }
		
        
    Vertices->v[i].value=value;
		

  }//for

	



}



//****************************************************************************80
graphe *mesh2graph(double *Atriangles)

  //****************************************************************************80
{
  int32_t numvert,numedge,numfaces;
  int32_t nummax;
  int32_t i,j,k;
  double value,v1,v2;
  graphe *g;

  numvert = Vertices->cur;
  printf("num vertices : %d\n", numvert);
  numedge = Edges->cur;
  printf("num edges : %d\n", numedge);
  numfaces = Faces->cur;
  printf("num faces : %d\n", numfaces);

  computeTrianglesArea(Atriangles);

  /*START BUILDING GRAPH*/

  /*
    graph with one node per face of the mesh, and an edge
    between faces if the faces share one side.
	
    each edge take the mean values of the  curvatures at their 
    vertices

  */


  nummax = numedge * 2;
  g = initgraphe(numfaces, nummax);
  printf("graph init done   building graph ....  \n");

  j = 0;

  for (i = 0; i < numedge ; i++){

    v1 = Vertices->v[Edges->e[i].v1].value;
    v2 = Vertices->v[Edges->e[i].v2].value;
    value=(v1+v2)/2;

		
		
    if((Edges->e[i].f1 == -1) || (Edges->e[i].f2 == -1)){
      printf("NOT PROPER FACE [%d]  f1: %d f2: %d  \n", i, Edges->e[i].f1, Edges->e[i].f2);
      exit(1);
    }else {
	
      addarete(g , Edges->e[i].f1, Edges->e[i].f2, value);
	 
    }		
  }//END FOR ....
	
	
  return g;

}


//****************************************************************************80
void LoadSaveCurvMesh(char *fileinname,char *curvOption,char *fileoutname) 

  //****************************************************************************80
{
  FILE *filein, *fileout;
  int32_t numvert,numedge,numfaces;

  filein = fopen(fileinname,"r");
  if (!filein)
  {
    fprintf(stderr, " cannot open file: %s\n",  fileinname);
    exit(0);
  }

  fileout = fopen(fileoutname,"w");
  if (!fileout)
  {
    fprintf(stderr, " cannot open file: %s\n",  fileoutname);
    exit(0);
  }
	
  /*LOAD MESH ########################### */ 

  LoadMeshPLY(filein);
  printf("Loading done ..\n");
  Edges = AllocEdges(1);
  printf("alloc done ..\n");
  ComputeEdges();
  printf("compute edges done ..\n");
  numvert = Vertices->cur;
  printf("num vertices : %d\n", numvert);
  numedge = Edges->cur;
  printf("num edges : %d\n", numedge);
  numfaces = Faces->cur;
  printf("num faces : %d\n", numfaces);

  computecurv(curvOption);
  SaveMeshPLYvalue(fileout);

  fclose(filein);
  fclose(fileout);
  TermineMesh();


}





//****************************************************************************80
void LoadSaveGraphMesh(char *graphfilename,char *meshfilename,char *mode,char *meshfileout)

  //****************************************************************************80
{
  int32_t numvert,numedge,numfaces;
  FILE *filein, *fileout;
  double EdgeValues[25*MAXADJFACES];
  double *Av;
  graphe *g;
  int32_t modus, nbedges;
  int32_t i,j,k,q;
  double value;
  int32_t value2;

  filein = fopen(meshfilename,"r");
  if (!filein)
  {
    fprintf(stderr, " cannot open file: %s\n", meshfilename);
    exit(0);
  }

  fileout = fopen(meshfileout,"w");
  if (!fileout)
  {
    fprintf(stderr, " cannot open file: %s\n",meshfileout);
    exit(0);
  }
	
  /*LOAD MESH ########################### */ 

  LoadMeshPLYvalue(filein);
  printf("Loading done ..\n");
  Edges = AllocEdges(1);
  printf("alloc done ..\n");
  ComputeEdges();
  printf("compute edges done ..\n");
  numvert = Vertices->cur;
  printf("num vertices : %d\n", numvert);
  numedge = Edges->cur;
  printf("num edges : %d\n", numedge);
  numfaces = Faces->cur;
  printf("num faces : %d\n", numfaces);
	
  /*READ GRAPH #######################*/	
  g= ReadGraphe(graphfilename,&Av);

	
  //cor=ReadCorrespondance(char *fileinname);

  if(!strcmp(mode, "min")){ modus=1;
  }else if(!strcmp(mode, "max")){modus=2;
  }else if(!strcmp(mode, "mean")){modus=3;
  }else {	printf("BAD FILTER MODE \n");
    usage();
    exit(1);
  }


  for (i = 0; i < numvert; i++){
    nbedges=Vertices->v[i].nedges;
    q=0;
    for(j=0;j<Vertices->v[i].nedges;j++){
      k=Vertices->v[i].edge[j];
      if((getweight(g,k)<1000)&&(getweight(g,k)>=0)){
	EdgeValues[q]=getweight(g,k);
	q=q+1;
      }else nbedges=nbedges-1;
    }
    if(nbedges!=0){
      if (modus==1) Vertices->v[i].value=minarray(EdgeValues,nbedges);
      else if (modus==2) Vertices->v[i].value=maxarray(EdgeValues,nbedges);
      else Vertices->v[i].value=meanarray(EdgeValues,nbedges);
		
      value= Vertices->v[i].value;
      //value = 256 * value;
      if (value < 0)printf("%g\n",value);
      value2 = (int32_t) value;
      if (value2 < 0) {value2=255;printf("%g\n",value);}
      if (value2 > 255) value2=255;
      Vertices->v[i].red = 255-value2;
      Vertices->v[i].green = 255-value2;
      Vertices->v[i].blue = 255-value2;
    }else{

      Vertices->v[i].value=FLT_MAX;
      Vertices->v[i].red = 255;
      Vertices->v[i].green = 1;
      Vertices->v[i].blue = 1;


    }	
  }

  SaveMeshPLYvalue(fileout); 

  fclose(fileout);
  fclose(filein);
  free(Av);
  terminegraphe(g);
  TermineMesh();


}

/* Returns the extremities of the edge which belongs to both faces i
   and j 
   returns 1 is such an edge exists
   returns 0 if such an edge does not exist.
*/
int getCommonEdge(int32_t i, int32_t j, int32_t *x, int32_t *y)
{
  int32_t k,z;
  u_int32_t first = 0;
  *x = -1;
  *y = -1;
  /* on parcourt tous les sommets de la face i */
  for(k = 0; k < 3; k++)
    for(z = 0; z <3; z++){
      printf("Vertices A = %d  and B = %d\n",Faces->f[i].vert[k], Faces->f[j].vert[z]);
      if( (Faces->f[i].vert[k] == Faces->f[j].vert[z]) && (Faces->f[i].vert[k] != (*x)) )
	if(*x == -1)
	  *x = Faces->f[i].vert[k]; 
	else
	  *y = Faces->f[i].vert[k];
    }
  if( ( *x != 0) && (*y != 0))
    return 1;
  else return 0;    
}



//****************************************************************************80
void LoadSaveGraphMeshJC(char *graphfilename,char *meshfilename,char *mode,char *meshfileout)
  // L'idée est de sortir 1) la surface de la statue et 2 un fichier polyline vtk de la LPE
  //****************************************************************************80
{
  int32_t numvert,numedge,numfaces;
  FILE *filein, *fileout, *filecurve;
  double EdgeValues[25*MAXADJFACES];
  double *Av;
  graphe *g;
  int32_t modus, nbedges;
  int32_t i,j,k,q;
  int32_t x, y, nbMarkedEdges=0;
  double value ;
  int32_t value2;
  u_int8_t * markEdge;
  char BUFFER[150];

  filein = fopen(meshfilename,"r");
  if (!filein)
  {
    fprintf(stderr, " cannot open file: %s\n", meshfilename);
    exit(0);
  }

  fileout = fopen(meshfileout,"w");
  if (!fileout)
  {
    fprintf(stderr, " cannot open file: %s\n",meshfileout);
    exit(0);
  }
	
  /*LOAD MESH ########################### */ 

  LoadMeshPLYvalue(filein);
  printf("Loading done ..\n");
  Edges = AllocEdges(1);
  printf("alloc done ..\n");
  ComputeEdges();
  printf("compute edges done ..\n");
  numvert = Vertices->cur;
  printf("num vertices : %d\n", numvert);
  numedge = Edges->cur;
  printf("num edges : %d\n", numedge);
  numfaces = Faces->cur;
  printf("num faces : %d\n", numfaces);

  markEdge = (u_int8_t *)malloc(sizeof(u_int8_t) * Edges->cur); 
  memset(markEdge,0,Edges->cur);

  /*READ GRAPH #######################*/	
  g= ReadGraphe(graphfilename,&Av);
	
  //cor=ReadCorrespondance(char *fileinname);

  if(!strcmp(mode, "min")){ modus=1;
  }else if(!strcmp(mode, "max")){modus=2;
  }else if(!strcmp(mode, "mean")){modus=3;
  }else {	printf("BAD FILTER MODE \n");
    usage();
    exit(1);
  }

  for (i = 0; i < numvert; i++){
    nbedges=Vertices->v[i].nedges;
    q=0;
    for(j=0;j<nbedges;j++){
      k=Vertices->v[i].edge[j];
      // value = poids de l'arete k = (i,j)
      value = getweight(g,k);
      if(value > 0){
	// l'arete appartient à la watershed
	if(markEdge[k] != 1){
	  markEdge[k] = 1;
	  nbMarkedEdges++;
	}
      }

      if((getweight(g,k)<1000)&&(getweight(g,k)>=0)){
	EdgeValues[q]=getweight(g,k);
	q=q+1;
      }else nbedges=nbedges-1;
    }
    if(nbedges!=0){
      if (modus==1) Vertices->v[i].value=minarray(EdgeValues,nbedges);
      else if (modus==2) Vertices->v[i].value=maxarray(EdgeValues,nbedges);
      else Vertices->v[i].value=meanarray(EdgeValues,nbedges);
		
      value= Vertices->v[i].value;
      //value = 256 * value;
      if (value < 0)printf("%g\n",value);
      value2 = (int32_t) value;
      if (value2 < 0) {value2=255;printf("%g\n",value);}
      if (value2 > 255) value2=255;
      Vertices->v[i].red = 255-value2;
      Vertices->v[i].green = 255-value2;
      Vertices->v[i].blue = 255-value2;
    }else{

      Vertices->v[i].value=FLT_MAX;
      Vertices->v[i].red = 255;
      Vertices->v[i].green = 1;
      Vertices->v[i].blue = 1;
    }	
  }
  
  printf("Nombre d'aretes marquees %d\n",nbMarkedEdges);

  /* SAVE MESH EN VTK */
  strcpy(BUFFER, meshfileout);
  strcat(BUFFER,"Mesh.vtk"); 

  fileout = fopen(BUFFER,"w");
  genheaderVTK(fileout,"graphToMeshConverter"); 
  SaveMeshVTK(fileout);
  fclose(fileout); 
  
  /* SAVE LINES IN VTK */
  strcpy(BUFFER, meshfileout);
  strcat(BUFFER,"Curve.vtk");
  filecurve = fopen(BUFFER,"w");
  genheaderVTK(filecurve,"graphToMeshConverterLINES");
  SaveWeightedCurveVTK(filecurve,markEdge,nbMarkedEdges,g->weight);
  fclose(filecurve);
  
  fileout = fopen(meshfileout,"w");
  SaveMeshPLYvalue(fileout); 
  fclose(fileout);

  //  fclose(fileout);
  fclose(filein);
  free(Av);
  terminegraphe(g);
  TermineMesh();


}
