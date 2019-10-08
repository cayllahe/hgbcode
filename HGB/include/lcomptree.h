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



// #ifndef lcomptree_h
// #define lcomptree_h

//#include <graphSegmentation.h>




typedef struct JCsoncell
{
  int32_t son;
  struct JCsoncell *next;
  struct JCsoncell *last;
} JCsoncell;

typedef struct
{
  int32_t father;            // value -1 indicates the root
  int32_t nbsons;            // value -1 indicates a deleted node
  JCsoncell *sonlist;
  JCsoncell *lastson;
  JCsoncell *nodeAsSon;
} JCctreenode;

typedef struct
{
  int32_t nbnodes;
  int32_t nbsoncells;
  int32_t root;
  JCctreenode *tabnodes; 
  JCsoncell *tabsoncells;
  int32_t *flags;
} JCctree;

int32_t saliencyTree(graphe *g, JCctree ** SaliencyTree, double *STaltitude);
int32_t componentTree(graphe *g, JCctree ** CompTree, double *CTaltitude);
int32_t lextinctionvalues(graphe *g, int32_t mode, double *Av, double **Bitmap);



int32_t BPTAO(/* INPUTS */
              graphe *g,
              /* OUTPUTS */
              JCctree ** SaliencyTree, /* Component tree - BPTAO  */
              double *STaltitude,  /* weights associated to the
                                    nodes of the BPTAO: Must
                                    be allocated elsewhere */
              MST* mst
              );

int CanonizeQBT(/* INPUTS */
                JCctree *ST, /* Component tree - QBT  */
                double *STaltitude,
                /* OUTPUTS */
                JCctree **CompTreeOut);

void componentTreeFree(JCctree * CT);
void componentTreePrint(JCctree *CT, double *Alt, int32_t root, int32_t *Att1);

//Edward
void componentTreePrintAttributes(JCctree *CT, double *Alt, int32_t root, int *STArea, double* STMaxWeight);
int getMST(graphe *g, MST** mst );

int componentTreeAttribute(JCctree *CT, double *Alt, int32_t root, int32_t *Area, double *maxWeight);
int componentTreeAttributeMin(JCctree *CT, double *Alt, int32_t root, int32_t *Area, double *minWeight);
int componentTreeAttributeMean(JCctree *CT, double *Alt, int32_t root, int32_t *Area, double * minWeight);


int getComponentAtLambda(JCctree *CT, double *Alt, const double lambda, const int vertex);

void computeQFZ_Attributes(graphe *g, JCctree **CompTree, double **STAltitudes , int** STArea , double** STMaxWeight);
void computeQFZ_AttributesMin(graphe *g, JCctree **CompTree, double **STAltitudes , int** STArea , double** STMaxWeight);
void computeQFZ_AttributesMean(graphe *g, JCctree **CompTree, double **STAltitudes , int** STArea , double** STMaxWeight);

int CanonizeQBTSegmentation(/* INPUTS */
                            JCctree *ST, /* Component tree - QBT  */
                            double *STaltitude,
                            
                            /* OUTPUTS */
                            JCctree **CompTreeOut,
                            double *STMaxWeight);

int CanonizeQBTSegmentationMin(/* INPUTS */
                            JCctree *ST, /* Component tree - QBT  */
                            double *STaltitude,
                            
                            /* OUTPUTS */
                            JCctree **CompTreeOut,
                            double *STMaxWeight);

int CanonizeQBTSegmentationMean(
                                JCctree *ST,
                                double *STaltitude,
                                
                                JCctree **CompTreeOut,
                                double *STMinWeight);


int32_t BPTAO_segmentation(/* INPUTS */
                           graphe *g,
                           /* OUTPUTS */
                           JCctree ** SaliencyTree, /* Component tree - BPTAO  */
                           double *STaltitude,  /* weights associated to the
                                                 nodes of the BPTAO: Must
                                                 be allocated elsewhere */
                           double *STOrigAltitude
                           );

int32_t BPTAO_segmentationMin(/* INPUTS */
                           graphe *g,
                           /* OUTPUTS */
                           JCctree ** SaliencyTree, /* Component tree - BPTAO  */
                           double *STaltitude,  /* weights associated to the
                                                 nodes of the BPTAO: Must
                                                 be allocated elsewhere */
                           double *STOrigAltitude
                           );



int treeUpdate( /* INPUTS */
               JCctree **CompTree, /* Component tree - QFZ */
               double  *STaltitudes,
               graphe *g,
               int edgeIdx 
);



JCctree * componentTreeAlloc(int32_t N);

int IncrementalQFZ(graphe *g, char* file);
int treeUpdateAltitude( /* INPUTS */
                       JCctree **CompTree, /* Component tree - QFZ */
                       double  *STaltitudes,
                       double  *STOrigAltitudes,
                       graphe *g,
                       int edgeIdx /* Edge index for {x,y}  */
/* OUTPUTS */
//       CompTree
);



double computeSM(JCctree *ST, graphe *g, double *STAltitudes);


int IncrementalQFZfromMST( graphe *g, char *file );


int initializeTree(/* INPUTS */
                   graphe *g,
                   /* OUTPUTS */
                   JCctree **CompTree,/* CompTree :  initalized component tree - QFZ */
                   double  *STaltitude
                   );

////// With attributes
int initializeTreeAttributes(/* INPUTS */
                             graphe *g,
                             /* OUTPUTS */
                             JCctree **CompTree,/* CompTree :  initalized component tree - QFZ */
                             double *STaltitude,
                             double * STAltitudesSegmentation,
                             int *STArea,
                             double *STMaxWeight
                             );

int treeUpdateAttributes( /* INPUTS */
                         JCctree **CompTree, /* Component tree - QFZ */
                         double  *STaltitudes,
                         double *STAltitudesSegmentation,
                         int* STArea,
                         double *STMaxWeight,
                         graphe *g,
                         int edgeIdx /* Edge index for {x,y}  */
);

void mergeAttributes(/* INPUTS */
                     JCctree **CompTree, /* Component tree - QFZ */
                     double  *STaltitudes,
                     int *STArea,
                     double *STMaxWeight,
                     int p1,
                     int p2,
                     int *np,
                     int *area1,
                     int *area2
);

void attachAttributes(/* INPUTS */
                      JCctree **CompTree, /* Component tree - QFZ */
                      double  *STaltitudes,
                      double *STAltitudesSegmentation,
                      int *STArea,
                      double *STMaxWeight,
                      const int idxP,
                      const int idxS, //  idxS is attached to idxP
                      int *area
);
void createNodeAttributes( /* INPUTS */
                          JCctree **CompTree, /* Component tree - QFZ */
                          double  *STaltitudes,
                          double *STAltitudesSegmentation,
                          int *STArea,
                          double *STMaxWeight,
                          graphe *g,
                          int edgeIdx, /* Edge index for {x,y}  */
                          int c1,
                          int c2,
                          int *n,int *area1, int *area2 );
int random_int(int min, int max);
void saveDOTfile(JCctree *CT, double *Alt, int32_t root, int *STArea, double* STMaxWeight , char * file, int nbsom);

// #endif
