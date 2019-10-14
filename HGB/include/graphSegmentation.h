/*
Copyright ESIEE (2019) 

edward.cayllahua@esiee.fr

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

/* authors : Edward Cayllahua - J. Cousty and Silvio Guimarães */

#ifndef graphSegmentation_h
#define graphSegmentation_h


#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <MST.h>
#include <mcweightgraph.h>
#include <mcunionfind.h>
#include <lcomptree.h>
#include <llca.h>
#include <mcsort.h>
#include <stdio.h>
#include <list.h>






//void clean_cache();

void naiveSegmentation(graphe *gs);
void segmentationRange(graphe *g);
void segmentationBranch(graphe *g);

void incrementalNaiveSegmentation(graphe *g);
void incrementalSegmentationRange(graphe* g);
void incrementalSegmentationBranch(graphe *g);

void incrementalSegmentationBranchCoarser(graphe *g);
void incrementalSegmentationInterval(graphe *g, int op);
void intervalSegmentation(graphe *g, int op, int th, double param);


void incrementalSegmentationBranchTimed(graphe *g);

double naiveMinimization(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx );
double rangeMinimization(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , NodeList**frange);
double branchMinimization(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx );

double branchMinimizationLambdaStars(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op);
double branchMinimizationLambdaStarsAll(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op);
double minimizationLambdasInterval(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op, int* numberInterval, int *max_numberInterval);


double minimizationLambdasIntervalVolume(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op, int* numberInterval, int *max_numberInterval);

double scaleSelection(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op, int th, double prank);


void reweightToMax(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes );

int computeLevels(graphe *g, double * sum);

// double computeD(JCctree *CT,double *STAltitudes , int *STArea, double* STMaxWeight , int *cx, int *cy , double Dif , double lambda );

double computeD( int AreaCx, int AreaCy, double STMaxWeightCx, double STMaxWeightCy, double Dif);

/**** Dual HGB method*/

void dualhgb(graphe *g);
double max_rule(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx);
double max_rule_parent(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx);
double max_rule_min(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx);
double max_rule_mean(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx);
double max_rule_trunc(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx);
void dualhgbInverse(graphe *g);
void dualhgbParent(graphe *g);
void dualhgbMin(graphe *g);
void dualhgbMean(graphe *g);
void dualhgbMeanNormal(graphe *g);
void dualhgbTruncated(graphe *g);


double computeDmin( int AreaCx, int AreaCy, double STMaxWeightCx, double STMaxWeightCy, double Dif);
double computeDmean( int AreaCx, int AreaCy, double STMinWeightCx, double STMinWeightCy, double Dif);

void dualhgbParentMin(graphe *g);

double max_rule_parent_min(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx);

#endif /* graphSegmentation_h */
