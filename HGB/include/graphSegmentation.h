//
//  graphSegmentation.h
//  SM_Edward
//
//  Created by Edward Cayllahua on 8/22/16.
//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

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
