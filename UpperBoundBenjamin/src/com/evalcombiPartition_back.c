
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <mccodimage.h>
#include <mcimage.h>
#include <mcweightgraph.h>
#include <lhierarchie.h>
#include <lcomptree.h>
#include <string.h>

#include "evalcombiPartition.h"


#define DEBUG
#ifdef DEBUG
#define printd(...) (printf(__VA_ARGS__))
#else
#define printd(...) {;}
#endif



inline int mini(int a, int b){
	return (a<=b)?a:b;
}

inline int maxi(int a, int b){
	return (a>=b)?a:b;
}

inline float minf(float a, float b){
	return (a<=b)?a:b;
}

inline float maxf(float a, float b){
	return (a>=b)?a:b;
}

inline double mind(double a, double b){
	return (a<=b)?a:b;
}

inline double maxd(double a, double b){
	return (a>=b)?a:b;
}

inline void * myCalloc(size_t num, size_t size)
{
	void * ref = calloc(num, size);
	if (ref==NULL)
	{ fprintf(stderr,"Cannot allocate enough memory\n"); exit(1); }
	return ref;
}

void usage(char *arg){
  fprintf(stderr,"#################################################################\n\n");
  fprintf(stderr,"USAGE: %s saliency_map.graph  ground_truth.pgm scoreFunction maxk\n",arg);
  fprintf(stderr,"\t scoreFunction in {'DHamming1', 'DHamming2', 'SHamming', 'DCovering1', 'DCovering2'; 'SCovering','BCE'}\n");
  fprintf(stderr,"#################################################################\n\n");
}

grayImage * readGrayImage(char * filename)
{

	struct xvimage *Im;
	//uint8_t *F;
	Im = readimage(filename);
	if(Im==NULL)
	{
		exit(1);
	}

	if(depth(Im) != 1) {
		fprintf(stderr,"%s: not yet implemented for this kind of data \n",filename);
		exit(1);
	} 
	//rs = rowsize(Im);
	//cs = colsize(Im);
	//N = rs*cs;

	switch(datatype(Im)){
	case 	VFF_TYP_1_BYTE: 
		//F = UCHARDATA(Im);
		break;
	default : 
		fprintf(stderr,"%s: not yet implemented for this kind of data \n",filename);
		exit(1);
	}
	
	return Im;
}






int GnbRegionsGT=0;

/*
Assume that there is no gap in region numbering, 
The segmentation labels will also be shifted such that label numbering start at 0.
*/
GroundTruth analyseGT(grayImage * image)
{
	uint8_t *seg = UCHARDATA(image);
	int32_t i, width=rowsize(image), heigth=colsize(image), N=width*heigth;
	printd("Ground thruth dimension:%dx%d\n",heigth, width);
	int min = 255;
	int max = 0;
	for(i=0;i<N;++i)
	{
		min = mini(min, seg[i]);
		max = maxi(max, seg[i]);
	}
	if(min>0)
	{
		for(i=0;i<N;++i)
		{
			seg[i]-=min;
		}
	}
	int nbElement =  max-min+1;
	int * area = myCalloc(nbElement, sizeof(int));
	for(i=0;i<N;++i)
		++area[seg[i]];
	printd("Nb regions in ground truth: %d\n",nbElement);
	GnbRegionsGT=nbElement;
	return (GroundTruth){.labelImage=image, .nbElement=nbElement, .area=area, .nbPixels=N};
}

void freeGT(GroundTruth gt)
{
	free(gt.area);
}

int * computeNodeArea(JCctree *CT)
{
	int * area = myCalloc(CT->nbnodes, sizeof(int));
	int i;
	JCctreenode *tabnodes = CT->tabnodes;
	for(i=0;i<CT->nbnodes;++i)
	{
		if(i<=CT->nbnodes/2)
		{
			area[i]=1;
		}else{
			JCsoncell * soncell;
			for( soncell = tabnodes[i].sonlist; soncell!=NULL; soncell=soncell->next)
			{
				int child = soncell->son;
				area[i]+=area[child];	
			}
		}
	}
	return area;
}



MatchingStat initMatchingStat(JCctree *CT, GroundTruth gt, int * nodeArea)
{
	int nbNodes = CT->nbnodes;
	int nbRegionInGT = gt.nbElement;
	int * regionGTSize = gt.area;
	int * regionTreeSize = nodeArea;
	uint8_t * bestMatch = myCalloc(nbNodes, sizeof(uint8_t));
	RegionMatching * bestMatchingScore = myCalloc(nbNodes, sizeof(RegionMatching));
	return (MatchingStat){.nbPixels=gt.nbPixels, .nbNodes = nbNodes, .nbRegionInGT = nbRegionInGT, .regionGTSize = regionGTSize, .regionTreeSize = regionTreeSize, .bestMatch = bestMatch, . bestMatchingScore = bestMatchingScore, .detailedMatchingScore=NULL};
}

void freeMatchingStat(MatchingStat s)
{
	free(s.bestMatch);
	free(s.bestMatchingScore);
	int i;
	for(i=0; i<s.nbNodes; ++i)
		free(s.detailedMatchingScore[i]);
	free(s.detailedMatchingScore);
}



MatchingStat computeMatchingStat(JCctree *CT, GroundTruth gt, int * nodeArea, RegionScoringFunction  regionScoringFunction, bool additif)
{
	int i,j;
	int nbE = gt.nbElement;
	int nbNodes=CT->nbnodes;
	printd("Init Matching struct\n");
	MatchingStat stat = initMatchingStat(CT, gt, nodeArea);
	RegionMatching ** matchingScore = myCalloc(CT->nbnodes, sizeof(RegionMatching *));
	for(i=0; i<nbNodes; ++i)
		matchingScore[i] = myCalloc(nbE, sizeof(RegionMatching ));
	
	JCctreenode *tabnodes = CT->tabnodes;
	uint8_t *seg = UCHARDATA(gt.labelImage);
	for(i=0;i<nbNodes;++i)
	{
		
		// Compute intersection with each region
		RegionMatching * c = matchingScore[i];
		if(i<=nbNodes/2)
		{
			++(c[seg[i]].cardInter);
		}else{
			JCsoncell * soncell;
			for( soncell = tabnodes[i].sonlist; soncell!=NULL; soncell=soncell->next)
			{
				int child = soncell->son;
				
				RegionMatching * cc = matchingScore[child];
				for(j=0;j<nbE;++j)
				{
					c[j].cardInter += cc[j].cardInter;
				}
			}
		}

		for(j=0;j<nbE;++j)
		{
			c[j].cardUnion = nodeArea[i] + stat.regionGTSize[j] - c[j].cardInter;
			c[j].score = regionScoringFunction(&(c[j]), nodeArea[i], stat.regionGTSize[j]);
		}

		if(additif)
		{
			float f=0;
			for(j=0;j<nbE;++j)
			{
				f += c[j].score;
			}
			stat.bestMatch[i] = 255; // invalid !
			stat.bestMatchingScore[i].cardInter = -1; // invalid !
			stat.bestMatchingScore[i].cardUnion = -1; // invalid !
			stat.bestMatchingScore[i].score = f; 
		}else{
		// compute  best match
			float maxScore = FLT_MIN;
			int indexMax=-1;
			for(j=0;j<nbE;++j)
			{
				
				if(c[j].score > maxScore )
				{
					maxScore = c[j].score;
					indexMax = j;
				}
			}
			stat.bestMatch[i] = indexMax;
			stat.bestMatchingScore[i] = c[indexMax]; 
		}
		//printd("Best match for %d(p:%d, area:%d) : %d(area:%d) score %d\n", i, CT->tabnodes[i].father, stat.regionTreeSize[i], indexMax, stat.regionGTSize[indexMax], c[indexMax].cardInter);
		
		//compute reverse scores
		for(j=0;j<nbE;++j)
			c[j].score = regionScoringFunction(&(c[j]), stat.regionGTSize[j], nodeArea[i]);
	}
	
	stat.detailedMatchingScore = matchingScore;
	return stat;
}









float scoreFusionReverseMatch(ReverseMatch m1, ReverseMatch m2, int nbR)
{

	int i =0;
	float score = 0;
	for(i=0;i<nbR;++i)
	{
		score += maxf(m1.scores[i],m2.scores[i]);
	}
	return score;
}



float computeFusionReverseMatch(ReverseMatch m1, ReverseMatch m2, ReverseMatch * res, int nbR)
{
	int i=0;
	float score = 0;

	if(res->scores == NULL)
		res->scores = myCalloc(nbR,sizeof(float));

	for(i=0;i<nbR;++i)
	{
		float v = maxf(m1.scores[i],m2.scores[i]);
		score += v;
		res->scores[i] = v;
	}
	return score;

}




/****************************************************************
Scoring functions

*****************************************************************/


///  META function

void initDynamicDataFromNodeAdditive(JCctree * CT, MatchingStat * matchingStat, double * level, int i, DynamicData *d, bool segmentationToGroundtruth, bool groundtruthToSegmentation)
{
	int p = CT->tabnodes[i].father;
	int refNode = (p!=-1)?p:i;
	
	d->minLevel = d->maxLevel = level[refNode];
	d->minRank = d->maxRank = refNode;



	int j;
	float reverseScore=0;
	float directScore=0;
	float k=0;
	if(segmentationToGroundtruth) 
	{
		RegionMatching * detailedMatch = matchingStat->detailedMatchingScore[i];
		int nbR =  matchingStat->nbRegionInGT;
		d->reverseMatchData.scores = myCalloc(nbR,sizeof(float));
		k+=1;
		for(j=0;j<nbR;++j)
		{
			d->reverseMatchData.scores[j] = detailedMatch[j].score;
			reverseScore += detailedMatch[j].score;
		}
	}
	
	if(groundtruthToSegmentation)
	{
		k+=1;
		directScore=matchingStat->bestMatchingScore[i].score;
	}
	
	
	d->directMatchScore = directScore;
	d->reverseMatchScore = reverseScore;
	d->score = 1.0 - (directScore + reverseScore)/ (k*matchingStat->nbPixels); 

}



float scoreFusionFunctionAdditive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl, bool segmentationToGroundtruth, bool groundtruthToSegmentation)
{
	float directMatchScore = 0; 
	float reverseMatchingScore = 0;
	float kk = 0;
	if(segmentationToGroundtruth)
	{
		reverseMatchingScore = scoreFusionReverseMatch(dk->reverseMatchData, dl->reverseMatchData, matchingStat->nbRegionInGT);
		kk+=1;
	}
	if(groundtruthToSegmentation)
	{
		directMatchScore = dk->directMatchScore + dl->directMatchScore;
		kk+=1;
	}
	
	return 1.0 - (directMatchScore + reverseMatchingScore)/(kk*matchingStat->nbPixels);
}

void fusionFunctionAdditive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl , DynamicData * res, float fusionScore, bool segmentationToGroundtruth, bool groundtruthToSegmentation)
{
	res->directMatchScore  = 0; 
	res->reverseMatchScore = 0;
	if(segmentationToGroundtruth)
	{
		res->reverseMatchScore = computeFusionReverseMatch(dk->reverseMatchData, dl->reverseMatchData, &(res->reverseMatchData), matchingStat->nbRegionInGT);
	}

	if(groundtruthToSegmentation)
	{
		res->directMatchScore = dk->directMatchScore + dl->directMatchScore;
	}

	res->minLevel = minf(dk->minLevel, dl->minLevel);
	res->maxLevel = maxd(dk->maxLevel, dl->maxLevel);
	res->minRank = mini(dk->minRank, dl->minRank);
	res->maxRank = maxi(dk->maxRank, dl->maxRank);
	res->score= fusionScore;
}


// specialization

float regionScoringFunctionIntersection(RegionMatching * r, int areaRegion1, int areaRegion2){
	return (float)r->cardInter;
}

float regionScoringFunctionIntersectionOverUnion(RegionMatching * r, int areaRegion1, int areaRegion2){
	return areaRegion1*r->cardInter/(float)r->cardUnion;	
}

float regionScoringFunctionBCE(RegionMatching * r, int areaRegion1, int areaRegion2){
	return (float)r->cardInter*minf((float)r->cardInter/(float)areaRegion1, (float)r->cardInter/(float)areaRegion2);	
}

// SegToGT

void initDynamicDataFromNodeDirectional1Additive(JCctree * CT, MatchingStat * matchingStat, double * level, int i, DynamicData *d)
{
	initDynamicDataFromNodeAdditive(CT, matchingStat, level, i, d, false, true);
}

float scoreFusionFunctionDirectional1Additive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl )
{
	return scoreFusionFunctionAdditive(CT, matchingStat, level, k, dk, l, dl, false, true);
}

void fusionFunctionDirectional1Additive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl , DynamicData * res, float fusionScore)
{
	fusionFunctionAdditive(CT, matchingStat, level, k, dk, l, dl, res, fusionScore, false, true);
}

// GTToSeg

void initDynamicDataFromNodeDirectional2Additive(JCctree * CT, MatchingStat * matchingStat, double * level, int i, DynamicData *d)
{
	initDynamicDataFromNodeAdditive(CT, matchingStat, level, i, d, true, false);
}

float scoreFusionFunctionDirectional2Additive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl )
{
	return scoreFusionFunctionAdditive(CT, matchingStat, level, k, dk, l, dl, true, false);
}

void fusionFunctionDirectional2Additive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl , DynamicData * res, float fusionScore)
{
	fusionFunctionAdditive(CT, matchingStat, level, k, dk, l, dl, res, fusionScore, true, false);
}

// Symmetric 

void initDynamicDataFromNodeSymmetricAdditive(JCctree * CT, MatchingStat * matchingStat, double * level, int i, DynamicData *d)
{
	initDynamicDataFromNodeAdditive(CT, matchingStat, level, i, d, true, true);
}

float scoreFusionFunctionSymmetricAdditive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl )
{
	return scoreFusionFunctionAdditive(CT, matchingStat, level, k, dk, l, dl, true, true);
}

void fusionFunctionSymmetricAdditive(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl , DynamicData * res, float fusionScore)
{
	fusionFunctionAdditive(CT, matchingStat, level, k, dk, l, dl, res, fusionScore, true, true);
}

///////////////




void freeDynamicData(DynamicData *d)
{	
	if(d->reverseMatchData.scores!=NULL)
		free(d->reverseMatchData.scores);
	free(d);
}




void printReverseMatchData(ReverseMatch r, int nbR)
{
	int i;
	printf("[");
	for(i=0; i< nbR; i++)
	{
		printf("%f", r.scores[i]);
		if(i!=nbR-1)
			printf(", ");
	}
	printf("]\n");
}


KCurve computeKScoreDynamique(JCctree * CT, double * level, grayImage * image, int maxK, int * nodeArea, Scoring scoring)
{

	int i= 0;
	int k,l;
	JCctreenode *tabnodes = CT->tabnodes;
	printd("Analysing GT\n");
	GroundTruth gt = analyseGT(image);
	printd("Computing pre-match\n");
	MatchingStat matchingStat = computeMatchingStat(CT, gt, nodeArea, scoring.regionScoringFunction, scoring.additif);	
	printd("Starting opt\n");
	KCurve * kCurves = myCalloc(CT->nbnodes, sizeof(KCurve));

	for(i=0; i <= CT->nbnodes/2; ++i)
	{
		DynamicData * d = myCalloc(1, sizeof(DynamicData));
		scoring.initFunction(CT, &matchingStat, level, i, d);
		kCurves[i].nbElement = 1;
		kCurves[i].data = d;
	}
	printd("Leaves processed\n");
	for(i=CT->nbnodes/2+1 ; i<CT->nbnodes; ++i)
	{	
		JCsoncell * soncell = tabnodes[i].sonlist;
		int child1 = soncell->son;
		int child2 = soncell->next->son;
		KCurve c1 = kCurves[child1];
		KCurve c2 = kCurves[child2];

		int nbSum = c1.nbElement + c2.nbElement; // unlimited combinations
		int nbElement = mini(nbSum,maxK);	
		
		DynamicData * d = myCalloc(nbElement, sizeof(DynamicData));
		kCurves[i].nbElement = nbElement;
		kCurves[i].data = d;
		for(k=0; k < nbElement; ++k)
			d[k].score = FLT_MAX ;

		// score for 1 region
		scoring.initFunction(CT, &matchingStat, level, i,&(d[0]));//just to be explicit
		/*if(i>=144)
		{
			printf("Node %d c1 %d c2 %d nbe %d score %f\n",i, child1, child2,nbElement, d[0].score);
			printReverseMatchData(d[0].reverseMatchData, matchingStat.nbRegionInGT);
		}*/
		int maxk = mini(c1.nbElement,maxK); // maximum number of regions usable for child 1
		for(k=0; k < maxk; ++k)
		{			
			DynamicData dk = c1.data[k];
			int maxl = mini(nbElement - k - 1,c2.nbElement); // maximum number of regions usable for child 2
			for( l=0; l < maxl; ++l)
			{			
				DynamicData dl = c2.data[l];	
				int ind=k+l+1; 
				float fusionScore = scoring.scoreFusionFunction(CT, &matchingStat, level, k, &dk, l, &dl);
				if(fusionScore < d[ind].score)
				{
					scoring.fusionFunction(CT, &matchingStat, level, k, &dk, l, &dl , &(d[ind]), fusionScore);
					/*if(i>=144)
						{
int nbR = matchingStat.nbRegionInGT;
printf("best scoring found %d %d %d %f(%f,%f)\n",ind+1,k+1,l+1,fusionScore,d[ind].directMatchScore,d[ind].reverseMatchScore);
printReverseMatchData(dk.reverseMatchData, nbR);
printReverseMatchData(dl.reverseMatchData, nbR);
printReverseMatchData(d[ind].reverseMatchData, nbR);
}*/
				}
			}
		}
		freeDynamicData(c1.data);
		freeDynamicData(c2.data);
	}

	KCurve result = kCurves[CT->root];


	free(kCurves);
	freeGT(gt);
	freeMatchingStat(matchingStat);
	return result;
}


int32_t main(int32_t argc, char ** argv) 
{
	graphe *g;
	grayImage * image;
	
	
	double *Av;
	int32_t i;

	if( (argc != 5)  ){
	usage(argv[0]);
	exit(1);
	}
	fprintf(stderr,"%s info: starting %s Measure %s\n", argv[0], argv[1], argv[3]);
	g = ReadGraphe(argv[1],&Av);
	if(g == NULL)
		exit(1);
	image = readGrayImage(argv[2]);
	
	Scoring scoreDirectional1Hamming = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeDirectional1Additive, 
		.scoreFusionFunction = &scoreFusionFunctionDirectional1Additive, 
		.fusionFunction = &fusionFunctionDirectional1Additive,
		.regionScoringFunction = &regionScoringFunctionIntersection,
		.additif = false
	};

	Scoring scoreDirectional2Hamming = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeDirectional2Additive, 
		.scoreFusionFunction = &scoreFusionFunctionDirectional2Additive, 
		.fusionFunction = &fusionFunctionDirectional2Additive,
		.regionScoringFunction = &regionScoringFunctionIntersection,
		.additif = false
	};

	Scoring scoreVanDongen = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeSymmetricAdditive,
		.scoreFusionFunction = &scoreFusionFunctionSymmetricAdditive,
		.fusionFunction = &fusionFunctionSymmetricAdditive,
		.regionScoringFunction = &regionScoringFunctionIntersection,
		.additif = false
	};

	Scoring scoreDirectional1Covering = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeDirectional1Additive, 
		.scoreFusionFunction = &scoreFusionFunctionDirectional1Additive, 
		.fusionFunction = &fusionFunctionDirectional1Additive,
		.regionScoringFunction = &regionScoringFunctionIntersectionOverUnion,
		.additif = false
	};

	Scoring scoreDirectional2Covering = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeDirectional2Additive, 
		.scoreFusionFunction = &scoreFusionFunctionDirectional2Additive, 
		.fusionFunction = &fusionFunctionDirectional2Additive,
		.regionScoringFunction = &regionScoringFunctionIntersectionOverUnion,
		.additif = false
	};

	Scoring scoreSymmetricCovering = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeSymmetricAdditive,
		.scoreFusionFunction = &scoreFusionFunctionSymmetricAdditive,
		.fusionFunction = &fusionFunctionSymmetricAdditive,
		.regionScoringFunction = &regionScoringFunctionIntersectionOverUnion,
		.additif = false
	};


	Scoring scoreBidirectionnalConsistency = (Scoring){ 
		.initFunction = &initDynamicDataFromNodeDirectional1Additive, 
		.scoreFusionFunction = &scoreFusionFunctionDirectional1Additive, 
		.fusionFunction = &fusionFunctionDirectional1Additive,
		.regionScoringFunction = &regionScoringFunctionBCE,
		.additif = true
	};

	
	Scoring scoring = scoreSymmetricCovering;//scoreVanDongen;//;

	
	if(strcmp(argv[3],"DHamming1")==0)
	{ scoring = scoreDirectional1Hamming; }
	else if(strcmp(argv[3],"DHamming2")==0)
	{ scoring = scoreDirectional2Hamming; }
	else if(strcmp(argv[3],"SHamming")==0)
	{ scoring = scoreVanDongen; }
	else if(strcmp(argv[3],"DCovering1")==0)
	{ scoring = scoreDirectional1Covering; }
	else if(strcmp(argv[3],"DCovering2")==0)
	{ scoring = scoreDirectional2Covering; }
	else if(strcmp(argv[3],"SCovering")==0)
	{ scoring = scoreSymmetricCovering; }
	else if(strcmp(argv[3],"BCE")==0)
	{ scoring = scoreBidirectionnalConsistency; }
	else{
		fprintf(stderr,"%s info: invalid error measure '%s'\n", argv[0], argv[3]);
		exit(1);
	}
	//'FMeasure','Jaccard','Error'



	char * stop;
	int maxK = (int)strtol(argv[4], &stop, 10);
	if(stop[0]!='\0' || maxK <0)
	{
		fprintf(stderr, "third argument must be a positive integer (base 10)\n");
		exit(1);
	}

	double *Val = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
  	saliencyTree(g,&CT,Val);
	int * nodeArea = computeNodeArea(CT);
	
	KCurve curve = computeKScoreDynamique(CT, Val, image, maxK, nodeArea, scoring);
	printf("#K\tscore\tminLevel\tmaxLevel\tminRank\tmaxRank\n");
	for(i=0; i<curve.nbElement; ++i)
	{
		printf("%d\t%f\t%f\t%f\t%d\t%d\n",i+1,curve.data[i].score,curve.data[i].minLevel,curve.data[i].maxLevel,curve.data[i].minRank,curve.data[i].maxRank);
	}
	freeDynamicData(curve.data);
	free(nodeArea);
	freeimage(image);
	terminegraphe(g);
	free(Av);
	free(Val); 
	componentTreeFree(CT);
	return 0;
}
