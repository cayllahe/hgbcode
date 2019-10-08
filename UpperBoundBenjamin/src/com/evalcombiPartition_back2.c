
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


//#define DEBUG
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
  fprintf(stderr,"USAGE: %s saliency_map.graph  ground_truth.pgm scoreFunction maxk [-threshold]\n",arg);
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


	switch(datatype(Im)){
	case 	VFF_TYP_1_BYTE: 
		break;
	default : 
		fprintf(stderr,"%s: not yet implemented for this kind of data \n",filename);
		exit(1);
	}
	
	return Im;
}






int GnbRegionsGT=0;

void remapGT(grayImage * image)
{
	int count[256];
	int i;

	for(i=0;i<256;++i)
		count[i]=0;

	uint8_t *seg = UCHARDATA(image);
	int width=rowsize(image), heigth=colsize(image), N=width*heigth;
	for(i=0;i<N;++i)
		++count[seg[i]];
	int num=0;
	for(i=0;i<256;++i)
	{
		if(count[i]!=0)
			count[i]=num++;
	}
	for(i=0;i<N;++i)
		seg[i] = count[seg[i]];
}


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
			if(indexMax == -1)
			{
				/*printf("nbe %d, i %d area %d\n", nbE, i, nodeArea[i]);
				for(j=0;j<nbE;++j)
 					printf("\ts %d : %f  inter %d union %d  area %d\n",j, c[j].score, c[j].cardInter, c[j].cardUnion, stat.regionGTSize[j] );*/
				printf("Failure : start crying now\n");
				exit(1);
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

/*
	1 region parition initialisation
*/

#define DeclareReverseInitDynamicDataFromNodeAdditive_true \
	RegionMatching * detailedMatch = matchingStat->detailedMatchingScore[i]; \
	int nbR =  matchingStat->nbRegionInGT; int j; \
	d->reverseMatchData.scores = myCalloc(nbR,sizeof(float)); \
	k+=1; \
	for(j=0;j<nbR;++j) \
	{ \
		d->reverseMatchData.scores[j] = detailedMatch[j].score; \
		reverseScore += detailedMatch[j].score; \
	}
#define DeclareReverseInitDynamicDataFromNodeAdditive_false 

#define DeclareDirectInitDynamicDataFromNodeAdditive_true \
	k+=1; \
	directScore=matchingStat->bestMatchingScore[i].score;
#define DeclareDirectInitDynamicDataFromNodeAdditive_false 

#define DeclareInitDynamicDataFromNodeAdditive(reverse,direct) void initDynamicDataFromNodeAdditive_##reverse##_##direct (JCctree * CT, MatchingStat * matchingStat, double * level, int i, DynamicData *d) \
{\
	int p = CT->tabnodes[i].father; int refNode = (p!=-1)?p:i; d->minRank = d->maxRank = refNode; \
	float reverseScore=0;	float directScore=0;	float k=0; \
	DeclareReverseInitDynamicDataFromNodeAdditive_##reverse \
	DeclareDirectInitDynamicDataFromNodeAdditive_##direct \
	d->directMatchScore = directScore; d->reverseMatchScore = reverseScore;	d->score = (directScore + reverseScore)/ (k*matchingStat->nbPixels); \
}

DeclareInitDynamicDataFromNodeAdditive(true,true)
DeclareInitDynamicDataFromNodeAdditive(true,false)
DeclareInitDynamicDataFromNodeAdditive(false,true)

#define CallInitDynamicDataFromNodeAdditive(reverse,direct) initDynamicDataFromNodeAdditive_##reverse##_##direct



/*
	Fusion scoring function
*/

#define DeclareReverseScoreFusionFunctionAdditive_true \
	reverseMatchingScore = scoreFusionReverseMatch(dk->reverseMatchData, dl->reverseMatchData, matchingStat->nbRegionInGT); \
	kk+=1;
#define DeclareReverseScoreFusionFunctionAdditive_false 

#define DeclareDirectScoreFusionFunctionAdditive_true \
	directMatchScore = dk->directMatchScore + dl->directMatchScore; \
	kk+=1;
#define DeclareDirectScoreFusionFunctionAdditive_false 

#define DeclareScoreFusionFunctionAdditive(reverse,direct) float scoreFusionFunctionAdditive_##reverse##_##direct (JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl) { \
	float directMatchScore = 0; float reverseMatchingScore = 0; float kk = 0; \
	DeclareReverseScoreFusionFunctionAdditive_##reverse \
	DeclareDirectScoreFusionFunctionAdditive_##direct \
	return (directMatchScore + reverseMatchingScore)/(kk*matchingStat->nbPixels); \
}



DeclareScoreFusionFunctionAdditive(true,false)
DeclareScoreFusionFunctionAdditive(false,true)
DeclareScoreFusionFunctionAdditive(true,true)

#define CallScoreFusionFunctionAdditive(reverse,direct) scoreFusionFunctionAdditive_##reverse##_##direct


/*
	Fusion function
*/

#define DeclareReverseFusionFunctionAdditive_true \
	res->reverseMatchScore = computeFusionReverseMatch(dk->reverseMatchData, dl->reverseMatchData, &(res->reverseMatchData), matchingStat->nbRegionInGT);
#define DeclareReverseFusionFunctionAdditive_false

#define DeclareDirectFusionFunctionAdditive_true \
	res->directMatchScore = dk->directMatchScore + dl->directMatchScore;
#define DeclareDirectFusionFunctionAdditive_false

#define DeclareFusionFunctionAdditive(reverse,direct) void fusionFunctionAdditive_##reverse##_##direct (JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl , DynamicData * res, float fusionScore) { \
	res->directMatchScore  = 0; 	res->reverseMatchScore = 0; \
	DeclareReverseFusionFunctionAdditive_##reverse \
	DeclareDirectFusionFunctionAdditive_##direct \
	res->minRank = mini(dk->minRank, dl->minRank);	res->maxRank = maxi(dk->maxRank, dl->maxRank);	res->score= fusionScore; \
}


DeclareFusionFunctionAdditive(true,true)
DeclareFusionFunctionAdditive(true,false)
DeclareFusionFunctionAdditive(false,true)

#define CallFusionFunctionAdditive(reverse,direct) fusionFunctionAdditive_##reverse##_##direct



// specialization

float regionScoringFunctionIntersection(RegionMatching * r, int areaRegion1, int areaRegion2){
	return (float)r->cardInter;
}

float regionScoringFunctionIntersectionOverUnion(RegionMatching * r, int areaRegion1, int areaRegion2){
	return (float)areaRegion1*(float)r->cardInter/(float)r->cardUnion;	
}

float regionScoringFunctionBCE(RegionMatching * r, int areaRegion1, int areaRegion2){
	return (float)r->cardInter*minf((float)r->cardInter/(float)areaRegion1, (float)r->cardInter/(float)areaRegion2);	
}

///////////////////////////////////////




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

/** -------------------------------------------------------------------------
Non horizontal cut optimization !
*/

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
			d[k].score = 0 ;

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
				if(fusionScore > d[ind].score)
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


/**
Horizontal cut optimization

*/

// find the highest node with a level strictly lower than curLevel (starting from curNode)
int findNextNode(JCctree *CT, double * level, int curNode, double curLevel)
{
	int nbPix = CT->nbnodes/2;
	while(level[curNode] >= curLevel && curNode > nbPix)
		--curNode;
	return curNode;
}


// init a Partition structure corresponding to the given threshold
int initPartition(JCctree *CT, double * level, double threshold, int curNode, Partition* partition,  int maxK)
{
	JCctreenode *tabnodes = CT->tabnodes;
	
	int count=0;
	int i;
	for(i=curNode+1; i < CT->nbnodes; ++i)
	{
		//int p=tabnodes[i].father;
		JCsoncell * soncell = tabnodes[i].sonlist;
		int child1 = soncell->son;
		int child2 = soncell->next->son;

		if(level[child1] < threshold)
		{
			if (count == maxK)
				return -1;
			partition->nodeNumbers[count++]=child1;
			
		}
		if(level[child2] < threshold)
		{
			if (count == maxK)
				return -1;
			partition->nodeNumbers[count++]=child2;
		}
		
	}

	partition->nbElement = count;
	return count;
}



#define DeclareReverseScorePartition_true \
	int j;	k+=1.0;	int nbRGT = matchingStat->nbRegionInGT;	RegionMatching ** reverseMatch = matchingStat->detailedMatchingScore;	float reverseScores[256]; \
	for(i=0;i<nbRGT;++i) \
		reverseScores[i]=0; \
	for(j=0; j<nbRP; ++j) 	{ \
		for(i=0;i<nbRGT;++i)	{ \
			if (reverseScores[i] < reverseMatch[partition->nodeNumbers[j]][i].score) \
				reverseScores[i]=reverseMatch[partition->nodeNumbers[j]][i].score; \
		}	}\
	for(i=0;i<nbRGT;++i) \
		reverseMatchingScore += reverseScores[i];
#define DeclareReverseScorePartition_false

#define DeclareDirectScorePartition_true \
	k+=1.0; \
	RegionMatching * directMatch = matchingStat->bestMatchingScore; \
	for(i=0; i<nbRP; ++i) \
		directMatchScore += directMatch[partition->nodeNumbers[i]].score;
#define DeclareDirectScorePartition_false

#define DeclareScorePartition(reverse,direct) float scorePartition_##reverse##_##direct (Partition * partition, MatchingStat * matchingStat) { \
	float directMatchScore = 0; float reverseMatchingScore = 0; float k = 0; \
	int nbRP = partition->nbElement, i;\
	DeclareReverseScorePartition_##reverse \
	DeclareDirectScorePartition_##direct \
	return (directMatchScore + reverseMatchingScore)/(k*matchingStat->nbPixels);  \
}


DeclareScorePartition(true,true)
DeclareScorePartition(true,false)
DeclareScorePartition(false,true)

#define CallScorePartition(reverse,direct) scorePartition_##reverse##_##direct


/*
float scorePartition(Partition * partition, MatchingStat * matchingStat)
{
	float directMatchScore = 0; float reverseMatchingScore = 0; float k = 0; 
	int nbRP = partition->nbElement, i;
	
	


	// direct match
	k+=1.0;
	RegionMatching * directMatch = matchingStat->bestMatchingScore;
	for(i=0; i<nbRP; ++i)
		directMatchScore += directMatch[partition->nodeNumbers[i]].score;



	int j;
	k+=1.0;
	int nbRGT = matchingStat->nbRegionInGT;
	RegionMatching ** reverseMatch = matchingStat->detailedMatchingScore;
	float reverseScores[256];
	for(i=0;i<nbRGT;++i)
		reverseScores[i]=0;
	for(j=0; j<nbRP; ++j)
	{
		for(i=0;i<nbRGT;++i)
		{
			if (reverseScores[i] < reverseMatch[partition->nodeNumbers[j]][i].score)
				reverseScores[i]=reverseMatch[partition->nodeNumbers[j]][i].score;
		}
	}
	for(i=0;i<nbRGT;++i)
		reverseMatchingScore += reverseScores[i];

	return (directMatchScore + reverseMatchingScore)/(k*matchingStat->nbPixels); 
}*/

ThresholdKCurve computeScoreCurve(JCctree *CT, double * level, grayImage * image, int maxK, int * nodeArea, Scoring  scoring)
{
	int curNode = CT->root;
	double curLevel = level[curNode];
	int nbPixels = CT->nbnodes/2;
	Partition partition;
	partition.nodeNumbers = myCalloc(maxK, sizeof(int));
	int i;
	

	ThresholdKCurve scores;
	scores.score = myCalloc(maxK, sizeof(float));
	scores.threshold = myCalloc(maxK, sizeof(double));
	scores.nbElement=maxK;
	for(i=0; i<maxK; ++i)
		scores.score[i] =-1 ;

	printd("Analysing GT\n");
	GroundTruth gt = analyseGT(image);
	printd("Computing pre-match\n");
	MatchingStat matchingStat = computeMatchingStat(CT, gt, nodeArea, scoring.regionScoringFunction, scoring.additif);	
	printd("Starting opt\n");

	int k=1;
	partition.nbElement = 1;
	partition.nodeNumbers[0] = CT->root;
	
	do {
		
		--k;
		//printf("-thrshold %f node %d/%d  k: %d\n", curLevel,curNode, CT->nbnodes, k);
		float score = scoring.scorePartitionFunction(&partition, &matchingStat);
		if (scores.score[k] < score)
		{
			scores.score[k] = score;	
			scores.threshold[k] = curLevel;
		}
		curNode = findNextNode(CT, level, curNode, curLevel);
		if(curNode <= nbPixels)
			break;
		
		k=initPartition(CT, level, curLevel, curNode, &partition, maxK);
		curLevel = level[curNode];
			
	} while( k!=-1 && curNode > nbPixels );

	for(i=1; i<maxK; ++i)
		if(scores.score[i]==-1)
		{
			scores.score[i] = scores.score[i-1];	
			scores.threshold[i] = scores.threshold[i-1];
		}

	free(partition.nodeNumbers);
	return scores;
}

/** -------------------------------------
Main 
*/


int32_t main(int32_t argc, char ** argv) 
{
	graphe *g;
	grayImage * image;
	
	
	double *Av;
	int32_t i;

	if( (argc != 5 && argc != 6)  ){
	usage(argv[0]);
	exit(1);
	}
	fprintf(stderr,"%s info: starting %s Measure %s\n", argv[0], argv[1], argv[3]);
	g = ReadGraphe(argv[1],&Av);
	if(g == NULL)
		exit(1);
	image = readGrayImage(argv[2]);
	
	Scoring scoreDirectional1Hamming = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(false,true), //initDynamicDataFromNodeDirectional1Additive, 
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(false,true) , //scoreFusionFunctionDirectional1Additive, 
		.fusionFunction = &CallFusionFunctionAdditive(false,true), //fusionFunctionDirectional1Additive,
		.regionScoringFunction = &regionScoringFunctionIntersection,
		.scorePartitionFunction = &CallScorePartition(false, true),
		.additif = false
	};

	Scoring scoreDirectional2Hamming = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(true,false), //initDynamicDataFromNodeDirectional2Additive, 
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(true,false) , //scoreFusionFunctionDirectional2Additive, 
		.fusionFunction = &CallFusionFunctionAdditive(true,false), //fusionFunctionDirectional2Additive,
		.regionScoringFunction = &regionScoringFunctionIntersection,
		.scorePartitionFunction = &CallScorePartition(true, false),
		.additif = false
	};

	Scoring scoreVanDongen = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(true,true), //initDynamicDataFromNodeSymmetricAdditive,
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(true,true) , //scoreFusionFunctionSymmetricAdditive,
		.fusionFunction = &CallFusionFunctionAdditive(true,true), //fusionFunctionSymmetricAdditive,
		.regionScoringFunction = &regionScoringFunctionIntersection,
		.scorePartitionFunction = &CallScorePartition(true, true),
		.additif = false
	};

	Scoring scoreDirectional1Covering = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(false,true), //initDynamicDataFromNodeDirectional1Additive, 
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(false,true) , //scoreFusionFunctionDirectional1Additive, 
		.fusionFunction = &CallFusionFunctionAdditive(false,true), //fusionFunctionDirectional1Additive,
		.regionScoringFunction = &regionScoringFunctionIntersectionOverUnion,
		.scorePartitionFunction = &CallScorePartition(false, true),
		.additif = false
	};

	Scoring scoreDirectional2Covering = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(true,false), //initDynamicDataFromNodeDirectional2Additive, 
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(true,false) , //scoreFusionFunctionDirectional2Additive, 
		.fusionFunction = &CallFusionFunctionAdditive(true,false), //fusionFunctionDirectional2Additive,
		.regionScoringFunction = &regionScoringFunctionIntersectionOverUnion,
		.scorePartitionFunction = &CallScorePartition(true, false),
		.additif = false
	};

	Scoring scoreSymmetricCovering = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(true,true), //initDynamicDataFromNodeSymmetricAdditive,
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(true,true) , //scoreFusionFunctionSymmetricAdditive,
		.fusionFunction = &CallFusionFunctionAdditive(true,true), //fusionFunctionSymmetricAdditive,
		.regionScoringFunction = &regionScoringFunctionIntersectionOverUnion,
		.scorePartitionFunction = &CallScorePartition(true, true),
		.additif = false
	};


	Scoring scoreBidirectionnalConsistency = (Scoring){ 
		.initFunction = &CallInitDynamicDataFromNodeAdditive(false,true), //initDynamicDataFromNodeDirectional1Additive, 
		.scoreFusionFunction = &CallScoreFusionFunctionAdditive(false,true) , //scoreFusionFunctionDirectional1Additive, 
		.fusionFunction = &CallFusionFunctionAdditive(false,true), //fusionFunctionDirectional1Additive,
		.regionScoringFunction = &regionScoringFunctionBCE,
		.scorePartitionFunction = &CallScorePartition(false, true),
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
		fprintf(stderr,"%s error: invalid measure '%s'\n", argv[0], argv[3]);
		exit(1);
	}
	//'FMeasure','Jaccard','Error'
	bool horizontalCut=false;
	if (argc==6)
	{
		if(strcmp(argv[5],"-threshold")==0)
		{ horizontalCut=true; }
		else{
			fprintf(stderr,"%s error: invalid option '%s'\n", argv[0], argv[5]);
			exit(1);
		}
	}
	char * stop;
	int maxK = (int)strtol(argv[4], &stop, 10);
	if(stop[0]!='\0' || maxK <0)
	{
		fprintf(stderr,"%s error: third argument must be a positive integer (base 10) %s\n", argv[0], argv[4]);
		exit(1);
	}

	double *level = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
  	saliencyTree(g,&CT,level);
	int * nodeArea = computeNodeArea(CT);
	remapGT(image);
	
	if(horizontalCut)
	{
		ThresholdKCurve curve = computeScoreCurve(CT, level, image, maxK, nodeArea, scoring);
		printf("#K\tscore\tthreshold\n");
		for(i=0; i<curve.nbElement; ++i)
		{
			printf("%d\t%f\t%f\n",i+1,curve.score[i],curve.threshold[i]);
		}
		free(curve.score);
		free(curve.threshold);
	}	
	else{
		KCurve curve = computeKScoreDynamique(CT, level, image, maxK, nodeArea, scoring);
		printf("#K\tscore\tminRank\tmaxRank\n");
		for(i=0; i<curve.nbElement; ++i)
		{
			printf("%d\t%f\t%d\t%d\n",i+1,curve.data[i].score,curve.data[i].minRank,curve.data[i].maxRank);
		}
		freeDynamicData(curve.data);
	}
	
	free(nodeArea);
	freeimage(image);
	terminegraphe(g);
	free(Av);
	free(level); 
	componentTreeFree(CT);
	return 0;
}
