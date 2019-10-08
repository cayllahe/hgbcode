
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "bptio.h"


#include <mccodimage.h>
#include <mcimage.h>
#include <lhierarchie.h>
#include <string.h>
#include "llca.h"
#include "evalcombiPartition.h"


//#define DEBUG
#ifdef DEBUG
#define printd(...) (printf(__VA_ARGS__))
#else
#define printd(...) {;}
#endif


/*
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
}*/

void usage(char *arg){
    fprintf(stderr,"#################################################################\n\n");
    fprintf(stderr,"USAGE: %s saliency_map.graph  -i inputFile.pgm -gt groundTruthFile -o outputFile\n",arg);
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"\t '-gt groundTruthFile.pgm' may be used several time\n");
    fprintf(stderr,"\t '-maxK N' with N a postitive integer: the maximum number of regions in the parition (default 100)\n");
    fprintf(stderr,"\t '-score scoreFunction' with scoreFunction in {'DHamming1', 'DHamming2', 'SHamming', 'DCovering1', 'DCovering2'; 'SCovering','BCE'} (default BCE)\n");
    fprintf(stderr,"\t '-action scoreFunction' with action in \
                    \n\t\t'upperBound' : compute the optimal non horizontal cuts (default);\
                    \n\t\t'threshold': compute the optimal horizontal cuts;\
                    \n\t\t'draw': draw the optimal non horizontal cuts from 2 to maxK regions in Khalimsky space (input graph must represent a 4-adjacency).\
                    \n\t\t'drawOne': draw one optimal non horizontal cut with maxK regions in Khalimsky space (input graph must represent a 4-adjacency)\n");   
    fprintf(stderr,"#################################################################\n\n");
}



// read an image, check if it is a byte 2d image
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
/*
insure that region labeling start at 0 and ends at N-1 (with N the number of regions)
*/
int remapGT(grayImage * image)
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
    return num;
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
	//printd("Nb regions in ground truth: %d\n",nbElement);
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


/*
Compute matching statistics between every tree node and gt regions
*/
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

		if(additif) // additif scoring
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
		}else{ // max scoring
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
			{ // simple integrity check
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

/** ****************************************************************************************
///  META function 
// some hugly macros for specialized init/scoring/fusion function generation
*/
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
	d->minRank = d->maxRank = i; \
    d->backTrackK = d->backTrackL = 0; \
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
    res->backTrackK = k+1; res->backTrackL = l+1; \
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
	{
        free(d->reverseMatchData.scores);
        d->reverseMatchData.scores=NULL;
    }
}

void freeKCurve( KCurve * curve, int nbNodes)
{
    int i;

    for(i=0; i<nbNodes;++i)
    {
        free(curve[i].data);
    }
    free(curve);
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

/** ****************************************************************************************
Non horizontal cut optimization !
*/

KCurve * computeKScoreDynamique(JCctree * CT, double * level, grayImage * image, int maxK, int * nodeArea, Scoring scoring)
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
				}
			}
		}

        freeDynamicData(c1.data);
        freeDynamicData(c2.data);

	}

	//KCurve result = kCurves[CT->root];

    freeDynamicData(kCurves[CT->root].data);
	//free(kCurves);
	freeGT(gt);
	freeMatchingStat(matchingStat);
	return kCurves;
}


/** ****************************************************************************************

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

    /*if (scores.score[maxK-1] == -1)
    {
        scores.score[maxK-1] = 0;	
	    scores.threshold[maxK-1] = 0;
    }
	for(i=maxK-2; i>=0; --i)
		if(scores.score[i]==-1)
		{
			scores.score[i] = scores.score[i+1];	
			scores.threshold[i] = scores.threshold[i+1];
		}
    */
	free(partition.nodeNumbers);
	return scores;
}

/** ****************************************************************************************
    Straightening 
*/

void backTrack2(JCctree * CT, KCurve * kCurves, int node, int numNodes, double kval, double * attr)
{

    JCctreenode *tabnodes = CT->tabnodes;
    int p=tabnodes[node].father;
    if(p!=-1) attr[p] += kval;
    //attr[node]+=1;//kval;
    /*if (node<=CT->nbnodes/2)
    {
        attr[node] = maxd(attr[node],kval);
        return;
    }*/
    if (node>CT->nbnodes/2)
    {
       
        JCsoncell * soncell = tabnodes[node].sonlist;
        int child1 = soncell->son;
        int child2 = soncell->next->son;
        KCurve curve = kCurves[node];
        DynamicData d = curve.data[numNodes-1];
        int btk = d.backTrackK;
        int btl = d.backTrackL;
        if(btk == 0 && btl == 0)
        {
           //attr[node] = maxd(attr[node],kval);
        }else{
        
            backTrack2(CT, kCurves, child1, btk, kval, attr);
            backTrack2(CT, kCurves, child2, btl, kval, attr);
        }
    }
}

void backTrack(JCctree * CT, KCurve * kCurves, int node, int numNodes, double kval, double * attr)
{
    if (node<=CT->nbnodes/2)
    {
        attr[node] = maxd(attr[node],kval);
    }
    else 
    {
        JCctreenode *tabnodes = CT->tabnodes;
        JCsoncell * soncell = tabnodes[node].sonlist;
        int child1 = soncell->son;
        int child2 = soncell->next->son;
        KCurve curve = kCurves[node];
        DynamicData d = curve.data[numNodes-1];
        int btk = d.backTrackK;
        int btl = d.backTrackL;
        if(btk == 0 && btl == 0)
        {
           attr[node] = maxd(attr[node],kval);
        }else{
            backTrack(CT, kCurves, child1, btk, kval, attr);
            backTrack(CT, kCurves, child2, btl, kval, attr);
        }
    }
}

bool testIncreasing(JCctree * CT, double * attr)
{
    int i;
    for(i=0;i<CT->root-1;++i)
    {
        if(attr[i] > attr[CT->tabnodes[i].father])
            return false;
    }
    return true;
}

double * straightenedAttribute(JCctree * CT, KCurve * kCurves, double * attr)
{
    //double * attr = myCalloc(CT->nbnodes, sizeof(double));
    KCurve opt = kCurves[CT->root];
    int maxK = opt.nbElement;
    int i;
    double prevScore = 0.0;
    for(i=1; i<=maxK; ++i)
    {
        double score =  kCurves[CT->root].data[i-1].score;
        //printf("backtracking i : %d\n",i);
        backTrack2(CT, kCurves, CT->root, i, maxd(score-prevScore,0), attr);
        prevScore=0;
    }
    return attr;
}

typedef struct SSaliencyCache {
    int32_t *Euler, *Depth, *Represent, *Number, **Minim;
} SaliencyCache;

void freeSaliencyCache(SaliencyCache ** saliencyCache)
{
    free((*saliencyCache)->Euler);
    free((*saliencyCache)->Represent);
    free((*saliencyCache)->Depth);
    free((*saliencyCache)->Number);
    free((*saliencyCache)->Minim[0]);
    free((*saliencyCache)->Minim);
    free(*saliencyCache);
    *saliencyCache = NULL;
}

void saliency(JCctree *ST, graphe *g, double * STAltitudes, SaliencyCache ** saliencyCache)
{
    int32_t logn, nbRepresent, u, x, y, c1;
    int32_t nbarcs = g->ind;
    int32_t *Euler, *Depth, *Represent, *Number, **Minim;
    if(*saliencyCache == NULL)
    {
        *saliencyCache = myCalloc(1, sizeof(SaliencyCache));
         /* Data structure for fast computation of LCA (least common ancestor) */
        (*saliencyCache)->Euler = Euler = myCalloc(2*ST->nbnodes-1, sizeof(int32_t));
        (*saliencyCache)->Represent  = Represent = myCalloc(2*ST->nbnodes-1, sizeof(int32_t));
        (*saliencyCache)->Depth = Depth = myCalloc(ST->nbnodes, sizeof(int32_t));
        (*saliencyCache)->Number = Number = myCalloc(ST->nbnodes, sizeof(int32_t));
        (*saliencyCache)->Minim = Minim = LCApreprocess(ST, Euler, Depth, Represent, Number, &nbRepresent, &logn);
    }else{
        Euler = (*saliencyCache)->Euler;
        Represent = (*saliencyCache)->Represent;
        Depth = (*saliencyCache)->Depth;
        Number = (*saliencyCache)->Number;
        Minim = (*saliencyCache)->Minim;
    }
 
    // For any edge of g
    for(u = 0; u < nbarcs; u++){
        x = g->tete[u]; y = g->queue[u];
        c1 = Represent[LowComAncFast(x, y, Euler, Number, Depth, Minim)];
        setweight(g,u,STAltitudes[c1]);
    }
}

void backTrackOptimKPart(JCctree * CT, KCurve * kCurves, int node, int numNodes, double * attr)
{
    JCctreenode *tabnodes = CT->tabnodes;
    int p=tabnodes[node].father;
    if(p!=-1) attr[p] = 255;

    if (node>CT->nbnodes/2)
    {
        JCsoncell * soncell = tabnodes[node].sonlist;
        int child1 = soncell->son;
        int child2 = soncell->next->son;
        KCurve curve = kCurves[node];
        DynamicData d = curve.data[numNodes-1];
        int btk = d.backTrackK;
        int btl = d.backTrackL;
        if(btk == 0 && btl == 0)
        {
          // attr[node] = 255;
        }else{
            backTrackOptimKPart(CT, kCurves, child1, btk, attr);
            backTrackOptimKPart(CT, kCurves, child2, btl, attr);
        }
    }
}

void graphToTimageKhalimski(graphe * g, grayImage * contour)
{
    int rs = g->rs;
    int cs = g->cs;
    int N = rs * cs;

    int rsIm = 2*rs - 1;
    //int csIm = 2*cs - 1;
    pcell p;
    uchar  * F = UCHARDATA(contour);
    
    int i, j, x, y;
    uchar tmp;

    for(x = 0; x < N; x++){
        for(p = g->listar[x]; p != NULL; p = p->next){
            y = g->queue[p->index];
            i = x%rs;
            j = x/rs;
            if(y == x) y = g->tete[p->index];
            if( (y - x) == 1)
            F[2*i+1 + (2*j)*rsIm] =  (uchar)getweight(g,p->index);
            if( (y - x) == rs)
            F[2*i + (2*j +1) * rsIm] = (uchar)getweight(g,p->index);
        }
    }

    for(i = 0; i < rs -1; i++)
        for(j = 0; j < cs - 1; j++)
        {
            tmp = maxi(F[2*i + (2*j+1)*rsIm], F[2*i + 1 + (2*j)*rsIm] );
            tmp = maxi(tmp,F[2*i+2 + (2*j+1)*rsIm]);
            F[2*i +1 + (2*j+1)*rsIm] = maxi(tmp, F[2*i+1 + (2*j + 2)*rsIm]);
        }
}

void extractOptimalKPartition(JCctree *ST, double * attr, KCurve * kCurves, int kVal, graphe * g, SaliencyCache ** saliencyCache, grayImage * contour)
{
    int i;
    for(i=0;i<ST->nbnodes;++i)
        attr[i]=0;
    backTrackOptimKPart(ST, kCurves, ST->root, kVal, attr);
    saliency(ST, g, attr,saliencyCache);
    graphToTimageKhalimski(g, contour);

}

void createPathName(char * buf, char * path, int i)
{
    strcpy(buf, path);
    sprintf(buf + strlen(buf),"_%d.pgm",i);

}

void generateKOptimalPartition(graphe * g, JCctree *ST, KCurve * kCurves, int kMax, char * path, bool allCuts){
    int i;
    SaliencyCache * saliencyCache = NULL;
    char  filename[1024];
    int rs = g->rs;
    int cs = g->cs;
    //int N = rs * cs;
 
    int rsIm = 2*rs - 1;
    int csIm = 2*cs - 1;
    grayImage * contour = allocimage(NULL, rsIm, csIm, 1, VFF_TYP_1_BYTE); 
    double * attr = myCalloc(ST->nbnodes, sizeof(double));
    if(allCuts)
    {
        for(i=2; i<=kMax; ++i)
        {
            extractOptimalKPartition(ST, attr, kCurves, i, g, &saliencyCache, contour);
            createPathName(filename, path, i);
            writeimage(contour,filename);
        }
    } else {
            extractOptimalKPartition(ST, attr, kCurves, kMax, g, &saliencyCache, contour);
            writeimage(contour,path);
    }
    
    free(attr);
    freeSaliencyCache(&saliencyCache);
    freeimage(contour);
    

}



 
/** ****************************************************************************************
Main 
*/

#define MAX_GT 32
typedef enum eOperation {UPPER_BOUND_OPTIM, THRESHOLD_OPTIM, DRAW_UPPER_BOUND_PARTITION, DRAW_ONE_UPPER_BOUND_PARTITION,  STRAIGHTENING} Operation;
int32_t main(int32_t argc, char ** argv) 
{
	graphe *g;
	grayImage * image;
	
	
	double *Av;
	int32_t i;
    

    
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

	
	Scoring scoring = scoreBidirectionnalConsistency;
    Operation operation = UPPER_BOUND_OPTIM;
    char * inFile=NULL;
    char * outFile=NULL;
    char * gtFiles[MAX_GT];
    int numGtFiles=0;
    int maxK=100;
    char * stop;


    int argNum;
    if(argc%2!=1) // options come in pairs
    {
        fprintf(stderr,"%s error: invalid number of arguments %d (cannot be odd)\n", argv[0], argc-1);
        usage(argv[0]);
        exit(1);
    }else{

        for(argNum=1; argNum<argc; argNum+=2)
        {
            
            if(strcmp(argv[argNum],"-i")==0)
            {
                inFile = argv[argNum+1];
            } else if(strcmp(argv[argNum],"-o")==0)
            {
                outFile = argv[argNum+1];
            } else if(strcmp(argv[argNum],"-gt")==0)
            {
                if (numGtFiles < MAX_GT )
                {
                    gtFiles[numGtFiles++] = argv[argNum+1];
                } else {
                    fprintf(stderr, "%s warning: too many ground truth provided (max=%d), remaining entries are ignored\n", argv[0], MAX_GT);
                }
            } else if(strcmp(argv[argNum],"-maxK")==0)
            {
                maxK = (int)strtol(argv[argNum+1], &stop, 10);
                if(stop[0]!='\0' || maxK <0)
                {
                    fprintf(stderr,"%s error: maxK argument must be a positive integer %s\n", argv[0], argv[argNum+1]);
                    usage(argv[0]);
                    exit(1);
                }
            } else if(strcmp(argv[argNum],"-score")==0)
            {
                char * argscore = argv[argNum+1];
                if(strcmp(argscore,"DHamming1")==0)
                { scoring = scoreDirectional1Hamming; }
                else if(strcmp(argscore,"DHamming2")==0)
                { scoring = scoreDirectional2Hamming; }
                else if(strcmp(argscore,"SHamming")==0)
                { scoring = scoreVanDongen; }
                else if(strcmp(argscore,"DCovering1")==0)
                { scoring = scoreDirectional1Covering; }
                else if(strcmp(argscore,"DCovering2")==0)
                { scoring = scoreDirectional2Covering; }
                else if(strcmp(argscore,"SCovering")==0)
                { scoring = scoreSymmetricCovering; }
                else if(strcmp(argscore,"BCE")==0)
                { scoring = scoreBidirectionnalConsistency; }
                else{
                    fprintf(stderr,"%s error: invalid measure '%s'\n", argv[0], argscore);
                    usage(argv[0]);
                    exit(1);
                }
            }
            else if(strcmp(argv[argNum],"-action")==0)
            {
                char * argscore = argv[argNum+1];
                if(strcmp(argscore,"upperBound")==0)
                { operation = UPPER_BOUND_OPTIM; }
                else if(strcmp(argscore,"threshold")==0)
                { operation = THRESHOLD_OPTIM; }
                else if(strcmp(argscore,"draw")==0)
                { operation = DRAW_UPPER_BOUND_PARTITION; }
                else if(strcmp(argscore,"drawOne")==0)
                { operation = DRAW_ONE_UPPER_BOUND_PARTITION; }
                else if(strcmp(argscore,"straightening")==0)
                { operation = STRAIGHTENING; }
                else{
                    fprintf(stderr,"%s error: invalid operation '%s'\n", argv[0], argscore);
                    usage(argv[0]);
                    exit(1);
                }
            }
        }
    }
	if(inFile == NULL){
	   fprintf(stderr,"%s error: missing input (-i) file argument\n", argv[0]);
            usage(argv[0]);
       exit(1);
	}
	
    if(outFile == NULL ){
	   fprintf(stderr,"%s error: missing output (-o) file argument\n", argv[0]);
            usage(argv[0]);
       exit(1);
	}

    if(numGtFiles == 0  ){
	   fprintf(stderr,"%s error: missing groundtruth (-gt) file argument\n", argv[0]);
            usage(argv[0]);
       exit(1);
	}

    // saliency tree
	g = ReadGraphe(inFile,&Av);
	if(g == NULL)
		exit(1);

	double *level = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
  	saliencyTree(g,&CT,level);
	int * nodeArea = computeNodeArea(CT);


    
	
    int igt;
	if(operation == THRESHOLD_OPTIM)
	{
        double * scoreL = myCalloc(maxK*numGtFiles, sizeof(double));
        double ** score = myCalloc(numGtFiles, sizeof(double *));
        int * nbRegions = myCalloc(numGtFiles, sizeof(int));
        for(igt=0; igt<numGtFiles; ++igt)
        {
            score[igt] = scoreL + igt*maxK;
        }
       
        double * thresholdLevel = myCalloc(maxK, sizeof(double));



        for(i=0; i<maxK*numGtFiles; ++i)
        {
            scoreL[i] = -1.0;
        }
        for(igt=0; igt<numGtFiles; ++igt)
        {
            image = readGrayImage(gtFiles[igt]);
            nbRegions[igt]  = remapGT(image);
            ThresholdKCurve curve = computeScoreCurve(CT, level, image, maxK, nodeArea, scoring);
            for(i=0; i<curve.nbElement; ++i)
            {
                if(curve.score[i]>=0)
                {
                    score[igt][i] = curve.score[i];
                    //sanity check!
                    if(igt==0){
                        thresholdLevel[i]=curve.threshold[i];
                    }
                    else{
                        if(thresholdLevel[i]!=curve.threshold[i])
                        {
                            fprintf(stderr, "Warning %s: inconsistency in thresholds levels detected ! %f != %f \n", argv[0],thresholdLevel[i], curve.threshold[i]);
                        }
                    }
                    
                }
            }
            free(curve.score);
		    free(curve.threshold);
            freeimage(image);
        }
        




        FILE *outfile = fopen(outFile,"w");
        if(outfile == NULL)
        {
            fprintf(stderr,"%s error: cannot open file %s for writing\n", argv[0], outFile);
            exit(1);
        }

        fprintf(outfile,"{\n");
        fprintf(outfile,"\"NUM_GT\":%d,\n", numGtFiles);
        fprintf(outfile,"\"MAXK\":%d,\n", maxK);
        fprintf(outfile,"\"RESULTS\":[\n");
        for(igt=0;igt<numGtFiles;++igt)
        {
            if(igt!=0)
                fprintf(outfile,",");
            fprintf(outfile,"{ \"FILE\":\"%s\", \"NUM_REGIONS\":%d, \"VALUES_HEADER\":[\"K\", \"score\", \"threshold\"], \"VALUES\":[",gtFiles[igt], nbRegions[igt]);
            for(i=0; i<maxK; ++i)
		    {
                if(score[igt][i]>=0)
                {
                    if(i!=0)
                        fprintf(outfile,",");
                    fprintf(outfile,"[%d,%f,%f]",i+1,score[igt][i],thresholdLevel[i]);
                }
		    }
            fprintf(outfile,"]}\n");
        }

		fprintf(outfile,"]}\n");
	
        fclose(outfile);




		
		

        free(score);
        free(scoreL);
        free(thresholdLevel);
	}else if(operation == STRAIGHTENING)
    {
        //if(numGtFiles>1)
         //   fprintf(stderr,"Warning %s: multiple ground truth files not supported for this operation, first provided GT will be used\n",argv[0]);
        double * attr = myCalloc(CT->nbnodes, sizeof(double));
        for(igt=0; igt<numGtFiles; ++igt)
        {
            image = readGrayImage(gtFiles[igt]);
            remapGT(image);

            KCurve * kCurves = computeKScoreDynamique(CT, level, image, maxK, nodeArea, scoring);
            straightenedAttribute(CT, kCurves, attr);
            freeKCurve(kCurves, CT->nbnodes);
            freeimage(image);
        }
		
		
        for(i=0;i<CT->nbnodes;++i)
            attr[i]/=(double)numGtFiles;
		
		
        if(!testIncreasing(CT, attr))
        {
            fprintf(stderr,"Error %s: straightening result is not increasing !\n",argv[0]);
            exit(1);
        }
        SaliencyCache * saliencyCache = NULL;
        saliency(CT, g, attr, &saliencyCache);
        freeSaliencyCache(&saliencyCache);
        SaveGraphe(g, outFile,Av);
        

        //////////////////////////////////////////////////////////////////////////////// bpt save
        int ls = strlen(outFile);
        outFile[ls-4]='b';
        outFile[ls-3]='p';
        outFile[ls-2]='t';
        outFile[ls-1]='\0';
        double ** attrs = myCalloc(2, sizeof(double *));
        attrs[0] = level;
        attrs[1] = attr;
        char ** names = myCalloc(2,sizeof(char*));
        names[0]="oldLevel";
        names[1]="level";
        int * parent = createParentArray(CT);
        if(!saveBPT(outFile, CT->nbnodes, parent, 2, attrs, names))
        {
            exit(1);
        }
        free(attrs);
        free(names);
        free(parent);
        //////////////////////////////////////////////////////////////////////////////// 
        

        free(attr);
    }	else if(operation == DRAW_UPPER_BOUND_PARTITION)
    {
        if(numGtFiles>1)
            fprintf(stderr,"Warning %s: multiple ground truth files not supported for this operation, first provided GT will be used\n",argv[0]);
        image = readGrayImage(gtFiles[0]);
        remapGT(image);
        KCurve * kCurves = computeKScoreDynamique(CT, level, image, maxK, nodeArea, scoring);
        generateKOptimalPartition(g, CT, kCurves, maxK, outFile, true);
        freeKCurve(kCurves, CT->nbnodes);
        freeimage(image);
    }	else if(operation == DRAW_ONE_UPPER_BOUND_PARTITION)
    {
        if(numGtFiles>1)
        {
            fprintf(stderr,"Warning %s: multiple ground truth files not supported for this operation, first provided GT will be used\n",argv[0]);
        }
        image = readGrayImage(gtFiles[0]);
        remapGT(image);
        KCurve * kCurves = computeKScoreDynamique(CT, level, image, maxK, nodeArea, scoring);
        generateKOptimalPartition(g, CT, kCurves, maxK, outFile, false);
        freeKCurve(kCurves, CT->nbnodes);
        freeimage(image);
    }else if(operation == UPPER_BOUND_OPTIM){
        double * scoreL = myCalloc(maxK*numGtFiles, sizeof(double));
        double ** score = myCalloc(numGtFiles, sizeof(double *));
        int * nbRegions = myCalloc(numGtFiles, sizeof(int));
        for(igt=0; igt<numGtFiles; ++igt)
        {
            score[igt] = scoreL + igt*maxK;
        }

        for(igt=0; igt<numGtFiles; ++igt)
        {
            image = readGrayImage(gtFiles[igt]);
            nbRegions[igt] = remapGT(image);
            KCurve * kCurves = computeKScoreDynamique(CT, level, image, maxK, nodeArea, scoring);
            KCurve curve =  kCurves[CT->root];
            for(i=0; i<curve.nbElement; ++i)
            {
                score[igt][i] = curve.data[i].score;
            }
            freeKCurve(kCurves, CT->nbnodes);
            freeimage(image);

        }


        FILE *outfile = fopen(outFile,"w");
        if(outfile == NULL)
        {
            fprintf(stderr,"%s error: cannot open file %s for writing\n", argv[0], outFile);
            exit(1);
        }
        fprintf(outfile,"{\n");
        fprintf(outfile,"\"NUM_GT\":%d,\n", numGtFiles);
        fprintf(outfile,"\"MAXK\":%d,\n", maxK);
        fprintf(outfile,"\"RESULTS\":[\n");
        for(igt=0;igt<numGtFiles;++igt)
        {
            if(igt!=0)
                fprintf(outfile,",");
            fprintf(outfile,"{ \"FILE\":\"%s\", \"NUM_REGIONS\":%d, \"VALUES_HEADER\":[\"K\", \"score\"], \"VALUES\":[",gtFiles[igt], nbRegions[igt]);
            for(i=0; i<maxK; ++i)
		    {
                if(i!=0)
                    fprintf(outfile,",");
                fprintf(outfile,"[%d,%f]",i+1,score[igt][i]);
		    }
            fprintf(outfile,"]}\n");
        }

		fprintf(outfile,"]}\n");
	
        fclose(outfile);


        free(score);
        free(scoreL);
        
	}
	
    
	free(nodeArea);

	terminegraphe(g);
	free(Av);
	free(level); 
	componentTreeFree(CT);
	return 0;
}
