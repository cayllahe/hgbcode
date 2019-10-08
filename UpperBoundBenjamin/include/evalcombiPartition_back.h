#ifndef __EVALCOMBIPARTITION_H__ 
#define __EVALCOMBIPARTITION_H__ 

typedef enum {GroundtruthToSegmentation, SegmentationToGroundtruth, Symmetric} Strategy;


typedef struct xvimage grayImage;
typedef uint8_t uchar;
typedef char bool;
#define true (1)
#define false (0)


typedef struct sGroundTruth{
	grayImage * labelImage;
	int nbPixels;
	int nbElement;
	int * area;
} GroundTruth;


typedef struct sRegionMatching{
	int cardInter;
	int cardUnion;
	float score;
} RegionMatching;

typedef struct sMatchingStat{
	int nbPixels;
	int nbRegionInGT;
	int nbNodes;
	int * regionGTSize;//
	int * regionTreeSize;
	uint8_t * bestMatch; // with i a region of the tree, bestMatch[i] is the region j of GT that maximises the score s(Ri,Rj)
	RegionMatching * bestMatchingScore; // with i a region of the tree, bestMatchingScore[i] is the matching stat of the couple (Ri,Rj) with j = bestMatch[i] 
	RegionMatching ** detailedMatchingScore; // with i a region of the tree, and j a region of GT, detailedMatchingScore[i][j] is the matching stat of the couple (Rj,Ri) 
	
} MatchingStat;



typedef struct sReverseMatch{
	float * scores;
} ReverseMatch;

typedef struct sDynamicData{
	float score;
	double minLevel;
	double maxLevel;
	int minRank;
	int maxRank;
	float directMatchScore;
	float reverseMatchScore;
	ReverseMatch reverseMatchData;
} DynamicData;

typedef struct sKCurve{
	int nbElement;
	DynamicData * data;
} KCurve;

typedef float (*RegionScoringFunction)(RegionMatching * r, int areaRegion1, int areaRegion2); // compute the score between region1 and region2
typedef void (*InitFromNodeFunction)(JCctree * CT, MatchingStat * matchingStat, double * level, int i, DynamicData *d); // initialise the dynamic data for the scoring a partition with a the single region i of the tree
typedef float (*ScoreFusionFunction)(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl ); // compute the scoring for the partition obtained from the fusion of the k-regions partition with dynamicData dk and the l-regions partition with dynamicData dl
typedef void (*FusionFunction)(JCctree * CT, MatchingStat * matchingStat, double * level, int k, DynamicData *dk, int l, DynamicData *dl , DynamicData * res, float fusionScore);// initialise the dynamic data for the scoring a partition obtained from the fusion of the k-regions partition with dynamicData dk and the l-regions partition with dynamicData dl

typedef struct sScoring{
	InitFromNodeFunction initFunction;
	ScoreFusionFunction scoreFusionFunction;
	FusionFunction fusionFunction;
	RegionScoringFunction  regionScoringFunction;
	bool additif;
} Scoring;


#endif