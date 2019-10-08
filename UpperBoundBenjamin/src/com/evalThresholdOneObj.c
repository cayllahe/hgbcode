

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
typedef struct xvimage grayImage;
typedef unsigned char uchar;
inline int mini(int a, int b){
	return (a<=b)?a:b;
}

inline int maxi(int a, int b){
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
  fprintf(stderr,"\t scoreFunction in {'FMeasure','Jaccard','Error'}\n");
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
		fprintf(stderr, "3D not yet implemented\n");
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
		fprintf(stderr,"Only byte image implemented \n");
		exit(1);
	}
	
	return Im;
}

typedef struct sstat{
	int TP;
	int FP;
} stat;





typedef float (*scoreFunction)(stat, int);

float scoreError(stat s, int nbPositif)
{
	return 1.0 - (float)(s.FP  + (nbPositif-s.TP))/(float)nbPositif;
}

float scoreFMeasure(stat s, int nbPositif)
{
	return 2.0*s.TP/(float)(s.TP + s.FP + nbPositif);
}

float scoreJaccard(stat s, int nbPositif)
{
	return s.TP/(float)(s.FP + nbPositif);
}

stat fusionStat(stat s1, stat s2)
{
	stat res = {.TP=s1.TP + s2.TP, .FP=s1.FP + s2.FP};
	return res;
}

stat* computeStat(JCctree *CT, grayImage * image, int * nbPoints)
{

	JCctreenode *tabnodes = CT->tabnodes;
	stat * attr = myCalloc(CT->nbnodes, sizeof(stat));
	
	
	int i;
	uint8_t *marker = UCHARDATA(image);

	uint8_t minM=255;
	for(i=0; i <= CT->nbnodes/2; ++i)
		if(marker[i]<minM)
			minM=marker[i];
	//int count=0;
	//for(i=0; i < CT->nbnodes/2; ++i)
	//	if(marker[i]>0)
	//		++count;
	
	
	for(i=0; i <= CT->nbnodes/2; ++i)
	{
		if(marker[i]==minM)
		{
			attr[i] = (stat){.TP=0, .FP=1};

		}else{
			attr[i] = (stat){.TP=1, .FP=0};
		}
	}
	
	for(i=CT->nbnodes/2+1 ; i<CT->nbnodes; ++i)
	{
		JCsoncell * soncell = tabnodes[i].sonlist;
		int c1 = soncell->son;
		int c2 = soncell->next->son;
		attr[i] = fusionStat(attr[c1], attr[c2]);
	}
	
	*nbPoints = attr[CT->root].TP;

	return attr;
	
}





int findNextNode(JCctree *CT, double * level, int curNode, double curLevel)
{
	int nbPix = CT->nbnodes/2;
	while(level[curNode] >= curLevel && curNode > nbPix)
		--curNode;
	return curNode;
}

#define MAXREGIONPARTITION 500

typedef struct {
	int nbElements;
	int nodeNumbers[MAXREGIONPARTITION];
	uchar available[MAXREGIONPARTITION];
} Partition ;



int initPartition(JCctree *CT, double * level, double threshold, Partition* partition, stat * statRegion)
{
	JCctreenode *tabnodes = CT->tabnodes;
	
	int count=0;
	int i;
	for(i=0; i< CT->nbnodes; ++i)
	{
		int p=tabnodes[i].father;
		if(level[i] < threshold && ( p==-1 || level[p]>=threshold) && statRegion[i].TP>0)
		{


			partition->nodeNumbers[count]=i;
			++count;
			if (count >= MAXREGIONPARTITION)
				return -1;
		}

		
		
	}
	partition->nbElements = count;
	return count;
}



typedef struct {
	int nbElements;
	float * score;
	double * thr;
} ScoreCurve;

void analysePartition(Partition * p, stat * regionStat, ScoreCurve * score, scoreFunction sF, int nbPositif, double threshold, double * level)
{
	int i;
	int k;
	int maxK = score->nbElements;
	int nbRegion = p->nbElements;
	for(i=0; i<nbRegion; ++i)
		p->available[i]=1;
	float best =0;
	stat curStat = {.TP = 0, .FP =0 };
	for (k=0; k < maxK; ++k)
	{
		float newBest = 0;
		int bestIndex=-1;
		for(i=0; i<nbRegion; ++i)
		{
			if(p->available[i])
			{
				stat fs = fusionStat(curStat, regionStat[p->nodeNumbers[i]]);
				float sc = sF(fs, nbPositif);
				if( sc > newBest )
				{
					newBest = sc;
					bestIndex = i;
				}
			}
		}

		if (newBest>best)
		{
			curStat = fusionStat(curStat, regionStat[p->nodeNumbers[bestIndex]]);
			best = newBest;
			p->available[bestIndex] = 0;
			if(fabs(threshold - 4.909893) < 0.0001)
				{
					printf("new best k=%d level:%f score %f, with node %d  (%d,%d), stat: (TP: %d,FP: %d, N: %d)\n", k, level[p->nodeNumbers[bestIndex]], best,  p->nodeNumbers[bestIndex],regionStat[p->nodeNumbers[bestIndex]].TP,regionStat[p->nodeNumbers[bestIndex]].FP, curStat.TP, curStat.FP, nbPositif);
				}
			if(score->score[k]<best)
			{
				/*if(fabs(threshold - 4.909893)<0.0001)
				{
					printf("new best k=%d score %f, with node %d  (%d,%d), stat: (TP: %d,FP: %d, N: %d)\n", k, best,  p->nodeNumbers[bestIndex],regionStat[p->nodeNumbers[bestIndex]].TP,regionStat[p->nodeNumbers[bestIndex]].FP, curStat.TP, curStat.FP, nbPositif);
				}*/
				score->score[k]=best;
				score->thr[k] =  threshold;
			}
		}else{return;}
		
	}
	

}

ScoreCurve computeScoreCurve(JCctree *CT, double * level, grayImage * image, int maxK, scoreFunction  sF)
{
	int curNode = CT->root;
	double curLevel = level[curNode]+1;
	int nbPixels = CT->nbnodes/2;
	Partition partition;
	int nbPositif;
	int i;
	stat * regionStat = computeStat(CT, image, &nbPositif);

	ScoreCurve scores;
	scores.score = myCalloc(maxK-1, sizeof(float));
	scores.thr = myCalloc(maxK-1, sizeof(double));
	scores.nbElements=maxK-1;
	for(i=0; i<maxK-1; ++i)
		scores.score[i] =0 ;


	while( curNode > nbPixels)
	{
		curNode = findNextNode(CT, level, curNode, curLevel);
		if(curNode <= nbPixels)
			break;
		curLevel = level[curNode];
		printf("-thrshold %f\n", curLevel);
		if(initPartition(CT, level, curLevel, &partition, regionStat)==-1)
			break;
		analysePartition(&partition, regionStat, &scores, sF, nbPositif, curLevel, level);
		
	}


	return scores;
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
	fprintf(stderr,"%s info: starting %s\n", argv[0], argv[1]);
	g = ReadGraphe(argv[1],&Av);
	if(g == NULL)
		exit(1);
	image = readGrayImage(argv[2]);
	
	scoreFunction scoreF;
	if(strcmp(argv[3],"FMeasure")==0)
	{ scoreF = &scoreFMeasure; }
	else if(strcmp(argv[3],"Jaccard")==0)
	{ scoreF = &scoreJaccard; }
	else if(strcmp(argv[3],"Error")==0)
	{ scoreF = &scoreError; }
	else{
		fprintf(stderr,"%s info: invalid error measure %s\n", argv[0], argv[3]);
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
	//printf("maxk = %d\n",maxK);

	double *Val = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
  	saliencyTree(g,&CT,Val);


	ScoreCurve curve = computeScoreCurve(CT, Val, image, maxK, scoreF);


	printf("#K\tscore\tthreshold\n");
	for(i=0; i<curve.nbElements; ++i)
	{
		printf("%d\t%f\t%f\n",i+1,curve.score[i], curve.thr[i]);
	}

	
	freeimage(image);
	terminegraphe(g);
	free(Av);
	free(Val); 
	componentTreeFree(CT);
	return 0;
}
