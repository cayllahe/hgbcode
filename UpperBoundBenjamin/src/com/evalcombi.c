

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

#define BLACK 0
#define WHITE 1

typedef struct sdynamicdata{
	float score;
	uchar color;
	int whiteMarker;
	int blackMarker;
	int minRank;
	int maxRank;
	int backTrackK;
	int backTrackL;
	stat s;
} dynamicdata;

typedef struct sKCurve{
	int nbElement;
	dynamicdata * data;
} KCurve;

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
	stat * attr = calloc(CT->nbnodes, sizeof(stat));
	if (attr==NULL)
		{ fprintf(stderr,"Cannot allocate enough memory\n"); exit(1); }
	
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

void backTrack(JCctree * CT, double * level, stat * attribute, KCurve * attributeDynamic, int node, int numNodes)
{
	
	JCctreenode *tabnodes = CT->tabnodes;
	JCsoncell * soncell = tabnodes[node].sonlist;
	int child1 = soncell->son;
	int child2 = soncell->next->son;
	KCurve curve = attributeDynamic[node];
	dynamicdata d = curve.data[numNodes];
	int btk = d.backTrackK;
	int btl = d.backTrackL;
	if (btk ==0 && btl==0)
	{
		printf("node %d (%d,%d)\n", node, attribute[node].TP, attribute[node].FP);
	}else{
		if(btk!=0)
			backTrack(CT, level, attribute, attributeDynamic, child1, btk);
		if(btl!=0)
			backTrack(CT, level, attribute, attributeDynamic, child2, btl);
	}
}

KCurve computeKScoreDynamique(JCctree * CT, double * level, grayImage * image, int maxK, int * totalPositif, scoreFunction scoreF)
{
	int nbPositif=0;
	int i= 0;
	int k,l;
	JCctreenode *tabnodes = CT->tabnodes;

//printf("computeStart entering\n");
	stat * attribute = computeStat(CT, image, &nbPositif);

	stat stat0={.TP=0, .FP=0};
	float score0 = scoreF(stat0, nbPositif);
	dynamicdata ddata0 = (dynamicdata){.score = score0, .color = WHITE, .whiteMarker = 1, .blackMarker = 0, .minRank = CT->root, .maxRank = 0, .s = stat0, .backTrackK=0, .backTrackL=0};
	*totalPositif = nbPositif;
//printf("computeStart exiting nbPositif %d\n", nbPositif);
	KCurve * attributeDynamic = calloc(CT->nbnodes, sizeof(KCurve));
	if (attributeDynamic==NULL)
		{ fprintf(stderr,"Cannot allocate enough memory\n"); exit(1); }
	for(i=0; i <= CT->nbnodes/2; ++i)
	{
		dynamicdata * d = calloc(2, sizeof(dynamicdata));
		if (d==NULL)
			{ fprintf(stderr,"Cannot allocate enough memory\n"); exit(1); }
		//min max should not be interpreted for k=0...
		d[0] = ddata0;
		
		stat statLeaf = attribute[i];
		int p = CT->tabnodes[i].father;
		int refNode = (p!=-1)?p:i;
	

		d[1] = (dynamicdata){.score = scoreF(statLeaf, nbPositif), .color = BLACK, .whiteMarker = 0, .blackMarker = 1, .minRank = refNode, .maxRank = refNode, .s = statLeaf, .backTrackK=0, .backTrackL=0};

		attributeDynamic[i] = (KCurve){.nbElement = 2, .data = d};
	}
//printf("Leaf init finsihed\n");
	for(i=CT->nbnodes/2+1 ; i<CT->nbnodes; ++i)
	{
		JCsoncell * soncell = tabnodes[i].sonlist;

		int child1 = soncell->son;
		int child2 = soncell->next->son;
		KCurve c1 = attributeDynamic[child1];
		KCurve c2 = attributeDynamic[child2];
		
		int nbSum = c1.nbElement + c2.nbElement -1;
		int nbElement = (nbSum <= maxK)? nbSum : maxK;
		
		int maxk = (c1.nbElement <= maxK) ? c1.nbElement : maxK;
		
		dynamicdata * d = calloc(nbElement, sizeof(dynamicdata));
		if (d==NULL)
			{ fprintf(stderr,"Cannot allocate enough memory\n"); exit(1); }
		attributeDynamic[i] = (KCurve){.data = d, .nbElement = nbElement};

		d[0] = ddata0;

		for(k=0; k < nbElement; ++k)
			d[k].score = FLT_MIN ;
		
		for(k=0; k < maxk; ++k)
		{
			dynamicdata dk = c1.data[k];
			int maxl = nbElement - k;
			maxl = (maxl > c2.nbElement)? c2.nbElement : maxl;
			for( l=0; l < maxl; ++l)
			{
				
				dynamicdata dl = c2.data[l];
				stat statFusion = fusionStat(dk.s, dl.s);
				float scoreFusion = scoreF(statFusion, nbPositif);
				if (scoreFusion > d[k+l].score)
				{
					uchar color;
					int wM=0;
					int bM=0;
					/*
					Hugly logic could certainly be simplified
					*/
					if (dk.color == BLACK && dl.color == BLACK)
					{
						color = BLACK;
						bM = 1;
					} else if ((dk.color == BLACK && dl.color == WHITE) || (dk.color == WHITE && dl.color == BLACK))
					{
						color = WHITE;
						bM = dk.blackMarker + dl.blackMarker;
						wM = dk.whiteMarker + dl.whiteMarker;
					} else { // white white
						color = WHITE;
						if(dk.blackMarker>0 &&  dl.blackMarker>0)
						{
							bM = dk.blackMarker + dl.blackMarker;
							wM = dk.whiteMarker + dl.whiteMarker;
						}else if(dk.blackMarker==0 &&  dl.blackMarker==0)
						{
							wM=1;
						}else{
							if(dk.blackMarker==0)
							{
								wM = dl.whiteMarker;
								bM = dl.blackMarker;
							}else{
								wM = dk.whiteMarker;
								bM = dk.blackMarker;
							}
						}
						
					}
					/*if(child2==308714 || child1==308714 || child1==308711 || child1==308716 || child1==308692 || child2==308711 || child2==308716 || child2==308692)
					{
						printf("c1:%d c2:%d best k:%d l:%d score:%f\n",child1, child2, k, l, scoreFusion);
					}*/
					d[k+l] = (dynamicdata){.score = scoreFusion, .color = color, .whiteMarker = wM, .blackMarker = bM, .minRank = mini(dk.minRank, dl.minRank), .maxRank = maxi(dk.maxRank, dl.maxRank), .s = statFusion, .backTrackK=k, .backTrackL=l};
				}
			}
		}

		stat statRegion = attribute[i];
		float scoreRegion = scoreF(statRegion, nbPositif);
		if (scoreRegion >= d[1].score)
		{
			/*if(child2==308714 || child1==308714 || child1==308711 || child1==308716 || child1==308692 || child2==308711 || child2==308716 || child2==308692)
			{printf("arg score i:%d (%d,%d) %f\n", i, statRegion.TP, statRegion.FP, scoreRegion);}

			if(i==308714  || i==308711 || i==308716 || i==308692 )
			{printf("plop score i:%d (%d,%d) %f\n", i, statRegion.TP, statRegion.FP, scoreRegion);}*/
			d[1] = (dynamicdata){.score = scoreRegion, .color = BLACK, .whiteMarker = 0, .blackMarker = 1, .minRank = i, .maxRank = i, .s = statRegion, .backTrackK=0, .backTrackL=0};
		}


	}

	backTrack(CT, level, attribute, attributeDynamic, CT->root, 4);


	KCurve result = attributeDynamic[CT->root];

	

	for(i=0 ; i<CT->root ; ++i)
		free(attributeDynamic[i].data);
	
	free(attributeDynamic);

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

	double *Val = calloc(2*g->nbsom, sizeof(double));
	if(Val == NULL){ 
		fprintf(stderr, "extinctionvalues: Cannot allocate memory\n"); 
		exit(1);
  	}
	JCctree *CT;
  	saliencyTree(g,&CT,Val);

	int nbPositif;
	KCurve curve = computeKScoreDynamique(CT, Val, image, maxK, &nbPositif, scoreF);
	printf("#K\tscore\tmarkerNegatif\tmarkerPositif\tminRank\tmaxRank\n");
	for(i=0; i<curve.nbElement; ++i)
	{
		printf("%d\t%f\t%d\t%d\t%d\t%d\n",i,curve.data[i].score,curve.data[i].whiteMarker,curve.data[i].blackMarker,curve.data[i].minRank,curve.data[i].maxRank);
	}

	
	freeimage(image);
	terminegraphe(g);
	free(Av);
	free(Val); 
	componentTreeFree(CT);
	return 0;
}
