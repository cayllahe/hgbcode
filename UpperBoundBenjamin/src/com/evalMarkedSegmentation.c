
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

typedef struct xvimage grayImage;
typedef uint8_t uchar;
typedef char bool;


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
  fprintf(stderr,"USAGE: %s saliency_map.graph  marker.pgm groundTruth.pgm [outFile]\n",arg);
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



uchar * markedSegmentation(JCctree *CT, grayImage * marker)
{
    int nbNodes = CT->nbnodes;
    uchar * ma = myCalloc(nbNodes, sizeof(uchar));

    int i;
    uchar *mark = UCHARDATA(marker);
	int width=rowsize(marker), heigth=colsize(marker), N=width*heigth;
    JCctreenode *tabnodes = CT->tabnodes;
    for(i=0;i<=nbNodes/2;++i)
    {
        ma[i] = mark[i]; // 1 for background and 2 for foreground  : nice bitmask !
    }
    for(i=0;i<nbNodes-1;++i)
    {
        int p = tabnodes[i].father;
        ma[p] = ma[p] | ma[i];
    }
    
    ma[nbNodes-1] = 1;


    for(i=nbNodes-2;i>=0;--i)
    {
        if(ma[i]==0)
        {
            ma[i] = ma[tabnodes[i].father];
        }else if(ma[i]==3)
        {
            ma[i]=1;
        }
    
    }

    uchar * result = myCalloc(N, sizeof(uchar));
    memcpy(result, ma, N*sizeof(uchar));
    free(ma);
    return result;
}

typedef struct sPrScore{
    double precision;
    double recall;
    double fm;
} PrScore;

PrScore evaluate(uchar * segmentation, grayImage * gtImage)
{
    uchar *gt = UCHARDATA(gtImage);
	int width=rowsize(gtImage), heigth=colsize(gtImage), nbPixel=width*heigth;
    uchar background = 0;
    uchar foreground = 0;

    int i;
    for(i=0;i<nbPixel;++i)
        if (gt[i]>foreground) foreground=gt[i];

    int tp=0;
    int fp=0;
    int fn=0;
    for(i=0;i<nbPixel;++i)
    {
        if(gt[i]==foreground && segmentation[i]==2)
        { ++tp; }
        else if(gt[i]==background && segmentation[i]==2)
        { ++fp; }
        else if(gt[i]==foreground && segmentation[i]==1)
        { ++fn; }
    }
    double precision = ((double)tp) /((double)tp+fp);
    double recall =  ((double)tp) /((double)tp+fn);
    PrScore result = {.precision = precision, .recall = recall, .fm=2.0*precision*recall/(precision+recall)};
    return result;
}


/** -------------------------------------
Main 
*/


int32_t main(int32_t argc, char ** argv) 
{
	graphe *g;
	grayImage * marker, * gt;
	
	
	double *Av;


	if( argc != 5 && argc != 4  ){
	usage(argv[0]);
	exit(1);
	}
	fprintf(stderr, "Marker evaluation: %s - %s\n",argv[1], argv[2]);
	g = ReadGraphe(argv[1],&Av);
	if(g == NULL)
		exit(1);
	marker = readGrayImage(argv[2]);
	gt = readGrayImage(argv[3]);
	

	double *level = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
  	saliencyTree(g,&CT,level);

	remapGT(marker);
	remapGT(gt);
	



    uchar * segmentation = markedSegmentation(CT, marker);
    PrScore score = evaluate(segmentation, gt);
    FILE *outfile;
    if(argc ==5)
    {
        char * outFile = argv[4];
        outfile = fopen(outFile,"w");
        if(outfile == NULL)
        {
            fprintf(stderr,"%s error: cannot open file %s for writing\n", argv[0], outFile);
            exit(1);
        }
    }else{
        outfile = stdout;
    }
   
    
    fprintf(outfile,"%0.5f\t%0.5f\t%0.5f\n", score.precision, score.recall, score.fm);
    if(argc==5)
        fclose(outfile);
/*
    int i;
    uchar *mark = UCHARDATA(marker);
    for(i=0;i<g->nbsom;++i)
        if (segmentation[i]==2)
            mark[i]=255;
        else if (segmentation[i]==1) 
            mark[i]=127;
        else mark[i]=0;

    writeimage(marker, "testsegmmark.pgm");*/
    
    free(segmentation);
	freeimage(gt);
    freeimage(marker);
	terminegraphe(g);
	free(Av);
	free(level); 
	componentTreeFree(CT);
	return 0;
}
