
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
#include "llca.h"

#include <errno.h>

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
  fprintf(stderr,"USAGE: %s input_saliency_map.graph  area output_saliency_map.graph\n",arg);
  fprintf(stderr,"#################################################################\n\n");
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



void saliency(JCctree *ST, graphe *g, double * STAltitudes)
{  
  int32_t logn, nbRepresent, u, x, y, c1;
  int32_t nbarcs = g->ind;

  /* Data structure for fast computation of LCA (least common ancestor) */
  int32_t *Euler, *Depth, *Represent, *Number, **Minim;

  Euler = myCalloc(2*ST->nbnodes-1, sizeof(int32_t));
  Represent = myCalloc(2*ST->nbnodes-1, sizeof(int32_t));
  Depth = myCalloc(ST->nbnodes, sizeof(int32_t));
  Number = myCalloc(ST->nbnodes, sizeof(int32_t));
  
  Minim = LCApreprocess(ST, Euler, Depth, Represent, Number, &nbRepresent, &logn);

  // For any edge of g
  for(u = 0; u < nbarcs; u++){
    x = g->tete[u]; y = g->queue[u];
    c1 = Represent[LowComAncFast(x, y, Euler, Number, Depth, Minim)];
    setweight(g,u,STAltitudes[c1]);
  }  
  free(Euler);
  free(Represent);
  free(Depth);
  free(Number);
  free(Minim[0]);
  free(Minim);

}


/** ****************************************************************************************
Main 
*/


int32_t main(int32_t argc, char ** argv) 
{
	graphe *g;

	
	
	double *Av;
	int i;

	if(argc != 4   ){
	usage(argv[0]);
	exit(1);
	}
	
	g = ReadGraphe(argv[1],&Av);
	if(g == NULL)
		exit(1);

        int nbPix = g->nbsom;

      

	char * stop;

	int minArea = (int)strtol(argv[2], &stop, 10);
	if(stop[0]!='\0' || minArea <=0)
	{
                double ratio = (double)strtod(argv[2], &stop);
                if(stop[0]!='\0' || ratio <0 || ratio >1)
                {
                fprintf(stderr,"%s error: second argument must be either a strictly positive integer or a floating point number in [0;1] - %s\n", argv[0], argv[2]);
		exit(1);
                }
                minArea = (int)(ratio * nbPix);
		
	}
        printf("Area %d\n",minArea);
	double *level = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
    int * MST;
  	saliencyTree2(g,&CT,level, &MST);
	int * nodeArea = computeNodeArea(CT);

    int nbLeaves = CT->nbnodes/2+1;
    JCctreenode *tabnodes = CT->tabnodes;
   
    for(i=nbLeaves ; i<CT->nbnodes; ++i)
    {
        JCsoncell * soncell = tabnodes[i].sonlist;
		int child1 = soncell->son;
		int child2 = soncell->next->son;
        if(nodeArea[i] < minArea || nodeArea[child1] < minArea || nodeArea[child2] < minArea)
            level[i]=0;
    }
    
    graphe * g2 = MSTToGraph(g, MST, level);
    free(MST);
    componentTreeFree(CT);

    saliencyTree2(g2,&CT,level, &MST);

    saliency(CT, g, level);
    SaveGraphe(g, argv[argc-1],Av);
    free(MST);
	free(nodeArea);
	terminegraphe(g);
    terminegraphe(g2);
	free(Av);
	free(level); 
	componentTreeFree(CT);
	return 0;
}
