
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "mat.h"
#include <vector>
#include <queue>
#include <set>
#include <map>
extern "C" { 
#include "base.h"
#include <mcweightgraph.h>
#include <lcomptree.h>
#include <lhierarchie.h>
}

//@TODO memory leaks eveywhere !

typedef unsigned short int uint16;





void usage(char *arg){
    fprintf(stderr,"#################################################################\n\n");
    fprintf(stderr,"USAGE: %s saliency_map.graph  outputFile.mat\n",arg);
    fprintf(stderr,"#################################################################\n\n");
}










/*
insure that region labeling start at 0 and ends at N-1 (with N the number of regions)
*/
void remap2int(int * image, int len)
{
	std::set<int> values;
    for(int i=0;i<len;++i)
        values.insert(image[i]);
    std::map<int,int> myMap;
    int c=1;
    for (std::set<int>::iterator it=values.begin(); it!=values.end(); ++it)
    {
         myMap[*it]=c++;
    }

    for(int i=0;i<len;++i)
        image[i] = myMap[image[i]];
}




// find the highest node with a level strictly lower than curLevel (starting from curNode)
int findNextNode(JCctree *CT, double * level, int curNode, double curLevel)
{
	int nbPix = CT->nbnodes/2;
	while(level[curNode] >= curLevel && curNode > nbPix)
		--curNode;
	return curNode;
}

void drawPixels(JCctree *CT, int node, int * image, int nbPixels)
{
    
    JCctreenode *tabnodes = CT->tabnodes;
    std::queue<int> q;
    q.push(node);

    if(node<nbPixels)
        image[node] = node;
    else
        while(!q.empty())
        {
            int n = q.front();
    
            q.pop();
            
            JCsoncell * soncell = tabnodes[n].sonlist;
            int child1 = soncell->son;
            int child2 = soncell->next->son;
    
            if(child1<nbPixels)
                image[child1] = node;
            else q.push(child1);
    
            if(child2<nbPixels)
                image[child2] = node;
            else q.push(child2);
        }
}

// init a Partition structure corresponding to the given threshold
int initPartition(JCctree *CT, double * level, double threshold, int curNode,  int maxK, int * image)
{
    int nbPixels = CT->nbnodes/2+1;
	JCctreenode *tabnodes = CT->tabnodes;
	
	int count=0;
	int i;
	for(i=curNode+1; i < CT->nbnodes; ++i)
	{
       //printf("%d/%d\n",i,CT->root);
		//int p=tabnodes[i].father;
		JCsoncell * soncell = tabnodes[i].sonlist;
		int child1 = soncell->son;
		int child2 = soncell->next->son;

		if(level[child1] < threshold)
		{
            
			if (count == maxK)
				return -1;
            drawPixels(CT, child1, image, nbPixels);
            ++count;
			
		}
		if(level[child2] < threshold)
		{
			if (count == maxK)
				return -1;
            drawPixels(CT, child2, image, nbPixels);
            ++count;
		}
		
	}

    remap2int(image,nbPixels);
	return count;
}


std::vector<int *> preparePartitions(JCctree *CT, double * level, int nbPartition, int scaledPartition=10,int maxRegions=3000)
{
    int curNode = CT->root;
	double curLevel = level[curNode];
    
	int nbPixels = CT->nbnodes/2+1;
    
    std::vector<int *> res;
    int k=1;
    int count=0;

    int nbPartition1 = nbPartition - scaledPartition;
    

    do {
		
		//--k;
		//printf("-thrshold %f node %d/%d  k: %d\n", curLevel,curNode, CT->nbnodes,k);

		curNode = findNextNode(CT, level, curNode, curLevel);
		if(curNode < nbPixels)
			break;
        int * image = new int[nbPixels];
		k = initPartition(CT, level, curLevel, curNode,  maxRegions, image);
        curLevel = level[curNode];
        if(k==-1)
            break;
		res.push_back(image);
		
        count ++;
			
	} while( k!=-1 && count<nbPartition1 );
    
    if(k!=-1)
    {
        double stepLevel = (curLevel-level[nbPixels])/scaledPartition;
        count = 0;
        do {
		
            //--k;
            //printf("-thrshold %f node %d/%d  k: %d\n", curLevel,curNode, CT->nbnodes,k);
    
            curNode = findNextNode(CT, level, curNode, curLevel-stepLevel);
            if(curNode < nbPixels)
                break;
            int * image = new int[nbPixels];
            k = initPartition(CT, level, curLevel-stepLevel, curNode,  maxRegions, image);
            curLevel = level[curNode];
            if(k==-1)
                break;
            res.push_back(image);
            
            count ++;
			
	   } while( k!=-1 && count<scaledPartition );
    }
    return res;
}

void saveAsMat(char * filename, std::vector<int *> imgs, int width, int  height, int nbe)
{
    
    int nbpix = width*height;
    mwSize dims[2];
    dims[0] = height;
    dims[1] = width;
    mxArray *cell_array_ptr;
    MATFile *pmat;
    cell_array_ptr = mxCreateCellMatrix(1,(mwSize)nbe);

    pmat = matOpen(filename, "w");
    if (pmat == NULL) {
        printf("Error creating file %s\n", filename);
        printf("(Do you have write permission in this directory?)\n");
        exit(EXIT_FAILURE);
    }
    for(unsigned int i=0; i<(mwIndex)nbe; i++){
        int * im = imgs[mini(i,imgs.size()-1)];
        mxArray *a = mxCreateNumericArray(2, dims, mxUINT16_CLASS, mxREAL);
        uint16 * data = (uint16 *)mxGetData(a);
        for(int j=0;j<nbpix;++j)
        {
            int r = j/width;
            int c = j%width;
            data[c*height+r]=im[j];
        }
        mxSetCell(cell_array_ptr,i,a);
    }
    int status = matPutVariable(pmat, "segs", cell_array_ptr);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    } 
    if (matClose(pmat) != 0) {
        printf("Error closing file %s\n",filename);
        exit(EXIT_FAILURE);
    }
}


/** ****************************************************************************************
Main 
*/


int32_t main(int32_t argc, char ** argv) 
{
    int nbpart = 50;

	graphe *g;
	

	double *Av;

    if(argc!=3)
    {
        usage(argv[0]);
        exit(1);
    }

    char * inFile=argv[1];
    char * outFile=argv[2];
    

    // saliency tree
	g = ReadGraphe(inFile,&Av);
	if(g == NULL)
		exit(1);
   
	double *level = (double *)myCalloc(2*g->nbsom, sizeof(double));

	JCctree *CT;
  	saliencyTree(g,&CT,level);
	 //printf("%d/%d", CT->nbnodes/2,g->rs*g->cs );
    std::vector<int *> imgs = preparePartitions(CT, level, nbpart);
    saveAsMat(outFile, imgs, g->rs, g->cs,nbpart);
	terminegraphe(g);
	free(Av);
	free(level); 
	componentTreeFree(CT);
	return 0;
}
