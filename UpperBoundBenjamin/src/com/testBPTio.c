#include "bptio.h"
#include <stdio.h>
#include <stdlib.h>
#include <lhierarchie.h>


double * computeNodeArea(JCctree *CT)
{
	double * area = myCalloc(CT->nbnodes, sizeof(double));
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

bool compareTabInt(int * t1, int * t2, int nbe)
{
    int i;
    for(i=0;i<nbe;++i)
        if(t1[i] != t2[i])
            return false;
    return true;
}

bool compareTabDouble(double * t1, double * t2, int nbe)
{
    int i;
    for(i=0;i<nbe;++i)
        if(t1[i] != t2[i])
            return false;
    return true;
}

int main(int argc, char ** argv)
{
    if(argc!=3) // options come in pairs
    {
        fprintf(stderr,"Error: usage %s in.graph out.bpt\n", argv[0]);
        exit(1);
    }
    char * inFile = argv[1];
    char * outFile = argv[2];
    double *Av;
    // saliency tree
	graphe * g = ReadGraphe(inFile,&Av);
	if(g == NULL)
		exit(1);

	double *level = myCalloc(2*g->nbsom, sizeof(double));

	JCctree *bpt;
  	saliencyTree(g,&bpt,level);
	double * nodeArea = computeNodeArea(bpt);

    double ** attrs = myCalloc(2, sizeof(double *));
    attrs[0] = level;
    attrs[1] = nodeArea;
    char ** names = myCalloc(2,sizeof(char*));
    names[0]="level";
    names[1]="nodeArea";
    int * parent = createParentArray(bpt);
    if(!saveBPT(outFile, bpt->nbnodes, parent, 2, attrs, names))
    {
        exit(1);
    }
    
    int * nparent;
    int nbNodes;
    int numAttr;
    double ** nattrs;
    char ** nnames;
    if(!readBPT(outFile, &nbNodes, &nparent, &numAttr, & nattrs, & nnames))
    {
        exit(1);
    }
    if(nbNodes != bpt->nbnodes)
    {
        printf("bug nbnodes\n");exit(1);
    }
     if(numAttr != 2)
    {
        printf("bug nbAttr\n");exit(1);
    }
    
    printf("Attr 1:%s\n",nnames[0]) ;
    printf("Attr 2:%s\n",nnames[1]) ;
//bool loadBPT(char * path, int * nbNodes, int ** parents, int * numAttr, double *** attrs, char *** attrNames); 
    if(!compareTabInt(parent, nparent, nbNodes))
    {
        printf("bug parent\n");exit(1);
    }
     if(!compareTabDouble(level, nattrs[0], nbNodes))
    {
        printf("bug attr 1\n");exit(1);
    }
    if(!compareTabDouble(nodeArea, nattrs[1], nbNodes))
    {
        printf("bug attr 2\n");exit(1);
    }
    return 0;
}