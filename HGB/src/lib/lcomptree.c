/*
 Copyright ESIEE (2009)
 
 m.couprie@esiee.fr
 
 This software is an image processing library whose purpose is to be
 used primarily for research and teaching.
 
 This software is governed by the CeCILL  license under French law and
 abiding by the rules of distribution of free software. You can  use,
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info".
 
 As a counterpart to the access to the source code and  rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability.
 
 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or
 data to be ensured and,  more generally, to use and operate it in the
 same conditions as regards security.
 
 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.
 */

/* authors : J. Cousty - L. Najman and M. Couprie */


#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mcweightgraph.h>
#include <mcunionfind.h>
#include <mcsort.h>
#include <lcomptree.h>
#include <MST.h>
#include <time.h>
#include <llca.h>




/* ==================================== */
JCctree * componentTreeAlloc(int32_t N)
/* ==================================== */
#undef F_NAME
#define F_NAME "componentTreeAlloc"
{
    JCctree *CT;
    CT = (JCctree *)malloc(sizeof(JCctree));
    CT->tabnodes = (JCctreenode *)malloc(N * sizeof(JCctreenode));
    CT->tabsoncells = (JCsoncell *)malloc(2*N * sizeof(JCsoncell));
    CT->flags = NULL;
    // CT->flags = (uint8_t *)calloc(N, sizeof(char));
    // If needed, must be allocated separatly
    //  memset(CT->flags, 0, N);
    if ((CT == NULL) || (CT->tabnodes == NULL) || (CT->tabsoncells == NULL))
    {
        fprintf(stderr, "%s : malloc failed\n", F_NAME);
        return NULL;
    }
    CT->nbnodes = N;
    CT->nbsoncells = 0;
    return CT;
} // componentTreeAlloc()

/* ==================================== */
void componentTreeFree(JCctree * CT)
/* ==================================== */
{
    free(CT->tabnodes);
    free(CT->tabsoncells);
    if(CT->flags != NULL)free(CT->flags);
    free(CT);
} // ComponentTreeFree()
/* ==================================== */
void componentTreeFreeQCT(JCctree * CT)
/* ==================================== */
{
    if(CT->tabnodes!=NULL)
        free(CT->tabnodes);
    if(CT->tabsoncells!=NULL)
        free(CT->tabsoncells);
    if(CT->flags != NULL)
        free(CT->flags);
    free(CT);
} // ComponentTreeFree()


void checkDegenerateNodes(JCctree* CT, double * STaltitudes,int32_t root ){
 JCsoncell *s;
    
    if(CT->tabnodes[root].nbsons<2 && root > CT->root ){
        printf("nodo degenerado %i at altitude %f\n",root, STaltitudes[root] );
    }
    
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
        
        checkDegenerateNodes(CT, STaltitudes, s->son );
    }

    
    
}



/* ==================================== */
void componentTreePrintAttributesDOT(JCctree *CT, double *Alt, int32_t root, int *STArea, double* STMaxWeight, FILE * fd)
/* ==================================== */
{
    //if(STArea[root]==1 && CT->tabnodes[root].nbsons<1)return;
    JCsoncell *s;
    fprintf(fd,"%i [label=\"%i\\n alt %f\\n nro %i\\n maxw %f\"] \n", root, root, Alt[root] , STArea[root], STMaxWeight[root] );
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
      //  if(STArea[s->son]==1 && CT->tabnodes[s->son].nbsons<1)continue;
        fprintf(fd,"%i->%i\n",root, s->son);
    }
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){

        componentTreePrintAttributesDOT(CT, Alt, s->son, STArea,STMaxWeight, fd );
    }
}

void saveDOTfile(JCctree *CT, double *Alt, int32_t root, int *STArea, double* STMaxWeight, char * filename, int nbsom ){

   // char * filename=  "/Users/edus/Dropbox/MatMorphology/codigo/SM_Edward/Graphs/dots/graph.dot";
    FILE * fd = NULL;
    
        fd = fopen(filename,"w");
        if (!fd)
        {
            fprintf(stderr, "%s: cannot open file: %s\n", F_NAME, filename);
            return;
        }
    
    fprintf(fd,"digraph{\n");
    
    componentTreePrintAttributesDOT(CT, Alt, root, STArea, STMaxWeight, fd);
    
    fprintf(fd,"{rank=same;");
    int i;
    for(i = 0; i< nbsom;i++){
       fprintf(fd," %i ", i);
    }
    
    
    fprintf(fd,"}\n}");
    fclose(fd);
    
}


/* ==================================== */
void componentTreePrintAttributes(JCctree *CT, double *Alt, int32_t root, int *STArea, double* STMaxWeight)
/* ==================================== */
{
    
    
    JCsoncell *s;
    char* c = NULL;
    if(root <= CT->nbnodes/2) c = "(leaf)"; else c = "(branch)";
    printf("---------------------\nnode %d %s\nAltitude: %lf\n", root, c, Alt[root]);
    printf("Area: %d\n",STArea[root]);
    //printf("MaxWeight: %d\n",STMaxWeight[root]);
    printf("children: ");
    
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
        if(s->son <= CT->nbnodes/2) c = "S"; else c = "C";
        printf("%c%d (altitude %lf) ",*c,s->son, Alt[s->son]);
    }
    printf("\n---------------------\n");
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next)
        componentTreePrintAttributes(CT, Alt, s->son, STArea,STMaxWeight );
}


/* ==================================== */
void componentTreePrint(JCctree *CT, double *Alt, int32_t root, int32_t *Att1)
/* ==================================== */
{
    JCsoncell *s;
    char* c = NULL;
    if(root <= CT->nbnodes/2) c = "(leaf)"; else c = "(branch)";
    printf("---------------------\nnode %d %s\nAltitude: %lf\n", root, c, Alt[root]);
    if( Att1 != NULL )
        printf("Attribut 1: %d\n",Att1[root]);
    
    printf("children: ");
    
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
        if(s->son <= CT->nbnodes/2) c = "S"; else c = "C";
        printf("%c%d (altitude %lf) ",*c,s->son, Alt[s->son]);
    }
    printf("\n---------------------\n");
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next)
        componentTreePrint(CT, Alt, s->son, Att1);
}

/* ==================================== */
void calculReversePointer(JCctree *CT, int32_t root)
/* ==================================== */
{
    JCsoncell *s;
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next)
    {
        calculReversePointer(CT, s->son);
        CT->tabnodes[s->son].father = root;
    }
}

/* The algorithm is an original adaptation/improvment to edge-weighted
 graph of the algorithm of Najman and Couprie published in IEEE IP
 in 2006 / Its description is provided in Cousty, Najman, Perret
 ISMM2013 and in Najman, Cousty, Perret ISMM 2013. */
int32_t BPTAO(/* INPUTS */
              graphe *g,
              /* OUTPUTS */
              JCctree ** SaliencyTree, /* Component tree - BPTAO  */
              double *STaltitude,  /* weights associated to the
                                    nodes of the BPTAO: Must
                                    be allocated elsewhere */
              MST* mst
              )
{
    int32_t i,x1,x2,n1,n2,z,k, STx, STy, nbsoncellsloc;
    JCctree *ST;
    int32_t *clefs;
    int32_t *STmap;
    JCsoncell * newsoncell1;
    JCsoncell * newsoncell2;
    Tarjan *T;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    
    if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de malloc\n");
    }
    
    if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
        fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
        exit(0);
    }
    
    setSizeMST(mst, narcs);
    
    
    for(i = 0; i < narcs; i++)
        clefs[i] = i;
    
    
    
    d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
    
  
    
    
    if( (ST = componentTreeAlloc((2*taille))) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
        exit(0);
    }
    ST->nbnodes = taille;
    ST->nbsoncells = 0;
    T = CreeTarjan(taille);
    for(i = 0; i < taille; i++)
    {
        STmap[i] = i;
        ST->tabnodes[i].nbsons = 0;
        ST->tabnodes[i].sonlist = NULL;
        TarjanMakeSet(T,i);
    }
    
    for(i = 0; i < narcs; i++){
        // for each edge of the graph or MST taken in increasing order of altitude
        n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
        if(n1 != n2) {
            /* If the two vertices do not belong to a same connected
             component at a lowest level */
            // Which component of ST n1 and n2 belongs to ?
            x1 = STmap[n1]; x2 = STmap[n2];
            // Create a new component
            z = ST->nbnodes; ST->nbnodes++;
            nbsoncellsloc = ST->nbsoncells;
            ST->tabnodes[z].nbsons = 2;
            // the altitude of the new component is the altitude of the edge
            // under consideration
            STaltitude[z] = getweight(g, clefs[i]);
            // add x1 and x2 to the lists of sons of the new component
            newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);
            newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
            ST->tabnodes[z].sonlist = newsoncell1;
            newsoncell1->son = x1;
            newsoncell1->next = newsoncell2;
            newsoncell2->son = x2;
            newsoncell2->next = NULL;
            ST->tabnodes[z].lastson = newsoncell2;
            ST->nbsoncells += 2
            ;
            
            // then link n1 and n2
            k = TarjanLink(T, n1, n2);
            STmap[k] = z;
            
            //printf("MST node z:%i edge x,y %i - %i , compl,compr %i,%i \n",z, g->tete[clefs[i]],g->queue[clefs[i]], n1, n2);
            mst->MSTedges[mst->size].x= g->tete[clefs[i]];
            mst->MSTedges[mst->size].y= g->queue[clefs[i]];
            mst->MSTedges[mst->size].weight = getweight(g, clefs[i]);//already calculated weight of edge
            
            //  printf("MST edge %c - %c :  %f \n", 97+g->tete[clefs[i]],97+g->queue[clefs[i]],mst->MSTedges[mst->size].weight);
            
            
            mst->size+=1;
            
        }//if(...)
        else{
            // printf("Discarded edge x,y %i - %i , components %i,%i\n", g->tete[clefs[i]],g->queue[clefs[i]],n1,n2 );
        }
    }//for(...)
    
    
    for(i = 0; i < taille; i++)
        STaltitude[i] = 0; /* altitudes of the leaves is 0 */
    
    /* Construct father relationship */
    ST->tabnodes[ST->nbnodes-1].father = -1;
    ST->root = ST->nbnodes - 1;
    calculReversePointer(ST, ST->root);
    
    
    
    
    
    (*SaliencyTree) = ST;
    
    /* liberation de la memoire */
    free(STmap);
    free(clefs);
    TarjanTermine(T);
    return ST->nbnodes - 1;
}


/* The algorithm is an original adaptation/improvment to edge-weighted
 graph of the algorithm of Najman and Couprie published in IEEE IP
 in 2006 / Its description is provided in Cousty, Najman, Perret
 ISMM2013 and in Najman, Cousty, Perret ISMM 2013. */
int32_t BPTMST(/* INPUTS */
              graphe *g,
              /* OUTPUTS */
              JCctree ** SaliencyTree, /* Component tree - BPTAO  */
              double *STaltitude,  /* weights associated to the
                                    nodes of the BPTAO: Must
                                    be allocated elsewhere */
              MST* mst
              )
{
    int32_t i,x1,x2,n1,n2,z,k, STx, STy, nbsoncellsloc;
    JCctree *ST;
    int32_t *clefs;
    int32_t *STmap;
    JCsoncell * newsoncell1;
    JCsoncell * newsoncell2;
    Tarjan *T;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    
    if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de malloc\n");
    }
    
    if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
        fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
        exit(0);
    }
    
    setSizeMST(mst, narcs);
    
    
    for(i = 0; i < narcs; i++)
        clefs[i] = i;
    
    
    
    d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
    
    
    
    
    if( (ST = componentTreeAlloc((2*taille))) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
        exit(0);
    }
    ST->nbnodes = taille;
    ST->nbsoncells = 0;
    T = CreeTarjan(taille);
    for(i = 0; i < taille; i++)
    {
        STmap[i] = i;
        ST->tabnodes[i].nbsons = 0;
        ST->tabnodes[i].sonlist = NULL;
        TarjanMakeSet(T,i);
    }
    
    for(i = 0; i < narcs; i++){
        // for each edge of the graph or MST taken in increasing order of altitude
        n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
        if(n1 != n2) {
            /* If the two vertices do not belong to a same connected
             component at a lowest level */
            // Which component of ST n1 and n2 belongs to ?
            x1 = STmap[n1]; x2 = STmap[n2];
            // Create a new component
            z = ST->nbnodes; ST->nbnodes++;
            nbsoncellsloc = ST->nbsoncells;
            ST->tabnodes[z].nbsons = 2;
            // the altitude of the new component is the altitude of the edge
            // under consideration
            STaltitude[z] = getweight(g, clefs[i]);
            // add x1 and x2 to the lists of sons of the new component
            newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);
            newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
            ST->tabnodes[z].sonlist = newsoncell1;
            newsoncell1->son = x1;
            newsoncell1->next = newsoncell2;
            newsoncell2->son = x2;
            newsoncell2->next = NULL;
            ST->tabnodes[z].lastson = newsoncell2;
            ST->nbsoncells += 2;
            
            // then link n1 and n2
            k = TarjanLink(T, n1, n2);
            STmap[k] = z;
            
            //printf("MST node z:%i edge x,y %i - %i , compl,compr %i,%i \n",z, g->tete[clefs[i]],g->queue[clefs[i]], n1, n2);
            mst->MSTedges[mst->size].x= g->tete[clefs[i]];
            mst->MSTedges[mst->size].y= g->queue[clefs[i]];
            mst->MSTedges[mst->size].weight = getweight(g, clefs[i]);//already calculated weight of edge
            mst->MSTedges[mst->size].idx = clefs[i];
            //printf("MST edge %i - %i :  %f \n", g->tete[clefs[i]], g->queue[clefs[i]],mst->MSTedges[mst->size].weight);
            
            
            mst->size+=1;
            
        }//if(...)
        else{
            // printf("Discarded edge x,y %i - %i , components %i,%i\n", g->tete[clefs[i]],g->queue[clefs[i]],n1,n2 );
        }
    }//for(...)
    
    

    
    (*SaliencyTree) = ST;
    
    /* liberation de la memoire */
    free(STmap);
    free(clefs);
    TarjanTermine(T);
    return ST->nbnodes - 1;
}



/** Computing binary partition tree, especially for segmentation.
 STOrigAltitude,
 */
int32_t BPTAO_segmentation(/* INPUTS */
                           graphe *g,
                           /* OUTPUTS */
                           JCctree ** SaliencyTree, /* Component tree - BPTAO  */
                           double *STaltitude,  /* weights associated to the
                                                 nodes of the BPTAO: Must
                                                 be allocated elsewhere */
                           double *STOrigAltitude
                           )
{
    int32_t i,x1,x2,n1,n2,z,k, nbsoncellsloc;
    JCctree *ST;
    int32_t *clefs;
    int32_t *STmap;
    JCsoncell * newsoncell1;
    JCsoncell * newsoncell2;
    Tarjan *T;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    
    if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de malloc\n");
    }
    
    if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
        fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
        exit(0);
    }
    
    
    
    for(i = 0; i < narcs; i++)
        clefs[i] = i;
    d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
    if( (ST = componentTreeAlloc((2*taille))) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
        exit(0);
    }
    ST->nbnodes = taille;
    ST->nbsoncells = 0;
    T = CreeTarjan(taille);
    for(i = 0; i < taille; i++)
    {
        STmap[i] = i;
        ST->tabnodes[i].nbsons = 0;
        ST->tabnodes[i].sonlist = NULL;
        TarjanMakeSet(T,i);
    }
    
    for(i = 0; i < narcs; i++){
        // for each edge of the graph or MST taken in increasing order of altitude
        n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
        if(n1 != n2) {
            /* If the two vertices do not belong to a same connected
             component at a lowest level */
            // Which component of ST n1 and n2 belongs to ?
            x1 = STmap[n1]; x2 = STmap[n2];
            // Create a new component
            z = ST->nbnodes; ST->nbnodes++;
            nbsoncellsloc = ST->nbsoncells;
            ST->tabnodes[z].nbsons = 2;
            // the altitude of the new component is the altitude of the edge
            // under consideration
            STaltitude[z] = getweight(g, clefs[i]);
            
            if(STaltitude[z] >= INFINITO)
                STOrigAltitude[z] = 0;
            else
                STOrigAltitude[z] = getweightOriginal(g, clefs[i]);
            // add x1 and x2 to the lists of sons of the new component
            newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);
            newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
            ST->tabnodes[z].sonlist = newsoncell1;
            newsoncell1->son = x1;
            newsoncell1->next = newsoncell2;
            newsoncell2->son = x2;
            newsoncell2->next = NULL;
            ST->tabnodes[z].lastson = newsoncell2;
            ST->nbsoncells += 2
            ;
            
            // then link n1 and n2
            k = TarjanLink(T, n1, n2);
            STmap[k] = z;
            
            
        }//if(...)
        else{
            // printf("Discarded edge x,y %i - %i , components %i,%i\n", g->tete[clefs[i]],g->queue[clefs[i]],n1,n2 );
        }
    }//for(...)
    
    
    for(i = 0; i < taille; i++)
        STaltitude[i] = 0; /* altitudes of the leaves is 0 */
    
    /* Construct father relationship */
    ST->tabnodes[ST->nbnodes-1].father = -1;
    ST->root = ST->nbnodes - 1;
    calculReversePointer(ST, ST->root);
    
    
    (*SaliencyTree) = ST;
    
    /* liberation de la memoire */
    free(STmap);
    free(clefs);
    TarjanTermine(T);
    return ST->nbnodes - 1;
}


int32_t BPTAO_segmentationMin(/* INPUTS */
                           graphe *g,
                           /* OUTPUTS */
                           JCctree ** SaliencyTree, /* Component tree - BPTAO  */
                           double *STaltitude,  /* weights associated to the
                                                 nodes of the BPTAO: Must
                                                 be allocated elsewhere */
                           double *STOrigAltitude
                           )
{
    int32_t i,x1,x2,n1,n2,z,k, nbsoncellsloc;
    JCctree *ST;
    int32_t *clefs;
    int32_t *STmap;
    JCsoncell * newsoncell1;
    JCsoncell * newsoncell2;
    Tarjan *T;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    
    if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de malloc\n");
    }
    
    if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
        fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
        exit(0);
    }
    
    
    
    for(i = 0; i < narcs; i++)
        clefs[i] = i;
    d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
    if( (ST = componentTreeAlloc((2*taille))) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
        exit(0);
    }
    ST->nbnodes = taille;
    ST->nbsoncells = 0;
    T = CreeTarjan(taille);
    for(i = 0; i < taille; i++)
    {
        STmap[i] = i;
        ST->tabnodes[i].nbsons = 0;
        ST->tabnodes[i].sonlist = NULL;
        TarjanMakeSet(T,i);
    }
    
    for(i = 0; i < narcs; i++){
        // for each edge of the graph or MST taken in increasing order of altitude
        n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
        if(n1 != n2) {
            /* If the two vertices do not belong to a same connected
             component at a lowest level */
            // Which component of ST n1 and n2 belongs to ?
            x1 = STmap[n1]; x2 = STmap[n2];
            // Create a new component
            z = ST->nbnodes; ST->nbnodes++;
            nbsoncellsloc = ST->nbsoncells;
            ST->tabnodes[z].nbsons = 2;
            // the altitude of the new component is the altitude of the edge
            // under consideration
            STaltitude[z] = getweight(g, clefs[i]);
            
            if(STaltitude[z] >= INFINITO)
                //STOrigAltitude[z] = INFINITO;
                STOrigAltitude[z] = 0;
            else
                STOrigAltitude[z] = getweightOriginal(g, clefs[i]);
            // add x1 and x2 to the lists of sons of the new component
            newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);
            newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
            ST->tabnodes[z].sonlist = newsoncell1;
            newsoncell1->son = x1;
            newsoncell1->next = newsoncell2;
            newsoncell2->son = x2;
            newsoncell2->next = NULL;
            ST->tabnodes[z].lastson = newsoncell2;
            ST->nbsoncells += 2
            ;
            
            // then link n1 and n2
            k = TarjanLink(T, n1, n2);
            STmap[k] = z;
            
            
        }//if(...)
        else{
            // printf("Discarded edge x,y %i - %i , components %i,%i\n", g->tete[clefs[i]],g->queue[clefs[i]],n1,n2 );
        }
    }//for(...)
    
    
    for(i = 0; i < taille; i++)
        STaltitude[i] = 0; /* altitudes of the leaves is 0 */
    
    /* Construct father relationship */
    ST->tabnodes[ST->nbnodes-1].father = -1;
    ST->root = ST->nbnodes - 1;
    calculReversePointer(ST, ST->root);
    
    
    (*SaliencyTree) = ST;
    
    /* liberation de la memoire */
    free(STmap);
    free(clefs);
    TarjanTermine(T);
    return ST->nbnodes - 1;
}



void compactRec(JCctree *CT, double *Val, int32_t root){
    JCsoncell *tmp;
    JCsoncell *f1 = CT->tabnodes[root].sonlist;
    JCsoncell *f2 = CT->tabnodes[root].lastson;
    
    if (f1 == NULL){ /* A leaf of the tree, that is a vertex of the initial graph */
        Val[root] = Val[CT->tabnodes[root].father];
    } else {
        compactRec(CT, Val, f1->son);
        compactRec(CT, Val, f2->son);
        
        if( (Val[root] == Val[f1->son]) && (f1->son > CT->nbnodes/2) ){
            /* Then f1 must be compacted */
            CT->tabnodes[root].sonlist = CT->tabnodes[f1->son].sonlist;
            tmp = CT->tabnodes[f1->son].lastson;
        } else tmp = f1;
        
        if( (Val[root] == Val[f2->son]) && (f2->son > CT->nbnodes/2) ){
            /* Then f1 must be compacted */
            tmp->next = CT->tabnodes[f2->son].sonlist;
            CT->tabnodes[root].lastson = CT->tabnodes[f2->son].lastson;
        } else  tmp->next = f2;
    }
}


int32_t QFZ(/* INPUTS */
            graphe *g,
            /* OUTPUTS */
            JCctree **CompTree, /* Component tree */
            double *CTaltitude  /* weights associated to the
                                 nodes of the component
                                 tree : Must be allocated
                                 elsewhere */
)
{
    //QFZ(g, CompTree, CTaltitude);
    //compactRec((*CompTree), CTaltitude, (*CompTree)->root);
    //calculReversePointer((*CompTree), (*CompTree)->root);
    return 1;
}







////////////////// Edward


/** Calculates MST, outputs JCctree St, MST
 */
int getMST(graphe *g, MST** mst ){
    
    double *STAltitudes;
    JCctree *ST;
    
    *mst = (MST*)malloc(sizeof(MST));
    
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr,"getQBT cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    
    // To obtain the binary partition tree by altitude ordering
    BPTMST(g, &ST, STAltitudes, *mst);
    free(STAltitudes);
    componentTreeFree(ST);
    return 0;
}



double computeSM(JCctree *ST, graphe *g, double *STAltitudes){

   // Compute Saliency map
    
    int32_t logn, nbRepresent, u, x, y, c1;
    /* Data structure for fast computation of LCA (least common ancestor) */
    int32_t *Euler, *Depth, *Represent, *Number, **Minim;
    //FILE *BitmapFile;
    Euler = (int32_t *)calloc(2*ST->nbnodes-1, sizeof(int32_t));
    Represent = (int32_t *)calloc(2*ST->nbnodes-1, sizeof(int32_t));
    Depth = (int32_t *)calloc(ST->nbnodes, sizeof(int32_t));
    Number = (int32_t *)calloc(ST->nbnodes, sizeof(int32_t));
    
    int32_t nbarcs = g->ind;
    
    if ((Euler == NULL) || (Represent == NULL)
        || (Depth == NULL) || (Number == NULL)) {
        fprintf(stderr, "%s : malloc failed\n", F_NAME);
        return -1;
    }
    
    double sum = 0 ;
    
    Minim = LCApreprocess(ST, Euler, Depth, Represent, Number, &nbRepresent, &logn);
    
    // For any edge of g
    for(u = 0; u < nbarcs; u++){
        x = g->tete[u]; y = g->queue[u];
        c1 = Represent[LowComAncFast(x, y, Euler, Number, Depth, Minim)];
        setweight(g,u,STAltitudes[c1]);
        
        sum = sum + STAltitudes[c1];
    }
    free(Euler);
    free(Represent);
    free(Depth);
    free(Number);
    free(Minim[0]);
    free(Minim);
    
    
    return sum;
   // componentTreeFree(ST);
   // free(STAltitudes);
    
    
    
    /*
    if(argc == 4){
        
        // Then the bitmap is loaded and applied to the edge weights
        if( (BitmapFile = fopen(argv[3],"r")) != NULL){
            fscanf(BitmapFile, "#NbMin: %d\n", &nbmin);
            Bitmap = calloc(nbmin+1, sizeof(double));
            if(Bitmap == NULL){
                fprintf(stderr, "%d cannod allocate (%d byte) memory\n", argv[0], nbmin);
                exit(1);
            }
            Bitmap[0] = 0;
            for(i = 0; i < nbmin; i++){
                fscanf(BitmapFile, "%lf\n", &(Bitmap[i+1]));
            }
            fclose(BitmapFile);
            
            for(u = 0; u < g->ind; u++){
                i = (int32_t)getweight(g, u);
                setweight(g,u, Bitmap[i]);
            }
            free(Bitmap);
        } else 
            fprintf(stderr,"%s cannot open %s\n", argv[0], argv[3]);
    }
     */
    
  

}


/** OUTPUTS: Edge of the corresponding n-th node */
int getEdge(/* INPUTS */
            JCctree **CompTree, /* Component tree */
            int n){
    
    JCctree * ST = *CompTree;
    
    if( ST->tabnodes[n].nbsons > 0){
        return n - ST->nbnodes;
    }
    
    
    return -1;
}

double weightNode(/* INPUTS */
               JCctree *ST, /* Component tree */
               double *STaltitude,
               int n){
    
    if( ST->tabnodes[n].nbsons > 0){
        return STaltitude[n];
    }
    return -1;
}


int CanonizeQBT(/* INPUTS */
                JCctree *ST, /* Component tree - QBT  */
                double *STaltitude,
                /* OUTPUTS */
                JCctree **CompTreeOut){
    
    
    JCctree* QCT;
    int i,p;
    if( (QCT = componentTreeAlloc((2* (ST->nbnodes)))) == NULL){
        fprintf(stderr, "QCT: error in Memory Alloc\n");
        exit(0);
    }
    QCT->nbnodes = ST->nbnodes;
    QCT->nbsoncells = 0;
    
    for(i = 0; i< QCT->nbnodes ; i++){
        QCT->tabnodes[i].nbsons = 0;
        QCT->tabnodes[i].sonlist = NULL;
        QCT->tabnodes[i].lastson = NULL;
        QCT->tabnodes[i].father = -99;
    }
    
    JCsoncell* c;
    
    for(i=0;i<ST->nbnodes;i++)
        QCT->tabnodes[i].father = ST->tabnodes[i].father;
    
    for(i = ST->nbnodes-1; i>=0; i--){
        if(ST->tabnodes[i].nbsons>0 && ST->tabnodes[i].father!= -1){
            p = QCT->tabnodes[i].father;            
            if( weightNode(ST,STaltitude,p) == weightNode(ST,STaltitude,i)) {
                
                c = ST->tabnodes[i].sonlist;
                while(c != NULL){
                    QCT->tabnodes[ c->son ].father = p;
                    c = c->next;
                }
                QCT->tabnodes[i].father  = i;
            }
        }
    }
    JCsoncell * newsoncell1;
    
    
    QCT->root = ST->root;
    
    // Building the children
    for(i = 0; i< QCT->nbnodes ; i++){
        p = QCT->tabnodes[i].father;
        if(p>=0 && p!=i ){
            newsoncell1 = &(QCT->tabsoncells[ QCT->nbsoncells]);
            newsoncell1->son = i;
            newsoncell1->next= NULL;
            
            if(QCT->tabnodes[p].lastson==NULL){
                QCT->tabnodes[p].sonlist = newsoncell1;
            }
            else{
                QCT->tabnodes[p].lastson->next = newsoncell1;
            }
            QCT->tabnodes[p].lastson = newsoncell1;
            QCT->tabnodes[p].nbsons+=1;
            QCT->nbsoncells += 1;
        }
    }
    
    *CompTreeOut = QCT;
    
    return 0;
}



int CanonizeQBTSegmentation(
                            JCctree *ST,
                            double *STaltitude,
 
                            JCctree **CompTreeOut,
                            double *STMaxWeight){
    
    
    JCctree* QCT;
    int i,p;
    if( (QCT = componentTreeAlloc((2* (ST->nbnodes)))) == NULL){
        fprintf(stderr, "QCT: error in Memory Alloc\n");
        exit(0);
    }
    QCT->nbnodes = ST->nbnodes;
    QCT->nbsoncells = 0;
    
    for(i = 0; i< QCT->nbnodes ; i++){
        QCT->tabnodes[i].nbsons = 0;
        QCT->tabnodes[i].sonlist = NULL;
        QCT->tabnodes[i].lastson = NULL;
        QCT->tabnodes[i].father = -99;
    }
    
    JCsoncell* c;
    
    for(i=0;i<ST->nbnodes;i++)
        QCT->tabnodes[i].father = ST->tabnodes[i].father;
    double maxTemp = 0;
    for(i = ST->nbnodes-1; i>=0; i--){
        if(ST->tabnodes[i].nbsons>0 && ST->tabnodes[i].father!= -1){
            p = QCT->tabnodes[i].father;
            if( weightNode(ST,STaltitude,p) == weightNode(ST,STaltitude,i)) {
                maxTemp = STMaxWeight[i];
                c = ST->tabnodes[i].sonlist;
                while(c != NULL){
                    QCT->tabnodes[ c->son ].father = p;
                    if(STMaxWeight[c->son]>maxTemp)maxTemp = STMaxWeight[c->son];
                    c = c->next;
                }
                QCT->tabnodes[i].father  = i;
                if(STMaxWeight[p] != 0 && maxTemp > STMaxWeight[p]  )
                    STMaxWeight[p]  = maxTemp;
            }
        }
    }
    JCsoncell * newsoncell1;
    
    
    QCT->root = ST->root;
    
    // Building the children
    for(i = 0; i< QCT->nbnodes ; i++){
        p = QCT->tabnodes[i].father;
        if(p>=0 && p!=i ){
            newsoncell1 = &(QCT->tabsoncells[ QCT->nbsoncells]);
            newsoncell1->son = i;
            newsoncell1->next= NULL;
            
            if(QCT->tabnodes[p].lastson==NULL){
                QCT->tabnodes[p].sonlist = newsoncell1;
            }
            else{
                QCT->tabnodes[p].lastson->next = newsoncell1;
            }
            QCT->tabnodes[p].lastson = newsoncell1;
            QCT->tabnodes[p].nbsons+=1;
            QCT->nbsoncells += 1;
        }
    }
    
    *CompTreeOut = QCT;
    
    return 0;

}


int CanonizeQBTSegmentationMin(
                            JCctree *ST,
                            double *STaltitude,
                            
                            JCctree **CompTreeOut,
                            double *STMinWeight){
    
    
    JCctree* QCT;
    int i,p;
    if( (QCT = componentTreeAlloc((2* (ST->nbnodes)))) == NULL){
        fprintf(stderr, "QCT: error in Memory Alloc\n");
        exit(0);
    }
    QCT->nbnodes = ST->nbnodes;
    QCT->nbsoncells = 0;
    
    for(i = 0; i< QCT->nbnodes ; i++){
        QCT->tabnodes[i].nbsons = 0;
        QCT->tabnodes[i].sonlist = NULL;
        QCT->tabnodes[i].lastson = NULL;
        QCT->tabnodes[i].father = -99;
    }
    
    JCsoncell* c;
    
    for(i=0;i<ST->nbnodes;i++)
        QCT->tabnodes[i].father = ST->tabnodes[i].father;
    double minTemp = 0;
    for(i = ST->nbnodes-1; i>=0; i--){
        if(ST->tabnodes[i].nbsons>0 && ST->tabnodes[i].father!= -1){
            p = QCT->tabnodes[i].father;
            if( weightNode(ST,STaltitude,p) == weightNode(ST,STaltitude,i)) {
                minTemp = STMinWeight[i];
                c = ST->tabnodes[i].sonlist;
                while(c != NULL){
                    QCT->tabnodes[ c->son ].father = p;
                    if(STMinWeight[c->son]!=0 && STMinWeight[c->son]<minTemp)minTemp = STMinWeight[c->son];
                   // if(ST->tabnodes[c->son].nbsons>0  && STMinWeight[c->son]<minTemp)minTemp = STMinWeight[c->son];
                    c = c->next;
                }
                QCT->tabnodes[i].father  = i;
                if(STMinWeight[p] != 0 && minTemp < STMinWeight[p]  )
                    STMinWeight[p]  = minTemp;
            }
        }
    }
    JCsoncell * newsoncell1;
    
    
    QCT->root = ST->root;
    
    // Building the children
    for(i = 0; i< QCT->nbnodes ; i++){
        p = QCT->tabnodes[i].father;
        if(p>=0 && p!=i ){
            newsoncell1 = &(QCT->tabsoncells[ QCT->nbsoncells]);
            newsoncell1->son = i;
            newsoncell1->next= NULL;
            
            if(QCT->tabnodes[p].lastson==NULL){
                QCT->tabnodes[p].sonlist = newsoncell1;
            }
            else{
                QCT->tabnodes[p].lastson->next = newsoncell1;
            }
            QCT->tabnodes[p].lastson = newsoncell1;
            QCT->tabnodes[p].nbsons+=1;
            QCT->nbsoncells += 1;
        }
    }
    
    *CompTreeOut = QCT;
    
    return 0;
    
}


int CanonizeQBTSegmentationMean(
                               JCctree *ST,
                               double *STaltitude,
                               
                               JCctree **CompTreeOut,
                               double *STMinWeight){
    
    
    JCctree* QCT;
    int i,p;
    if( (QCT = componentTreeAlloc((2* (ST->nbnodes)))) == NULL){
        fprintf(stderr, "QCT: error in Memory Alloc\n");
        exit(0);
    }
    QCT->nbnodes = ST->nbnodes;
    QCT->nbsoncells = 0;
    
    for(i = 0; i< QCT->nbnodes ; i++){
        QCT->tabnodes[i].nbsons = 0;
        QCT->tabnodes[i].sonlist = NULL;
        QCT->tabnodes[i].lastson = NULL;
        QCT->tabnodes[i].father = -99;
    }
    
    JCsoncell* c;
    
    for(i=0;i<ST->nbnodes;i++)
        QCT->tabnodes[i].father = ST->tabnodes[i].father;
    
    double minTemp = 0;
    double tmp = 0;
    int count = 0;
    for(i = ST->nbnodes-1; i>=0; i--){
        if(ST->tabnodes[i].nbsons>0 && ST->tabnodes[i].father!= -1){
            count = 0;
            p = QCT->tabnodes[i].father;
            if( weightNode(ST,STaltitude,p) == weightNode(ST,STaltitude,i)) {
                count+=1;
                tmp+= STMinWeight[i];
                
                c = ST->tabnodes[i].sonlist;
                while(c != NULL){
                    count+=1;
                    
                    QCT->tabnodes[ c->son ].father = p;
                    tmp+= STMinWeight[c->son];
                    //if(STMinWeight[c->son]!=0 && STMinWeight[c->son]<minTemp)minTemp = STMinWeight[c->son];
                    // if(ST->tabnodes[c->son].nbsons>0  && STMinWeight[c->son]<minTemp)minTemp = STMinWeight[c->son];
                    c = c->next;
                }
                QCT->tabnodes[i].father  = i;
                
                tmp = tmp/count;
                STMinWeight[p]  = (tmp+STMinWeight[p])/2;
            }
        }
    }
    JCsoncell * newsoncell1;
    
    
    QCT->root = ST->root;
    
    // Building the children
    for(i = 0; i< QCT->nbnodes ; i++){
        p = QCT->tabnodes[i].father;
        if(p>=0 && p!=i ){
            newsoncell1 = &(QCT->tabsoncells[ QCT->nbsoncells]);
            newsoncell1->son = i;
            newsoncell1->next= NULL;
            
            if(QCT->tabnodes[p].lastson==NULL){
                QCT->tabnodes[p].sonlist = newsoncell1;
            }
            else{
                QCT->tabnodes[p].lastson->next = newsoncell1;
            }
            QCT->tabnodes[p].lastson = newsoncell1;
            QCT->tabnodes[p].nbsons+=1;
            QCT->nbsoncells += 1;
        }
    }
    
    *CompTreeOut = QCT;
    
    return 0;
    
}


/* ==================================== */
int componentTreeAttribute(JCctree *CT, double *Alt, int32_t root, int32_t *Area, double * maxWeight)
/* ==================================== */
{
    int nvertices = 0;
    JCsoncell *s;
    if(CT->tabnodes[root].sonlist == NULL){
        Area[root] =1;
        maxWeight[root] = 0;
        return 1; //leaf node, vertex
    }
    double max = maxWeight[root];
        
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
        nvertices+=componentTreeAttribute(CT, Alt, s->son, Area, maxWeight);
        
        if(max!= 0 && maxWeight[s->son] > max){
            max = maxWeight[s->son];
        }
    }
    
    
    
    Area[root] = nvertices;
    maxWeight[root] = max;
    
    return nvertices;
}

/* ==================================== */
int componentTreeAttributeMin(JCctree *CT, double *Alt, int32_t root, int32_t *Area, double * minWeight)
/* ==================================== */
{
    int nvertices = 0;
    JCsoncell *s;
    if(CT->tabnodes[root].sonlist == NULL){
        Area[root] =1;
        minWeight[root] = 0;
        return 1; //leaf node, vertex
    }
    double min = minWeight[root];
    
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
        nvertices+=componentTreeAttributeMin(CT, Alt, s->son, Area, minWeight);
        
         if(min!= 0 && minWeight[s->son] !=0 && minWeight[s->son] < min){
         //if(min!= 0 &&  CT->tabnodes[s->son].nbsons >0  && minWeight[s->son] < min){
            min = minWeight[s->son];
        }
    }
    
    
    
    Area[root] = nvertices;
    minWeight[root] = min;
    
    return nvertices;
}


/* ==================================== */
int componentTreeAttributeMean(JCctree *CT, double *Alt, int32_t root, int32_t *Area, double * minWeight)
/* ==================================== */
{
    int nvertices = 0;
    JCsoncell *s;
    if(CT->tabnodes[root].sonlist == NULL){
        Area[root] =1;
        minWeight[root] = 0;
        return 1; //leaf node, vertex
    }
  
    
    for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
        nvertices+=componentTreeAttributeMean(CT, Alt, s->son, Area, minWeight);
    }

    Area[root] = nvertices;
    
    return nvertices;
}






int getComponentAtLambda(JCctree *CT, double *Alt, const double lambda, const int vertex){
  
    /*
    FILE *log = NULL;
    log = fopen("/tmp/jumps.txt", "a");
    
    */
    
    int x = vertex;
    int par_x = CT->tabnodes[vertex].father;
    
  //  fprintf(log, " %i: ", x);
    
    while( par_x!= -1 && !(Alt[par_x] >= lambda )  ){
        x = par_x;
//        fprintf(log, " %i-> ", x);
        par_x = CT->tabnodes[x].father;
    
    }
    
    //fprintf(log, " \n ");
    //fclose(log);
    
    return x;
}


void computeQFZ_Attributes(graphe *g , JCctree **CompTree, double **STAltitudes , int** STArea , double** STMaxWeight)
{
    JCctree *QCT = *CompTree;
    
    if( QCT != NULL){
        componentTreeFree(  QCT );
        free(*STAltitudes);
        free(*STMaxWeight);
        free(*STArea);
    }
    
    *STAltitudes = NULL;
    *STMaxWeight = NULL;
    *STArea = NULL;
    JCctree *ST = NULL;
    
    *STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    *STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    *STArea = (int *)calloc(2*g->nbsom, sizeof(int));
    
    if (STAltitudes == NULL){
        fprintf(stderr,"computeQFZ_Attributes cannot allocate memory for STAltitudes\n");
        exit(0);
    } // To obtain the binary partition tree by altitude ordering
    
    BPTAO_segmentation(g, &ST, *STAltitudes, *STMaxWeight );
    // Computing the Quase Flat Zones
    CanonizeQBTSegmentation(ST, *STAltitudes, &QCT,*STMaxWeight );
    //CanonizeQBT(ST, STAltitudes, &QCT );
    
    // Computing the cardinality of the components
    componentTreeAttribute(QCT, *STAltitudes, QCT->root, *STArea, *STMaxWeight);
    
    *CompTree = QCT;
    
    componentTreeFree(ST);

}

void computeQFZ_AttributesMin(graphe *g , JCctree **CompTree, double **STAltitudes , int** STArea , double** STMinWeight)
{
    JCctree *QCT = *CompTree;
    
    if( QCT != NULL){
        componentTreeFree(  QCT );
        free(*STAltitudes);
        free(*STMinWeight);
        free(*STArea);
    }
    
    *STAltitudes = NULL;
    *STMinWeight = NULL;
    *STArea = NULL;
    JCctree *ST = NULL;
    
    *STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    *STMinWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    *STArea = (int *)calloc(2*g->nbsom, sizeof(int));
    
    if (STAltitudes == NULL){
        fprintf(stderr,"computeQFZ_Attributes cannot allocate memory for STAltitudes\n");
        exit(0);
    } // To obtain the binary partition tree by altitude ordering
    
    BPTAO_segmentationMin(g, &ST, *STAltitudes, *STMinWeight );
    // Computing the Quase Flat Zones
    CanonizeQBTSegmentationMin(ST, *STAltitudes, &QCT,*STMinWeight );
    //CanonizeQBT(ST, STAltitudes, &QCT );
    
    // Computing the cardinality of the components
    componentTreeAttributeMin(QCT, *STAltitudes, QCT->root, *STArea, *STMinWeight);
    
    *CompTree = QCT;
    
    componentTreeFree(ST);
    
}


void computeQFZ_AttributesMean(graphe *g , JCctree **CompTree, double **STAltitudes , int** STArea , double** STMinWeight)
{
    JCctree *QCT = *CompTree;
    
    if( QCT != NULL){
        componentTreeFree(  QCT );
        free(*STAltitudes);
        free(*STMinWeight);
        free(*STArea);
    }
    
    *STAltitudes = NULL;
    *STMinWeight = NULL;
    *STArea = NULL;
    JCctree *ST = NULL;
    
    *STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    *STMinWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    *STArea = (int *)calloc(2*g->nbsom, sizeof(int));
    
    if (STAltitudes == NULL){
        fprintf(stderr,"computeQFZ_Attributes cannot allocate memory for STAltitudes\n");
        exit(0);
    } // To obtain the binary partition tree by altitude ordering
    
    BPTAO_segmentationMin(g, &ST, *STAltitudes, *STMinWeight );
    // Computing the Quase Flat Zones
    CanonizeQBTSegmentationMean(ST, *STAltitudes, &QCT,*STMinWeight );
    //CanonizeQBT(ST, STAltitudes, &QCT );
    
    // Computing the cardinality of the components
    componentTreeAttributeMean(QCT, *STAltitudes, QCT->root, *STArea, *STMinWeight);
    
    *CompTree = QCT;
    
    componentTreeFree(ST);
    
}



//////////// Incremental

int initializeTree(/* INPUTS */
                   graphe *g,
                   /* OUTPUTS */
                   JCctree **CompTree,/* CompTree :  initalized component tree - QFZ */
                   double  *STaltitude
                   ){
    
    int i ,nbsoncellsloc;
    JCsoncell * newsoncell1, *lastchild;
    JCctree *CT = *CompTree; //component tree
    int32_t taille = g->nbsom;
    
    CT->nbnodes = taille;
    CT->nbsoncells = 0;
    
    for(i = 0; i < taille; i++){
        CT->tabnodes[i].nbsons = 0;
        CT->tabnodes[i].sonlist = NULL;
        CT->tabnodes[i].lastson = NULL;
        CT->tabnodes[i].nodeAsSon = NULL;
    }
    
    for(i = 0; i < taille; i++)
        STaltitude[i] = 0; /* altitudes of the leaves is 0 */
    
    int z;
    
    z = CT->nbnodes; CT->nbnodes++;
    STaltitude[z] = 999999;
    CT->tabnodes[z].nbsons = taille;
    CT->tabnodes[z].nodeAsSon = NULL;
    
//// Create root node of the tree z and all vertices as their initial children.
    CT->root = z;
    nbsoncellsloc = CT->nbsoncells;
    newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
    newsoncell1->son = 0;
    newsoncell1->next = NULL;
    newsoncell1->last = NULL;
    
    CT->tabnodes[z].sonlist = newsoncell1;
    CT->tabnodes[z].lastson = newsoncell1;
    CT->tabnodes[z].father = -1;
    CT->tabnodes[0].father = z;
    CT->tabnodes[0].nodeAsSon = newsoncell1;
    CT->nbsoncells += 1;
    
    lastchild =newsoncell1;     
    for(i = 1; i < taille; i++){
        nbsoncellsloc = CT->nbsoncells;
        newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
        newsoncell1->son = i;
        newsoncell1->next = NULL;
        CT->tabnodes[i].father = CT->root;
        CT->nbsoncells += 1;
        CT->tabnodes[i].nodeAsSon = newsoncell1;
        newsoncell1->last = lastchild;
        lastchild->next =newsoncell1;
        lastchild = newsoncell1;
    }
    CT->tabnodes[z].lastson = lastchild;
    
    return 0;
}



void findTransition(/* INPUTS */
                    JCctree *QFZ,
                    graphe *g,
                    double  *STaltitude,
                    int leafIdx, /* leaf index for Edge{x,y}*/
                    int edgeIdx,
                    /* OUTPUTS */
                    int *c,
                    int *p
                    ){
    
    int idxSon = leafIdx;
    int idxFather =  QFZ->tabnodes[idxSon].father;
    
   // if(STaltitude[idxSon] <= g->weight[edgeIdx]+1 && !(STaltitude[idxFather]<= g->weight[edgeIdx]+1)){
    double lambda = g->weight[edgeIdx]+1;
    //double lambda = g->weight[edgeIdx];
     if(STaltitude[idxSon] <= lambda && !(STaltitude[idxFather]<= lambda)){
        *c = idxSon; *p = idxFather;
        return;
    }
    else{
        findTransition(QFZ, g,STaltitude, idxFather , edgeIdx , c,p );
    }
    
}


void detach(JCctree **CompTree, /* Component tree - QFZ */
            double  *STaltitudes,
            const int idxP,   // idxS is dettached from idxS
            const int idxS ){
    
    JCctree * CT = *CompTree;
    JCsoncell *node = CT->tabnodes[idxS].nodeAsSon;
    
    if(CT->tabnodes[idxP].sonlist == node && CT->tabnodes[idxP].lastson == node ){
        // We are dettaching the last son of idxP, we have to mark it to delete
        CT->tabnodes[idxP].sonlist = NULL;
        CT->tabnodes[idxP].lastson = NULL;
        CT->tabnodes[idxP].nbsons = 0;
       // printf("idxP %i se quedo sin hijos \n", idxP);
        return; //detach
    }
    if(CT->tabnodes[idxP].sonlist == node && CT->tabnodes[idxP].lastson != node ){
        // There are still nodes left
        CT->tabnodes[idxP].sonlist = node->next;
        node->next->last = NULL;
        CT->tabnodes[idxP].nbsons-=1;
        return;
    }
    else{
        // it is the last;
        if(CT->tabnodes[idxP].lastson == node){
            CT->tabnodes[idxP].lastson = node->last;
            node->last->next = NULL;
            CT->tabnodes[idxP].nbsons-=1;
        }// it is in the middle;
        else{
            node->last->next = node->next;
            node->next->last = node->last;
            CT->tabnodes[idxP].nbsons-=1;
        }
    }
}


void attach(/* INPUTS */
            JCctree **CompTree, /* Component tree - QFZ */
            double  *STaltitudes,
            const int idxP,
            const int idxS //  idxS is attached to idxP
){
    JCctree * CT = *CompTree;
    //int  nbsoncellsloc;
    JCsoncell * newsoncell1 ;
    
    if(CT->tabnodes[idxS].father == idxP)return;//no need to attach as it is already son of idxP
    if(idxS == idxP){return;}//no attaching to its self}

    CT->tabnodes[idxP].nbsons += 1;
    
    detach(CompTree, STaltitudes , CT->tabnodes[idxS].father,  idxS);
    
    newsoncell1 = CT->tabnodes[idxS].nodeAsSon;
    newsoncell1->next = NULL;
    newsoncell1->last = CT->tabnodes[idxP].lastson;
    
    if(CT->tabnodes[idxP].lastson != NULL){
        CT->tabnodes[idxP].lastson->next = newsoncell1 ;
        CT->tabnodes[idxP].lastson = newsoncell1;
    }
    else{
        CT->tabnodes[idxP].lastson = newsoncell1;
        CT->tabnodes[idxP].sonlist = newsoncell1;
    }
    CT->tabnodes[idxS].father = idxP;
    
}




void merge(/* INPUTS */
           JCctree **CompTree, /* Component tree - QFZ */
           double  *STaltitudes,
           int p1,
           int p2,
           int *np
           ){
    
    if(p1==p2){*np = p1; return;}
    
    JCctree * CT = *CompTree;
    int larger, smaller, i ;
    if(CT->tabnodes[p1].nbsons > CT->tabnodes[p2].nbsons){
        larger = p1;
        smaller = p2;
    }
    else{
        larger = p2;
        smaller = p1;
    }
    // Children of the smaller one will go into the larger one.
    JCsoncell * child = CT->tabnodes[smaller].sonlist;
    
    if(CT->tabnodes[larger].lastson){
        CT->tabnodes[larger].lastson->next = child ;
        child->last = CT->tabnodes[larger].lastson;
    }
    else{
        // in case larger does not have a last son
        CT->tabnodes[larger].sonlist = child ;
        printf("larger %i has  %i sons, merging with %i with %i sons \n", larger, CT->tabnodes[larger].nbsons,smaller , CT->tabnodes[smaller].nbsons );
    }
    
    CT->tabnodes[larger].lastson = CT->tabnodes[smaller].lastson;
    
    for(i=0; i< CT->tabnodes[smaller].nbsons; i++ ){
        CT->tabnodes[child->son].father = larger;
        CT->tabnodes[larger].nbsons += 1;
        child = child->next;
    }
    
    detach(CompTree, STaltitudes , CT->tabnodes[smaller].father,  smaller);
    
    CT->tabnodes[smaller].sonlist = NULL;
    CT->tabnodes[smaller].lastson = NULL;
    CT->tabnodes[smaller].father = smaller; //marked for delete
    CT->tabnodes[smaller].nbsons = 0;
    
    *np = larger;
    
}




void createNode( /* INPUTS */
                JCctree **CompTree, /* Component tree - QFZ */
                double  *STaltitudes,
                graphe *g,
                int edgeIdx, /* Edge index for {x,y}  */
                int c1,
                int c2,
                /* OUTPUTS */
                int *n ){
    JCctree * CT = *CompTree;
    int z, nbsoncellsloc;
    double lambda;
    lambda = g->weight[edgeIdx]+1;
    
    if(STaltitudes[c1] ==  STaltitudes[c2] &&STaltitudes[c1] == lambda  && CT->tabnodes[c1].nbsons > 0 && CT->tabnodes[c2].nbsons >0 ){
        //merge
        merge(CompTree , STaltitudes, c1,c2, n);
        return ;
    }
    if(STaltitudes[c1] == lambda && STaltitudes[c2] < lambda ){
        attach(CompTree , STaltitudes, c1 , c2);
        *n = c1;
        return ;
    }
    if(STaltitudes[c2] == lambda && STaltitudes[c1] < lambda ){
        attach(CompTree , STaltitudes, c2 , c1);
        *n = c2;
        return ;
    }
    
    JCsoncell *newsoncell1;
    z = CT->nbnodes; CT->nbnodes++;
    
    nbsoncellsloc = CT->nbsoncells;
    
    CT->tabnodes[z].lastson = NULL;
    CT->tabnodes[z].sonlist = NULL;
    CT->tabnodes[z].nbsons = 0;
    
    STaltitudes[z] = g->weight[edgeIdx]+1;
    
    attach(CompTree, STaltitudes , z , c1);
    attach(CompTree, STaltitudes , z , c2);
    
    newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
    newsoncell1->son = z;
    newsoncell1->next =NULL;
    newsoncell1->last =NULL;
    
    CT->tabnodes[z].father = CT->root;
    CT->tabnodes[CT->tabnodes[z].father].nbsons += 1;
    CT->tabnodes[z].nodeAsSon = newsoncell1;
    
    if(CT->tabnodes[CT->root].lastson != NULL){
        CT->tabnodes[CT->tabnodes[z].father].lastson->next = newsoncell1 ;
        newsoncell1->last = CT->tabnodes[CT->tabnodes[z].father].lastson ;
        CT->tabnodes[CT->tabnodes[z].father].lastson = newsoncell1;
    }
    else{ // root was left with no children
        CT->tabnodes[CT->root].lastson = newsoncell1 ;
        CT->tabnodes[CT->root].sonlist = newsoncell1 ;
    }
    
    *n = z;
    CT->nbsoncells+=1;
}


int treeUpdate( /* INPUTS */
               JCctree **CompTree, /* Component tree - QFZ */
               double  *STaltitudes,
               graphe *g,
               int edgeIdx /* Edge index for {x,y}  */
/* OUTPUTS */
//       CompTree
){
    
    int c1,p1,c2,p2, n, tmp;
    int p1p, p2p, np;
    //printf("{%i,%i} and w=%f \n" , g->tete[edgeIdx], g->queue[edgeIdx], g->weight[edgeIdx]);
    findTransition(*CompTree, g,STaltitudes , g->tete[edgeIdx] ,  edgeIdx  , &c1,&p1);
    findTransition(*CompTree, g,STaltitudes , g->queue[edgeIdx] ,  edgeIdx  , &c2,&p2);
    
    createNode(CompTree, STaltitudes, g, edgeIdx , c1,c2 , &n);
    
    while(p1 != p2){
        if(STaltitudes[p2] < STaltitudes[p1]){
            // SWAP
            tmp = p1;p1= p2; p2 = tmp;
        }
        if(STaltitudes[p1] < STaltitudes[p2]){
            p1p = (*CompTree)->tabnodes[ p1 ].father ;
            attach(CompTree ,STaltitudes, p1 , n );
            n = p1;
            p1 = p1p;                                    
        }
        else{
            p1p = (*CompTree)->tabnodes[ p1 ].father;  p2p = (*CompTree)->tabnodes[ p2 ].father;
            //Merging
            merge(CompTree , STaltitudes, p1,p2, &np);
            attach(CompTree , STaltitudes , np , n );
            n = np;
            p1 = p1p;
            p2 = p2p;
        }
    }
    
    
    return 0;
}



////////----------------------------- With Attributes ---------------


int initializeTreeAttributes(/* INPUTS */
                   graphe *g,
                   /* OUTPUTS */
                   JCctree **CompTree,/* CompTree :  initalized component tree - QFZ */
                   double  *STaltitude,
                   double * STAltitudesSegmentation,
                   int *STArea,
                   double *STMaxWeight
                   ){
    
    int i ,nbsoncellsloc;
    JCsoncell * newsoncell1, *lastchild;
    JCctree *CT = *CompTree; //component tree
    int32_t taille = g->nbsom;
    
    CT->nbnodes = taille;
    CT->nbsoncells = 0;
    
    for(i = 0; i < taille; i++){
        CT->tabnodes[i].nbsons = 0;
        CT->tabnodes[i].sonlist = NULL;
        CT->tabnodes[i].lastson = NULL;
        CT->tabnodes[i].nodeAsSon = NULL;
    }
    
    for(i = 0; i < taille; i++){
        STaltitude[i] = 0; /* altitudes of the leaves is 0 */
        STArea[i] = 1; /* altitudes of the leaves is 0 */
        STMaxWeight = 0;
    }

    int z;
    z = CT->nbnodes; CT->nbnodes++;
    STaltitude[z] = 9999999999;
    CT->tabnodes[z].nbsons = taille;
    CT->tabnodes[z].nodeAsSon = NULL;
    
    //// Create root node z of the tree and attach all vertices as their initial children.
    CT->root = z;
    nbsoncellsloc = CT->nbsoncells;
    newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
    newsoncell1->son = 0;
    newsoncell1->next = NULL;
    newsoncell1->last = NULL;
    
    CT->tabnodes[z].sonlist = newsoncell1;
    CT->tabnodes[z].lastson = newsoncell1;
    CT->tabnodes[z].father = -1;
    CT->tabnodes[0].father = z;
    CT->tabnodes[0].nodeAsSon = newsoncell1;
    CT->nbsoncells += 1;
    
    STArea[z]=STArea[z]+ STArea[0];
    
    lastchild =newsoncell1;
    for(i = 1; i < taille; i++){
        nbsoncellsloc = CT->nbsoncells;
        newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
        newsoncell1->son = i;
        newsoncell1->next = NULL;
        STArea[z]=STArea[z]+ STArea[i];
        CT->tabnodes[i].father = CT->root;
        CT->nbsoncells += 1;
        CT->tabnodes[i].nodeAsSon = newsoncell1;
        newsoncell1->last = lastchild;
        lastchild->next =newsoncell1;
        lastchild = newsoncell1;
    }
    CT->tabnodes[z].lastson = lastchild;
    
    
    for(i =g->nbsom ; i< 2*g->nbsom; i++ ){
        STAltitudesSegmentation[i] = STaltitude[i] ;
    }
    
    
    return 0;
}


void attachAttributes(/* INPUTS */
            JCctree **CompTree, /* Component tree - QFZ */
            double  *STaltitudes,
            double *STAltitudesSegmentation,
            int *STArea,
            double *STMaxWeight,
            const int idxP,
            const int idxS, //  idxS is attached to idxP
            int *area
){
    JCctree * CT = *CompTree;
    //int  nbsoncellsloc;
    JCsoncell * newsoncell1 ;
    
    if(idxS == idxP){return;}//no attaching to its self}

    //////// Attributes  //////
    
    int mp = STArea[idxP] ;
    STArea[idxP] = STArea[idxP] - *area + STArea[idxS];
    *area = mp;
    
    if(STMaxWeight[idxS]>STMaxWeight[idxP])
        STMaxWeight[idxP] = STMaxWeight[idxS];
    ////////////////////////////
    
    if(CT->tabnodes[idxS].father == idxP){
        return;//no need to attach as it is already son of idxP
    }

    
    
    CT->tabnodes[idxP].nbsons += 1;

    
    detach(CompTree, STaltitudes , CT->tabnodes[idxS].father,  idxS);
    
    newsoncell1 = CT->tabnodes[idxS].nodeAsSon;
    newsoncell1->next = NULL;
    newsoncell1->last = CT->tabnodes[idxP].lastson;
    
    if(CT->tabnodes[idxP].lastson != NULL){
        CT->tabnodes[idxP].lastson->next = newsoncell1 ;
        CT->tabnodes[idxP].lastson = newsoncell1;
    }
    else{
        CT->tabnodes[idxP].lastson = newsoncell1;
        CT->tabnodes[idxP].sonlist = newsoncell1;
    }
    CT->tabnodes[idxS].father = idxP;
    
}

void mergeAttributes(/* INPUTS */
           JCctree **CompTree, /* Component tree - QFZ */
           double  *STaltitudes,
           int *STArea,
           double *STMaxWeight,
           int p1,
           int p2,
           int *np,
           int *area1,
           int *area2
           ){
    
    if(p1==p2){*np = p1; return;}
    
    JCctree * CT = *CompTree;
    int larger, smaller, i ;
    if(CT->tabnodes[p1].nbsons > CT->tabnodes[p2].nbsons){
        larger = p1;
        smaller = p2;
        
        //////// Attributes  //////
        int mp = STArea[larger] ;
        STArea[larger] = STArea[larger] + STArea[smaller] - *area1 - *area2;
        *area1 = mp;
        mp = STArea[smaller] ;
        STArea[smaller] = 0;
        *area2 = mp;
        
        if(STMaxWeight[smaller] > STMaxWeight[larger])
            STMaxWeight[larger] = STMaxWeight[smaller];
        ////////////////////////////
        
    }
    else{
        larger = p2;
        smaller = p1;
        
        //////// Attributes  //////
        int mp = STArea[larger] ;
        STArea[larger] = STArea[larger] + STArea[smaller] - *area1 - *area2;
        *area2 = mp;
        mp = STArea[smaller] ;
        STArea[smaller] = 0;
        *area1 = mp;
        
        if(STMaxWeight[smaller] > STMaxWeight[larger])
            STMaxWeight[larger] = STMaxWeight[smaller];
        
        ////////////////////////////
        
    }
    // Children of the smaller one will go into the larger one.
    JCsoncell * child = CT->tabnodes[smaller].sonlist;
    
    if(CT->tabnodes[larger].lastson){
        CT->tabnodes[larger].lastson->next = child ;
        child->last = CT->tabnodes[larger].lastson;
    }
    else{
        // in case larger does not have a last son
        CT->tabnodes[larger].sonlist = child ;
        printf("larger %i has  %i sons, merging with %i with %i sons \n", larger, CT->tabnodes[larger].nbsons,smaller , CT->tabnodes[smaller].nbsons );
    }
    
    CT->tabnodes[larger].lastson = CT->tabnodes[smaller].lastson;
    
    for(i=0; i< CT->tabnodes[smaller].nbsons; i++ ){
        CT->tabnodes[child->son].father = larger;
        CT->tabnodes[larger].nbsons += 1;
        child = child->next;
    }
    
    detach(CompTree, STaltitudes , CT->tabnodes[smaller].father,  smaller);
    
    

    
    CT->tabnodes[smaller].sonlist = NULL;
    CT->tabnodes[smaller].lastson = NULL;
    CT->tabnodes[smaller].father = smaller; //marked for delete
    CT->tabnodes[smaller].nbsons = 0;
    
    *np = larger;
    
}

void createNodeAttributes( /* INPUTS */
                JCctree **CompTree, /* Component tree - QFZ */
                double  *STaltitudes,
                double *STAltitudesSegmentation,
                int *STArea,
                double *STMaxWeight,
                graphe *g,
                int edgeIdx, /* Edge index for {x,y}  */
                int c1,
                int c2,
                /* OUTPUTS */
                int *n,
                int *area1, int *area2){
    JCctree * CT = *CompTree;
    int z, nbsoncellsloc;
    double lambda;
    lambda = g->weight[edgeIdx]+1;
   //  lambda = g->weight[edgeIdx];
    
    if(STaltitudes[c1] ==  STaltitudes[c2] &&STaltitudes[c1] == lambda  && CT->tabnodes[c1].nbsons > 0 && CT->tabnodes[c2].nbsons >0 ){
        mergeAttributes(CompTree , STaltitudes, STArea, STMaxWeight, c1,c2, n, area1, area2);
        if( g->originalWeight[edgeIdx] > STMaxWeight[*n]){
            STMaxWeight[*n] = g->originalWeight[edgeIdx];
        }
        
        return ;
    }
    if(STaltitudes[c1] == lambda && STaltitudes[c2] < lambda ){
        attachAttributes(CompTree , STaltitudes, STAltitudesSegmentation, STArea, STMaxWeight, c1 , c2, area1);
        
        *area2 = STArea[c2];
        // Max weight
        if( g->originalWeight[edgeIdx] > STMaxWeight[c1]){
            STMaxWeight[c1] = g->originalWeight[edgeIdx];
        }
        *n = c1;
        return ;
    }
    if(STaltitudes[c2] == lambda && STaltitudes[c1] < lambda ){
        attachAttributes(CompTree , STaltitudes,STAltitudesSegmentation ,  STArea, STMaxWeight, c2 , c1, area2);
        *area1 = STArea[c1];
        // Max weight
        if( g->originalWeight[edgeIdx] > STMaxWeight[c2]){
            STMaxWeight[c2] = g->originalWeight[edgeIdx];
        }
        
        *n = c2;
        return ;
    }
    
    JCsoncell *newsoncell1;
    z = CT->nbnodes; CT->nbnodes++;
    
    nbsoncellsloc = CT->nbsoncells;
    
    CT->tabnodes[z].lastson = NULL;
    CT->tabnodes[z].sonlist = NULL;
    CT->tabnodes[z].nbsons = 0;
    
    STaltitudes[z] = lambda;
    STAltitudesSegmentation[z] = lambda-1; // For compatibilty with normal QFZ 
    
    attachAttributes(CompTree, STaltitudes, STAltitudesSegmentation, STArea, STMaxWeight, z, c1, area1);
    attachAttributes(CompTree, STaltitudes, STAltitudesSegmentation, STArea,STMaxWeight, z, c2, area2);
    /////////////// Attributes /////////////
    //Area
    *area1 = STArea[c1];
    *area2 = STArea[c2];
    //// get the maximum weight
    STMaxWeight[z] = g->originalWeight[edgeIdx];
    if(STMaxWeight[c1] > STMaxWeight[z]){
        STMaxWeight[z] = STMaxWeight[c1];
    }
    if(STMaxWeight[c2] > STMaxWeight[z]){
        STMaxWeight[z] = STMaxWeight[c2];
    }
    ///////////////////////////////////////
    
    
    newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
    newsoncell1->son = z;
    newsoncell1->next =NULL;
    newsoncell1->last =NULL;
    
    CT->tabnodes[z].father = CT->root;
    CT->tabnodes[CT->tabnodes[z].father].nbsons += 1;
    CT->tabnodes[z].nodeAsSon = newsoncell1;
    
    if(CT->tabnodes[CT->root].lastson != NULL){
        CT->tabnodes[CT->tabnodes[z].father].lastson->next = newsoncell1 ;
        newsoncell1->last = CT->tabnodes[CT->tabnodes[z].father].lastson ;
        CT->tabnodes[CT->tabnodes[z].father].lastson = newsoncell1;
    }
    else{ // root was left with no children
        CT->tabnodes[CT->root].lastson = newsoncell1 ;
        CT->tabnodes[CT->root].sonlist = newsoncell1 ;
    }
    *n = z;
    CT->nbsoncells+=1;
    
    
    
}






int treeUpdateAttributes( /* INPUTS */
               JCctree **CompTree, /* Component tree - QFZ */
               double  *STaltitudes,
               double *STAltitudesSegmentation,
               int* STArea,
               double *STMaxWeight,
               graphe *g,
               int edgeIdx /* Edge index for {x,y}  */
/* OUTPUTS */
//       CompTree
){
    
    int c1,p1,c2,p2, n, tmp;
    int p1p, p2p, np;
    int area1, area2;
    area1 = area2 = 0;
    //printf("{%i,%i} and w=%f \n" , g->tete[edgeIdx], g->queue[edgeIdx], g->weight[edgeIdx]);
    findTransition(*CompTree, g,STaltitudes , g->tete[edgeIdx] ,  edgeIdx  , &c1,&p1);
    findTransition(*CompTree, g,STaltitudes , g->queue[edgeIdx] ,  edgeIdx  , &c2,&p2);
    
    
    //printf("c1 %i c2 %i \n" ,  c1, c2);
    
    createNodeAttributes(CompTree, STaltitudes, STAltitudesSegmentation, STArea, STMaxWeight, g, edgeIdx , c1,c2 , &n, &area1, &area2);
    
    while(p1 != p2){
        if(STaltitudes[p2] < STaltitudes[p1]){
            // SWAP
            tmp = p1;p1= p2; p2 = tmp;
            tmp = area2;
            area2 = area1;
            area1 = tmp;
        }
        if(STaltitudes[p1] < STaltitudes[p2]){
            p1p = (*CompTree)->tabnodes[ p1 ].father ;
            attachAttributes(CompTree ,STaltitudes,STAltitudesSegmentation, STArea, STMaxWeight, p1 , n , &area1);
            n = p1;
            p1 = p1p;
        }
        else{
            p1p = (*CompTree)->tabnodes[ p1 ].father;  p2p = (*CompTree)->tabnodes[ p2 ].father;
            //Merging
            mergeAttributes(CompTree , STaltitudes, STArea, STMaxWeight, p1,p2, &np, &area1, &area2);
            int dummy = 0;
            attachAttributes(CompTree , STaltitudes,STAltitudesSegmentation, STArea, STMaxWeight, np , n , &dummy);
            n = np;
            p1 = p1p;
            p2 = p2p;
        }
    }
    
   /*
    for(int i =g->nbsom ; i< 2*g->nbsom; i++ ){
        if(STaltitudes[i]!= INFINITO) STAltitudesSegmentation[i] = STaltitudes[i] - 1;        
    }*/
    
    return 0;
}





//------------------------------------------------------------------------------------------------------------------



int IncrementalQFZ( /* INPUTS */
                   graphe *g,
                   char* file
                   ){
    
    double  *STaltitudes;
    int narcs, i;
    JCctree *CompTree; /* Component tree - QFZ */
    
    
    int32_t taille = g->nbsom;
    
    if( (CompTree = componentTreeAlloc((3*taille))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STaltitudes = (double *)calloc(3*g->nbsom, sizeof(double));
    if (STaltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    
    
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    
    
    initializeTree(g,&CompTree,STaltitudes);
    
    narcs = g->ind;
    
    for(i = 0; i < narcs; i++){
        treeUpdate(&CompTree, STaltitudes , g, i);
    }
    
    end = clock();
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("total update took %f seconds to execute \n", cpu_time_used);
   printf("File\t %s \t Number of edges \t %i \t time \t %f \t \n",file, narcs, cpu_time_used );

    /*
    printf("nro nodos %i, hijos %i \n",(CompTree)->nbnodes, (CompTree)->nbsoncells);
    
     taille = (CompTree)->nbnodes;
     for(i = 0; i < taille; i++){
     printf("node i %i has parent %i node at altitude %f \n", i ,  (CompTree)->tabnodes[i].father , STaltitudes[i]);
     }
    */
    
   // printf("Number of edges %i\n",narcs);
    
    
    //componentTreePrint(CompTree, STaltitudes, CompTree->root, NULL);
    
    
    
    return 0;
}




int random_int(int min, int max){
    return min + rand() % (max+1 - min);
}

int IncrementalQFZfromMST( /* INPUTS */
                          graphe *g,
                          char *file
                          ){
    double  *STaltitudes;
    int narcs, i;
    JCctree *CompTree; /* Component tree - QFZ */
    
    int32_t taille = g->nbsom;
    
    if( (CompTree = componentTreeAlloc((3*taille))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STaltitudes = (double *)calloc(3*g->nbsom, sizeof(double));
    if (STaltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    
    MST *mst;
    graphe *newG;
    getMST(g, &mst);
    graphFromMST2(&newG, mst,g->nbsom ,mst->size);
    
//////////////// Shuffle values ////////////////////////
    int* randomkeys = (int *)calloc(mst->size, sizeof(int));
    int* created = (int *)calloc(mst->size, sizeof(int));
    int num;
    for(i=0;i<mst->size;i++){
        
        while(1){
            num = random_int(0,mst->size-1);
                if(created[num]==0){
                    randomkeys[i] = num;
                    created[num] = 1;
                    break;
                }
            }
    }
    
    clock_t start, end;
    double cpu_time_used;
    start = clock();
////////////////////////////////////
    
    narcs = newG->ind;
    
    initializeTree(newG,&CompTree,STaltitudes);

    for(i = 0; i < narcs; i++){
        treeUpdate(&CompTree, STaltitudes , newG, randomkeys[i]);
        //componentTreePrint(CompTree, STaltitudes, CompTree->root, NULL);
    }
    
    end = clock();
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("total update took %f seconds to execute \n", cpu_time_used);
    //checkDegenerateNodes(CompTree, STaltitudes, CompTree->root );
    //componentTreePrint(CompTree, STaltitudes, CompTree->root , NULL);
    
    if(CompTree->tabnodes[CompTree->root].nbsons > 1){
        printf("not connected graph \n");
    }
    
    CompTree->root = CompTree->tabnodes[CompTree->root].sonlist->son;
    CompTree->tabnodes[CompTree->root].father = -1;
    
    taille = CompTree->nbnodes;
    for(i = newG->nbsom; i < taille; i++)
        STaltitudes[i] -= 1;
    
    double sm ;
    sm = computeSM(CompTree, g, STaltitudes);
    printf("sum out: %f\n", sm);
    
    /*
    taille = (CompTree)->nbnodes;
     for(i = 0; i < taille; i++){
     printf("node i %i has parent %i node at altitude %f \n", i ,  (CompTree)->tabnodes[i].father , STaltitudes[i]);
     }
    */
    free(randomkeys);
    free(created);
    
    return 0;
}











