//
//  graphSegmentation.c
//  SM_Edward
//
//  Created by Edward Cayllahua on 8/22/16.
//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//




#include<graphSegmentation.h>

#include <time.h>

clock_t diss_ini, diss_fin;
double diss_sum=0;


double min(double v1, double v2){
    if(v1<=v2)
        return v1;
    else return v2;
}
double max(double v1, double v2){
    if(v1>=v2)
        return v1;
    else return v2;
}


void naiveSegmentation(graphe *g){
    
    
    clock_t ini_all, fin_all, ini_QFZ, fin_QFZ, ini_mini, fin_mini, ini_set, fin_set;
    double sum_QFZ=0, sum_mini = 0 , sum_all=0, sum_loop=0, sum_set=0;
    
    int i;
    double  newf;
    MST *mst;
    graphe *newG;
    //////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    double *STAltitudes = NULL;
    
    // The MST is computed and stored in mst
    getMST(g, &mst);
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO );
    /// Set all edges to large for later SM computation ////
    computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    ini_all = clock();
    
    for(i=0;i<mst->size;i++){
        
        
        ini_mini = clock();
        newf = naiveMinimization(newG, mst,  &CompTree ,STAltitudes , STArea, STMaxWeight, i ) -1 ;
        fin_mini = clock();
        sum_mini += ((double) (fin_mini - ini_mini)) / CLOCKS_PER_SEC;
        printf("New value for edge: %i,%i = %f , edge %i of %i \n",mst->MSTedges[i].x, mst->MSTedges[i].y,newf , i ,mst->size  );
        ini_set = clock();
        setweight(newG,i,newf);
        fin_set = clock();
        sum_set += ((double) (fin_set - ini_set)) / CLOCKS_PER_SEC;
        
        ini_QFZ = clock();
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        fin_QFZ = clock();
        sum_QFZ += ((double) (fin_QFZ - ini_QFZ)) / CLOCKS_PER_SEC;

    }
    
    fin_all = clock();
    sum_all += ((double) (fin_all - ini_all)) / CLOCKS_PER_SEC;
    sum_loop = sum_all - sum_QFZ - sum_mini - sum_set;
    
    //printf("1\t all\t %f \t loop: \t %f \t QFZ: \t %f \t MINI: \t %f \t DIS \t %f \t SETS\t %f \n ", sum_all , sum_loop, sum_QFZ , sum_mini, diss_sum, sum_set);
    
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);

    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    componentTreeFree(CompTree);
}


void segmentationRange(graphe *g){

    clock_t ini_all, fin_all, ini_QFZ, fin_QFZ, ini_mini, fin_mini, ini_set, fin_set;
    double sum_QFZ=0, sum_mini = 0 , sum_all=0, sum_loop=0, sum_set=0;
    
    int i;
    double  newf;
    MST *mst;
    graphe *newG;
    //////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    double *STAltitudes = NULL;
    
    // The MST is computed and stored in mst
    getMST(g, &mst);
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO);
    /// Set all edges to large value for later SM computation ////
    NodeList* frange;
    ListInit(&frange);
    computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    ini_all = clock();
    
    for(i=0;i<mst->size;i++){
        
        
        ini_mini = clock();
        newf =  rangeMinimization(newG, mst,  &CompTree ,STAltitudes , STArea, STMaxWeight, i, &frange ) - 1 ;
        fin_mini = clock();
        sum_mini += ((double) (fin_mini - ini_mini)) / CLOCKS_PER_SEC;
        
        //printf("New value for edge: %i,%i = %f , edge %i of %i \n",mst->MSTedges[i].x, mst->MSTedges[i].y,newf , i ,mst->size  );
        ini_set = clock();
        setweight(newG,i,newf);
        fin_set = clock();
        sum_set += ((double) (fin_set - ini_set)) / CLOCKS_PER_SEC;
        
        ini_QFZ = clock();
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        fin_QFZ = clock();
        sum_QFZ += ((double) (fin_QFZ - ini_QFZ)) / CLOCKS_PER_SEC;
        

        
    }
    
    
    fin_all = clock();
    sum_all += ((double) (fin_all - ini_all)) / CLOCKS_PER_SEC;
    sum_loop = sum_all - sum_QFZ - sum_mini - sum_set;
    
    printf("2\t all\t %f \t loop: \t %f \t QFZ: \t %f \t MINI: \t %f \t DIS \t %f \t SETS\t %f \n ", sum_all , sum_loop, sum_QFZ , sum_mini, diss_sum, sum_set);
    
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
    NodeList* tmp = frange, *tmp2 = frange;
    while(tmp!=NULL){
        tmp2= tmp->next;
        free(tmp);
        tmp = tmp2;
    }
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    componentTreeFree(CompTree);
}

void segmentationBranch(graphe *g){
    
    clock_t ini_all, fin_all, ini_QFZ, fin_QFZ, ini_mini, fin_mini, ini_set, fin_set;
    double sum_QFZ=0, sum_mini = 0 , sum_all=0, sum_loop=0, sum_set=0;
    
    int i;
    double  newf;
    MST *mst;
    graphe *newG;
    //////
    JCctree* CompTree = NULL; // Component tree - QFZ
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    double *STAltitudes = NULL;
    // The MST is computed and stored in mst
    
    getMST(g, &mst);
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO);
    /// Set all edges to large value for later SM computation ////
    
    computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    
    ini_all = clock();
    
    for(i=0;i<mst->size;i++){
        
       
        
        ini_mini = clock();
        newf = branchMinimization(newG, mst,  &CompTree ,STAltitudes , STArea, STMaxWeight, i ) - 1;
        fin_mini = clock();
        sum_mini += ((double) (fin_mini - ini_mini)) / CLOCKS_PER_SEC;
    
        //  printf("New value for edge: %i,%i = %f , edge %i of %i \n",mst->MSTedges[i].x, mst->MSTedges[i].y,newf , i ,mst->size  );
        ini_set = clock();
        setweight(newG,i,newf);
        fin_set = clock();
        sum_set += ((double) (fin_set - ini_set)) / CLOCKS_PER_SEC;
        
        
        
        ini_QFZ = clock();
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        fin_QFZ = clock();
        sum_QFZ += ((double) (fin_QFZ - ini_QFZ)) / CLOCKS_PER_SEC;
    
    }

    
    fin_all = clock();
    sum_all += ((double) (fin_all - ini_all)) / CLOCKS_PER_SEC;
    sum_loop = sum_all - sum_QFZ - sum_mini - sum_set;
    printf("3\t all\t %f \t loop: \t %f \t QFZ: \t %f \t MINI: \t %f \t DIS \t %f \t SETS\t %f \n ", sum_all , sum_loop, sum_QFZ , sum_mini, diss_sum, sum_set);
    
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
  //  char * filename=  "/tmp/normalGraph.dot";
  //  saveDOTfile(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight,  filename, g->nbsom);
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    componentTreeFree(CompTree);
    
    
}






double branchMinimization(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx )
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    double lambda_star = -9999;
    
    cx = x;
    cy = y;
    
    lambda = STAltitudes[cx];
    double prevLambda=-1;
   
    /*
    FILE *log = NULL;
    log = fopen("/tmp/jumps.txt", "a");
    fprintf(log,"for edge %i - %i Diss: " , x, y );
    */
   
    // componentTreePrintAttributes(QCT, STAltitudes, QCT->root, STArea, STMaxWeight );
    
    //D = computeD(CT, STAltitudes, STArea, STMaxWeight, &cx, &cy , Dif , lambda );
    D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
    
   //  fprintf(log,"at lambda: %f D:  %f \n" , lambda, D );
    
   
   // int steps = 1;
    
    while( D >lambda ){
      //  steps+=1;
        
        px = CT->tabnodes[cx].father;
        py = CT->tabnodes[cy].father;
        
        prevLambda = lambda;
        
        lambda = min(STAltitudes[px ] , STAltitudes[py] );
        
        diss_ini = clock();
            //D = computeD(CT, STAltitudes, STArea, STMaxWeight, &cx, &cy , Dif , lambda );
        
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
        diss_fin = clock();
        
        diss_sum+=((double) (diss_fin - diss_ini)) / CLOCKS_PER_SEC;
        //   fprintf(log,"at lambda: %f D:  %f \n" , lambda, D );
        
        if(STAltitudes[px ] == lambda )
            cx =px;
        if(STAltitudes[py ] == lambda )
            cy =py;
    }
   
    if(D <= prevLambda) lambda_star = prevLambda + 1;
    else lambda_star = D;
    
    
   //  char * filename=  "/tmp/BranchGraph.dot";
   //  saveDOTfile(CT, STAltitudes, CT->root, STArea, STMaxWeight,  filename);
  
    //fprintf(log,"%i \n" , steps );
   // fclose(log);
    
    return lambda_star;
}



double branchMinimizationLambdaStars(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    NodeList* lambdaStars;
    ListInitNormal(&lambdaStars);
    
    double lambda_star = -9999;
    
    cx = x;
    cy = y;
    
    lambda = STAltitudes[cx];
    double prevLambda=-1;
    
    
    D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );

    while( D >lambda ){
        px = CT->tabnodes[cx].father;
        py = CT->tabnodes[cy].father;
        prevLambda = lambda;
        lambda = min(STAltitudes[px ] , STAltitudes[py] );
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
        
        if(STAltitudes[px ] == lambda )
            cx =px;
        if(STAltitudes[py ] == lambda )
            cy =py;
    }
    if(D <= prevLambda) lambda_star = prevLambda + 1;
    else lambda_star = D;
    
    NodeList* LastNode = lambdaStars;
    LastNode->value = lambda_star ; // The first Lambda Star

    double prevD = D;
    
    while(px!=py && D<=lambda  ){
        prevLambda = lambda;
        lambda = min(STAltitudes[px ] , STAltitudes[py] );
        prevD = D;
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
        
        px = CT->tabnodes[cx].father;
        py = CT->tabnodes[cy].father;
        
        if(STAltitudes[px ] == lambda )
            cx =px;
        if(STAltitudes[py ] == lambda )
            cy =py;
        
        if(D < prevLambda )
            lambda_star = prevLambda + 1;
        else lambda_star = prevD;
        
        addValueNormal(&lambdaStars, &LastNode, lambda_star);
    }
    if(D < prevLambda )
        lambda_star = prevLambda + 1;
    else lambda_star = prevD;

    addValueNormal(&lambdaStars, &LastNode, lambda_star);
    
    
    switch(op){
        case 1: //Min: first value on the list
            lambda_star  = lambdaStars->value;
            break;
        case 2: // Max, last value on the list
            lambda_star  = LastNode->value;
            break;
        case 3: // Mean
            lambda_star = listMeanValue(lambdaStars);
            break;
        case 4: // Median
            lambda_star = listMedianValue(lambdaStars);
            break;
    }
    
    
    
    while(lambdaStars!= NULL){
        LastNode = lambdaStars;
        lambdaStars = lambdaStars->next;
        free(LastNode);
    }
    
    return lambda_star;
}

double branchMinimizationLambdaStarsAll(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    NodeList* lambdaStars;
    lambdaStars = NULL;
    NodeList* LastNode = lambdaStars;
    double lambda_star = -9999;
    
    NodeList* lambdaStarsLastInterval = NULL;
    NodeList* LastNodeLI = lambdaStarsLastInterval;
    
    cx = x;
    cy = y;
    
    lambda = STAltitudes[cx];
    double prevLambda=-1;
    
    
    do{
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        
        while( D >lambda ){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            
            if(STAltitudes[px ] == lambda )
                cx =px;
            if(STAltitudes[py ] == lambda )
                cy =py;
        }
        if(D <= prevLambda) lambda_star = prevLambda + 1;
        else lambda_star = D;
        
        
        addValueNormal(&lambdaStars, &LastNode, lambda_star);
        cleanList(&lambdaStarsLastInterval);
        addValueNormal(&lambdaStarsLastInterval, &LastNodeLI, lambda_star);
        
        double prevD = D;
        
        while(px!=py && D<=lambda  ){
            prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            prevD = D;
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            
            if(STAltitudes[px ] == lambda )
                cx =px;
            if(STAltitudes[py ] == lambda )
                cy =py;
            
            if(D < prevLambda )
                lambda_star = prevLambda + 1;
            else lambda_star = prevD;
            
            addValueNormal(&lambdaStars, &LastNode, lambda_star);
            addValueNormal(&lambdaStarsLastInterval, &LastNodeLI, lambda_star);
        }
        //if(D < prevLambda )
        //    lambda_star = prevLambda + 1;
        //else lambda_star = prevD;
        
       // addValueNormal(&lambdaStars, &LastNode, lambda_star);
    }
    while(px!=py);
    
    switch(op){
        case 1: //Min: first value on the list
            lambda_star  = lambdaStars->value;
            break;
        case 2: // Max, last value on the list
            lambda_star  = LastNode->value;
            break;
        case 3: // Mean
            lambda_star = listMeanValue(lambdaStars);
            break;
        case 4: // Median
            lambda_star = listMedianValue(lambdaStars);
            break;
       //  values on last interval
        case 5:// Min last interval
            lambda_star  = lambdaStarsLastInterval->value;
            break;
        case 6: // Max last interval
            lambda_star  = LastNodeLI->value;
            break;
        case 7: // Mean last interval
            lambda_star = listMeanValue(lambdaStarsLastInterval);
            break;
        case 8:// Median last interval
            lambda_star = listMedianValue(lambdaStarsLastInterval);
            break;
    }
    
    
    
    while(lambdaStars!= NULL){
        LastNode = lambdaStars;
        lambdaStars = lambdaStars->next;
        free(LastNode);
    }
    
    return lambda_star;
}

int  getNegativeIntervals(NodeList* positiveInterval, NodeList** negativeInterval, NodeList** lastNinterval, double *total_length){
    
    int counter = 0;
    //double lbound = -999999999999, ubound;
    double lbound = 0, ubound;
    double tlength = 0;
    
    NodeList* tmp = positiveInterval;
    NodeList* lastNode = NULL;
    while(tmp != NULL){
        ubound = tmp->value-0.001;
        addValueNormal( negativeInterval, &lastNode, lbound );
        addValueNormal( negativeInterval, &lastNode, ubound );
        
        tlength+=fabs(ubound-lbound);
        
        if((tmp->next)->next!= NULL)
            lbound = (tmp->next)->value + 0.001;
        tmp=(tmp->next)->next;
        counter++;
    }
    *lastNinterval = lastNode;
    *total_length = tlength;
    
    
    
    return counter;
}
void printInterval(int x,int y, NodeList* list){
    NodeList* tmp = list;
    
    
    
    while(tmp!=NULL){
      //  if(  (tmp->next->value - tmp->value) < 9000000000.000)
      //  printf("%.3f\t%.3f\t%.3f\n", tmp->value, tmp->next->value , (tmp->next->value - tmp->value));
        printf("\t%.3f\t%.3f \n", tmp->value, tmp->next->value  );
        
        tmp = tmp->next->next;
    }
}


void printIntervalDis(int x,int y, NodeList* list){
    NodeList* tmp = list;
    printf("Edge %i %i\n", x, y);
    while(tmp!=NULL){
        //  if(  (tmp->next->value - tmp->value) < 9000000000.000)
        //  printf("%.3f\t%.3f\t%.3f\n", tmp->value, tmp->next->value , (tmp->next->value - tmp->value));
       if(tmp->next->bound ==1)
           printf("\t %.8f\t%.8f D=%.8f >> BOUND \n", tmp->value, tmp->next->value , tmp->dis );
       else
           printf("\t %.8f\t%.8f D=%.8f \n", tmp->value, tmp->next->value , tmp->dis );
        
        tmp = tmp->next->next;
    }
}




void filterInterval( NodeList** list, int threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list;
    ite = tmp;
    
    while(ite!=NULL){
        ite = tmp->next->next;
        if( fabs(tmp->next->value - tmp->value) < threshold){
            if( tmp->last!=NULL ){
                tmp->last->next = tmp->next->next;
                if(ite!=NULL){
                    ite->last = tmp->last;
                }
                
                free(tmp->next);
                free(tmp); 
            }
            else{
                ite->last = NULL;
                *list = ite;
                free(tmp->next);
                free(tmp);
            }
        }
        tmp = ite;
    }
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
}

void filterIntervalLength( NodeList** list, double threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
//    if( (*lastNode)->value - (*list)->value <= threshold){
//        printf("size of interval less than threshold \n");
//        return;
//    }
    
    double sumLenght = 0 ;
    
    while(ite!=NULL){
        
        firstp = tmp;
        sumLenght = 0;
        while(1){
            
            sumLenght+= fabs(tmp->next->value - tmp->value);
            if(tmp->next->bound!=1)
                tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
            }
        }
        
        ite = lastp->next;
        
        
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( sumLenght < threshold){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
}





void filterIntervalArea( NodeList** list, double threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
    double a, b, c;
    double sumVol = 0 ;
    
    while(ite!=NULL){
        
        firstp = tmp;
        sumVol = 0;
        while(1){
            c = tmp->next->value - tmp->value;
            a = fabs(tmp->dis - tmp->value);
            b = fabs(tmp->dis - tmp->next->value);
            sumVol+= c*(a+b)/2;
            if(tmp->next->bound!=1)
               tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
                }
        }

        ite = lastp->next;
       // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( sumVol < threshold){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
   // printf("scale %.8f\n", (*lastNode)->value  );
}


void filterIntervalVolume( NodeList** list, double threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
    double a, b, c;
    double sumVol = 0 ;
    
    while(ite!=NULL){
        
        firstp = tmp;
        sumVol = 0;
        while(1){
            c = tmp->next->value - tmp->value;
            a = fabs(tmp->dis - tmp->value);
            b = fabs(tmp->dis - tmp->next->value);
            if( fabs(tmp->dis) > 0 )
                sumVol+= (( c*(a+b)/2 )*fabs(tmp->dis));
            else
                sumVol+= (c*(a+b)/2 );
            
            if(tmp->next->bound!=1)
                tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
            }
        }
        
        ite = lastp->next;
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( sumVol < threshold){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
    // printf("scale %.8f\n", (*lastNode)->value  );
}



void filterIntervalVolumeAdapt( NodeList** list, NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
    double a, b, c, count=0;
    double sumVol = 0 , avgVol = 0;
    
    
  //  printf("para filtering \n");
   NodeList* fp = *list;
    
    while(ite!=NULL){
        
        sumVol = 0;
        while(1){
            //printf("   valores  %.3f address %i \n", tmp->value ,  tmp);
            c = tmp->next->value - tmp->value;
            a = fabs(tmp->dis - tmp->value);
            b = fabs(tmp->dis - tmp->next->value);
            
            sumVol+= (( c*(a+b)/2 ));
            
            
            if(tmp->next->bound!=1)
                tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
            }
        }
        ite = lastp->next;
        tmp = ite ;
        
        if( fp != *list ){
            count++;
            avgVol+=sumVol;
        }
        else
            fp = ite;
        
    }
    
    //*lastNode = (*list)->next;
    
    
    
    
    //return;
    
    if(count>0){
        avgVol=(avgVol/count);
      //  printf("average %.5f\n",avgVol);
    }
        //avgVol=(avgVol/count)/700;
    else
        avgVol = 9999999999;
    
    tmp = *list;
    ite = tmp;
    lastp = *list;
    
    while(ite!=NULL){
        
        firstp = tmp;
        sumVol = 0;
        while(1){
            c = tmp->next->value - tmp->value;
            a = fabs(tmp->dis - tmp->value);
            b = fabs(tmp->dis - tmp->next->value);
            
            sumVol+= ( c*(a+b)/2 );
            
            
            if(tmp->next->bound!=1)
                tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
            }
        }
        
        ite = lastp->next;
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( sumVol < avgVol){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
    // printf("scale %.8f\n", (*lastNode)->value  );
}


void filterIntervalDepths( NodeList** list, double threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
    double a, b, c;
    double maxDepth = -999999999999;
    
    while(ite!=NULL){
        
        firstp = tmp;
        maxDepth = -999999999999;
        while(1){
            c = tmp->next->value - tmp->value;
            a = fabs(tmp->dis - tmp->value);
            b = fabs(tmp->dis - tmp->next->value);
           
            maxDepth = max(maxDepth,max(a,b));
            
            if(tmp->next->bound!=1)
                tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
            }
        }
        
        ite = lastp->next;
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( maxDepth < threshold){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
    // printf("scale %.8f\n", (*lastNode)->value  );
}




void computeDepths( NodeList** list,  NodeList** lastNode, double* min_d, double* max_d){
    
    NodeList* ite;
    //NodeList* tmp = *list, *tmp2;
    
    ite = *list;
    
    double a, b, c;
    double min, max, lb_max, lb_min, ub_max, ub_min;
    min=999999999;
    max=-100;
    lb_max = ub_max = lb_min = ub_min = 0;
    
    while(ite!=NULL){
        
        c = ite->next->value - ite->value;
        
        a = fabs(ite->dis - ite->value);
        b = fabs(ite->dis - ite->next->value);
        
        if( ((a+b)/2) >= max){
            lb_max = ite->value;
            ub_max = ite->next->value;
            max = (a+b)/2;
        }
        if( ((a+b)/2) <= min ){
            lb_min = ite->value;
            ub_min = ite->next->value;
            min = (a+b)/2;
        }
        ite = ite->next->next;
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
    }
    
 
      *min_d = lb_min +( (ub_min - lb_min)/2.0) ;
      *max_d = lb_max + ((ub_max - lb_max )/2.0);
    
    
}


void computeDepthsPosBound( NodeList** list,  NodeList** lastNode, double* min_d, double* max_d){
    
    NodeList* ite;
    //NodeList* tmp = *list, *tmp2;
    
    ite = *list;
    
    double a, b, c;
    double min, max, lb_max, lb_min, ub_max, ub_min;
    min=999999999;
    max=-100;
    lb_max = ub_max = lb_min = ub_min = 0;
    
    while(ite!=NULL){
        
        c = ite->next->value - ite->value;
        
        a = fabs(ite->dis - ite->value);
        b = fabs(ite->dis - ite->next->value);
        
        if( ((a+b)/2) >= max){
            lb_max = ite->value;
            ub_max = ite->next->value;
            max = (a+b)/2;
        }
        if( ((a+b)/2) <= min ){
            lb_min = ite->value;
            ub_min = ite->next->value;
            min = (a+b)/2;
        }
        ite = ite->next->next;
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
    }
    
    //    *min_d = (ub_min - lb_min)/2;
    //    *max_d = (ub_max - lb_max )/2;
    
    
    *min_d = lb_min ;
    *max_d = lb_max ;
    
    
}

void computeDepthsNegBound( NodeList** list,  NodeList** lastNode, double* min_d, double* max_d){
    
    NodeList* ite;
    //NodeList* tmp = *list, *tmp2;
    
    ite = *list;
    
    double a, b, c;
    double min, max, lb_max, lb_min, ub_max, ub_min;
    min=999999999;
    max=-100;
    lb_max = ub_max = lb_min = ub_min = 0;
    
    while(ite!=NULL){
        
        c = ite->next->value - ite->value;
        
        a = fabs(ite->dis - ite->value);
        b = fabs(ite->dis - ite->next->value);
        
        if( ((a+b)/2) >= max){
            lb_max = ite->value;
            ub_max = ite->next->value;
            max = (a+b)/2;
        }
        if( ((a+b)/2) <= min ){
            lb_min = ite->value;
            ub_min = ite->next->value;
            min = (a+b)/2;
        }
        ite = ite->next->next;
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
    }
    
    //    *min_d = (ub_min - lb_min)/2;
    //    *max_d = (ub_max - lb_max )/2;
    
    
    *min_d = ub_min ;
    *max_d = ub_max ;
    //  *min_d = lb_min +( (ub_min - lb_min)/2.0) ;
    // *max_d = lb_max + ((ub_max - lb_max )/2.0);
    
    
}




double getRankNegative(NodeList* list, NodeList* lastNode, double obs){
    
    double rank_val = -1111;
    
    NodeList* ite;
    NodeList* tmp = lastNode;
    ite = tmp;
    
    while(ite!=NULL){
        if( ( (ite->value - ite->last->value) - obs ) >  0 ){
            rank_val = ite->value - obs;
            break;
        }
        else{
            obs-=(ite->value - ite->last->value);
            ite = ite->last->last;
        }
    }
    return rank_val;
    
    
    
}

double getRankNegativeVolume(NodeList* list, NodeList* lastNode, double obs){
    
    double rank_val = -1111;
    
    NodeList* ite;
    NodeList* tmp;
    
    ite = list;
    ite->value = 0 ;
    double length = 0;
    
    while(ite!=NULL){
        length+= fabs( ite->next->value - ite->value );
        ite=ite->next->next;
    }
    
    //obs = floor( length*obs );
    obs = length*obs ;
    
    tmp = lastNode;
    ite = tmp;
    
    while(ite!=NULL){
        if( ( (ite->value - ite->last->value) - obs ) >  0 ){
            rank_val = ite->value - obs;
            break;
        }
        else{
            obs-=(ite->value - ite->last->value);
            ite = ite->last->last;
        }
    }
    return rank_val;
}


double getRankPositive(NodeList* list, double p){
    
    double rank_val = -1111;
    NodeList* ite = list;
    double total_l = 0, p_obs;
    
    while(ite!=NULL){
        total_l+=fabs(  ite->next->value - ite->value );
        
        ite = ite->next->next;
    }
    //p_obs = floor( total_l*p );
    p_obs = floor( total_l*p );
    
    ite = list;
    while(ite!=NULL){
        if( (  ite->next->value - ite->value ) - p_obs  >=  0 ){
            rank_val = ite->value+p_obs;
            break;
        }
        else{
            p_obs-=(  ite->next->value - ite->value );
            ite = ite->next->next;
        }
        
    }
    return rank_val;
}





double getRankPositiveBound(NodeList* list, double p){
    
    double rank_val = -1111;
    NodeList* ite = list;
    double total_l = 0, p_obs;
    
    while(ite!=NULL){
        total_l+=(  ite->next->value - ite->value );
        
        ite = ite->next->next;
    }
    //p_obs = floor( total_l*p );
    p_obs = ( total_l*p );
    
    ite = list;
    while(ite!=NULL){
        if( (  ite->next->value - ite->value ) - p_obs  >=  0 ){
            rank_val = ite->value;
            break;
        }
        else{
            p_obs-=(  ite->next->value - ite->value );
            ite = ite->next->next;
        }
        
    }
    return rank_val;
}



double getRankNegativeBound(NodeList* list, NodeList* lastNode, double obs){
    
    double rank_val = -1111;
    
    NodeList* ite;
    NodeList* tmp = lastNode;
    ite = tmp;
    
    while(ite!=NULL){
        if( ( (ite->value - ite->last->value) - obs ) >  0 ){
            
            if(ite->last->value!=0)
                //rank_val = ite->last->value ;
                rank_val = ite->value ;
            else
                rank_val = ite->value;
            break;
        }
        else{
            obs-=(ite->value - ite->last->value);
            ite = ite->last->last;
        }
    }
    return rank_val;
}



double minimizationLambdasInterval(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op, int* numberInterval, int *max_numberInterval)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    
    double lambdaUpper = -9999;
    double lambdaLower = -9999;
    double total_l;
    NodeList* positiveInterval = NULL;
    NodeList* negativeInterval = NULL;
    NodeList* LastNode = positiveInterval;
    NodeList* lastNinterval = NULL;
    int cn_intervals = -1;
    
    cx = x;
    cy = y;
    int T = -1;
    double p_obs=0;
    lambda = STAltitudes[cx];
    double prevLambda=-1, lambda_star =-9999, l_1,l_2;
    int cinterval=0;
    
      //  printf("Edge %i,%i \n",x,y );
    do{
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        
        while( D >lambda ){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
          //  printf( "%.2f,%i,%i,%.2f,%.2f,%.2f\n " ,  lambda ,STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy]  ,  D );
            if(STAltitudes[px ] == lambda )
                cx =px;
            if(STAltitudes[py ] == lambda )
                cy =py;
        }
        if(D <= prevLambda) lambdaLower = prevLambda + 0.0001;
        else
            lambdaLower = D;
        
        double prevD = D;
        
        while( D<=lambda  &&( px!=-1&&py!=-1)){
            prevLambda = lambda;
            
            if(px!=-1&&py!=-1)
                lambda = min(STAltitudes[px ] , STAltitudes[py] );
            else{
               if(py!=-1&& px == -1 )
                   lambda = STAltitudes[py];
                else
                    if(px!=-1&& py == -1 )
                        lambda = STAltitudes[px];
                    else
                        lambda =STAltitudes[cx]; // arrived to root
            }
            
            prevD = D;
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            
            if(px!=-1 && STAltitudes[px ] == lambda )
                cx =px;
            if(py!=-1 && STAltitudes[py ] == lambda )
                cy =py;
            
        }
        if(prevD < prevLambda )
            lambdaUpper = prevLambda + 0.0001;
        else lambdaUpper = prevLambda;
        
        addValueNormal(&positiveInterval, &LastNode, lambdaLower);
        addValueNormal(&positiveInterval, &LastNode, lambdaUpper);
       // if(cinterval>=1)
       //    printf("  bounds <%.3f , %.3f > \n", lambdaLower, lambdaUpper );
        *numberInterval+=1; // counting number of intervals
        cinterval+=1;
    }
    while(px!=-1 && py!=-1);
    
    //if(cinterval>=1){
        //printInterval(x,y,positiveInterval);
        //printf(" intervals %i\n", cinterval);
        //}
    switch(op){
        case 1: //Min: min value on the first observation interval
            lambda_star  = positiveInterval->value;
            break;
        case 2: // Max, first value on the last observation interval
            lambda_star  = LastNode->last->value;
            
                            //printInterval(x,y,positiveInterval);
                            //printf("lambda %.3f\n", lambda_star);
            
            break;
        case 3:  // min lowet bound of filtered positive intervals
            T=10;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 4:
            T=100;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 5:
            T=500;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            
            break;
        case 6:
            T=10;
            
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;            
            break;
        case 7:
            T=100;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 8:
            T=500;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 9:
            T=700;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 10:
            T=900;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
            
        case 11:
            //P = 0.3
            //P = 0.1
            //P = 0.01
            //P = 0.05
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.01 );
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
//            if( cinterval >3){
//            printInterval(x,y,negativeInterval);
//            printf("lambda %.3f\n",lambda_star );
//            }
            
            break;
        case 12:
            //P = 0.3;
            //P = 0.1
            //P=0.01
            //P = 0.005
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.005 );
           
            //for boundaries
            //lambda_star = getRankNegativeBound(negativeInterval, lastNinterval,p_obs);
            //For rank value only
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
           // printf("lambda %.3f \n", lambda_star);
            
            break;
            
        case 13:
            // 0.01 on positive
            //p_obs = floor( total_l*0.01 );
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.9);
            break;
        case 14:
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.99);
            break;
            
        case 15:
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.01);
            break;
        
        case 16:
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.001);
//            if(cinterval>1){
//                printInterval(x,y,positiveInterval);
//                printf("lambda %.3f\n", lambda_star);
//            }
            break;
        
        case 17:
            
            
            LastNode->value = LastNode->last->value + 10;
            lambda_star = getRankPositiveBound(positiveInterval, 0.99999);
            
            
                            printInterval(x,y,positiveInterval);
                            printf("lambda %.3f\n", lambda_star);
            
            
            break;
            
            
        case 20:
            
        //Silvio's :
            //identify negative intervals
           // getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            //lambda_star = lastNinterval->value;
           //lambda_star  = positiveInterval->value;
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.01 );
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
            
            //if(cinterval>1)
               // printInterval(x,y,negativeInterval);
         //  printf("lambda %.3f \n", lambda_star);
            break;
        case 21:
            
            lambda_star  = LastNode->last->value;
            break;
    }
    
    if(cinterval>  *max_numberInterval)
    {
        *max_numberInterval = cinterval;
        
    }
   // printf("  ->  edge %i,%i   lambda %.3f \n",x,y, lambda_star);
    
    while(positiveInterval!= NULL){
        LastNode = positiveInterval;
        positiveInterval = positiveInterval->next;
        free(LastNode);
    }
    
    return lambda_star;
}




double minimizationLambdasIntervalVolume(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op, int* numberInterval, int *max_numberInterval)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    
    double lambdaUpper = -9999;
    double lambdaLower = -9999;
    double total_l;
    NodeList* positiveInterval = NULL;
    NodeList* positiveIntervalVol = NULL;
    NodeList* negativeIntervalVol = NULL;
    
    
    
    NodeList* negativeInterval = NULL;
    NodeList* LastNode = positiveInterval;
    NodeList* ln_neg = negativeIntervalVol;
    NodeList* ln_pos = positiveIntervalVol;
    NodeList* lastNinterval = NULL;
    int cn_intervals = -1;
    
    cx = x;
    cy = y;
    int T = -1;
    double p_obs=0;
    
    
    lambda = STAltitudes[cx];
    double prevLambda=-1, lambda_star =-9999;
    int cinterval=0;
    double prevD = -1;
    double min_d, max_d ;
    //  printf("Edge %i,%i \n",x,y );
    do{
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        prevD = -1;
        prevLambda = lambda;
        while( D >lambda ){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            //prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            if(D>lambda && D!= prevD ){
                addValueDis(&negativeIntervalVol, &ln_neg, prevLambda, D);
                addValueDis(&negativeIntervalVol, &ln_neg, lambda, D);
                prevD = D;
                prevLambda = lambda;
            }
            else{
                if(D<=lambda){
//                    addValueDis(&positiveIntervalVol, &ln_pos, prevLambda, D);
//                    addValueDis(&positiveIntervalVol, &ln_pos, lambda, D);
//                    ln_pos->bound = 1;
                    break;
                    }
                }
            //  printf( "%.2f,%i,%i,%.2f,%.2f,%.2f\n " ,  lambda ,STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy]  ,  D );
            if(STAltitudes[px ] == lambda )
                cx =px;
            if(STAltitudes[py ] == lambda )
                cy =py;
        }
        
        if(D <= prevLambda) {
            lambdaLower = prevLambda + 0.0001;
        }
        else{
            lambdaLower = D;
            addValueDis(&negativeIntervalVol, &ln_neg, prevLambda, D);
            addValueDis(&negativeIntervalVol, &ln_neg, lambdaLower, D);
        }
        
        ln_neg->bound = 1;
        prevLambda = lambdaLower;
        prevD = -1;
        
        while( D<=lambda  &&( px!=-1&&py!=-1)){
            
            
            if(px!=-1&&py!=-1)
                lambda = min(STAltitudes[px ] , STAltitudes[py] );
            else{
                if(py!=-1&& px == -1 )
                    lambda = STAltitudes[py];
                else
                    if(px!=-1&& py == -1 )
                        lambda = STAltitudes[px];
                    else{
                        //lambda = STAltitudes[cx]; // arrived to root
                        break;
                    }
                
            }
            
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            if(D<=lambda ){
                addValueDis(&positiveIntervalVol, &ln_pos, prevLambda, D);
                addValueDis(&positiveIntervalVol, &ln_pos, lambda, D);
                prevD = D;
                 prevLambda = lambda;
            }
          
            
            
            if(px!=-1 && STAltitudes[px ] == lambda )
                cx =px;
            if(py!=-1 && STAltitudes[py ] == lambda )
                cy =py;
            
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;

            
//            if(px!=-1 && STAltitudes[px ] == lambda )
//                cx =px;
//            if(py!=-1 && STAltitudes[py ] == lambda )
//                cy =py;
          //  prevLambda = lambda;
        }
        if(prevD < prevLambda )
            lambdaUpper = prevLambda + 0.0001;
        else lambdaUpper = prevLambda;
        
        ln_pos->bound = 1;
        
        
        addValueNormal(&positiveInterval, &LastNode, lambdaLower);
        addValueNormal(&positiveInterval, &LastNode, lambdaUpper);
        // if(cinterval>=1)
        //    printf("  bounds <%.3f , %.3f > \n", lambdaLower, lambdaUpper );
        *numberInterval+=1; // counting number of intervals
        cinterval+=1;
    }
    while(px!=-1 && py!=-1);
    
    //if(cinterval>=1){
    //printInterval(x,y,positiveInterval);
    //printf(" intervals %i\n", cinterval);
    //}
    switch(op){
        case 1: //Min: min value on the first observation interval
            lambda_star  = positiveInterval->value;
            break;
        case 2: // Max, first value on the last observation interval
            lambda_star  = LastNode->last->value;
            break;
        case 3:  // min lowet bound of filtered positive intervals
            T=10;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 4:
            T=100;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 5:
            T=500;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 6:
            T=10;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 7:
            T=100;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 8:
            T=500;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 9:
            T=700;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 10:
            T=900;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
            
        case 11:
            
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.01 );
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
            //            if( cinterval >3){
            //            printInterval(x,y,negativeInterval);
            //            printf("lambda %.3f\n",lambda_star );
            //            }
            
            break;
        case 12:
            
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.005 );
            
            //for boundaries
            //lambda_star = getRankNegativeBound(negativeInterval, lastNinterval,p_obs);
            //For rank value only
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
            // printf("lambda %.3f \n", lambda_star);
            
            break;
            
        //volume;
        case 15:
            T = 1650000;
            //printIntervalDis(x,y,positiveIntervalVol );
            filterIntervalArea(&positiveIntervalVol, T, &LastNode);
           // printInterval(x,y,positiveInterval );
            
            //////
//            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
//            p_obs = floor( total_l*0.005 );
//            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
//            printf("lambda %.3f\n",lambda_star );
//            printf("New lambda %.3f\n",positiveIntervalVol->value);
            lambda_star =positiveIntervalVol->value;
            
            break;
        case 16:
            
            T = 50000000;
            //printIntervalDis(x,y,positiveIntervalVol );
            negativeIntervalVol->value = -99999999;
            
            filterIntervalArea(&negativeIntervalVol, T, &LastNode);
            
            //////
//            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
//            p_obs = floor( total_l*0.005 );
//            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
//            printf("lambda %.3f\n",lambda_star );
//            printf("New lambda %.3f\n",LastNode->value);
            lambda_star = LastNode->value;
            break;
        case 17:
            
            T = 15;
            
            negativeIntervalVol->value = -99999999;
            
            filterIntervalArea(&negativeIntervalVol, T, &LastNode);
            
            p_obs = 0.005;
            
            lambda_star = getRankNegativeVolume(negativeIntervalVol, LastNode,p_obs);
            
            // printInterval(x,y,positiveInterval );
            
            //////
            //            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            //            p_obs = floor( total_l*0.005 );
            //            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
            //            printf("lambda %.3f\n",lambda_star );
            //            printf("New lambda %.3f\n",LastNode->value);
           // lambda_star = LastNode->value;
            break;
            
        case 18:
            // Using depth negative, lambda, min
            negativeIntervalVol->value = 0;
            computeDepths(&negativeIntervalVol, &LastNode, &min_d, &max_d);
            lambda_star = min_d;
            break;
        case 19:
            // Using depth negative, lambda, max
            negativeIntervalVol->value = 0;
            computeDepths(&negativeIntervalVol, &LastNode, &min_d, &max_d);
            lambda_star = max_d;
            break;
        case 20:
            // Using depth negative, bound min
            negativeIntervalVol->value = 0;
            computeDepthsNegBound(&negativeIntervalVol, &LastNode, &min_d, &max_d);
            lambda_star = min_d;
            break;
        case 21:
            // Using depth negative, bound min
            negativeIntervalVol->value = 0;
            computeDepthsNegBound(&negativeIntervalVol, &LastNode, &min_d, &max_d);
            lambda_star = max_d;
            break;
        case 22:
            // Using depth positive, lambda, min
            ln_pos->value = ln_pos->last->value+1;
            computeDepths(&positiveIntervalVol, &ln_pos, &min_d, &max_d);
            lambda_star = min_d;
            break;
        case 23:
            // Using depth positive, lambda, max
            ln_pos->value = ln_pos->last->value+1;
            computeDepths(&positiveIntervalVol, &ln_pos, &min_d, &max_d);
            lambda_star = max_d;
            break;
        case 24:
            ln_pos->value = ln_pos->last->value+1;
            computeDepthsPosBound(&positiveIntervalVol, &ln_pos, &min_d, &max_d);
            lambda_star = min_d;
            break;
        case 25:
            ln_pos->value = ln_pos->last->value+1;
            computeDepthsPosBound(&positiveIntervalVol, &ln_pos, &min_d, &max_d);
            lambda_star = max_d;
            break;
            
            
            
    }
    
    if(cinterval>  *max_numberInterval)
    {
        *max_numberInterval = cinterval;
        
    }
    // printf("  ->  edge %i,%i   lambda %.3f \n",x,y, lambda_star);
    
    while(positiveInterval!= NULL){
        LastNode = positiveInterval;
        positiveInterval = positiveInterval->next;
        free(LastNode);
    }
    
    return lambda_star;
}



double scaleSelection(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , int op, int th, double prank)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    
    double lambdaUpper = -9999;
    double lambdaLower = -9999;
   // double total_l;
   // NodeList* positiveInterval = NULL;
    NodeList* positiveIntervalVol = NULL;
    NodeList* negativeIntervalVol = NULL;
    
    //NodeList* negativeInterval = NULL;
    //NodeList* LastNode = positiveInterval;
    NodeList* ln_neg = negativeIntervalVol;
    NodeList* ln_pos = positiveIntervalVol;
  //  NodeList* lastNinterval = NULL;
  //  int cn_intervals = -1;
    
    cx = x;
    cy = y;
    int T = -1;
    double p_obs=0;
    
    
    lambda = STAltitudes[cx];
    double prevLambda=-99999999, lambda_star =-9999;
    int cinterval=0;
    double prevD = -1;
    
    
   // printf("Edge %i,%i   intervals \n",x,y );
    do{
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        prevD = -1;
        prevLambda = lambda;
        while( D >lambda ){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            if(D>lambda  ){
                addValueDis(&negativeIntervalVol, &ln_neg, prevLambda, D);
                addValueDis(&negativeIntervalVol, &ln_neg, lambda, D);
                prevD = D;
                //prevLambda = lambda;
                if(STAltitudes[px ] == lambda )
                    cx =px;
                if(STAltitudes[py ] == lambda )
                    cy =py;
            }
        }
        if(D <= prevLambda) {
            lambdaLower = prevLambda + 0.009;
            //lambdaLower = prevLambda ;
        }
        else{
            lambdaLower = D ;
            addValueDis(&negativeIntervalVol, &ln_neg, prevLambda, D);
            addValueDis(&negativeIntervalVol, &ln_neg, lambdaLower, D);
        }
        
        ln_neg->bound = 1;
        prevLambda = lambdaLower;
        prevD = -99999999;
        
//        if(lambda == prevLambda ){
//            if(px!=-1 && STAltitudes[px ] == lambda )
//                cx =px;
//            if(py!=-1 && STAltitudes[py ] == lambda )
//                cy =py;
//        }
        
        
        while( D<=lambda  &&( px!=-1&&py!=-1)){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            
            if(px!=-1&&py!=-1)
                lambda = min(STAltitudes[px ] , STAltitudes[py] );
            else{
                if(py!=-1&& px == -1 )
                    lambda = STAltitudes[py];
                else
                    if(px!=-1&& py == -1 )
                        lambda = STAltitudes[px];
                    else{
                         break;
                    }
            }
            
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            if(D<=lambda  ){
                addValueDis(&positiveIntervalVol, &ln_pos, prevLambda, D);
                addValueDis(&positiveIntervalVol, &ln_pos, lambda, D);
                prevD = D;
                prevLambda = lambda;
                if(px!=-1 && STAltitudes[px ] == lambda )
                    cx =px;
                if(py!=-1 && STAltitudes[py ] == lambda )
                     cy =py;
            }
//            else if(prevLambda == lambda ){
//                    if(px!=-1 && STAltitudes[px ] == lambda )
//                        cx =px;
//                    if(py!=-1 && STAltitudes[py ] == lambda )
//                        cy =py;
//                 }
            

        }
        if(prevD < prevLambda )
            lambdaUpper = prevLambda + 0.009;
            //lambdaUpper = prevLambda ;
        else lambdaUpper = prevLambda;
        
        ln_pos->bound = 1;
        
        //addValueNormal(&positiveInterval, &LastNode, lambdaLower);
        //addValueNormal(&positiveInterval, &LastNode, lambdaUpper);
        
        cinterval+=1;
    }
    while(px!=-1 && py!=-1);

    switch(op){
        case 1: //Min: min value on positive observation intervals
            lambda_star  = positiveIntervalVol->value;
            break;
        case 2: // Max,  last upper bound on negative intervals
            
//            printf("\t Neg \n");
//            printIntervalDis(x,y,negativeIntervalVol);
//            printf("\t Pos \n");
//            printIntervalDis(x,y,positiveIntervalVol);
//            printf("\t All positive intervals \n" );
//            printInterval(x,y,positiveInterval );
//
//            printf("lambda last note, %.3f \t %.3f  \n ", LastNode->last->value , ln_neg->value);
            lambda_star = ln_neg->value;
            
            break;
        case 3:  // On positive intervals, apply length thresholding and  min-rule
            //T=10;
            T = th;
            filterIntervalLength(&positiveIntervalVol,T,&ln_pos  );
            lambda_star  = positiveIntervalVol->value;
            
           // printf("lambda star note, %.3f   \n ", lambda_star);
            break;
            
        case 4:
            // On negative intervals, apply length thresholding and  max-rule
            T=th;            
           // printf("\t Neg \n");
            //printIntervalDis(x,y,negativeIntervalVol);
            
            negativeIntervalVol->value = -9999999999;  // force negative to have a large interval by the right (to deal with case of only one negative interval)
            
            filterIntervalLength(&negativeIntervalVol,T,&ln_neg  );
            
            //printf("\t Neg Filtered\n");
            //printIntervalDis(x,y,negativeIntervalVol);
            
            lambda_star  = ln_neg->value;
            break;
        case 5:
            // On positive intervals, apply rank and  min-rule
            ln_pos->value =ln_pos->last->value+1; // limit the upper bound 
            lambda_star = getRankPositive(positiveIntervalVol, prank);
            break;
        case 6:
            // On negative intervals, apply rank and  max-rule
            lambda_star = getRankNegativeVolume(negativeIntervalVol, ln_neg,prank);
            break;
            
        case 7:
            
            // On positive intervals, apply area and  min-rule
            T = th;
            filterIntervalArea(&positiveIntervalVol, T, &ln_pos);
            lambda_star =positiveIntervalVol->value;
            
            break;
        case 8:
            // On negative intervals, apply area and  max-rule
            T = th;
            negativeIntervalVol->value = -99999999;
            filterIntervalArea(&negativeIntervalVol, T, &ln_neg);
            lambda_star = ln_neg->value;
            break;
        case 9:
             // On negative intervals, apply area, then ranking and  max-rule
            T = th;
            negativeIntervalVol->value = -99999999;
            filterIntervalArea(&negativeIntervalVol, T, &ln_neg);
            p_obs = prank;
            lambda_star = getRankNegativeVolume(negativeIntervalVol, ln_neg,p_obs);
            break;      
        case 10 ://On positive intervals, apply volume filter and min-rule    
            filterIntervalVolume(&positiveIntervalVol, th, &ln_pos);
            lambda_star = positiveIntervalVol->value;
            break;
        case 11:// On negative intervals, apply volume filter and  max-rule
            negativeIntervalVol->value = -99999999;
            filterIntervalVolume(&negativeIntervalVol, th, &ln_neg);
            lambda_star = ln_neg->value;
            break;

        case 12:// On positive intervals, apply depth filter and  min-rule            
            filterIntervalDepths(&positiveIntervalVol, th, &ln_pos);
            lambda_star = positiveIntervalVol->value;
            break;
        case 13:// On negative intervals, apply depth filter and max-rule
            negativeIntervalVol->value = -99999999;
            filterIntervalDepths(&negativeIntervalVol, th, &ln_neg);
            lambda_star = ln_neg->value;
            break;
                
            
    }
    
    
//    while(positiveInterval!= NULL){
//        LastNode = positiveInterval;
//        positiveInterval = positiveInterval->next;
//        free(LastNode);
//    }

    while(positiveIntervalVol!= NULL){
        ln_pos = positiveIntervalVol;
        positiveIntervalVol = positiveIntervalVol->next;
        free(ln_pos);
    }
    
    while(negativeIntervalVol!= NULL){
        ln_neg = negativeIntervalVol;
        negativeIntervalVol = negativeIntervalVol->next;
        free(ln_neg);
    }
    
    return lambda_star;
}




void reweightToMax(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes)
{
    JCctree* CT = *CompTree;
    
    
    double   lambda;
    int i;
    
     for(i=0;i<mst->size;i++){
        
         
        int cx, cy, px, py;
        int x = mst->MSTedges[i].x;
        int y = mst->MSTedges[i].y;
         
         cx = x;
         cy = y;
        lambda = STAltitudes[cx];
         double newf = -9999;
            while(  cx!=cy ){ // search common parent
                px = CT->tabnodes[cx].father;
                py = CT->tabnodes[cy].father;
                
                if(px!=-1&&py!=-1)
                    lambda = min(STAltitudes[px ] , STAltitudes[py] );
                else{
                    if(py!=-1&& px == -1 )
                        lambda = STAltitudes[py];
                    else  if(px!=-1&& py == -1 )
                            lambda = STAltitudes[px];
                          else
                            lambda =STAltitudes[cx]; // arrived to root
                }
                px = CT->tabnodes[cx].father;
                py = CT->tabnodes[cy].father;
                if(px!=-1 && STAltitudes[px ] == lambda )
                    cx =px;
                if(py!=-1 && STAltitudes[py ] == lambda )
                    cy =py;
            }
         
         px = CT->tabnodes[cx].father;
         
         if(px == CT->root || px==-1 ){
             newf= STAltitudes[cx]+1;
         }
         else{
             
             newf= STAltitudes[px];
         }
         
         setweight(g,i,newf);
         
           //  printf("edge %i,%i   lambda %.3f \n",x,y, newf);
         
         
     }
    
    
    
}



/*
void clean_cache(){

    const int size = 120*1024*1024; // Allocate 20M. Set much larger then L2
    char *c = (char *)malloc(size);
    for (int i = 0; i < 0xffff; i++)
        for (int j = 0; j < size; j++)
            c[j] = i*j;

}*/


void incrementalSegmentationBranchTimed(graphe *g){
    clock_t ini_all, fin_all, ini_QFZ, fin_QFZ, ini_mini, fin_mini, ini_set, fin_set;
    double sum_QFZ=0, sum_mini = 0 , sum_all=0, sum_loop=0, sum_set=0;
    
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
    
    //////////////////  Incremental Segmentation  //////////////
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMaxWeight );
    
    
    /// Set all edges to large values for later SM computation ////
    
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    
    ini_all = clock();
    
    for(i=0;i<mst->size;i++){
        
        ini_mini = clock();
        newf = branchMinimization(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ) - 1;
        fin_mini = clock();
        sum_mini += ((double) (fin_mini - ini_mini)) / CLOCKS_PER_SEC;
        
        //  printf("New value for edge: %i,%i = %f , former %f\n",mst->MSTedges[id].x, mst->MSTedges[id].y,newf ,mst->MSTedges[id].weight );
       // ini_set = clock();
        setweight(newG,i,newf);
       // fin_set = clock();
       // sum_set += ((double) (fin_set - ini_set)) / CLOCKS_PER_SEC;
        
        ini_QFZ = clock();
        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation,  STArea, STMaxWeight, newG, i);
        fin_QFZ = clock();
        sum_QFZ += ((double) (fin_QFZ - ini_QFZ)) / CLOCKS_PER_SEC;
    }
    
    //componentTreePrintAttributes(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight );
    //  char * filename=  "/tmp/InGraph.dot";
    //  saveDOTfile(CompTree, STAltitudesSegmentation, CompTree->root, STArea, STMaxWeight,  filename, g->nbsom);
    //////////////////////////////////////////////////////////////
    
    fin_all = clock();
    sum_all += ((double) (fin_all - ini_all)) / CLOCKS_PER_SEC;
    sum_loop = sum_all - sum_QFZ - sum_mini - sum_set;
    
    printf("4\t all\t %f \t loop: \t %f \t QFZ: \t %f \t MINI: \t %f \t DIS \t %f \t SETS\t %f \n ", sum_all , sum_loop, sum_QFZ , sum_mini, diss_sum, sum_set);
    
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}

void incrementalSegmentationBranch(graphe *g){
    
    
    clock_t startAll, endALL; double cpu_time_usedALL;
   
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
/////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));

    
    //////////////////  Incremental Segmentation  //////////////
    
    
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMaxWeight );

    
/// Set all edges to large values for later SM computation ////
    
    
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    startAll = clock();
    for(i=0;i<mst->size;i++){
        newf = branchMinimization(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ) - 1;
        printf("New value for edge: %i,%i = %f , former %f\n",mst->MSTedges[i].x, mst->MSTedges[i].y,newf ,mst->MSTedges[i].weight );
        setweight(newG,i,newf);
        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation,  STArea, STMaxWeight, newG, i);
    }
    
    endALL = clock();    
    cpu_time_usedALL = ((double) (endALL - startAll)) / CLOCKS_PER_SEC;
    
    printf("\t total\t %.5f\n",cpu_time_usedALL );
    //componentTreePrintAttributes(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight );
  //  char * filename=  "/tmp/InGraph.dot";
  //  saveDOTfile(CompTree, STAltitudesSegmentation, CompTree->root, STArea, STMaxWeight,  filename, g->nbsom);
//////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}




void incrementalSegmentationInterval(graphe *g, int op){
    
    
  //  clock_t startAll, endALL; double cpu_time_usedALL;
    
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
    
    
    //////////////////  Incremental Segmentation  //////////////
    
    
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMaxWeight );
    
    
    /// Set all edges to large values for later SM computation ////
    
    
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    //startAll = clock();
    
    int numberInterval=0;
    int max_numberInterval=0;
    
    for(i=0;i<mst->size;i++){
        NodeList* lambdaStarInterval;
        ListInit(&lambdaStarInterval);
        
        
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 1 ;
    //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 0.0001 ;
    newf =  minimizationLambdasIntervalVolume(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 0.0001 ;
        
        
        
        
        setweight(newG,i,newf);
        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation,  STArea, STMaxWeight, newG, i);
    }
    
    //endALL = clock();
    //cpu_time_usedALL = ((double) (endALL - startAll)) / CLOCKS_PER_SEC;
    
    printf("\t mean number of intervals %i, max number of intervals  %i\n", numberInterval/mst->size, max_numberInterval );
   
//    if(op>=20){
//        //Set values to Max on last interval
//        reweightToMax(newG, mst,  &CompTree ,STAltitudesSegmentation);
//    }
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}





void intervalSegmentation(graphe *g, int op, int th, double prank){
    
    
    //  clock_t startAll, endALL; double cpu_time_usedALL;
    
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
    
    
    //////////////////  Incremental Segmentation  //////////////
    
    
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMaxWeight );
    
    
    /// Set all edges to large values for later SM computation ////
    
    
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    //startAll = clock();
    
    int numberInterval=0;
    int max_numberInterval=0;
    
    for(i=0;i<mst->size;i++){
        //NodeList* lambdaStarInterval;
        //ListInit(&lambdaStarInterval);
        
        
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 1 ;
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 0.0001 ;
        //newf =  scaleSelection(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, th, prank) - 0.000001 ;
        newf =  scaleSelection(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, th, prank) - 0.000001  ;
        
        
        
        
        setweight(newG,i,newf);
        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation,  STArea, STMaxWeight, newG, i);
    }
    
    //endALL = clock();
    //cpu_time_usedALL = ((double) (endALL - startAll)) / CLOCKS_PER_SEC;
    
    //printf("\t mean number of intervals %i, max number of intervals  %i\n", numberInterval/mst->size, max_numberInterval );
    
    //    if(op>=20){
    //        //Set values to Max on last interval
    //        reweightToMax(newG, mst,  &CompTree ,STAltitudesSegmentation);
    //    }
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}




double naiveMinimization(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx )
{
    
    JCctree* CT = *CompTree;
    int cx, cy;
    double   D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
        
    cx = x;
    cy = y;
    
   // char * filename=  "/Users/edus/Dropbox/MatMorphology/codigo/SM_Edward/Graphs/dots/graph.dot";
   // saveDOTfile(CT, STAltitudes, CT->root, STArea, STMaxWeight, filename );
    lambda = 0;
    //D = computeD(CT, STAltitudes, STArea, STMaxWeight, &cx, &cy , Dif , lambda );
    
    D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
    
    while(D > lambda){
        lambda+=1;
      //  D = computeD(CT, STAltitudes, STArea, STMaxWeight, &cx, &cy , Dif , lambda );
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        diss_ini = clock();
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        diss_fin = clock();
        diss_sum+=((double) (diss_fin - diss_ini)) / CLOCKS_PER_SEC;
    }
    
    return lambda;
}


void incrementalNaiveSegmentation(graphe *g){
    clock_t ini_all, fin_all, ini_QFZ, fin_QFZ, ini_mini, fin_mini, ini_set, fin_set;
    double sum_QFZ=0, sum_mini = 0 , sum_all=0, sum_loop=0, sum_set=0;
    
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double* STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    
    //////////////////  Incremental Segmentation  //////////////
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation , STArea, STMaxWeight );
    
    ini_all = clock();
    
    for(i=0;i<mst->size;i++){
        ini_mini = clock();
        newf =  naiveMinimization(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ) -1 ;
        fin_mini = clock();
        sum_mini += ((double) (fin_mini - ini_mini)) / CLOCKS_PER_SEC;
      //  printf("New value for edge: %i,%i = %f , former %f\n",mst->MSTedges[id].x, mst->MSTedges[id].y,newf ,mst->MSTedges[id].weight );
        ini_set = clock();
        setweight(newG,i,newf);
        fin_set = clock();
        sum_set += ((double) (fin_set - ini_set)) / CLOCKS_PER_SEC;
        
        ini_QFZ = clock();
        treeUpdateAttributes(&CompTree, STAltitudes, STAltitudesSegmentation, STArea, STMaxWeight, newG, i);
        fin_QFZ = clock();
        sum_QFZ += ((double) (fin_QFZ - ini_QFZ)) / CLOCKS_PER_SEC;
        
    }

    //componentTreePrintAttributes(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight );
    //char * filename=  "/Users/edus/Dropbox/MatMorphology/codigo/SM_Edward/Graphs/dots/graph.dot";
    // saveDOTfile(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight,  filename);
    //////////////////////////////////////////////////////////////
    
    fin_all = clock();
    sum_all += ((double) (fin_all - ini_all)) / CLOCKS_PER_SEC;
    sum_loop = sum_all - sum_QFZ - sum_mini - sum_set;
    
    printf("5\t all\t %f \t loop: \t %f \t QFZ: \t %f \t MINI: \t %f \t DIS \t %f \t SETS\t %f \n ", sum_all , sum_loop, sum_QFZ , sum_mini, diss_sum, sum_set);
    
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);

    
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STArea);
    free(STAltitudesSegmentation);
    free(mst->MSTedges);
    free(mst);
    componentTreeFree(CompTree);
    
}


/*
double computeD(JCctree *CT,double *STAltitudes , int *STArea, double* STMaxWeight , int *cx, int *cy , double Dif , double lambda ){
    
    double  sc1_c2, sc2_c1, D;
    
    *cx = getComponentAtLambda(CT,STAltitudes, lambda, *cx );
    *cy = getComponentAtLambda(CT,STAltitudes, lambda, *cy );
    sc1_c2 = (Dif -  STMaxWeight[*cx] )* STArea[*cx];
    sc2_c1 = (Dif -  STMaxWeight[*cy] )* STArea[*cy];
    
    D = max(sc1_c2 , sc2_c1);
    return D;

}
*/
 // D = computeD(AreaCx, AreaCy, STMaxWeightCx, STMaxWeightCy, Dif);

double computeD( int AreaCx, int AreaCy, double STMaxWeightCx, double STMaxWeightCy, double Dif){
    
    double  sc1_c2, sc2_c1, D;
    
    sc1_c2 = (Dif -  STMaxWeightCx )* AreaCx;
    sc2_c1 = (Dif -  STMaxWeightCy )* AreaCy;
    
    D = max(sc1_c2 , sc2_c1);
    
    if(D<0)D=0; // thersholding dissimilarity measure
    
    return D;
}



double rangeMinimization(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx , NodeList**frange)
{
    JCctree* CT = *CompTree;
    int cx, cy;
    double  lambda,D;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    double lambda_star = -9999;
    
    cx = x;
    cy = y;
    
    //char * filename=  "/Users/edus/Dropbox/MatMorphology/codigo/SM_Edward/Graphs/dots/graph.dot";
    //saveDOTfile(*CompTree, STAltitudes, (*CompTree)->root, STArea, STMaxWeight,  filename);
    
    NodeList* node = *frange;
    lambda = node->value;
    double prevLambda =  -1;
    
    //D = computeD(CT, STAltitudes, STArea, STMaxWeight, &cx, &cy , Dif , lambda );
    D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
    
    while( D > lambda ){
        prevLambda = lambda; node= node->next;
        lambda = node->value;
        
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        //D = computeD(CT, STAltitudes, STArea, STMaxWeight, &cx, &cy , Dif , lambda );
        diss_ini = clock();
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        diss_fin = clock();
        diss_sum+=((double) (diss_fin - diss_ini)) / CLOCKS_PER_SEC;
    }
    
    
    if(D <= prevLambda) lambda_star = prevLambda + 1;
    else lambda_star = D;
    if(lambda_star - 1 != prevLambda)
        addValue(frange, node, lambda_star-1);
            
    return lambda_star;
}



void incrementalSegmentationRange(graphe *g){

    clock_t ini_all, fin_all, ini_QFZ, fin_QFZ, ini_mini, fin_mini, ini_set, fin_set;
    double sum_QFZ=0, sum_mini = 0 , sum_all=0, sum_loop=0, sum_set=0;
    
    
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double* STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
  
    //////////////////  Incremental Segmentation  //////////////
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMaxWeight );
    
    
    
    NodeList* frange;
    ListInit(&frange);
   /// Set all edges to large value for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    
    ini_all = clock();
    
   
    for(i=0;i<mst->size;i++){
        
        ini_mini = clock();
        newf =  rangeMinimization(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i, &frange ) - 1 ;
        fin_mini = clock();
        sum_mini += ((double) (fin_mini - ini_mini)) / CLOCKS_PER_SEC;
        
      //  printf("New value for edge: %i,%i = %f , former %f\n",mst->MSTedges[id].x, mst->MSTedges[id].y,newf ,mst->MSTedges[id].weight );
        ini_set = clock();
        setweight(newG,i,newf);
        fin_set = clock();
        sum_set += ((double) (fin_set - ini_set)) / CLOCKS_PER_SEC;
        
        
        ini_QFZ = clock();
        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation ,  STArea, STMaxWeight, newG, i);
        fin_QFZ = clock();
        sum_QFZ += ((double) (fin_QFZ - ini_QFZ)) / CLOCKS_PER_SEC;
    }
    //componentTreePrintAttributes(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight );
    //char * filename=  "/Users/edus/Dropbox/MatMorphology/codigo/SM_Edward/Graphs/dots/graph.dot";
    // saveDOTfile(CompTree, STAltitudes, CompTree->root, STArea, STMaxWeight,  filename);
    //////////////////////////////////////////////////////////////
    
    fin_all = clock();
    sum_all += ((double) (fin_all - ini_all)) / CLOCKS_PER_SEC;
    sum_loop = sum_all - sum_QFZ - sum_mini - sum_set;
    
    printf("6 \t all\t %f \t loop: \t %f \t QFZ: \t %f \t MINI: \t %f \t DIS \t %f \t SETS\t %f \n ", sum_all , sum_loop, sum_QFZ , sum_mini, diss_sum, sum_set);
    
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    
    
    
    NodeList* tmp = frange, *tmp2 = frange;
    while(tmp!=NULL){
        tmp2= tmp->next;
        free(tmp);
        tmp = tmp2;
    }
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STArea);
    free(STAltitudesSegmentation);
    free(mst->MSTedges);
    free(mst);
    componentTreeFree(CompTree);
    
    
}



int computeLevels(graphe *g, double *sum){
    
    *sum = 0;
    TriRapideStochastiqueDouble(g->weight, 0, g->nbmaxar-1);
    double last=-INFINITO;
    int counter = 0;
    int i;
    for(i = 0; i< g->nbmaxar ; i++ ){
        
        if(g->weight[i] != last )
            counter++;
        last =g->weight[i];
        // printf("weight %f \n", g->weight[i]);
        if(g->weight[i]!= INFINITO){

            *sum +=g->weight[i];
        }
    }
    return counter;
    
    
}




///////// Dual HGB method
void dualhgb(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        newf =  max_rule(newG, mst,  &CompTree  , STAltitudes , STArea, STMaxWeight, i)  ;
        setweight(newG,i,newf);
    }
    
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
   // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}


void dualhgbInverse(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    double maxVal = -1;
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        newf =  max_rule(newG, mst,  &CompTree  , STAltitudes , STArea, STMaxWeight, i)  ;
        setweight(newG,i,newf);
        if(newf > maxVal )
            maxVal = newf;
    }
    //////////////////////////////////////////////////////////////
    for(i=mst->size-1;i>=0;i--){
        printf("Edge %i,%i:  %.3f \n",mst->MSTedges[i].x ,mst->MSTedges[i].y, getweight(newG,i));
    }
    printf("nuevo \n");
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i, fabs( maxVal - getweight(newG,i)));
        printf("Edge %i,%i:  %.3f \n",mst->MSTedges[i].x ,mst->MSTedges[i].y, getweight(newG,i));
    }
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}


double max_rule(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx){
    
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double   D=0, lambda, max_lambda=-9999;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    cx = x; px = cx;
    cy = y; py = cy;
    
    lambda = 0;
    
    do{
        lambda+=1;
        
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        if( D >lambda ){
            max_lambda = lambda;
        }
       // px = CT->tabnodes[cx].father;
       // py = CT->tabnodes[cy].father;
    }
    //while(  px != py &&  CT->tabnodes[cx].father !=-1);
    while(D>(lambda ) );
    
//    if(D > lambda){
//       max_lambda = D;
//    }
    if(max_lambda<0){max_lambda = 0;}
        
    printf("Edge %i,%i:  %.3f \n",x,y,max_lambda);
        
        
        return max_lambda;
//        px = CT->tabnodes[cx].father;
//        py = CT->tabnodes[cy].father;
//
//
//    }
//    while(D>(lambda ) );
//
//    if(max_lambda<0){max_lambda = 0;}
//
//   // printf("Edge %i,%i:  %.3f \n",x,y,max_lambda);
//
//
//    return max_lambda;
    
}



void dualhgbParent(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        newf =  max_rule_parent(newG, mst,  &CompTree  , STAltitudes , STArea, STMaxWeight, i)  ;
        setweight(newG,i,newf);
       
    }

    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}



double max_rule_parent(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx){
    
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double   D=0, lambda, max_lambda=-9999;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    cx = x; px = cx;
    cy = y; py = cy;
    double al_cx = -99 , al_cy = -99;
    
    lambda = 1;
    do{
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        D  = computeD( STArea[cx], STArea[cy], STMaxWeight[cx], STMaxWeight[cy], Dif );
        if(al_cx < 0 && STAltitudes [ cx ] >0)
            al_cx = STAltitudes [ cx ] ;
        if(al_cy < 0 && STAltitudes [ cy ] >0)
            al_cy = STAltitudes [ cy ] ;
        
        if( D>lambda ){
            max_lambda = lambda;
        }
        
        lambda+=1;
        
    }
    while(D>(lambda ) );
    
    double th;
    if(max(al_cx, al_cy) >0){
        th = max(al_cx, al_cy);
        if ( max_lambda > th ) max_lambda = th;
    }
    if(max_lambda<0){max_lambda = 0;}
    // printf("Edge %i,%i:  %.3f possible %.3f \n",x,y,max_lambda, min( th , max_lambda)  );
     printf("Edge %i,%i:  %.3f possible %.3f \n",x,y,max_lambda,  max_lambda);
    
    return max_lambda;
    
}



void dualhgbMin(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMinWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        computeQFZ_AttributesMin(newG, &CompTree, &STAltitudes, &STArea , &STMinWeight);
        newf =  max_rule_min(newG, mst,  &CompTree  , STAltitudes , STArea, STMinWeight, i)  ;
        setweight(newG,i,newf);
    }
    
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMinWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}



double max_rule_min(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMinWeight, int edgeIdx){
    
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double   D=0, lambda, max_lambda=-9999;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    cx = x; px = cx;
    cy = y; py = cy;
    
    lambda = 0;
    
    
    do{
        lambda+=1;
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        
        D = computeDmin(STArea[cx], STArea[cy] ,STMinWeight[cx], STMinWeight[cy], Dif );
        if( D>lambda ){
            max_lambda = lambda;
        }
        px = CT->tabnodes[cx].father;
        py = CT->tabnodes[cy].father;
    }
    //while(  px != py &&  CT->tabnodes[cx].father !=-1);
    while(D>(lambda ) );
    
   // if(D > lambda){
   //    max_lambda = D;
   // }
    
    if(max_lambda<0){max_lambda = 0;}
    
    printf("Edge %i,%i:  %.3f \n",x,y,max_lambda);
    
    
    return max_lambda;
    
}

double max_rule_mean(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMinWeight, int edgeIdx){
    
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double   D=0, lambda, max_lambda=-9999;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    cx = x; px = cx;
    cy = y; py = cy;
    
    lambda = 0;
    
    
    do{
        lambda+=1;
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        
        D = computeDmean(STArea[cx], STArea[cy] ,STMinWeight[cx], STMinWeight[cy], Dif );
        if( D>lambda ){
            max_lambda = lambda;
        }
        px = CT->tabnodes[cx].father;
        py = CT->tabnodes[cy].father;
    }
    //while(  px != py &&  CT->tabnodes[cx].father !=-1);
     while(D>(lambda ) );
    
    //if(D > lambda){
    //    max_lambda = D;
   // }
    
    if(max_lambda<0){max_lambda = 0;}
    
   //2 printf("Edge %i,%i:  %.3f \n",x,y,max_lambda);
    
    
    return max_lambda;
    
}
double computeDmean( int AreaCx, int AreaCy, double STMinWeightCx, double STMinWeightCy, double Dif){
    
    double  sc1_c2, sc2_c1, D;
    
    sc1_c2 = (Dif -  STMinWeightCx )* AreaCx;
    sc2_c1 = (Dif -  STMinWeightCy )* AreaCy;
    
    D = max(sc1_c2 , sc2_c1);
    
    if(D<0)D=0; // thersholding dissimilarity measure
    
    return D;
}

double computeDmin( int AreaCx, int AreaCy, double STMinWeightCx, double STMinWeightCy, double Dif){
    
    double  sc1_c2, sc2_c1, D;
    
    sc1_c2 = (Dif -  STMinWeightCx )* AreaCx;
    sc2_c1 = (Dif -  STMinWeightCy )* AreaCy;
    
    D = min(sc1_c2 , sc2_c1);
    
    if(D<0)D=0; // thersholding dissimilarity measure
    
    return D;
}



void dualhgbParentMin(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        computeQFZ_AttributesMin(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        newf =  max_rule_parent_min(newG, mst,  &CompTree  , STAltitudes , STArea, STMaxWeight, i)  ;
        setweight(newG,i,newf);
        
    }
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}



double max_rule_parent_min(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx){
    
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double   D=0, lambda, max_lambda=-9999;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    cx = x; px = cx;
    cy = y; py = cy;
    double al_cx = -99 , al_cy = -99;
    
    lambda = 0;
    do{
        lambda+=1;
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        D  = computeDmin( STArea[cx], STArea[cy], STMaxWeight[cx], STMaxWeight[cy], Dif );
        if(al_cx < 0 && STAltitudes [ cx ] >0)
            al_cx = STAltitudes [ cx ] ;
        if(al_cy < 0 && STAltitudes [ cy ] >0)
            al_cy = STAltitudes [ cy ] ;
        
        if( D>lambda ){
            max_lambda = lambda;
        }
        
    }
    while(D>(lambda ) );
    
    double th;
    if(max(al_cx, al_cy) >0){
        th = max(al_cx, al_cy);
        if ( max_lambda > th ) max_lambda = th;
    }
    if(max_lambda<0){max_lambda = 0;}
    // printf("Edge %i,%i:  %.3f possible %.3f \n",x,y,max_lambda, min( th , max_lambda)  );
   // printf("Edge %i,%i:  %.3f possible %.3f \n",x,y,max_lambda,  max_lambda);
    
    return max_lambda;
    
}



/////////////////////////

void dualhgbMean(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    int* STArea = NULL;
    double * STMinWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
       // computeQFZ_AttributesMin(newG, &CompTree, &STAltitudes, &STArea , &STMinWeight);
        computeQFZ_AttributesMean(newG, &CompTree, &STAltitudes, &STArea , &STMinWeight);
        newf =  max_rule_min(newG, mst,  &CompTree  , STAltitudes , STArea, STMinWeight, i)  ;
        setweight(newG,i,newf);
    }
    
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STMinWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}


void dualhgbMeanNormal(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    int* STArea = NULL;
    double * STMinWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        // computeQFZ_AttributesMin(newG, &CompTree, &STAltitudes, &STArea , &STMinWeight);
        computeQFZ_AttributesMean(newG, &CompTree, &STAltitudes, &STArea , &STMinWeight);
        newf =  max_rule_mean(newG, mst,  &CompTree  , STAltitudes , STArea, STMinWeight, i)  ;
        setweight(newG,i,newf);
    }
    
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STMinWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}



void dualhgbTruncated(graphe *g){
    
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, 0); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double * STMaxWeight = NULL;
    
    /// Set all edges to large values for later SM computation ////
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    
    for(i=mst->size-1;i>=0;i--){
        setweight(newG,i,INFINITO);
        computeQFZ_Attributes(newG, &CompTree, &STAltitudes, &STArea , &STMaxWeight);
        newf =  max_rule_trunc(newG, mst,  &CompTree  , STAltitudes , STArea, STMaxWeight, i)  ;
        setweight(newG,i,newf);
    }
    
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) );
    
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    // free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}



double max_rule_trunc(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double* STMaxWeight, int edgeIdx){
    
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    double   D=0, lambda, max_lambda=-9999;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    cx = x; px = cx;
    cy = y; py = cy;
    
    lambda = 0;
    
    do{
        lambda+=1;
        
        cx = getComponentAtLambda(CT,STAltitudes, lambda, cx );
        cy = getComponentAtLambda(CT,STAltitudes, lambda, cy );
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        if(D < Dif){
            D = Dif;// trunc D to always be more than the interdiference between regions
        }
        if( D >lambda ){
            max_lambda = lambda;
        }
        // px = CT->tabnodes[cx].father;
        // py = CT->tabnodes[cy].father;
    }
    //while(  px != py &&  CT->tabnodes[cx].father !=-1);
    while(D>(lambda ) );
    
    //    if(D > lambda){
    //       max_lambda = D;
    //    }
    if(max_lambda<0){max_lambda = 0;}
    
   // printf("Edge %i,%i:  %.3f \n",x,y,max_lambda);
    
    
    return max_lambda;
    //        px = CT->tabnodes[cx].father;
    //        py = CT->tabnodes[cy].father;
    //
    //
    //    }
    //    while(D>(lambda ) );
    //
    //    if(max_lambda<0){max_lambda = 0;}
    //
    //   // printf("Edge %i,%i:  %.3f \n",x,y,max_lambda);
    //
    //
    //    return max_lambda;
    
}
