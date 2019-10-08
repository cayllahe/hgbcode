//
//  MST.h
//  SM_Edward
//
//  Created by Edward Cayllahua on 8/22/16.
//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

#ifndef MST_h
#define MST_h

#include <stdlib.h>
#include <sys/types.h>
#include <stdio.h>

#define INFINITO 9999999999
#define MSTConstant 0

typedef struct{
    int x, y;
    int idx;
    double weight;
}edge;


typedef edge* edges;

typedef struct{
    edges MSTedges;
    int size;
    
}MST;



void setSizeMST(MST* mst, int n);

#endif /* MST_h */
