//
//  MST.c
//  SM_Edward
//
//  Created by Edward Cayllahua on 8/22/16.
//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

#include<MST.h>

void setSizeMST(MST* mst, int n){
    
    
    if(( mst->MSTedges = (edge*)malloc(sizeof(edge) * n)) == NULL){
        fprintf(stderr,"setSizeMST: allocation error\n");
        exit(0);
    }
    mst->size=0;

}
