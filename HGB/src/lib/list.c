//
//  list.c
//  SM_Edward
//
//  Created by Edward Cayllahua on 5/26/17.
//  Copyright © 2017 Edward Cayllahua. All rights reserved.
//

#include <list.h>


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#define INFINITO 9999999999

/**
    ListInit and addValue are used for range minimization,
    ListInit initializes with a first node with infinite value (the first value in  the range of values.)
    addValue  adds a node with value val before a node n
 */

void addValue(NodeList** lista, NodeList* insertBefore, double val){
    
    if(insertBefore == (*lista) ){
        NodeList* node = (NodeList*)malloc(sizeof(NodeList));
        node->value = val;
        insertBefore->last = node;
        node->last = NULL;
        node->next = insertBefore;
        *lista = node;
    }
    else{
        NodeList* node = (NodeList*)malloc(sizeof(NodeList));
        node->value = val;
        insertBefore->last->next = node;
        node->last = insertBefore->last;
        node->next = insertBefore;
        insertBefore->last = node;        
    }
}



void ListInit(NodeList** lista ){
    
    NodeList* node1 = (NodeList*)malloc(sizeof(NodeList));
    
    node1->value = INFINITO;
    node1->last = NULL;
    node1->next = NULL;
        
    (* lista) = node1;
    
}


/**
    ListInitNormal This initializes a regular list with a first node set to INFINITY
    addValueNormal add a value to a regular list
 */

void ListInitNormal(NodeList** lista ){
    
    NodeList* node1 = (NodeList*)malloc(sizeof(NodeList));
    
    node1->value = INFINITO;
    node1->last = NULL;
    node1->next = NULL;
    (* lista) = node1;
}

void addValueNormal(NodeList** lista, NodeList** lastNode, double val){
    
   // printf("ver %i - %i\n", *lastNode, *lista);
    if( (*lista) == NULL){
        //(*lista)->value = val; // only updates for first value of the list
        ListInitNormal(lista);
        (*lista)->value = val;
        *lastNode = *lista;
        return; 
    }
    
    //if( (*lastNode)->value == val){
    //    return; // does not allow same values on last node
    //}
    
    NodeList* node = (NodeList*)malloc(sizeof(NodeList));
    node->value = val;
    
    node->last = (*lastNode);
    node->next = NULL;
    node->last->next = node;
    *lastNode = node;
}


void addValueDis(NodeList** lista, NodeList** lastNode, double val, double dis){
    
    // printf("ver %i - %i\n", *lastNode, *lista);
    if( (*lista) == NULL){
        //(*lista)->value = val; // only updates for first value of the list
        ListInitNormal(lista);
        (*lista)->value = val;
        (*lista)->dis = dis;
        (*lista)->bound =0;
        *lastNode = *lista;
        return;
    }
    
    //if( (*lastNode)->value == val){
    //    return; // does not allow same values on last node
    //}
    
    NodeList* node = (NodeList*)malloc(sizeof(NodeList));
    node->value = val;
    node->dis = dis;
    node->bound = 0;
    node->last = (*lastNode);
    node->next = NULL;
    node->last->next = node;
    *lastNode = node;
    
}

double listMeanValue(NodeList* lista){
    double sum = 0;
    int count = 0;
    NodeList* node = lista;
    while(node!=NULL){
        sum = sum + node->value;
        count++;
        node=node->next;
    }
    
    return sum/count;
    
}

double listMedianValue(NodeList* lista){
    
    double count = 0;
    int i;
    NodeList* node = lista;
    while(node!=NULL){
        count++;
        node=node->next;
    }
    node = lista;
    count = ceil(count/2.0);
    for(i = 1 ; i<count; i++){
        node= node->next;
    }
    
    return node->value;
    
}

void cleanList(NodeList** list){
    
    NodeList* tmp;
    
    while(*list!= NULL){
        tmp = *list;
        (*list) = (*list)->next;
        free(tmp);
    }
    
    
}
