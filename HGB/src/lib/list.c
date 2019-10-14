/*
Copyright ESIEE (2019) 

edward.cayllahua@esiee.fr

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

/* authors : Edward Cayllahua - J. Cousty and Silvio Guimarães */
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
