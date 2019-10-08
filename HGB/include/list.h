//
//  list.h
//  SM_Edward
//
//  Created by Edward Cayllahua on 5/26/17.
//  Copyright © 2017 Edward Cayllahua. All rights reserved.
//

#ifndef list_h
#define list_h

#include <stdio.h>


typedef struct Node
{
    double value;
    double dis;
    int bound;
    struct Node *next;
    struct Node *last;
} NodeList;

void addValue(NodeList** lista, NodeList* before, double val);
void ListInit(NodeList** lista );

void ListInitNormal(NodeList** lista );
void addValueNormal(NodeList** lista, NodeList** lastNode, double val);
void addValueDis(NodeList** lista, NodeList** lastNode, double val, double dis);

double listMeanValue(NodeList* lista);
double listMedianValue(NodeList* lista);
void cleanList(NodeList** list);


#endif /* list_h */
