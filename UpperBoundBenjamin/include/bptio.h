#ifndef __BPTIO_H__
#define __BPTIO_H__

#include "base.h"
//#include <mccodimage.h>
//#include <mcimage.h>
#include <mcweightgraph.h>
//#include <lhierarchie.h>
#include <lcomptree.h>

#define BPTIO_NBNODES_KEY "NBNODES"
#define BPTIO_NBATTRIBUTES_KEY "NBATTR"
#define BPTIO_HEADEREND_KEY "END"
#define BPTIO_NAME_KEY "NAME"

bool saveBPT(char * path, int nbnodes, int * parents, int numAttr, double ** attrs, char ** attrNames); 
bool readBPT(char * path, int * nbnodes, int ** parents, int * numAttr, double *** attrs, char *** attrNames); 
int * createParentArray(JCctree * bpt);


#endif 
