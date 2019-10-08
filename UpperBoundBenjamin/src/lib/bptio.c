#include "bptio.h"
#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>

bool starts_with(const char *string, const char *prefix)
{
    while(*prefix)
    {
        if(*prefix++ != *string++)
            return 0;
    }

    return 1;
}

// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char * trimws(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

bool readCharKeyVal(char * line, char ** v)
{
    char * pe = strchr(line, '=');
    if (pe==NULL)
    {
        return false;
    }
    ++pe;
    pe= trimws(pe);
    char * res = myCalloc(sizeof(char), strlen(pe)+1);
    strcpy(res, pe);
    *v = res;
    return true;
}

bool readIntKeyVal(char * line, int * v)
{
    char * pe = strchr(line, '=');
    if (pe==NULL)
    {
        return false;
    }
    ++pe;

    char * stop;
    int vt = strtol(pe, &stop, 10);
    if(*pe=='\0' || stop==pe)
    {
        return false;
    }
    *v=vt;
    
    return true;
}

int * createParentArray(JCctree * bpt)
{
    int * parent = myCalloc(sizeof(int), bpt->nbnodes);
    int i;
    JCctreenode *tabnodes = bpt->tabnodes;
    for(i=0; i<bpt->nbnodes; i++)
    {
        parent[i] = tabnodes[i].father;
    }
    return parent;
}


bool saveBPT(char * path, int nbnodes, int * parents, int numAttr, double ** attrs, char ** attrNames)
{
    int nbn = nbnodes;
    FILE * f = fopen(path, "w");
    if(f==NULL)
    {
        perror("saveBPT: ");
        return false;
    }
    fprintf(f, "%s=%d\n",BPTIO_NBNODES_KEY, nbn);
    fprintf(f, "%s=%d\n",BPTIO_NBATTRIBUTES_KEY, numAttr);
    fprintf(f, "%s\n",BPTIO_HEADEREND_KEY);

    
    int nbw = fwrite(parents, sizeof(int), nbn, f);
    if(nbw != nbn)
    {
        fprintf(stderr, "saveBPT: error writting file %s\n", path);
        fclose(f);
    }
    
    int i;
    for(i=0; i<numAttr; ++i)
    {
        fprintf(f, "%s=%s\n",BPTIO_NAME_KEY, attrNames[i]);
        fprintf(f, "%s\n",BPTIO_HEADEREND_KEY);
        nbw=fwrite(attrs[i], sizeof(double), nbn, f);
        if(nbw != nbn)
        {
            fprintf(stderr, "saveBPT: error writting file %s\n", path);
            fclose(f);
            return false;
        }
    }
    fclose(f);
    return true;
}

#define BUFSIZE 4096

/**
    lots of memory leaks in case of errors...
*/
bool readBPT(char * path, int * nbnodes, int ** parent, int * numAttr, double *** attrs, char *** attrNames)
{
    FILE * f = fopen(path, "r");
    if(f==NULL)
    {
        perror("loadBPT: ");
        return false;
    }
    char buf[BUFSIZE] ; 

    bool flag=true;
    int nbNodes=-1;
    int nbAttr=0;
    do{
        char * res = fgets(buf, BUFSIZE, f);
        //printf("Line %s\n",res);
        if(res==NULL)
        {
            fprintf(stderr, "loadBPT: error reading file %s\n", path);
            fclose(f);
            return false;
        }
        char * b = trimws(buf);
        //printf("Trimmed %s\n",b);
        if(strcmp(b,BPTIO_HEADEREND_KEY)==0)
        {
            flag=false;
        }else if(starts_with(b, BPTIO_NBNODES_KEY))
        {
            if(!readIntKeyVal(b, &nbNodes))
            {
                fprintf(stderr, "loadBPT: error reading file %s\n", path);
                fprintf(stderr, "\t invalid key/value on line %s\n", res);
                fclose(f);
                return false;
            }
        }else if(starts_with(b, BPTIO_NBATTRIBUTES_KEY))
        {
            if(!readIntKeyVal(b, &nbAttr))
            {
                fprintf(stderr, "loadBPT: error reading file %s\n", path);
                fprintf(stderr, "\t invalid key/value on line %s\n", res);
                fclose(f);
                return false;
            }
        }else{
            
        }
    
    }while(flag);

    if(nbNodes<=0)
    {
        fprintf(stderr, "loadBPT: error reading file %s\n", path);
        fprintf(stderr, "\t invalid or missing mandatory node numbers (key: %s)\n",BPTIO_NBNODES_KEY);
        fclose(f);
        return false;

    }

    if(nbAttr<=0)
    {
        fprintf(stderr, "loadBPT: error reading file %s\n", path);
        fprintf(stderr, "\t  attributes number must be positive (key: %s)\n",BPTIO_NBATTRIBUTES_KEY);
        fclose(f);
        return false;

    }

    int * parentt = myCalloc(nbNodes, sizeof(int));
    int nbr = fread(parentt, sizeof(int), nbNodes, f);
    if(nbr!=nbNodes)
    {
        fprintf(stderr, "loadBPT: error reading file %s\n", path);
        fprintf(stderr, "\t Error reading parent array\n");
        free(parentt);
        fclose(f);
        return false;
    }
    
    
    if(nbAttr >0)
    {
        double ** attrt = myCalloc(nbAttr, sizeof(double *));
        char ** names = myCalloc(nbAttr ,sizeof(char*));
        int i,j;
        for(i=0;i<nbAttr;++i)
        {

            
            flag=true;
            do{
                char * res = fgets(buf, BUFSIZE, f);
                if(res==NULL)
                {
                    fprintf(stderr, "loadBPT: error reading file %s\n", path);
                    fclose(f);
                    return false;
                }
                char * b = trimws(buf);
        
                if(strcmp(b,BPTIO_HEADEREND_KEY)==0)
                {
                    flag=false;
                }else if(starts_with(b, BPTIO_NAME_KEY))
                {
                    if(!readCharKeyVal(b, &(names[i])))
                    {
                        fprintf(stderr, "loadBPT: error reading file %s\n", path);
                        fprintf(stderr, "\t invalid key/value on line %s\n", res);
                        fclose(f);
                        return false;
                    }
                }
                
            
            }while(flag);

            attrt[i] = myCalloc(nbNodes,sizeof(double));
            nbr = fread(attrt[i], sizeof(double), nbNodes, f);
            if(nbr!=nbNodes)
            {
                fprintf(stderr, "loadBPT: error reading file %s\n", path);
                fprintf(stderr, "\t Error reading attribute array %d, %s\n", i+1, names[i]);
                fclose(f);
                return false;
            }

        }
       
        
         *attrs = attrt;
        *attrNames = names;
    }

    long curPos = ftell(f);
    fseek(f, 0L, SEEK_END);
    long endPos = ftell(f);
    if(curPos!=endPos)
    {
        fprintf(stderr, "loadBPT: warning, end of file not reach but all declared data were read! %s\n", path);
        fprintf(stderr,"%d byte left\n",endPos-curPos);
    }
    *parent=parentt;
    *numAttr=nbAttr;
   
    *nbnodes = nbNodes;
    return true;
}


