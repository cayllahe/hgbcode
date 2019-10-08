//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

#include <stdio.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
//#include <mcweightgraph.h>
//#include <lhierarchie.h>
#include <dealpng.h>
#include <graphSegmentation.h>
#include <time.h>

void usage(char *arg){
    printf("#################################################################\n\n");
    printf("USAGE: %s IN.png  OUT.graph  \n",arg);
    printf("#################################################################\n\n");
}

char *get_filename_ext( char *filename) {
    char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}


void loadGraph(char* argv,  graphe **g, double **Av ){
// Load a graph from the file argv
// if file is png, it saves it as a graph then loads it into g
// if its .graph it loads it directly into g
    
    char * extension1, *png="png" ;//, tempFile[200];
    int fd;
    char tempFile[] = "/tmp/fileGraphXXXXXXX"; //creates temp file with this format
    extension1 = get_filename_ext(argv);
    int y;
    if(strcmp( extension1, png) == 0 ){
        //Load png into a graph and saves it to tempFile
        fd = mkstemp(tempFile);
       // strcat(tempFile, argv); // Temporary file to save graph
       // strcat(tempFile, ".graph");
        PNG_container* pngfile;
        pngfile = (PNG_container *)calloc(1, sizeof(PNG_container));
        read_png_file(pngfile, argv);
        process_file(pngfile,tempFile );
        
        *g = ReadGraphe(tempFile, Av);
        
        for (y=0; y<pngfile->height; y++)
            free(pngfile->row_pointers[y]);
            
        free(pngfile->row_pointers);
        png_destroy_read_struct(&pngfile->png_ptr, &pngfile->info_ptr, &pngfile->end_info);
        
        
        
        free(pngfile);
        
    }
    else{
        // I assume  is a .graph file format
        *g = ReadGraphe(argv , Av);
    }
    

}

int main(int argc,char ** argv){
    
    clock_t startAll, endALL; double cpu_time_usedALL;
    
    graphe *g;
    double *Av;
 
    if(argc!=6){
        usage(argv[0]);
        exit(1);
    }
    
    loadGraph(argv[1], &g , &Av);
    
    startAll = clock();
    
    switch (atoi(argv[3]) ) {
            
        case 1:
            //naive minimization
            // Non incremental QFZ
            naiveSegmentation(g);
            break;
            
        case 2:
            //range minimization
            // Non incremental QFZ
            segmentationRange(g);
            break;
        case 3:
            // branch minimization
            // Non incremental QFZ
            segmentationBranch(g);
            break;
        case 4:
            // efficient scale
            // incremental QFZ
           // incrementalSegmentationBranch(g);
             incrementalSegmentationBranchTimed(g);
            break;
        case 5:
            //naive minimization
            // incremental QFZ
            incrementalNaiveSegmentation(g);
            break;
        case 6:
            // range minimization
            // incremental QFZ
            incrementalSegmentationRange(g);
            break;
        default:
            printf("Bad argument: type of segmentation \n");
            terminegraphe(g);
            free(Av);
            return 1;            
    }
    
    
    
    endALL = clock();
    
     SaveGraphe(g, argv[2],Av) ;
    
   // cpu_time_usedALL = ((double) (endALL - startAll)) / CLOCKS_PER_SEC;
   // printf("%i \t  %s \t %s \t  took \t %f \t seconds \t levels \t %i \n",  atoi(argv[3]), argv[4] , argv[5] , cpu_time_usedALL);
    
   // double sum;
   // int levels = computeLevels(g, &sum);
  //  printf("%i \t  %s \t %s \t  took \t %f \t seconds \t levels \t %i \t sum \t %f\n",  atoi(argv[3]), argv[4] , argv[5] , cpu_time_usedALL, levels, sum );
    
   // printf("%i \t  %s \t %s \t  took \t %f \t seconds \t levels \t %i \t sum \t %f\n",  atoi(argv[3]), argv[4] , argv[5] , cpu_time_usedALL );
    
    
    
   
    terminegraphe(g);
    free(Av);
    return 1;
}
