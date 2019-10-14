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

#include <stdio.h>
#include <dealpng.h>
#include <math.h>

#include <mcweightgraph.h>

void abort_(const char * s, ...)
{
    va_list args;
    va_start(args, s);
    vfprintf(stderr, s, args);
    fprintf(stderr, "\n");
    va_end(args);
    abort();
}


void read_png_file(PNG_container* pngfile, char* file_name)
{
    char header[8];    // 8 is the maximum size that can be checked
    
    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
        abort_("[read_png_file] File %s could not be opened for reading", file_name);
    fread(header, 1, 8, fp);
    if (png_sig_cmp(header, 0, 8))
        abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);
    
    
    /* initialize stuff */
    pngfile->png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!pngfile->png_ptr)
        abort_("[read_png_file] png_create_read_struct failed");
    
    pngfile->info_ptr = png_create_info_struct(pngfile->png_ptr);
    if (!pngfile->info_ptr)
        abort_("[read_png_file] png_create_info_struct failed");
    
    if (setjmp(png_jmpbuf(pngfile->png_ptr)))
        abort_("[read_png_file] Error during init_io");
    
    png_init_io(pngfile->png_ptr, fp);
    png_set_sig_bytes(pngfile->png_ptr, 8);
    
    png_read_info(pngfile->png_ptr, pngfile->info_ptr);
    
    pngfile->width = png_get_image_width(pngfile->png_ptr, pngfile->info_ptr);
    pngfile->height = png_get_image_height(pngfile->png_ptr, pngfile->info_ptr);
    pngfile->color_type = png_get_color_type(pngfile->png_ptr, pngfile->info_ptr);
    pngfile->bit_depth = png_get_bit_depth(pngfile->png_ptr, pngfile->info_ptr);
    
    pngfile->number_of_passes = png_set_interlace_handling(pngfile->png_ptr);
    
    
    if(pngfile->bit_depth == 16)
        png_set_strip_16(pngfile->png_ptr);
    
    if(pngfile->color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(pngfile->png_ptr);
    
    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    if(pngfile->color_type == PNG_COLOR_TYPE_GRAY && pngfile->bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(pngfile->png_ptr);
    
    
   
    
    
    if(pngfile->color_type == PNG_COLOR_TYPE_GRAY ||
       pngfile->color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(pngfile->png_ptr);
    
    
    png_read_update_info(pngfile->png_ptr, pngfile->info_ptr);
    
    
    /* read file */
    if (setjmp(png_jmpbuf(pngfile->png_ptr)))
        abort_("[read_png_file] Error during read_image");
    int y;
    pngfile->row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * pngfile->height);
    for (y=0; y<pngfile->height; y++)
        pngfile->row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(pngfile->png_ptr,pngfile->info_ptr));
    
    png_read_image(pngfile->png_ptr, pngfile->row_pointers);
    
    fclose(fp);
}




void write_png_file(PNG_container* pngfile,char* file_name)
{
    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
        abort_("[write_png_file] File %s could not be opened for writing", file_name);
    
    
    /* initialize stuff */
    pngfile->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!pngfile->png_ptr)
        abort_("[write_png_file] png_create_write_struct failed");
    
    pngfile->info_ptr = png_create_info_struct(pngfile->png_ptr);
    if (!pngfile->info_ptr)
        abort_("[write_png_file] png_create_info_struct failed");
    
    pngfile->end_info = png_create_info_struct(pngfile->png_ptr);  
    if (!pngfile->end_info) {
        png_destroy_read_struct(&pngfile->png_ptr, &pngfile->info_ptr, NULL);
        fclose(fp);
        abort_("[write_png_file] png_create_info_struct end failed");
    }
    
    
    if (setjmp(png_jmpbuf(pngfile->png_ptr)))
        abort_("[write_png_file] Error during init_io");
    
    png_init_io(pngfile->png_ptr, fp);
    
    
    /* write header */
    if (setjmp(png_jmpbuf(pngfile->png_ptr)))
        abort_("[write_png_file] Error during writing header");
    
    
    pngfile->color_type = PNG_COLOR_TYPE_RGB;
    
    
    png_set_IHDR(pngfile->png_ptr, pngfile->info_ptr, pngfile->width, pngfile->height,
                 pngfile->bit_depth, pngfile->color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
    png_write_info(pngfile->png_ptr, pngfile->info_ptr);
    
    
    /* write bytes */
    if (setjmp(png_jmpbuf(pngfile->png_ptr)))
        abort_("[write_png_file] Error during writing bytes");
    
    png_write_image(pngfile->png_ptr, pngfile->row_pointers);
    
    
    /* end write */
    if (setjmp(png_jmpbuf(pngfile->png_ptr)))
        abort_("[write_png_file] Error during end of write");
    
    png_write_end(pngfile->png_ptr, NULL);
    int y;
    /* cleanup heap allocation */
    for (y=0; y<pngfile->height; y++)
        free(pngfile->row_pointers[y]);
    free(pngfile->row_pointers);
    
    fclose(fp);
}






void addEdgeRGB(graphe *g, int i , int j,png_byte* iPtr, png_byte* jPtr ){

    // Computing the Euclidian distance betweeen RGB values
    
    
  //  double value =  sqrt( pow(iPtr[0]-jPtr[0],2) + pow(iPtr[1]-jPtr[1],2) +  pow(iPtr[2]-jPtr[2],2) );
    
    // L1 -norm
     //   double value =  abs(iPtr[0]-jPtr[0])  + abs(iPtr[1]-jPtr[1]) +  abs(iPtr[2]-jPtr[2]) ;
    //INTEGERS
    double value =  pow(iPtr[0]-jPtr[0],2.) + pow(iPtr[1]-jPtr[1],2.) +  pow(iPtr[2]-jPtr[2],2.) + 2 ;
    
    addarete(g, i, j, value);

}

void createEdgeFromPNG_4_Connectivity(graphe *g, PNG_container* pngfile,  int x, int y , int rs, int cs){
    int i, j;
    i = x * cs + y;
    
    
    png_byte* CentralRow = pngfile->row_pointers[x];
    png_byte* centralPtr = &(CentralRow[y*3]);
    
    png_byte* tmpRow, *tmpPtr;
    
    if (x-1 >= 0){
        j = (x-1) * cs + y;
        tmpRow = pngfile->row_pointers[x-1];
        tmpPtr = &(tmpRow[y*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
        
    }
    if(y+1 <cs){
        j = (x) * cs + (y+1);
        tmpRow = pngfile->row_pointers[x];
        tmpPtr = &(tmpRow[(y+1)*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
    }
    
    /*
    if(y-1 >= 0){
        j = (x) * cs + (y-1);
        tmpRow = pngfile->row_pointers[x];
        tmpPtr = &(tmpRow[(y-1)*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
    }
     if (x+1 < rs){
     j = (x+1) * cs + y;
     tmpRow = pngfile->row_pointers[x+1];
     tmpPtr = &(tmpRow[y*3]);
     addEdgeRGB(g, i,j , centralPtr, tmpPtr );
     }
    */
    
    /*
     x-1, y-1
     x-1, y
     x-1, y+1
     
     x, y-1
     x, y+1
     
     x+1, y-1
     x+1,y
     x+1,y+1
     */


}


void createEdgeFromPNG_8_Connectivity(graphe *g, PNG_container* pngfile,  int x, int y , int rs, int cs){
    int i, j;
    i = x * cs + y;
    
    
    png_byte* CentralRow = pngfile->row_pointers[x];
    png_byte* centralPtr = &(CentralRow[y*3]);
    
    png_byte* tmpRow, *tmpPtr;
    
    if (x-1 >= 0){
        j = (x-1) * cs + y;
        tmpRow = pngfile->row_pointers[x-1];
        tmpPtr = &(tmpRow[y*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
        
    }
    if(y+1 <cs){
        j = (x) * cs + (y+1);
        tmpRow = pngfile->row_pointers[x];
        tmpPtr = &(tmpRow[(y+1)*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
    }
    
    if( (x-1)>=0 &&  (y+1)<cs){
        j = (x-1) * cs + (y+1);
        tmpRow = pngfile->row_pointers[x-1];
        tmpPtr = &(tmpRow[(y+1)*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
    }
    if( (x+1)<rs &&  (y+1)<cs){
        j = (x+1) * cs + (y+1);
        tmpRow = pngfile->row_pointers[x+1];
        tmpPtr = &(tmpRow[(y+1)*3]);
        addEdgeRGB(g, i,j , centralPtr, tmpPtr );
    }
    
    /*
     if(y-1 >= 0){
     j = (x) * cs + (y-1);
     tmpRow = pngfile->row_pointers[x];
     tmpPtr = &(tmpRow[(y-1)*3]);
     addEdgeRGB(g, i,j , centralPtr, tmpPtr );
     }
     if (x+1 < rs){
     j = (x+1) * cs + y;
     tmpRow = pngfile->row_pointers[x+1];
     tmpPtr = &(tmpRow[y*3]);
     addEdgeRGB(g, i,j , centralPtr, tmpPtr );
     }
     */
    
    /*
     x-1, y-1
     x-1, y
     x-1, y+1
     
     x, y-1
     x, y+1
     
     x+1, y-1
     x+1,y
     x+1,y+1
     */
    
    
}


void process_file(PNG_container* pngfile, char *outputfile)
{        
    graphe *g;
    int cs = pngfile->width;
    int rs = pngfile->height;
    int N = rs*cs;
    double *Fv;
    
    g = initgraphe(N, 2*(4*N-rs-cs));  // 4 connectivity
    
    //setSize(g,rs,cs);
    setSize(g,cs,rs);
    
    
    int x,y;
    double grayValue;
    
    int i, j; // i->j edge
    
    for (x=0; x<pngfile->height; x++){
        // A la fila
        png_byte* row = pngfile->row_pointers[x];
        for (y=0; y<pngfile->width; y++) {
            // por columna
            png_byte* ptr = &(row[y*3]);
            
            //printf("Pixel at position [ %d - %d ] has RGBA values: %d - %d - %d - %d\n",
            //       x, y, ptr[0], ptr[1], ptr[2], ptr[3]);
            
           // grayValue = (ptr[0] + ptr[1] + ptr[2])/3;
            //printf("Pixel at position [ %d - %d ] has RGB values: %d - %d - %d \n", x, y, ptr[0], ptr[1], ptr[2]);
            
            createEdgeFromPNG_8_Connectivity(g,pngfile, x,y,rs,cs);
                       
        }
    }
    
    Fv = calloc(N, sizeof(double));
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    //SaveGraphe(g, "/Users/edus/Documents/MorphoGraph/SM_Edward/Graphs/pngGraph.graph", Fv);
    SaveGraphe(g, outputfile, Fv);
    terminegraphe(g);
    free(Fv);
    
    
}

