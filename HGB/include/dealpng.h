//
//  dealpng.h
//  SM_Edward
//
//  Created by Edward Cayllahua on 8/30/16.
//  Copyright © 2016 Edward Cayllahua. All rights reserved.
//

#ifndef dealpng_h
#define dealpng_h


#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#define PNG_DEBUG 3
#include <png.h>



typedef struct{

int width, height;
png_byte color_type;
png_byte bit_depth;

png_structp png_ptr;
png_infop info_ptr;
png_infop end_info;
int number_of_passes;
png_bytep * row_pointers;
    
}PNG_container;



void read_png_file(PNG_container* pngfile, char* file_name);
void abort_(const char * s, ...);
void write_png_file(PNG_container* pngfile,char* file_name);
void process_file(PNG_container* pngfile,char *outputfile);

#endif /* dealpng_h */
