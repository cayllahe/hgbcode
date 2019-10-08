/*
 Copyright ESIEE (2009)
 
 m.couprie@esiee.fr
 
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

/* authors : J. Cousty - L. Najman and M. Couprie */

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mccodimage.h>
#include <mcimage.h>
#include <mcweightgraph.h>

#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
//#define new_min(x,y) ((x) <= (y)) ? (x) : (y)


void usage(char *arg){
    printf("#################################################################\n\n");
    printf("USAGE: %s IN.ppm mode OUT.graph \n",arg);
    printf("When mode = 1 ->  4 connected,  Euclidean distance in Lab space is used as gradient\n");
    printf("When mode = 2 -> 8 connected,  Euclidean distance in Lab space is used as gradient\n");
    printf("When mode = 3 -> 4 connected,  Euclidean distance in RGB space is used as gradient\n");
    printf("When mode = 4 -> 8 connected, Euclidean distance in RGB space is used as gradient\n");
    printf("When mode = 5 -> 4 connected, gradient distance in Max \n");
    printf("When mode = 6 -> 8 connected, gradient distance in Max\n");
    
    printf("#################################################################\n\n");
}

int32_t max(int32_t a, int32_t b){
    if (a<b) return b;
    else return a;
}

double decompand(double V){
    if (V <= 0.04045) return V/12.92;
    else return pow( (V+0.055) / 1.055, 2.4);
}

#define KAPPA 903.3
#define EPSILON  0.008856

double simplified(double v, double ref){
    double vv;
    vv = v / ref;
    
    if(vv > EPSILON) vv = pow(vv, 1.0/3.0);
    else vv = (KAPPA * vv + 16) / 116;
    
    return vv;
}

void lab(uint8_t R, uint8_t G, uint8_t B, double *l, double* a, double *b){
    double dr, dg, db, x, y, z;
    
    dr = decompand( ((double) R)/255 );
    dg = decompand( ((double) G)/255 );
    db = decompand( ((double) B)/255 );
    
    x = 0.4124564*dr +  0.3575761*dg  + 0.1804375*db;
    y = 0.2126729*dr +  0.7151522*dg  + 0.0721750*db;
    z = 0.0193339*dr +  0.1191920*dg  + 0.9503041*db;
    
    x = simplified(x, 0.95047);
    y = simplified(y, 1.0);
    z = simplified(z, 1.08883);
    
    
    
    (*l) = 116 * y - 16;
    (*a) = 500 * (x - y);
    (*b) = 200 * (y - z);
}

double labgrad(uint8_t Rx, uint8_t Gx, uint8_t Bx, uint8_t Ry, uint8_t Gy, uint8_t By){
    double lx, ax, bx, ly, ay, by;
    lab(Rx, Gx, Bx, &lx, &ax, &bx);
    lab(Ry, Gy, By, &ly, &ay, &by);
    
    return sqrt( (lx-ly) * (lx-ly) + (ax-ay) * (ax-ay) + (bx-by)*(bx-by) );
}

//RGB
double euclidian(uint8_t Rx, uint8_t Gx, uint8_t Bx, uint8_t Ry, uint8_t Gy, uint8_t By){
    
    double value = sqrt( pow(Rx - Ry,2.) + pow(Gx - Gy ,2.) +  pow(Bx - By,2.) )   ;
     //double value =  pow(Rx - Ry,2.) + pow(Gx - Gy ,2.) +  pow(Bx - By,2.)    +  2;
    
    
    return value;
}
// Gray
double euclidianGray(uint8_t GrayX, uint8_t GrayY){
    double value = sqrt ( pow(GrayX - GrayY,2.) + 2 ) ;
    
    return value;
}

/** Case of a 4-connected graph where each edge is weighted by the
 absolute difference of intensity between its extremity
 pixels */
void LAB4connectedRGB(uint8_t* fR, uint8_t* fG, uint8_t*  fB, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(2*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}

/** Case of a 8-connected graph where each edge is weighted by the
 absolute difference of intensity between its extremity
 pixels */
void LAB8connectedRGB(uint8_t* fR, uint8_t* fG, uint8_t*  fB, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(4*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1  ){
                x = j * rs + i;
                y = (j+1) * rs + i +1;
                addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1 && j>0  ){
                x = j * rs + i;
                y = (j-1) * rs + i +1;
                addarete(g, x, y, labgrad(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
            
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}



/** Case of a 4-connected graph where each edge is weighted by the
 absolute difference of intensity between its extremity
 pixels */
void euclidian4connectedRGB(uint8_t* fR, uint8_t* fG, uint8_t*  fB, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(2*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}


/** Case of a 8-connected graph where each edge is weighted by the
 absolute difference of intensity between its extremity
 pixels */
void euclidian8connectedRGB(uint8_t* fR, uint8_t* fG, uint8_t*  fB, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(4*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1  ){
                x = j * rs + i;
                y = (j+1) * rs + i +1;
                addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
            if (  (i < rs - 1)  && j>0  ){
                x = j * rs + i;
                y = (j-1) * rs + i +1;
                addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
            }
            
            
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}



/** Case of a 4-connected graph in Gray image where each edge is weighted by the
 absolute difference of intensity between its extremity
 pixels */
void euclidian4connectedGray(uint8_t* fGray, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(2*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, euclidianGray(fGray[x], fGray[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, euclidianGray(fGray[x], fGray[y]) );
            }
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}

/** Case of a 8-connected graph where each edge is weighted by the
 absolute difference of intensity between its extremity
 pixels */
void euclidian8connectedGray(uint8_t* fGray, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(4*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, euclidianGray(fGray[x], fGray[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, euclidianGray(fGray[x], fGray[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1  ){
                x = j * rs + i;
                y = (j+1) * rs + i +1;
                addarete(g, x, y, euclidianGray(fGray[x], fGray[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1 && j>0  ){
                x = j * rs + i;
                y = (j-1) * rs + i +1;
                addarete(g, x, y, euclidianGray(fGray[x], fGray[y]) );
            }
            
            
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}


/** Case of a 4-connected graph where on gradients */
void gradient4connectedMax(uint8_t* fR, uint8_t* fG, uint8_t*  fB, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(2*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                //addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
                addarete(g, x, y,  new_max(fR[x], fR[y] ));
                
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                //addarete(g, x, y, euclidian(fR[x], fG[x], fB[x], fR[y], fG[y], fB[y]) );
                addarete(g, x, y,  new_max(fR[x], fR[y] ));
            }
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}
void gradient4connectedMaxGray(uint8_t* fGray, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(2*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, new_max(fGray[x], fGray[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, new_max(fGray[x], fGray[y]) );
            }
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}


void gradient8connectedMax(uint8_t* fR, uint8_t* fG, uint8_t*  fB, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(4*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, new_max( fR[x], fR[y] ) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, new_max( fR[x], fR[y] ) );
            }
            
            if (i < rs - 1 && j < cs - 1  ){
                x = j * rs + i;
                y = (j+1) * rs + i +1;
                addarete(g, x, y, new_max( fR[x], fR[y] ) );
            }
            
            if (  (i < rs - 1)  && j>0  ){
                x = j * rs + i;
                y = (j-1) * rs + i +1;
                addarete(g, x, y, new_max( fR[x], fR[y] ) );
            }
            
            
            
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}



void gradient8connectedMaxGray(uint8_t* fGray, int rs, int cs, char * file){
    
    graphe *g;
    int N, i , j, x, y;
    N = rs*cs;
    
    double *Fv;
    
    g = initgraphe(N, 2*(4*N-rs-cs));
    setSize(g,rs,cs);
    for (j = 0; j < cs; j++){
        for (i = 0; i < rs; i++){
            
            if (i < rs - 1){
                x = j * rs + i;
                y = j * rs + i + 1;
                addarete(g, x, y, new_max(fGray[x], fGray[y]) );
            }
            if (j < cs - 1){
                x = j * rs + i;
                y = (j+1) * rs + i;
                addarete(g, x, y, new_max(fGray[x], fGray[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1  ){
                x = j * rs + i;
                y = (j+1) * rs + i +1;
                addarete(g, x, y, new_max(fGray[x], fGray[y]) );
            }
            
            if (i < rs - 1 && j < cs - 1 && j>0  ){
                x = j * rs + i;
                y = (j-1) * rs + i +1;
                addarete(g, x, y, new_max(fGray[x], fGray[y]) );
            }
        }
    }
    
    Fv = calloc(N, sizeof(double));
    if(Fv == NULL){
        fprintf(stderr," cannot allocate (%lu bytes) memory\n" , N*sizeof(double)  );
        terminegraphe(g);
        exit(1);
    }
    /* The area of any pixel is 1 */
    
    for(i = 0; i < N; i++){
        Fv[i] = 1.0;
    }
    
    SaveGraphe(g, file, Fv);
    terminegraphe(g);
    free(Fv);
}


int checkFormat(char* filename){
    int BUFFERSIZE = 10000;
    char buffer[BUFFERSIZE];
    FILE *fd = NULL;
    char *read;
    
    fd = fopen(filename,"r");

    if (!fd)
    {
        fprintf(stderr, " file not found: %s\n", filename);
        exit(1);
    }
    
    read = fgets(buffer, BUFFERSIZE, fd); /* P5: raw int32_t bw  ; P2: ascii bw */
    /* P6: raw int32_t rgb ; P3: ascii rgb */
    if (!read)
    {
        fprintf(stderr, "IMage: fgets returned without reading\n");
        return 0;
    }
    if (buffer[0] != 'P')
    {   fprintf(stderr,"image : invalid image format\n");
        return 0;
    }
    int type = -1;
    
    
    switch (buffer[1])
    {
        case '3':
        case '6': type = 1; //flag as RGB
                  break;
        case '5':
                  type = 0;  //flag as gray image
                  break;
       //Edward
        case '2': type = 2;   // Gray image on different format
                  break;
        default:
            fprintf(stderr,"Invalid image format\n");
            exit(1);
    } /* switch */
    
    
    return type;
}



int32_t main(argc, argv)
/* =============================================================== */
int32_t argc; char **argv;
{
    struct xvimage *R, *G, *B, *Gray;
    
    int32_t rs, cs, mode;
    uint8_t *fR, *fG, *fB ,  *fGray;
    
    double l, b;
    
    if(argc != 4){
        usage(argv[0]);
        exit(1);
    }
    
    mode = atoi(argv[2]);
    int imagetype = checkFormat(argv[1]);
    if(  imagetype == 1 ){
        // RGB image
        readrgbimage(argv[1], &R, &G, &B);
        if(depth(R) != 1) {
            fprintf(stderr, "%s: 3D not yet implemented\n", argv[0]);
            exit(1);
        }
        rs = rowsize(R);
        cs = colsize(R);
        
        fR = UCHARDATA(R);
        fG = UCHARDATA(G);
        fB = UCHARDATA(B);
        Gray = 0;
        fGray = 0;
    }
    else{
        if (imagetype == 0) {
            // Gray image
            //readimage
            Gray = readimage(argv[1]);
            fGray = UCHARDATA(Gray);
            rs = rowsize(Gray);
            cs = colsize(Gray);
            fR = fG = fB = 0;
            R = G = B= 0;
        }
        if (imagetype == 2) {
            // Gray image different format
            Gray = readimage2(argv[1]);
            fGray = UCHARDATA(Gray);
            rs = rowsize(Gray);
            cs = colsize(Gray);
            fR = fG = fB = 0;
            R = G = B= 0;
        }
        else{
            printf("format not supported by image2graph\n" );
            exit(1);
        }
    }
    switch(mode){
        case 1:
            LAB4connectedRGB( fR, fG, fB , rs, cs, argv[3]);
            break;
        case 2:
            LAB8connectedRGB( fR, fG, fB , rs, cs, argv[3]);
            break;
            
        case 3:
            if ( imagetype == 1)
                euclidian4connectedRGB( fR, fG, fB , rs, cs, argv[3]);
            else
                euclidian4connectedGray(fGray,rs, cs, argv[3] );
            break;
        case 4:
            if ( imagetype == 1)
                euclidian8connectedRGB( fR, fG, fB , rs, cs, argv[3]);
            else
               euclidian8connectedGray(fGray,rs, cs, argv[3] );
            break;
        
        case 5:
            if ( imagetype == 1)
                gradient4connectedMax( fR, fG, fB , rs, cs, argv[3]);
            else
                gradient4connectedMaxGray(fGray,rs, cs, argv[3] );
            break;
        case 6:
            if ( imagetype == 1)
                gradient8connectedMax( fR, fG, fB , rs, cs, argv[3]);
            else
                gradient8connectedMaxGray(fGray,rs, cs, argv[3] );
            break;
            
        default :
            fprintf(stderr, "%s: mode %d not available\n",argv[0], mode);
            exit(0);
    }
    if(R)freeimage(R);
    if(G)freeimage(G);
    if(B)freeimage(B);
    if(Gray)freeimage(Gray);
        
    return 0;
}

