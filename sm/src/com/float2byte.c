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
/*! \file float2byte.c

\brief converts a "float" image to a "byte" image

<B>Usage:</B> float2byte in.pgm [mode] out.pgm

<B>Description:</B>

Depending on the value given for the (optional) parameter <B>mode</B>:
\li   <B>mode</B> = 0 (default) : for all x, out[x] = min(255, arrondi(in[x])).
\li   <B>mode</B> = 1 : for all x, out[x] = arrondi(in[x]) modulo 256.
\li   <B>mode</B> = 2 : scales values in the range 0-255.
\li   <B>mode</B> = 4 : truncation of the square root in the range 0-255.
\li   <B>mode</B> = 5 : truncation of the log in the range 0-255.

<B>Types supported:</B> float 2d, float 3d

<B>Category:</B> convert
\ingroup  convert

\author Michel Couprie
*/

/*
%TEST float2byte %IMAGES/2dfloat/f2fish1.pgm 0 %RESULTS/float2byte_f2fish1_0.pgm
%TEST float2byte %IMAGES/2dfloat/f2fish1.pgm 1 %RESULTS/float2byte_f2fish1_1.pgm
%TEST float2byte %IMAGES/2dfloat/f2fish1.pgm 2 %RESULTS/float2byte_f2fish1_2.pgm
%TEST float2byte %IMAGES/2dfloat/f2fish1.pgm 4 %RESULTS/float2byte_f2fish1_4.pgm
%TEST float2byte %IMAGES/2dfloat/f2fish1.pgm 5 %RESULTS/float2byte_f2fish1_5.pgm
*/

/*
   Michel Couprie - mai 1998

   Modif : decembre 1999 - mode 3 (trunchisto)
 */

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mcimage.h>
#include <mccodimage.h>
#include <mcutil.h>
#include <lhisto.h>
#ifdef HP
#define _INCLUDE_XOPEN_SOURCE
#endif
#include <math.h>

/* =============================================================== */
int main(int argc, char **argv)
/* =============================================================== */
{
  struct xvimage * imagefloat;
  struct xvimage * imagebyte;
  float *L;
  uint8_t *B;
  int32_t tmp, mode = 0;
  float Min, Max, t;
  double T;
  index_t x, rs, cs, ds, N;

  if ((argc < 3) || (argc > 4))
  {
    fprintf(stderr, "usage: %s in1.pgm [mode] out.pgm \n", argv[0]);
    fprintf(stderr, "mode = 0 (trunc) | 1 (modulo) | 2 (scale) | \n");
    fprintf(stderr, "       4 (square root) | 5 (log)\n");
    exit(1);
  }

  imagefloat = readimage(argv[1]);  
  if (imagefloat == NULL)
  {
    fprintf(stderr, "%s: readimage failed\n", argv[0]);
    exit(1);
  }

  if (datatype(imagefloat) != VFF_TYP_FLOAT)
  {
    fprintf(stderr, "%s: image type must be float\n", argv[0]);
    exit(1);
  }

  if (argc > 3) mode = atoi(argv[2]);

  rs = rowsize(imagefloat);
  cs = colsize(imagefloat);
  ds = depth(imagefloat);
  N = rs * cs * ds;
  L = FLOATDATA(imagefloat);
  
  imagebyte = allocimage(NULL, rs, cs, ds, VFF_TYP_1_BYTE);
  if (imagebyte == NULL)
  {
    fprintf(stderr, "%s: allocimage failed\n", argv[0]);
    exit(1);
  }
  B = UCHARDATA(imagebyte);
  imagebyte->xdim = imagefloat->xdim;
  imagebyte->ydim = imagefloat->ydim;
  imagebyte->zdim = imagefloat->zdim;

  switch(mode)
  {
    case 0:
      for (x = 0; x < N; x++)
      {
        tmp = arrondi(L[x]);
        B[x] = (uint8_t)mcmin(tmp,255);
      }
      break;
    case 1:
      for (x = 0; x < N; x++)
      {
        tmp = arrondi(L[x]);
        B[x] = (uint8_t)(tmp % 256);
      }
      break;
    case 2:
      Min = Max = L[0];
      for (x = 0; x < N; x++) 
        if (L[x] > Max) Max = L[x]; else if (L[x] < Min) Min = L[x];
      for (x = 0; x < N; x++) 
      {
        t = ((L[x]-Min) * 255.0) / (float)(Max-Min);
        tmp = arrondi(t);
        B[x] = (uint8_t)mcmin(255,tmp);
      }
      break;
    case 4:
      for (x = 0; x < N; x++)
      {
        T = sqrt((double)(L[x]));
        tmp = arrondi(T);
        tmp = mcmin(255,tmp);
        tmp = mcmax(0,tmp);
        B[x] = (uint8_t)tmp;
      }
      break;
    case 5:
      for (x = 0; x < N; x++)
      {
        T = log((double)(L[x]));
        tmp = arrondi(T);
        tmp = mcmin(255,tmp);
        tmp = mcmax(0,tmp);
        B[x] = (uint8_t)tmp;
      }
      break;
    default:
      fprintf(stderr, "usage: %s in1.pgm [mode] out.pgm \n", argv[0]);
      fprintf(stderr, "mode = 0 (trunc) | 1 (modulo) | 2 (scale) | \n");
      fprintf(stderr, "       4 (square root) | 5 (log)\n");
      exit(1);
  } /* switch(mode) */

  writeimage(imagebyte, argv[argc-1]);
  freeimage(imagefloat);
  freeimage(imagebyte);

  return 0;
} /* main */
