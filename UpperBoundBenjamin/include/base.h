#ifndef __BASE_H__
#define __BASE_H__

#include <sys/types.h>
#include <stdint.h>

typedef uint8_t uchar;

#ifndef __cplusplus
    typedef char bool;
#endif

#define true (1)
#define false (0)

int mini(int a, int b);

int maxi(int a, int b);

float minf(float a, float b);

float maxf(float a, float b);

double mind(double a, double b);

double maxd(double a, double b);

void * myCalloc(size_t num, size_t size);

#endif