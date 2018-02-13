#ifndef DEMO_H
#define DEMO_H

#include <time.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void set_rank(int rank);

double random_number();
void random_seed();

void empty_array(int n,double **x);
void ones_array(int n,double **x);
void random_array(int n, double **array);
double sum_array(int n, double *x);

int pow2(int p);


void delete_array(double **x);

void print_global(const char* format, ... );
void print_debug(const char* format, ... );


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
