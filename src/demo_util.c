#include "demo_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <stdarg.h>   /* for va_list, va_start, va_end, vprintf */
#include <string.h>   /* For atoi */

static int s_rank;

/* --------------------------------------------------------------------
    Arrays : 

    Routines that allocate memory and possibly assign
    values.  Call these as follows : 

        double *x;    // No memory allocated yet
        int n = 10;

        empty_array(n,&x);  // x is now ready for use 
        x[0] = 1;
        x[1] = 2;
        // ...

        delete_array(&x);  // clean up when you are done!

    Array routine : 
        empty_array     : Allocate memory, but don't assign any values.
        ones_array      : Create an array of ones
        linspace_array  : Analogous to Matlab's 'linspace' 
        random_array    : Array of random numbers in [0,1]
        delete_array    : Deletes memory for x
   ------------------------------------------------------------------ */ 

void empty_array(int n,double **x)
{
    *x = malloc(n*sizeof(double));
}

void ones_array(int n,double **x)
{
    int i;
    *x = malloc(n*sizeof(double));
    for(i = 0; i < n; i++)
    {
        (*x)[i] = 1;
    }
}

void linspace_array(double a,double b,int n,double **x)
{
    double h = (b-a)/(n-1);
    int i;
    empty_array(n,x);
    for(i = 0; i < n; i++)
    {
        (*x)[i] = a + i*h;
    }
}

void random_array(int n, double **x)
{
    *x = malloc(n*sizeof(double));    
    random_seed();   
    int i;

    for(i=0; i < n; i++)
    {
        (*x)[i] = random_number();
    }
}

void delete_array(double **x)
{
    free(*x);
}


/* --------------------------------------------------------------------
    Operations on arrays

    Routines : 
        sum_array  : sum entries and return scalar.
    ----------------------------------------------------------------- */

double sum_array(int n, double *x)
{
    int i;
    double s;

    s = 0;
    for(i = 0; i < n; i++)
    {
        s += x[i];
    }
    return s;
}

/* --------------------------------------------------------------------
    Input/Output

    Read routines from the command line; 
    Print values either from node 0 only or from each processor. 

        print_global  : Print only from node 0
        print_debug   : print from each processor

    Example : 

        print_global("hello!\n");

        returns : 
        Processor [0] : hello!

        print_debug("hello!\n");

        returns : 
        Processor [0] : hello!
        Processor [1] : hello!
        Processor [2] : hello!
        Processor [3] : hello!
        
    ------------------------------------------------------------------ */

void read_int(int argc, char** argv, char arg[], int* value,int *err)
{
    *err = 1;  /* Nothing found yet */
    int arg_index = 1;     /* Skip first argument */
    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], arg) == 0)
        {
            arg_index++;
            *value = atoi(argv[arg_index++]);   
            *err = 0;         
            return;
        }
        else
        {
            arg_index++;
        }
    }
}


void print_global(const char* format, ... )
{
    /* Only print if on processor 0 */
    if (s_rank == 0)
    {
        va_list arglist;
        printf( "Processor [0] : " );
        va_start( arglist, format );
        vprintf( format, arglist );
        va_end( arglist );
    }
}

void print_debug(const char* format, ... )
{
    /* Include rank number in print statement */
    va_list arglist;
    printf( "Processor [%d] : ",s_rank);
    va_start( arglist, format );
    vprintf( format, arglist );
    va_end( arglist );
}

/* -------------------------------------------------
    Miscellaneous routines
   ----------------------------------------------- */ 
void set_rank(int  rank)
{
    s_rank = rank;
}

void sleep(double t_total)
{
    double t0, t1;
    t0 = clock();
    t1 = t0;
    while ((t1-t0)/CLOCKS_PER_SEC < t_total)
    {
        t1 = clock();
    }
}

double random_number()
{
  return (double) rand() / (double) RAND_MAX ;
}

void random_seed()
{
    srand(time(NULL));
}

int pow2(int p)
{
    /* Compute n = 2^p */
    int n,i;

    n = 1;
    for(i = 0; i < p; i++)
    {
        n *= 2;
    }
    return n;
}

