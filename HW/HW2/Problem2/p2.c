#include "demo7.h"
#include <demo_util.h>
#include <mpi.h>
#include <math.h>

#define ROOT_PROCESSOR 0

/* Function Prototype */
double trapz(double low, double high, int n, double width);

/* Integrate this over [-1,1] */
double f(double x)
{
    double fx;
	fx = (x-1)*(x-1)*exp(-1*x*x);
    return fx;
}

/* Trapz function for integration */
double trapz(double low, double high, int n, double width)
{
	int i;
	double lowLimit = low;
    double integral = 0; // hold the value of integration
	integral = (f(low) + f(high)) / 2.0;
	
	for (i = 1; i < n; i++)
	{
		lowLimit = lowLimit + width;
		integral += f(lowLimit);
	}
	return integral*width;
}

void main(int argc, char** argv)
{
    int n_global;
	double range[2];
	double exact_ans = 1.87259295726583875460;
    /* MPI variables */
    int my_rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    set_rank(my_rank);  /* Used in printing */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (my_rank == 0)
    {        
        int p0,err;
        read_int(argc,argv, "-p",&p0,&err);
        if (err > 0)
        {
            print_global("Command line argument '-p' not found.\n");
            exit(0);
        }
        n_global = pow2(p0);     /* Number of sub-intervals used for integration */

        /* Hardwire values */
        double a,b;
        a = -1;
        b = 1;  

        /* Your Node P=0 work goes here */
        /* Send sub-interval to other processors */
        double w = (b-a)/nprocs;
		int p;
        for(p = 1; p < nprocs; p++)
        {
            /* pass two values to processor p */
            range[0] = a + p*w;
            range[1] = range[0] + w;
            int tag = 0;
            int dest = p;
            MPI_Send((void*) range,2,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
        }
		range[0] = a;
		range[1] = range[0] + w;
    }
    else
    {
        MPI_Status status;
        int count;
        /* Receive range values */
        int source = 0;
        int tag = 0;
        MPI_Recv((void*) range,2,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&status);   
        MPI_Get_count(&status,MPI_DOUBLE,&count);         
    }

    /* Broadcast value of N to all processors;  compute number of panels in each 
    subinterval */
	MPI_Bcast(&n_global, 1, MPI_INT, (int)ROOT_PROCESSOR, MPI_COMM_WORLD);
	int n_local = n_global/nprocs;
	double w_local = (range[1]-range[0])/n_local;

    /* Every processor now knows its range and number of panels */
    /* Apply trapezoidal rule (for loop) */
    double integral_local = 0;
	integral_local = trapz(range[0], range[1], n_local, w_local);
	
	/* Call MPI_Reduce to get final integral */
	double integral_global = 0;
	MPI_Reduce((void *)&integral_local, (void *)&integral_global, 1, MPI_DOUBLE, MPI_SUM, (int)ROOT_PROCESSOR, MPI_COMM_WORLD);
	
	/* Node 0 prints the results - print_global() */
	print_global("%10d   %.16f   %.16f   %.16f\n", n_global, integral_global, exact_ans, (integral_global-exact_ans));
    MPI_Finalize();
}
