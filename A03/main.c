#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sieve_original.h"


int main(int argc, char **argv){
	int id;	// My Id
	int p;	// Number of processes
	unsigned int n;
	int exp;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);

	if (argc <2){
		printf("Please provide maximum number n\n");
		MPI_Finalize();
		exit(1);
	}

	exp = atoi(argv[1]);
	n = pow(10,(double)exp);


	sieve_original(id,n,p);


	MPI_Finalize();
}
