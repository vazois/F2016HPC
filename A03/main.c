#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sieve_original.h"
#include "sieve_odd.h"


void validate(unsigned int a, unsigned int b, char *msg){

	char *check = a == b ? "(PASS)" : "(FAIL)";
	printf("Comparison between %s,%s\n",msg,check);

}

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


	unsigned int res_orig=sieve_original(id,n,p);
	MPI_Barrier(MPI_COMM_WORLD);
	unsigned int res_odd=sieve_odd(id,n,p);
	MPI_Barrier(MPI_COMM_WORLD);

	if(id==0) validate(res_orig,res_odd,"original to odd version");



	MPI_Finalize();
}
