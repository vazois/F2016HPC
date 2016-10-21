#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "sieve_original.h"
#include "sieve_odd.h"
#include "sieve_local.h"
#include "sieve_cache.h"


void validate(uint64_t a, uint64_t b, char *msg){

	char *check = a == b ? "(PASS)" : "(FAIL)";
	printf("Comparison between %s,%s\n",msg,check);

}

int main(int argc, char **argv){
	int id;	// My Id
	int p;	// Number of processes
	uint64_t n;
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

	/*uint64_t res_orig=sieve_original(id,n,p);
	MPI_Barrier(MPI_COMM_WORLD);
	uint64_t res_odd=sieve_odd(id,n,p);
	MPI_Barrier(MPI_COMM_WORLD);
	uint64_t res_local_odd=sieve_local_odd(id,n,p);


	if(id==0) validate(res_orig,res_odd,"original to odd version");
	if(id==0) validate(res_orig,res_local_odd,"original to local odd version");*/


	/*if(id==0){
		unsigned int *sieve;
		unsigned int sqrt_n = (unsigned int)sqrt((double)n);
		unsigned int size = localSieve(sqrt_n,&sieve);

		//unsigned int i=0;
		//for(i = 0;i<size;i++){
		//	printf("%d ",sieve[i]);
		//}
		//printf("\n");
	}*/

	uint64_t res_local_odd=sieve_local_cache(id,n,p);


	MPI_Finalize();
}
