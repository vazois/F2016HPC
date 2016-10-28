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

	if (argc <3){
		printf("Please provide maximum number n and BSIZE for cached version\n");
		MPI_Finalize();
		exit(1);
	}

	exp = atoi(argv[1]);
	n = pow(10,(double)exp);
	BSIZE = atoi(argv[2]);


	uint64_t rlc0,rlc1;
	//uint64_t res_orig=sieve_original(id,n,p);
	//MPI_Barrier(MPI_COMM_WORLD);
	//uint64_t res_odd=sieve_odd(id,n,p);
	//MPI_Barrier(MPI_COMM_WORLD);
	uint64_t res_local_odd=sieve_local_odd(id,n,p,&rlc0);
	MPI_Barrier(MPI_COMM_WORLD);
	uint64_t res_local_cache=sieve_local_cache(id,n,p,&rlc1);

	//if(rlc0!=rlc1) printf("{%d},%"PRIu64",%"PRIu64"\n",id,rlc0,rlc1);

	//if(id==0) validate(res_orig,res_odd,"original to odd version");
	//if(id==0) validate(res_orig,res_local_odd,"original to local odd version");
	//if(id==0) validate(res_orig,res_local_cache,"original to local cache aware version");
	if(id==0) validate(res_local_odd,res_local_cache,"local odd to local cache aware version");


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}
