#ifndef SIEVE_ORIGINAL
#define SIEVE_ORIGINAL


#include <math.h>
#include <stdint.h>
#include <inttypes.h>

uint64_t sieve_original(int id, uint64_t n,uint64_t p){
	if(id == 0) printf("Executing original sieve\n");
	uint64_t low_value = 2 + ((uint64_t)id)*(n-1)/p;
	uint64_t  high_value = 1 + ((uint64_t)(id+1))*(n-1)/p;
	uint64_t size = high_value - low_value + 1;

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();
	//printf("([%d],%"PRIu64",%"PRIu64",%"PRIu64")\n",id,low_value,high_value,size);

	uint64_t proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (uint64_t) sqrt((double) n)) {
		if (id == 0) printf ("Too many processes\n");
		MPI_Finalize();
		exit (1);
	}

	char *marked = (char *) malloc (size);
	if (marked == NULL) {
		printf ("[%d] Cannot allocate enough memory\n",id);
		MPI_Finalize();
		exit (1);
	}

	uint64_t i;
	for (i = 0; i < size; i++) marked[i] = 0;
	uint64_t index;
	if (id == 0) index = 0;
	uint64_t prime = 2;
	uint64_t first;

	do{
		if(prime * prime > low_value){
			first = prime * prime - low_value;
		}else{
			 if ((low_value % prime)== 0){//if multiple of current prime start from it// 26 % 2 == 0 -> 26,28,30,32 -> 0,2,4,6,8
				 first = 0;
			 }else{//find next multiple of current prime// 51%2==1 -> 51+1=52,54,56,58 -> 1,3,5,7
				 first = prime - (low_value % prime);
			 }
		}
		for (i = first; i < size; i += prime) marked[i] = 1;
		if(id == 0){
			while (marked[++index] == 1);
			prime = index + 2;
		}
		if( p > 1 ) MPI_Bcast (&prime,  1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	}while(prime * prime <= n);



	uint64_t count = 0;
	uint64_t global_count;
	for(i = 0; i < size; i++){
		if (marked[i] == 0) count++;
	}
	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();

	if(!id){
		printf ("There are {%"PRIu64"} primes less than or equal to %"PRIu64"\n",global_count, (uint64_t)n);
		printf ("Elapsed time of <Original Sieve> for (%"PRIu64") processes %10.6f\n", p, elapsed_time);
	}
	return global_count;
}

#endif
