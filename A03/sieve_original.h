#ifndef SIEVE_ORIGINAL
#define SIEVE_ORIGINAL

#include <math.h>

void sieve_original(int id, unsigned int n,unsigned int p){
	if(id == 0) printf("Executing original sieve\n");
	unsigned int low_value = 2 + id*(n-1)/p;
	unsigned int high_value = 1 + (id+1)*(n-1)/p;
	unsigned int size = high_value - low_value + 1;

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();
	printf("([%d],%d,%d,%d)\n",id,low_value,high_value,size);

	unsigned int proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (int) sqrt((double) n)) {
		if (!id) printf ("Too many processes\n");
		MPI_Finalize();
		exit (1);
	}

	char *marked = (char *) malloc (size);
	if (marked == NULL) {
		printf ("[%d] Cannot allocate enough memory\n",id);
		MPI_Finalize();
		exit (1);
	}

	unsigned int i;
	for (i = 0; i < size; i++) marked[i] = 0;
	unsigned int index;
	if (!id) index = 0;
	unsigned int prime = 2;
	unsigned int first;

	do{
		if(prime * prime > low_value){
			first = prime * prime - low_value;
		}else{
			 if (!(low_value % prime)){
				 first = 0;
			 }else{
				 first = prime - (low_value % prime);
			 }
		}
		for (i = first; i < size; i += prime) marked[i] = 1;
		if(!id){
			while (marked[++index]);
			prime = index + 2;
		}
		if (p > 1) MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
	}while(prime * prime <= n);

	unsigned int count = 0;
	unsigned int global_count;
	for(i = 0; i < size; i++){
		if (!marked[i]) count++;
	}
	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_UNSIGNED, MPI_SUM,0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();

	if(!id){
		printf ("There are %d primes less than or equal to %d\n",global_count, n);
		printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}


}

#endif
