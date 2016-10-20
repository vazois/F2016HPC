#ifndef SIEVE_LOCAL
#define SIEVE_LOCAL

#include <math.h>

#define ODD_INDEX(x) (((x-3)/2))
#define INDEX_ODD(x) (2*x + 3)

uint64_t sieve_local_odd(int id, uint64_t n,uint64_t p){
	if(id == 0) printf("Executing odd numbers sieve\n");
	uint64_t low_value = (id==0) ? 3 : 2 + ((uint64_t)id)*((n-1))/p;
	if(!(id==0) ) low_value = (low_value % 2 == 0) ? low_value+1 : low_value;//Start with odd
	uint64_t high_value = 1 + ((uint64_t)(id+1))*((n-1))/p;
	uint64_t size = (high_value - low_value)/2 + 1;

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();
	uint64_t proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (unsigned int) sqrt((double) n)) {
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

	uint64_t count =0;
	uint64_t index = 0;
	uint64_t prime = 3;
	uint64_t i,j;
	if (id == 0){
		//printf("([%d],%"PRIu64",%"PRIu64",%"PRIu64")\n",id,low_value,high_value,size);
		for (i = 0; i < size; i++) marked[i] = 0;
		while(prime * prime <= high_value){
			j = prime * prime;
			//printf("{%"PRIu64"}\n",j);
			for( ; j <= high_value;j+=(prime<<1)) marked[ODD_INDEX(j)]=1;
			while (marked[++index]);
			prime = INDEX_ODD(index);
		}

		for(i = 0;i<size;i++){
			if (marked[i] == 0){
				//printf("[%"PRIu64"],{%"PRIu64"}\n",i,INDEX_ODD(i));
				count++;
			}
		}
		count++;//count 2
		free(marked);
		printf ("There are {%"PRIu64"} primes less than or equal to %"PRIu64"\n",count, high_value);
	}else{
		uint64_t sqrt_n = (uint64_t)ceil(sqrt((double)n));
		uint64_t first;

		char *sieve = (char *) malloc (sqrt_n);
		if(id == 1) printf("sqrt_n:%"PRIu64"\n",sqrt_n);

		for (i = 0; i < size; i++) marked[i] = 0;
		for (i = 0; i < sqrt_n; i++) sieve[i] = 0;

		index = 0;
		prime = 3;
		while(prime*prime <=n){
			//if(id==1) printf("%"PRIu64"\n",prime);
			if(prime * prime > low_value){
					first = prime * prime;
			}else{
				if ((low_value % prime)== 0){//if multiple of current prime start from 0 local index
					first = low_value;
				}else{//find next multiple of current prime//convert it to local odd index//
					first = (low_value/prime)*prime;
					first = ( first % 2 == 0 ) ? first+prime : first + (prime << 1);//
				}
			}
			uint64_t i = first;
			uint64_t offset = ODD_INDEX(low_value);
			while( i <= high_value ){
				marked[ODD_INDEX(i) - offset]=1;
				i+=(prime <<1);
			}

			//Find next prime to test
			j = prime * prime;
			//printf("{%"PRIu64"}\n",j);
			for( ; j <= sqrt_n;j+=(prime<<1)){ sieve[ODD_INDEX(j)]=1; }
			while (sieve[++index]);
			prime = INDEX_ODD(index);
		}

		for(i = 0; i < size; i++) if (marked[i] == 0) count++;
		free(marked);
	}

	uint64_t global_count=0;
	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();
	if(!id){
		printf ("There are {%"PRIu64"} primes less than or equal to %"PRIu64"\n",global_count, (uint64_t)n);
		printf ("SIEVE LOCAL ODD (%"PRIu64") %10.6f\n", p, elapsed_time);
	}

	return global_count;
}


#endif
