#ifndef SIEVE_CACHE
#define SIEVE_CACHE

#define ODD_INDEX(x) (((x-3)/2))
#define INDEX_ODD(x) (2*x + 3)

unsigned int localSieve(unsigned int n, unsigned int **sieve){
	unsigned int size = (n/2) - 1;
	char *marked = (char *) malloc (size);
	unsigned int i,j;

	for(i = 0;i<size;i++) marked[i] = 1;

	unsigned int prime = 3;
	unsigned int pp = prime * prime;

	while(pp < n){
		for(j = pp;j < n;j+=(prime<<1)){
			marked[ODD_INDEX(j)] = 0;
		}

		for(j = prime+2; j < n; j+=2){
			if(marked[ODD_INDEX(j)]==1){
				prime=j; break;
			}
		}
		pp = prime * prime;
	}

	unsigned int count = 1;
	for(i = 0; i <size; i++ ) count+=marked[i];
	//printf ("In the local sieve there are {%d} primes less than or equal to %d\n",count, n);
	*sieve = (unsigned int *) malloc(sizeof(unsigned int) * (count-1));

	j=0;
	for(i = 0; i < size; i++){
		if(marked[i] == 1){
			(*sieve)[j] = INDEX_ODD(i);
			j++;
		}
	}

	return count - 1;
}

uint64_t sieve_local_cache(int id, uint64_t n,uint64_t p){
	if(id == 0) printf("Executing odd numbers sieve with cache awareness \n");
	uint64_t low_value = (id==0) ? 3 : 2 + ((uint64_t)id)*((n-1))/p;
	if(!(id==0) ) low_value = (low_value % 2 == 0) ? low_value+1 : low_value;//Start with odd
	uint64_t high_value = 1 + ((uint64_t)(id+1))*((n-1))/p;
	uint64_t size = (high_value - low_value)/2 + 1;

	printf("([%d],%"PRIu64",%"PRIu64",%"PRIu64")\n",id,low_value,high_value,size);

	//MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();
	uint64_t proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (unsigned int) sqrt((double) n)) {
		if (id == 0) printf ("Too many processes\n");
		MPI_Finalize();
		exit (1);
	}

	unsigned int *sieve;
	unsigned int sqrt_n = (unsigned int) sqrt((double)n);
	unsigned int psize = localSieve(sqrt_n,&sieve);

	uint64_t i,j,k;
	unsigned int bsize = 5000;
	uint64_t csize = 0;
	char *marked = (char*)malloc(bsize/2);
	uint64_t prime;
	uint64_t first;
	uint64_t c = 0;

	i = low_value;
	while(i <= high_value){
		csize = (i+bsize < high_value) ? bsize : (high_value - i);//chunk size
		uint64_t lo = i;
		uint64_t hi = i+csize;
		for(j = 0 ; j < csize;j++) marked[j] = 1;

		//if(id==1) printf("lo,hi,csize: %"PRIu64",%"PRIu64",%"PRIu64"\n",lo,hi,csize);
		//printf("low_value,high_value,csize: %"PRIu64",%"PRIu64",%d\n",i,high_value,bsize);

		for(j=0;j<psize;j++){
			prime = sieve[j];
			first = lo;

			if(prime * prime > lo){
				first = prime * prime;
			}else{
				if ((lo % prime)== 0){
					first = lo;
				}else{
					first = (lo/prime)*prime;
					first = first + (prime << ( first & 1 ));
				}
			}
			//if(id==1) printf("p:%"PRIu64",lo:%"PRIu64",f:%"PRIu64",hi:%"PRIu64"\n",prime,lo,first,hi);

			uint64_t offset = ODD_INDEX(lo);
			for(k=first;k<=hi;k+=(prime<<1)){
				//if(ODD_INDEX(k) - offset >= bsize/2){
				//	printf("Error out of bounds %"PRIu64",%"PRIu64"\n",ODD_INDEX(k),offset);
				//	return 0;

				//if(id==1) printf("}%"PRIu64",%"PRIu64",%"PRIu64",%d\n",k,ODD_INDEX(k),offset,bsize/2);
				//if(id==1) printf("}%"PRIu64"\n",k);
				marked[ODD_INDEX(k) - offset] = 0;
			}
			//break;

		}
		//break;
		for(j = 0; j < csize/2; j++){
			//if(id==1 && marked[j]==1){
			//	printf("%"PRIu64", ",lo + 2*j);
			//}
			c = c + marked[j];
		}
		//if(id==1) printf("\n");

		i+=bsize;
	}

	if(id == 0) c++;//Do not forget to count 2 also
	printf("[%d],%"PRIu64"\n",id,c);
	uint64_t global_count=0;
	if (p > 1) MPI_Reduce (&c, &global_count, 1, MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();

	if(id==0){
		printf("There are {%"PRIu64"} primes less than or equal to %"PRIu64"\n",global_count, (uint64_t)n);
		printf ("Elapsed time of <Sieve with cache awareness> for (%"PRIu64") processes %10.6f\n", p, elapsed_time);
	}

	return global_count;
}

#endif
