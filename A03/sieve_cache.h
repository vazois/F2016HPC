#ifndef SIEVE_CACHE
#define SIEVE_CACHE

//#define BSIZE 5000
#define ODD_INDEXC(x) (((x-3)/2))
#define INDEX_ODDC(x) (2*x + 3)
//#define ODD_INDEXC(x) (((x-3)>>1))
//#define INDEX_ODDC(x) (x<<1 + 3)

uint64_t BSIZE = 500000;

unsigned int localSieve(unsigned int n, unsigned int **sieve){
	unsigned int size = (n/2) - 1;
	char *marked = (char *) malloc (size);
	unsigned int i,j;

	for(i = 0;i<size;i++) marked[i] = 1;

	unsigned int prime = 3;
	unsigned int pp = prime * prime;

	while(pp < n){
		for(j = pp;j < n;j+=(prime<<1)){
			marked[ODD_INDEXC(j)] = 0;
		}

		for(j = prime+2; j < n; j+=2){
			if(marked[ODD_INDEXC(j)]==1){
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
			(*sieve)[j] = INDEX_ODDC(i);
			j++;
		}
	}

	return count - 1;
}

void mark(char **marked,uint64_t p1,uint64_t offset,uint64_t lo,uint64_t odds){
	uint64_t first = p1 * p1;
	if (first <= lo){
		if ((lo % p1) == 0){
			first = lo;
		}else{
			first = (lo/p1)*p1;
			first = first + (p1 << ( first & 1 ));
		}
	}
	uint64_t k;
	for(k=ODD_INDEX(first) - offset;k<=odds;k+=(p1)){ (*marked)[k] = 0; }

}

uint64_t sieve_local_cache(int id, uint64_t n,uint64_t p, uint64_t *lcount){
	if(id == 0) printf("Executing odd numbers sieve with cache awareness \n");
	uint64_t low_value = (id==0) ? 3 : 2 + ((uint64_t)id)*((n-1))/p;
	if(!(id==0) ) low_value = (low_value % 2 == 0) ? low_value+1 : low_value;//Start with odd
	uint64_t high_value = 1 + ((uint64_t)(id+1))*((n-1))/p;
	uint64_t size = (high_value - low_value)/2 + 1;

	//printf("([%d],%"PRIu64",%"PRIu64",%"PRIu64")\n",id,low_value,high_value,size);
	if(id == 0) printf("Cached Array Size:%"PRIu64"\n",BSIZE);
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

	uint64_t i,j;
	uint64_t bsize = BSIZE;
	uint64_t csize = 0;
	char *marked = (char*)malloc(bsize/2 + 1);
	uint64_t prime;
	//uint64_t first;
	uint64_t c = 0;
	uint64_t offset;
	//uint64_t k=0;

	i = low_value;
	while(i <= high_value){
		csize = (i+bsize <= high_value) ? bsize : (high_value - i);//chunk size
		uint64_t lo = i;
		//uint64_t hi = i+csize;
		uint64_t odds = (i+csize - lo)/2;
		uint64_t offset = ODD_INDEX(lo);

		for(j = 0 ; j < bsize/2;j++) marked[j] = 1;

		//if(id==0) printf("lo,hi,csize: %"PRIu64",%"PRIu64",%"PRIu64"\n",lo,hi,csize);
		//if(id==0) printf("low_value,high_value,csize: %"PRIu64",%"PRIu64",%d\n",low_value,high_value,bsize);

		/*for(j=0;j<7;j+=3){
			prime = sieve[j];
			if(prime == 0 && j+1 >=psize){
				printf("%d,%"PRIu64",%d\n",id,j,psize);
			}
			first = prime * prime;
			if (first <= lo){
				if ((lo % prime) == 0){
					first = lo;
				}else{
					first = (lo/prime)*prime;
					first = first + (prime << ( first & 1 ));
				}
			}
			for(k=ODD_INDEX(first) - offset;k<=odds;k+=(prime)){ marked[k] = 0; }
		}*/

		uint64_t p1 = sieve[0];
		uint64_t p2 = sieve[1];
		uint64_t p3 = sieve[2];
		uint64_t p4;

		mark(&marked,p1,offset,lo,odds);
		mark(&marked,p2,offset,lo,odds);
		mark(&marked,p3,offset,lo,odds);

		for(j=3;j<psize;j+=4){
			p1 = sieve[j];
			p2 = sieve[j+1];
			p3 = sieve[j+2];
			p4 = sieve[j+3];

			mark(&marked,p1,offset,lo,odds);
			mark(&marked,p2,offset,lo,odds);
			mark(&marked,p3,offset,lo,odds);
			mark(&marked,p4,offset,lo,odds);
		}

		for(j = 0; j <= odds; j++){ c = c + marked[j]; }
		i+=bsize;
	}

	if(id == 0) c++;//Do not forget to count 2 also
	uint64_t global_count=0;
	if (p > 1) MPI_Reduce (&c, &global_count, 1, MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();

	if(id==0){
		printf("There are {%"PRIu64"} primes less than or equal to %"PRIu64"\n",global_count, (uint64_t)n);
		printf ("Elapsed time of <Sieve with cache awareness> for (%"PRIu64") processes %10.6f\n", p, elapsed_time);
	}

	*lcount = c;
	return global_count;
}

#endif
