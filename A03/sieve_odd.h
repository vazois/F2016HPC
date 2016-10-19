#ifndef SIEVE_ODD
#define SIEVE_ODD

#include <math.h>

#define ODD_INDEX(x) (((x-3)/2))
#define INDEX_ODD(x) (2*x + 3)

unsigned int sieve_odd(int id, unsigned int n,unsigned int p){
	if(id == 0) printf("Executing odd numbers sieve\n");
	unsigned int low_value = (id==0) ? 3 : 2 + id*((n-1))/p;
	if(!(id==0) ) low_value = (low_value % 2 == 0) ? low_value+1 : low_value;//Start with odd
	unsigned int high_value = 1 + (id+1)*((n-1))/p;
	unsigned int size = (high_value - low_value)/2 + 1;

	/*int k;
	int c=-1;
	if(id==c){
	for(k=low_value; k<=high_value;k+=2){
		printf("%d{}%d,",(ODD_INDEX(k) - ODD_INDEX(low_value)),k);
	}
	printf("\n");
	}*/

	//MPI_Barrier(MPI_COMM_WORLD);
	//if(id==c)printf("([%d],%d,%d,%d)\n",id,low_value,high_value,size);
	MPI_Barrier(MPI_COMM_WORLD);

	double elapsed_time = -MPI_Wtime();

	unsigned int proc0_size = (n-1)/p;

	if ((2 + proc0_size) < (int) sqrt((double) n)) {
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

	unsigned int i;
	for (i = 0; i < size; i++) marked[i] = 0;
	unsigned int index;
	if (id == 0) index = 0;
	unsigned int prime = 3;
	unsigned int first;

	do
	{
		//if(id == 0) printf("current prime: %d\n",prime);
		if(prime * prime > low_value){
			//first = ODD_INDEX(prime * prime) - ODD_INDEX(low_value);//Find relative index for odd number// prime * prime will be odd always
			first = prime * prime;
		}else{
			 if ((low_value % prime)== 0){//if multiple of current prime start from 0 local index
				 //first = ODD_INDEX(low_value);
				 //first = 0;
				 first = low_value;
			 }else{//find next multiple of current prime//convert it to local odd index//
				 //first = low_value/prime;
				 //first = ( first % 2 == 0 ) ? first+prime : first + (prime << 1);//
				 //first = ODD_INDEX(first) - ODD_INDEX(low_value);

				 first = (low_value/prime)*prime;
				 first = ( first % 2 == 0 ) ? first+prime : first + (prime << 1);//
			 }
		}

		i=first;
		unsigned int offset = ODD_INDEX(low_value);
		while( i <= high_value){
			marked[ODD_INDEX(i) - offset]=1;
			i+=(prime <<1);
		}

		if(id == 0){
			while (marked[++index]);
			prime = INDEX_ODD(index);
		}
		if( p > 1 ) MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
	}while(prime * prime <= n);

	/*if(id==c){
	for(i = low_value; i < high_value; i+=2){
		if(marked[ODD_INDEX(i)-ODD_INDEX(low_value)]==0){
			printf("[%d],%d\n",ODD_INDEX(i)-ODD_INDEX(low_value) ,i);
		}
	}
	}*/

	unsigned int count = (id==0) ? 1 : 0;
	unsigned int global_count=0;
	for(i = 0; i < size; i++){
		if (marked[i] == 0) count++;
	}
	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_UNSIGNED, MPI_SUM,0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();
	if(!id){
		printf ("There are %d primes less than or equal to %d\n",global_count, n);
		printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}
	return global_count;
}

#endif
