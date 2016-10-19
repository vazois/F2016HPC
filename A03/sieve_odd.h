#ifndef SIEVE_ODD
#define SIEVE_ODD

#include <math.h>

#define ODD_INDEX(x) (((x-3)/2))
#define INDEX_ODD(x) (2*x + 3)

void sieve_odd(int id, unsigned int n,unsigned int p){
	if(id == 0) printf("Executing odd numbers sieve\n");
	unsigned int low_value = (id==0) ? 3 : 2 + id*((n-1))/p;
	if(!(id==0) ) low_value = (low_value % 2 == 0) ? low_value+1 : low_value;//Start with odd
	unsigned int high_value = 1 + (id+1)*((n-1))/p;
	unsigned int size = (high_value - low_value)/2 + 1;

	int k;
	if(id==1){
	for(k=low_value; k<=high_value;k+=2){
		printf("%d{}%d,",(ODD_INDEX(k) - ODD_INDEX(low_value)),k);
	}
	printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();
	printf("([%d],%d,%d,%d)\n",id,low_value,high_value,size);

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


	do{
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

				 first = low_value/prime;
				 first = ( first % 2 == 0 ) ? first+prime : first + (prime << 1);//
			 }
		}

		//for (i = first; i < size; i += prime) marked[i] = 1;

		i=first;
		unsigned int offset = ODD_INDEX(low_value);
		while( i < high_value){
			marked[ODD_INDEX(i) - offset]=1;
			i+=(prime <<1);
		}

		if(id == 0){
			while (marked[++index] == 1)
			prime = INDEX_ODD(index);
			printf("id == 0 odd\n");
		}
		break;
		if( p > 1 ) MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
	}while(prime * prime <= n);

	unsigned int count = 0;
	unsigned int global_count;
	for(i = 0; i < size; i++){
		if (marked[i] == 0) count++;
	}
	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_UNSIGNED, MPI_SUM,0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();
	if(!id){
		printf ("There are %d primes less than or equal to %d\n",global_count, n);
		printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
	}

}

#endif
