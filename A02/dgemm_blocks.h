#ifndef DGEMM_BLOCKS_H
#define DGEMM_BLOCKS_H

void dgemm_blocks_ijk(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){

	for(unsigned int i = 0; i < N; i+=Bz){
		for(unsigned int j = 0; j < N; j+=Bz){
			for(unsigned int k = 0; k < N; k+=Bz){

				//Blocked Matrix
				for(unsigned int ii = i; ii < i + Bz; ii++){
					for(unsigned int jj = j; jj < j + Bz; jj++){
						double rC = C[ii * N + jj];
						for(unsigned int kk = k; kk < k + Bz; kk++){
								rC+= A[ii * N + kk] * B[kk*N + jj];
						}
						C[ii * N + jj] = rC;
					}
				}

			}
		}
	}
}

#endif
