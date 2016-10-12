#ifndef DGEMM_BLOCKS_H
#define DGEMM_BLOCKS_H

//1
void dgemm_blocks_ijk(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int i = 0; i < N; i+=Bz){
		for(unsigned int j = 0; j < N; j+=Bz){
			for(unsigned int k = 0; k < N; k+=Bz){

				//Blocked Matrix
				for(unsigned int ii = i; ii < i + Bz; ii++){
					for(unsigned int jj = j; jj < j + Bz; jj++){
						register double rC = C[ii * N + jj];
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

//2
void dgemm_blocks_jik(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int j = 0; j < N; j+=Bz){
		for(unsigned int i = 0; i < N; i+=Bz){
			for(unsigned int k = 0; k < N; k+=Bz){

				//Blocked Matrix
				for(unsigned int jj = j; jj < j + Bz; jj++){
					for(unsigned int ii = i; ii < i + Bz; ii++){
						register double rC = C[ii * N + jj];
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

//3
void dgemm_blocks_ikj(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int i = 0; i < N; i+=Bz){
		for(unsigned int k = 0; k < N; k+=Bz){
			for(unsigned int j = 0; j < N; j+=Bz){

				//Blocked Matrix
				for(unsigned int ii = i; ii < i + Bz; ii++){
					for(unsigned int kk = k; kk < k + Bz; kk++){
						double rA = A[ii * N + kk];
						for(unsigned int jj = j; jj < j + Bz; jj++){
							C[ii * N + jj]+= rA * B[kk*N + jj];
						}
					}
				}
				//
			}
		}
	}
}

//4
void dgemm_blocks_kij(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int k = 0; k < N; k+=Bz){
		for(unsigned int i = 0; i < N; i+=Bz){
			for(unsigned int j = 0; j < N; j+=Bz){

				//Blocked Matrix
				for(unsigned int kk = k; kk < k + Bz; kk++){
					for(unsigned int ii = i; ii < i + Bz; ii++){
						double rA = A[ii * N + kk];
						for(unsigned int jj = j; jj < j + Bz; jj++){
							C[ii * N + jj]+= rA * B[kk*N + jj];
						}
					}
				}
				//
			}
		}
	}
}

//5
void dgemm_blocks_jki(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int j = 0; j < N; j+=Bz){
		for(unsigned int k = 0; k < N; k+=Bz){
			for(unsigned int i = 0; i < N; i+=Bz){

				//Blocked Matrix
				for(unsigned int jj = j; jj < j + Bz; jj++){
					for(unsigned int kk = k; kk < k + Bz; kk++){
						register double rB = B[kk*N + jj];
						for(unsigned int ii = i; ii < i + Bz; ii++){
							C[ii * N + jj]+= A[ii * N + kk] * rB;
						}
					}
				}
				//
			}
		}
	}
}

//6
void dgemm_blocks_kji(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int k = 0; k < N; k+=Bz){
		for(unsigned int j = 0; j < N; j+=Bz){
			for(unsigned int i = 0; i < N; i+=Bz){

				//Blocked Matrix
				for(unsigned int kk = k; kk < k + Bz; kk++){
					for(unsigned int jj = j; jj < j + Bz; jj++){
						register double rB = B[kk*N + jj];
						for(unsigned int ii = i; ii < i + Bz; ii++){
							C[ii * N + jj]+= A[ii * N + kk] * rB;
						}
					}
				}
				//
			}
		}
	}
}

#endif
