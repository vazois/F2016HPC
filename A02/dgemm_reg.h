#ifndef DGEMM_REG_H
#define DGEMM_REG_H

//1
void dgemm_reg_ijk(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int i = 0; i < N ;i++){
		for(unsigned int j = 0; j < N ;j++){
			register double rC = 0;
			for(unsigned int k = 0; k < N ;k++){
				rC += A[ i * N + k ] * B[ k * N + j ];
			}
			C[ i * N + j ] = rC;
		}
	}
}

//2
void dgemm_reg_jik(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int j = 0; j < N ;j++){
		for(unsigned int i = 0; i < N ;i++){
			register double rC = 0;
			for(unsigned int k = 0; k < N ;k++){
				rC += A[ i * N + k ] * B[ k * N + j ];
			}
			C[ i * N + j ] = rC;
		}
	}
}

//3
void dgemm_reg_kij(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int k = 0 ; k < N ; k++){
		for(unsigned int i=0; i < N ; i++){
			double rA = A[i * N + k];
			for(unsigned int j=0; j < N ; j++){
				C[i * N + j]+= rA * B[k*N + j];
			}
		}
	}
}

//4
void dgemm_reg_ikj(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int i=0; i < N ; i++){
		for(unsigned int k = 0 ; k < N ; k++){
			double rA = A[i * N + k];
			for(unsigned int j=0; j < N ; j++){
				C[i * N + j]+= rA * B[k*N + j];
			}
		}
	}
}

//5
void dgemm_reg_jki(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int j=0; j < N ; j++){
		for(unsigned int k=0; k < N ; k++){
			double rB = B[k * N + j];
			for(unsigned int i=0; i < N ; i++){
				C[i * N + j] += A[i * N + k] * rB;
			}
		}
	}
}

//6
void dgemm_reg_kji(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int k=0; k < N ; k++){
		for(unsigned int j=0; j < N ; j++){
			double rB = B[k * N + j];
			for(unsigned int i=0; i < N ; i++){
				C[i * N + j] += A[i * N + k] * rB;
			}
		}
	}
}





#endif
