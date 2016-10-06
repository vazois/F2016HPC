#ifndef DGEMM_REG_H
#define DGEMM_REG_H

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





#endif
