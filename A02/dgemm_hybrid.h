#ifndef DGEMM_HYBRID_H
#define DGEMM_HYBRID_H

void dgemm_hybrid(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int i = 0; i < N; i+=Bz){
		for(unsigned int j = 0; j < N; j+=Bz){
			for(unsigned int k = 0; k < N; k+=Bz){

				//Blocked Matrix
				for(unsigned int ii = i; ii < i + Bz; ii+=2){
					for(unsigned int jj = j; jj < j + Bz; jj+=2){
						unsigned int iC = ii * N + jj; unsigned int iiC = iC + N;
						register double rC0 = C[iC];
						register double rC1 = C[iC+1];
						register double rC2 = C[iiC];
						register double rC3 = C[iiC+1];


						for(unsigned int kk = k; kk < k + Bz; kk+=2){
							register int iA = ii * N + kk; register int iiA = iA + N;
							register int iB = kk * N + jj; register int iiB = iB + N;

							register double rA0 = A[iA]; register double rA1 = A[iiA];
							register double rB0 = B[iB]; register double rB1 = B[iB + 1];

							rC0+=rA0 * rB0; rC1+=rA0 * rB1; rC2+=rA1 * rB0; rC3+=rA1 * rB1;

							rA0 = A[iA+1]; rA1 = A[iiA+1]; rB0 = B[iiB]; rB1 = B[iiB + 1];

							rC0+=rA0 * rB0; rC1+=rA0 * rB1; rC2+=rA1 * rB0; rC3+=rA1 * rB1;

						}
						C[iC] = rC0;
						C[iC+1] = rC1;
						C[iiC] = rC2;
						C[iiC+1] = rC3;
					}
				}

			}
		}
	}
}


#endif
