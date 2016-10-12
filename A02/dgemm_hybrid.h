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

void dgemm_hybrid2(const double *A, const double *B, double *&C, unsigned int N, unsigned int Bz){
	for(unsigned int i = 0; i < N; i+=Bz){
		for(unsigned int j = 0; j < N; j+=Bz){
			for(unsigned int k = 0; k < N; k+=Bz){

				//Blocked Matrix
				for(unsigned int ii = i; ii < i + Bz; ii+=4){
					for(unsigned int jj = j; jj < j + Bz; jj+=4){
						unsigned int iC = ii * N + jj;
						unsigned int iiC = iC + N;
						unsigned int iiiC = iiC + N;
						unsigned int iiiiC = iiiC + N;

						register double rC0 = C[iC]; register double rC1 = C[iC+1]; register double rC2 = C[iC+2]; register double rC3 = C[iC+3];
						register double rC4 = C[iiC]; register double rC5 = C[iiC+1]; register double rC6 = C[iiC+2]; register double rC7 = C[iiC+3];
						register double rC8 = C[iiiC]; register double rC9 = C[iiiC+1]; register double rC10 = C[iiiC+2]; register double rC11 = C[iiiC+3];
						register double rC12 = C[iiiiC]; register double rC13 = C[iiiiC+1]; register double rC14 = C[iiiiC+2]; register double rC15 = C[iiiiC+3];


						for(unsigned int kk = k; kk < k + Bz; kk+=4){
							register int iA = ii * N + kk;
							register int iiA = iA + N;
							register int iiiA = iiA + N;
							register int iiiiA = iiiA + N;

							register int iB = kk * N + jj;
							register int iiB = iB + N;
							register int iiiB = iiB + N;
							register int iiiiB = iiiB + N;

							//LOAD(1)
							register double rA0 = A[iA]; register double rA1 = A[iiA]; register double rA2 = A[iiiA]; register double rA3 = A[iiiiA];
							register double rB0 = B[iB]; register double rB1 = B[iB + 1]; register double rB2 = B[iB + 2]; register double rB3 = B[iB + 3];
							//COMPUTE
							rC0+=rA0 * rB0; rC1+=rA0 * rB1; rC2+=rA0 * rB2; rC3+=rA0 * rB3;
							rC4+=rA1 * rB0; rC5+=rA1 * rB1; rC6+=rA1 * rB2; rC7+=rA1 * rB3;
							rC8+=rA2 * rB0; rC9+=rA2 * rB1; rC10+=rA2 * rB2; rC11+=rA2 * rB3;
							rC12+=rA3 * rB0; rC13+=rA3 * rB1; rC14+=rA3 * rB2; rC15+=rA3 * rB3;

							//LOAD(2)
							rA0 = A[iA+1]; rA1 = A[iiA+1]; rA2 = A[iiiA+1]; rA3 = A[iiiiA+1];
							rB0 = B[iiB]; rB1 = B[iiB + 1]; rB2 = B[iiB + 2]; rB1 = B[iiB + 3];
							//COMPUTE
							rC0+=rA0 * rB0; rC1+=rA0 * rB1; rC2+=rA0 * rB2; rC3+=rA0 * rB3;
							rC4+=rA1 * rB0; rC5+=rA1 * rB1; rC6+=rA1 * rB2; rC7+=rA1 * rB3;
							rC8+=rA2 * rB0; rC9+=rA2 * rB1; rC10+=rA2 * rB2; rC11+=rA2 * rB3;
							rC12+=rA3 * rB0; rC13+=rA3 * rB1; rC14+=rA3 * rB2; rC15+=rA3 * rB3;

							//LOAD(3)
							rA0 = A[iA+2]; rA1 = A[iiA+2]; rA2 = A[iiiA+2]; rA3 = A[iiiiA+2];
							rB0 = B[iiiB]; rB1 = B[iiiB + 1]; rB2 = B[iiiB + 2]; rB1 = B[iiiB + 3];
							//COMPUTE
							rC0+=rA0 * rB0; rC1+=rA0 * rB1; rC2+=rA0 * rB2; rC3+=rA0 * rB3;
							rC4+=rA1 * rB0; rC5+=rA1 * rB1; rC6+=rA1 * rB2; rC7+=rA1 * rB3;
							rC8+=rA2 * rB0; rC9+=rA2 * rB1; rC10+=rA2 * rB2; rC11+=rA2 * rB3;
							rC12+=rA3 * rB0; rC13+=rA3 * rB1; rC14+=rA3 * rB2; rC15+=rA3 * rB3;

							//LOAD(4)
							rA0 = A[iA+3]; rA1 = A[iiA+3]; rA2 = A[iiiA+3]; rA3 = A[iiiiA+3];
							rB0 = B[iiiiB]; rB1 = B[iiiiB + 1]; rB2 = B[iiiiB + 2]; rB1 = B[iiiiB + 3];
							//COMPUTE
							rC0+=rA0 * rB0; rC1+=rA0 * rB1; rC2+=rA0 * rB2; rC3+=rA0 * rB3;
							rC4+=rA1 * rB0; rC5+=rA1 * rB1; rC6+=rA1 * rB2; rC7+=rA1 * rB3;
							rC8+=rA2 * rB0; rC9+=rA2 * rB1; rC10+=rA2 * rB2; rC11+=rA2 * rB3;
							rC12+=rA3 * rB0; rC13+=rA3 * rB1; rC14+=rA3 * rB2; rC15+=rA3 * rB3;

						}
						C[iC] = rC0; C[iC+1] = rC1; C[iC+2] = rC2; C[iC+3] = rC3;
						C[iiC] = rC4; C[iiC+1] = rC5; C[iiC+2] = rC6; C[iiC+3] = rC7;
						C[iiiC] = rC8; C[iiiC+1] = rC9; C[iiiC+2] = rC10; C[iiiC+3] = rC11;
						C[iiiiC] = rC12; C[iiiiC+1] = rC13; C[iiiiC+2] = rC14; C[iiiiC+3] = rC15;
					}
				}

			}
		}
	}
}


#endif
