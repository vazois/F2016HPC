#include "time/Time.h"
#include "tools/Utils.h"
#include "tools/ArgParser.h"

#include <cmath>

void initMatrix(double *&A, double *&B, unsigned int N){
	Utils<double> u;
	for(unsigned int i = 0 ; i < N*N;i++){
		A[i] = u.uni(0,1);
		B[i] = u.uni(0,1);
	}
}

void zeros(double *&C, unsigned int N){
	for(unsigned int i = 0 ; i < N*N;i++){
		C[i] = 0;
	}
}

void dgemm0(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int i = 0; i < N ;i++){
		for(unsigned int j = 0; j < N ;j++){
			for(unsigned int k = 0; k < N ;k++){
				C[i * N + j] += A[ i * N + k ] * B[ k * N + j ];
			}
		}
	}

}

void dgemm1(const double *A, const double *B, double *&C, unsigned int N){
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

void dgemm2(const double *A, const double *B, double *&C, unsigned int N){

	for(unsigned int i = 0; i < N ;i+=2){
		for(unsigned int j = 0; j < N ;j+=2){
			register int iC = i * N + j; register int iiC = iC + N;
			register double rC00 = 0.0f;
			register double rC01 = 0.0f;
			register double rC10 = 0.0f;
			register double rC11 = 0.0f;

			for(unsigned int k = 0; k < N ;k+=2){
				register int iA = i * N + k; register int iiA = iA + N;
				register int iB = k * N + j; register int iiB = iB + N;

				register double rA00 = A[ iA ]; register double rA01 = A[ iA + 1 ];
				register double rA10 = A[ iiA ]; register double rA11 = A[ iiA + 1 ];
				register double rB00 = B[ iB ]; register double rB01 = B[ iB + 1 ];
				register double rB10 = B[ iiB ]; register double rB11 = B[ iiB + 1 ];

				rC00 += rA00 * rB00; rC00+=rA01 * rB10;

				rC01 += rA00 * rB01; rC01+=rA01 * rB11;
				rC10 += rA10 * rB00; rC10+=rA11 * rB10;
				rC11 += rA10 * rB01; rC11+=rA11 * rB11;
			}
			C[ iC ] = rC00;
			C[ iC + 1 ] = rC01;
			C[ iiC ] = rC10;
			C[ iiC + 1 ] = rC11;
		}
	}
}

void dgemm3(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int i = 0; i < N ;i+=2){
		for(unsigned int j = 0; j < N ;j+=4){
			register int iC = i * N + j; register int iiC = iC + N;

			register double rC00 = 0.0f;
			register double rC01 = 0.0f;
			register double rC10 = 0.0f;
			register double rC11 = 0.0f;
			register double rC02 = 0.0f;
			register double rC12 = 0.0f;
			register double rC03 = 0.0f;
			register double rC13 = 0.0f;

			for(unsigned int k = 0; k < N ;k+=2){
				register int iA = i * N + k; register int iiA = iA + N;
				register int iB = k * N + j; register int iiB = iB + N;

				register double rA0 = A[iA]; register double rA1 = A[iiA];
				register double rB0 = B[iB]; register double rB1 = B[iB + 1]; register double rB2 = B[iB + 2]; register double rB3 = B[iB + 3];

				rC00 += rA0 * rB0;
				rC01 += rA0 * rB1;
				rC02 += rA0 * rB2;
				rC03 += rA0 * rB3;

				rC10 += rA1 * rB0;
				rC11 += rA1 * rB1;
				rC12 += rA1 * rB2;
				rC13 += rA1 * rB3;

				rA0 = A[iA+1]; rA1 = A[iiA+1];
				rB0 = B[iiB]; rB1 = B[iiB + 1]; rB2 = B[iiB + 2]; rB3 = B[iiB + 3];

				rC00 += rA0 * rB0;
				rC01 += rA0 * rB1;
				rC02 += rA0 * rB2;
				rC03 += rA0 * rB3;

				rC10 += rA1 * rB0;
				rC11 += rA1 * rB1;
				rC12 += rA1 * rB2;
				rC13 += rA1 * rB3;
			}

			C[ iC ] = rC00;
			C[ iC + 1 ] = rC01;
			C[ iC + 2 ] = rC02;
			C[ iC + 3 ] = rC03;

			C[ iiC ] = rC10;
			C[ iiC + 1 ] = rC11;
			C[ iiC + 2 ] = rC12;
			C[ iiC + 3 ] = rC13;
		}
	}
}

void dgemm3b(const double *A, const double *B, double *&C, unsigned int N){
	for(unsigned int i = 0; i < N ;i+=4){
		for(unsigned int j = 0; j < N ;j+=2){
			register int iC = i * N + j; register int iiC = iC + N;

			register double rC00 = 0.0f;
			register double rC01 = 0.0f;
			register double rC10 = 0.0f;
			register double rC11 = 0.0f;
			register double rC02 = 0.0f;
			register double rC12 = 0.0f;
			register double rC03 = 0.0f;
			register double rC13 = 0.0f;

			for(unsigned int k = 0; k < N ;k+=2){
				register int iA = i * N + k; register int iiA = iA + N;
				register int iB = k * N + j; register int iiB = iB + N;

				register double rA0 = A[iA]; register double rA1 = A[iiA]; register double rA2 = A[iiA + N]; register double rA3 = A[iiA + (N<<1)];
				register double rB0 = B[iB]; register double rB1 = B[iB + 1];

				rC00 += rA0 * rB0;
				rC01 += rA1 * rB0;
				rC02 += rA2 * rB0;
				rC03 += rA3 * rB0;

				rC10 += rA0 * rB1;
				rC11 += rA1 * rB1;
				rC12 += rA2 * rB1;
				rC13 += rA3 * rB1;

				rA0 = A[iA+1]; rA1 = A[iiA+1]; rA2 = A[iiA+1+N]; rA3 = A[iiA+1+(N<<1)];
				rB0 = B[iiB]; rB1 = B[iiB + 1];

				rC00 += rA0 * rB0;
				rC01 += rA1 * rB0;
				rC02 += rA2 * rB0;
				rC03 += rA3 * rB0;

				rC10 += rA0 * rB1;
				rC11 += rA1 * rB1;
				rC12 += rA2 * rB1;
				rC13 += rA3 * rB1;
			}

			C[ iC ] = rC00;
			C[ iiC ] = rC01;
			C[ iiC + N ] = rC02;
			C[ iiC + (N<<1) ] = rC03;

			C[ iC + 1 ] = rC10;
			C[ iiC + 1 ] = rC11;
			C[ iiC + 1 + N ] = rC12;
			C[ iiC + 1 + (N<<1) ] = rC13;
		}
	}
}

void cmpResults(double *A,double *B, double *C, double *D, unsigned int N, std::string a, std::string b){
	double diff = std::abs(C[0] - D[0]);
	double maxA = std::abs(C[0]);
	double maxB = std::abs(D[0]);

	for(unsigned int i = 0; i <N *N; i++){
		if(std::abs(C[i] - D[i]) > diff) diff = std::abs(C[i] - D[i]);
		if(std::abs(C[i]) > maxA) maxA = std::abs(A[i]);
		if(std::abs(D[i]) > maxB) maxB = std::abs(B[i]);
	}
	diff/=maxA*maxB;
	std::cout<<"maximum difference between "<<a << " and "<<b  <<" is " << diff <<std::endl;
}

int main(int argc, char **argv){
	ArgParser ap;
	ap.parseArgs(argc,argv);

	double *A,*B,*C, *D;
	int N = 0 ;


	if ( ap.count() == 0 ){
		std::cout<< "Please provide matrix size argument: dgemm -n=64" << std::endl;
		return 1;
	}else{
		if( !ap.exists("-n") ){
			std::cout<< "Execute the program using -n for the matrix size: ./dgemm -n=64" << std::endl;
			return 1;
		}
		N = ap.getInt("-n");
	}

	std::cout << "Matrix dimensions (" << N << "x" << N << ")" << std::endl;

	A = new double[N * N];
	B = new double[N * N];
	C = new double[N * N];
	D = new double[N * N];

	initMatrix(A,B,N);
	zeros(C,N);
	zeros(D,N);
	Time<secs> t;

	t.start();
	dgemm0(A,B,C,N);
	double d0 = t.lap("Elapsed time for dgemm0 in secs");

	t.start();
	dgemm1(A,B,D,N);
	double d1 = t.lap("Elapsed time for dgemm1 in secs");
	cmpResults(A,B,C,D,N,"dgemm0","dgemm1");

	zeros(D,N);
	t.start();
	dgemm2(A,B,D,N);
	double d2 = t.lap("Elapsed time for dgemm2 in secs");
	cmpResults(A,B,C,D,N,"dgemm0","dgemm2");

	zeros(D,N);
	t.start();
	dgemm3(A,B,D,N);
	double d3 = t.lap("Elapsed time for dgemm3 in secs");
	cmpResults(A,B,C,D,N,"dgemm0","dgemm3");

	zeros(D,N);
	t.start();
	dgemm3b(A,B,D,N);
	double d3b = t.lap("Elapsed time for dgemm3b in secs");
	cmpResults(A,B,C,D,N,"dgemm0","dgemm3b");


	uint64_t fop = N;
	fop *=fop*fop*2;
	std::cout << "Total number of double precision floating point operations: " << fop << std::endl;
	d0 = (fop / d0)/1000000000;
	std::cout << "Dgemm0 GFLOPS: " << d0 << std::endl;
	d1 = (fop / d1)/1000000000;
	std::cout << "Dgemm1 GFLOPS: " << d1 << std::endl;
	d2 = (fop / d2)/1000000000;
	std::cout << "Dgemm2 GFLOPS: " << d2 << std::endl;
	d3 = (fop / d3)/1000000000;
	std::cout << "Dgemm3 GFLOPS: " << d3 << std::endl;
	d3 = (fop / d3b)/1000000000;
	std::cout << "Dgemm3b GFLOPS: " << d3b << std::endl;

	delete A;
	delete B;
	delete C;
	delete D;


	return 0;

}
