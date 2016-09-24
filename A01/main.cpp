#include "time/Time.h"
#include "tools/Utils.h"
#include "tools/ArgParser.h"

void initMatrix(double *&A, double *&B, double *&C, unsigned int N){
	Utils<double> u;
	for(unsigned int i = 0 ; i < N;i++){
		A[i] = u.uni(-2,5);
		B[i] = u.uni(-2,5);
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
			double rC = 0.0f;
			for(unsigned int k = 0; k < N ;k++){
				rC += A[ i * N + k ] * B[ k * N + j ];
			}
			C[ i * N + j ] = 0;
		}
	}

}

void dgemm2(const double *A, const double *B, double *&C, unsigned int N){

	for(unsigned int i = 0; i < N ;i+=2){
		for(unsigned int j = 0; j < N ;j+=2){
			double rC0 = 0.0f;
			double rC1 = 0.0f;
			double rC2 = 0.0f;
			double rC3 = 0.0f;

			for(unsigned int k = 0; k < N ;k+=2){
				rC += A[ i * N + k ] * B[ k * N + j ];
			}
			C[ i * N + j ] = 0;
		}
	}
}

int main(int argc, char **argv){
	ArgParser ap;
	ap.parseArgs(argc,argv);

	double *A,*B,*C;
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

	initMatrix(A,B,C,N);
	Time<secs> t;

	t.start();
	dgemm0(A,B,C,N);
	t.lap("Elapsed time for dgemm0 in secs");

	t.start();
	dgemm1(A,B,C,N);
	t.lap("Elapsed time for dgemm1 in secs");

	delete A;
	delete B;
	delete C;


	return 0;

}
