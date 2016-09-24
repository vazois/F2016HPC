#include "time/Time.h"
#include "tools/Utils.h"
#include "tools/ArgParser.h"

void initMatrix(double *&A, double *&B, unsigned int N){
	Utils<double> u;
	for(unsigned int i = 0 ; i < N;i++){
		A[i] = u.uni(-2,5);
		B[i] = u.uni(-2,5);
	}
}

void dgemm0(const double *A, const double *B, double *&B, unsigned int N){


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

	initMatrix(A,B,N);

	for(int i =0;i<10;i++){
		std::cout << A[i] <<"x" << B[i] << std::endl;
	}

	delete A;
	delete B;
	delete C;


	return 0;

}
