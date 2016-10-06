#include "time/Time.h"
#include "tools/Utils.h"
#include "tools/ArgParser.h"
#include "dgemm_reg.h"
#include "dgemm_blocks.h"
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

void cmpResults(double *A,double *B, double *&C, double *D, unsigned int N, std::string a, std::string b){
	double diff = std::abs(C[0] - D[0]);
	double maxA = std::abs(C[0]);
	double maxB = std::abs(D[0]);

	for(unsigned int i = 0; i <N *N; i++){
		if(std::abs(C[i] - D[i]) > diff) diff = std::abs(C[i] - D[i]);
		if(std::abs(C[i]) > maxA) maxA = std::abs(A[i]);
		if(std::abs(D[i]) > maxB) maxB = std::abs(B[i]);
		D[i]=0.0;
	}
	diff/=maxA*maxB;
	std::string CMP = diff == 0 ? "(PASSED)" : "(FAILED)";
	std::cout<<"maximum difference between "<<a << " and "<<b  <<" is " << diff <<" "<<CMP <<std::endl;
}



int main(int argc, char **argv){
	ArgParser ap;
	ap.parseArgs(argc,argv);

	double *A,*B,*C, *D;
	int N = 0 ;
	int Bz = 0 ;

	if ( ap.count() == 0 ){
		std::cout<< "Please provide matrix size and block size: example dgemm -n=256 -b=16 " << std::endl;
		return 1;
	}else{
		if( !ap.exists("-n") ){
			std::cout<< "Execute the program using -n for the matrix size: ./dgemm -n=256 -b=16" << std::endl;
			return 1;
		}

		if( !ap.exists("-b") ){
			std::cout<< "Execute the program using -b for the block size: ./dgemm -n=256 -b=16" << std::endl;
			return 1;
		}

		N = ap.getInt("-n");
		Bz = ap.getInt("-b");
	}

	std::cout << "Matrix dimensions (" << N << "x" << N << ")" << std::endl;

	A = new double[N * N];
	B = new double[N * N];
	C = new double[N * N];
	D = new double[N * N];

	uint64_t fop = N;
	fop *=fop*fop*2;
	std::cout << "Total number of double precision floating point operations: " << fop << std::endl;

	initMatrix(A,B,N);
	zeros(C,N);
	Time<secs> t;

	//Register Reuse MMUL
	t.start();//(1)
	dgemm_reg_ijk(A,B,C,N);
	double dr_ijk = t.lap("Elapsed time for dgemm_reg_ijk in secs");

	//zeros(D,N);
	t.start();//(2)
	dgemm_reg_jik(A,B,D,N);
	double dr_jik = t.lap("Elapsed time for dgemm_reg_jik in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_ijk","dgemm_reg_jik");

	//GFLOPS
	dr_ijk = (fop / dr_ijk)/1000000000;
	std::cout << "Dgemm_reg_ijk GFLOPS: " << dr_ijk << std::endl;
	dr_jik = (fop / dr_jik)/1000000000;
	std::cout << "Dgemm_reg_jik GFLOPS: " << dr_jik << std::endl;
////////////////////////////////////////////////////////////////////////////////////////////

	//Blocking MMUL
	//zeros(D,N);
	t.start();//(1)
	dgemm_blocks_ijk(A,B,D,N,Bz);
	double db_ijk = t.lap("Elapsed time for dgemm_blocks_ijk in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_ijk","dgemm_blocks_ijk");

	//GFLOPS
	db_ijk = (fop / db_ijk)/1000000000;
	std::cout << "Dgemm_blocks_ijk GFLOPS: " << db_ijk << std::endl;




	//cmpResults(A,B,C,D,N,"dgemm","dgemm3");






	delete A;
	delete B;
	delete C;
	delete D;


	return 0;

}
