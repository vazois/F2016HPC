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
	std::cout<<"maximum difference between "<<a << " and "<<b  <<" is: " << diff <<"! "<<CMP <<std::endl;
}


void test_reg_approach(unsigned int N){
	double *A,*B,*C, *D;
	A = new double[N * N];
	B = new double[N * N];
	C = new double[N * N];
	D = new double[N * N];

	std::cout << "Testing register version..." <<std::endl;
	std::cout << "Matrix dimensions (" << N << "x" << N << ")" << std::endl;

	uint64_t fop = N;
	fop *=fop*fop*2;
	std::cout << "Total number of double precision floating point operations: " << fop << std::endl;

	initMatrix(A,B,N);
	zeros(C,N);
	Time<secs> t;

	//Start Benchmark
	t.start();//(1)
	dgemm_reg_ijk(A,B,C,N);
	double dr_ijk = t.lap("Elapsed time for dgemm_reg_ijk in secs");

	t.start();//(2)
	dgemm_reg_jik(A,B,D,N);
	double dr_jik = t.lap("Elapsed time for dgemm_reg_jik in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_ijk","dgemm_reg_jik");

	t.start();//(3)
	dgemm_reg_kij(A,B,D,N);
	double dr_kij = t.lap("Elapsed time for dgemm_reg_kij in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_ijk","dgemm_reg_kij");

	t.start();//(4)
	dgemm_reg_ikj(A,B,D,N);
	double dr_ikj = t.lap("Elapsed time for dgemm_reg_ikj in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_ikj","dgemm_reg_ikj");

	t.start();//(5)
	dgemm_reg_jki(A,B,D,N);
	double dr_jki = t.lap("Elapsed time for dgemm_reg_jki in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_jki","dgemm_reg_jki");

	t.start();//(6)
	dgemm_reg_kji(A,B,D,N);
	double dr_kji = t.lap("Elapsed time for dgemm_reg_kji in secs");
	cmpResults(A,B,C,D,N,"dgemm_reg_jki","dgemm_reg_kji");

	//GFLOPS
	dr_ijk = (fop / dr_ijk)/1000000000;
	std::cout << "Dgemm_reg_ijk GFLOPS: " << dr_ijk << std::endl;
	dr_jik = (fop / dr_jik)/1000000000;
	std::cout << "Dgemm_reg_jik GFLOPS: " << dr_jik << std::endl;
	dr_kij = (fop / dr_kij)/1000000000;
	std::cout << "Dgemm_reg_kij GFLOPS: " << dr_kij << std::endl;
	dr_ikj = (fop / dr_ikj)/1000000000;
	std::cout << "Dgemm_reg_ikj GFLOPS: " << dr_ikj << std::endl;
	dr_jki = (fop / dr_jki)/1000000000;
	std::cout << "Dgemm_reg_jki GFLOPS: " << dr_jki << std::endl;
	dr_kji = (fop / dr_kji)/1000000000;
	std::cout << "Dgemm_reg_kji GFLOPS: " << dr_kji << std::endl;

	delete A;
	delete B;
	delete C;
	delete D;
}


void test_blocked_approach(unsigned int N, unsigned int Bz){
	double *A,*B,*C, *D;
	A = new double[N * N];
	B = new double[N * N];
	C = new double[N * N];
	D = new double[N * N];

	std::cout << "Testing blocked version..." << std::endl;
	std::cout << "Matrix dimensions (" << N << "x" << N << ")" << std::endl;
	std::cout << "Block Size: (" << Bz <<"x" << Bz <<")"<<std::endl;

	uint64_t fop = N;
	fop *=fop*fop*2;
	std::cout << "Total number of double precision floating point operations: " << fop << std::endl;

	initMatrix(A,B,N);
	zeros(C,N);
	Time<secs> t;

	//Start Benchmark
	t.start();//(1)
	dgemm_blocks_ijk(A,B,C,N,Bz);
	double db_ijk = t.lap("Elapsed time for dgemm_blocks_ijk in secs");

	//GFLOPS
	db_ijk = (fop / db_ijk)/1000000000;
	std::cout << "Dgemm_blocks_ijk GFLOPS: " << db_ijk << std::endl;


	delete A;
	delete B;
	delete C;
	delete D;
}


int main(int argc, char **argv){
	ArgParser ap;
	ap.parseArgs(argc,argv);

	double *A,*B,*C, *D;
	int N = 0 ;
	int Bz = 0 ;
	int MD = 0;

	if ( ap.count() == 0 ){
		std::cout<< "Please provide matrix size and block size: example ./dgemm -md=0 -n=256 -b=16 " << std::endl;
		std::cout<< "-md : mode of execution, 0 register version, 1 blocked version, 2 hybrid version." << std::endl;
		std::cout<< "-n : matrix dimension " << std::endl;
		std::cout<< "-b : block size, required only for blocked and hybrid version" << std::endl;
		return 1;
	}else{

		if( !ap.exists("-md") ){
			std::cout << "Please provide mode of execution using -md. Run without arguments to get menu options!!!" << std::endl;
			return 1;
		}

		if( !ap.exists("-n") ){
			std::cout << "Please provide matrix dimension using -n. Run without arguments to get menu options!!!" << std::endl;
			return 1;
		}

		MD = ap.getInt("-md");
		N = ap.getInt("-n");

		if( (MD == 1  || MD == 2) &&  !ap.exists("-b") ){
			std::cout << "Please provide block size using -b. Run without arguments to get menu options!!!" << std::endl;
			return 1;
		}else if((MD == 1  || MD == 2) && ap.exists("-b")){
			Bz = ap.getInt("-b");
		}
	}

////////////////////////////////////////////////////////////////////////////////////////////

	if(MD==0){
		test_reg_approach(N);
	}else if(MD==1){
		test_blocked_approach(N,Bz);
	}


	return 0;

}
