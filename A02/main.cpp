//#include "tools/Utils.h"
#include "tools/ArgParser.h"
#include "time/CTime.h"
#include "dgemm_reg.h"
#include "dgemm_blocks.h"
#include "dgemm_hybrid.h"
#include <cmath>
#include <inttypes.h>


double randDouble(){
	srand(time(NULL));
	double X=((double)rand()/(double)RAND_MAX);

	return X;
}

void initMatrix(double *&A, double *&B, unsigned int N){
	for(unsigned int i = 0 ; i < N*N;i++){
		//if(i == 0) std::cout << randDouble() << std::endl;
		A[i] = randDouble();
		B[i] = randDouble();
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

	std::string cmp_result = diff == 0 ? "(PASSED)" : "(FAILED)";
	std::cout<<"maximum difference between "<<a << " and "<<b  <<" is: " << diff <<"! "<<cmp_result <<std::endl;
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

	//Start Benchmark
	start_clock();//(1)
	dgemm_reg_ijk(A,B,C,N);
	stop_clock();
	double dr_ijk = secf();

	start_clock();//(2)
	dgemm_reg_jik(A,B,D,N);
	stop_clock();
	double dr_jik = secf();
	cmpResults(A,B,C,D,N,"dgemm_reg_ijk","dgemm_reg_jik");

	start_clock();//(3)
	dgemm_reg_kij(A,B,D,N);
	stop_clock();
	double dr_kij = secf();
	cmpResults(A,B,C,D,N,"dgemm_reg_ijk","dgemm_reg_kij");

	start_clock();
	dgemm_reg_ikj(A,B,D,N);
	stop_clock();
	double dr_ikj = secf();
	cmpResults(A,B,C,D,N,"dgemm_reg_ikj","dgemm_reg_ikj");

	start_clock();//(5)
	dgemm_reg_jki(A,B,D,N);
	stop_clock();
	double dr_jki = secf();
	cmpResults(A,B,C,D,N,"dgemm_reg_jki","dgemm_reg_jki");

	start_clock();//(6)
	dgemm_reg_kji(A,B,D,N);
	stop_clock();
	double dr_kji = secf();
	cmpResults(A,B,C,D,N,"dgemm_reg_jki","dgemm_reg_kji");

	//Timing
	printTime(dr_ijk,"Elapsed time of dgemm_reg_ijk in secs: ");
	printTime(dr_jik,"Elapsed time of dgemm_reg_jik in secs: ");
	printTime(dr_kij,"Elapsed time of dgemm_reg_kij in secs: ");
	printTime(dr_ikj,"Elapsed time of dgemm_reg_ikj in secs: ");
	printTime(dr_jki,"Elapsed time of dgemm_reg_kji in secs: ");
	printTime(dr_kji,"Elapsed time of dgemm_reg_jki in secs: ");


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
	zeros(D,N);

	//////////////////////////////////////////////////////////////////
	//Start Benchmark for blocking versions
	start_clock();//(1)
	dgemm_blocks_ijk(A,B,C,N,Bz);
	stop_clock();
	double db_ijk = secf();

	start_clock();//(2)
	dgemm_blocks_jik(A,B,D,N,Bz);
	stop_clock();
	double db_jik = secf();
	cmpResults(A,B,C,D,N,"dgemm_blocks_ijk","dgemm_blocks_jik");

	start_clock();//(3)
	dgemm_blocks_ikj(A,B,D,N,Bz);
	stop_clock();
	double db_ikj = secf();
	cmpResults(A,B,C,D,N,"dgemm_blocks_ijk","dgemm_blocks_ikj");

	start_clock();//(4)
	dgemm_blocks_kij(A,B,D,N,Bz);
	stop_clock();
	double db_kij = secf();
	cmpResults(A,B,C,D,N,"dgemm_blocks_ijk","dgemm_blocks_kij");

	start_clock();//(5)
	dgemm_blocks_jki(A,B,D,N,Bz);
	stop_clock();
	double db_jki = secf();
	cmpResults(A,B,C,D,N,"dgemm_blocks_ijk","dgemm_blocks_jki");

	start_clock();//(5)
	dgemm_blocks_kji(A,B,D,N,Bz);
	stop_clock();
	double db_kji = secf();
	cmpResults(A,B,C,D,N,"dgemm_blocks_ijk","dgemm_blocks_kji");

	//Time for blocking versions
	printTime(db_ijk,"Elapsed time of dgemm_blocks_ijk in secs: ");
	printTime(db_jik,"Elapsed time of dgemm_blocks_jik in secs: ");
	printTime(db_ikj,"Elapsed time of dgemm_blocks_ikj in secs: ");
	printTime(db_kij,"Elapsed time of dgemm_blocks_kij in secs: ");
	printTime(db_jki,"Elapsed time of dgemm_blocks_jki in secs: ");
	printTime(db_kji,"Elapsed time of dgemm_blocks_kji in secs: ");

	//GFLOPS for blocking version
	db_ijk = (fop / db_ijk)/1000000000;
	std::cout << "Dgemm_blocks_ijk GFLOPS: " << db_ijk << std::endl;
	db_jik = (fop / db_jik)/1000000000;
	std::cout << "Dgemm_blocks_jik GFLOPS: " << db_jik << std::endl;
	db_ikj = (fop / db_ikj)/1000000000;
	std::cout << "Dgemm_blocks_ikj GFLOPS: " << db_ikj << std::endl;
	db_kij = (fop / db_kij)/1000000000;
	std::cout << "Dgemm_blocks_kij GFLOPS: " << db_kij << std::endl;
	db_jki = (fop / db_jki)/1000000000;
	std::cout << "Dgemm_blocks_jki GFLOPS: " << db_jki << std::endl;
	db_kji = (fop / db_kji)/1000000000;
	std::cout << "Dgemm_blocks_kji GFLOPS: " << db_kji << std::endl;

}

void test_hybrid_approach(unsigned int N, unsigned int Bz){
	double *A,*B,*C, *D;
	A = new double[N * N];
	B = new double[N * N];
	C = new double[N * N];
	D = new double[N * N];

	std::cout << "Testing hybrid version..." << std::endl;
	std::cout << "Matrix dimensions (" << N << "x" << N << ")" << std::endl;
	std::cout << "Block Size: (" << Bz <<"x" << Bz <<")"<<std::endl;

	uint64_t fop = N;
	fop *=fop*fop*2;
	std::cout << "Total number of double precision floating point operations: " << fop << std::endl;

	initMatrix(A,B,N);
	zeros(C,N);
	zeros(D,N);

	/////////////////////////////////////////////////////////////////
	//Hybrid
	start_clock();
	dgemm_hybrid(A,B,C,N,Bz);
	stop_clock();
	double db_hybrid = secf();

	start_clock();
	dgemm_hybrid2(A,B,D,N,Bz);
	stop_clock();
	double db_hybrid2 = secf();
	cmpResults(A,B,C,D,N,"dgemm_hybrid","dgemm_hybrid2");

	start_clock();
	dgemm_hybrid3(A,B,D,N,Bz);
	stop_clock();
	double db_hybrid3 = secf();
	cmpResults(A,B,C,D,N,"dgemm_hybrid","dgemm_hybrid3");

	//Timing for hybrid version
	printTime(db_hybrid,"Elapsed time of dgemm_hybrid in secs: ");
	printTime(db_hybrid2,"Elapsed time of dgemm_hybrid2 in secs: ");
	printTime(db_hybrid3,"Elapsed time of dgemm_hybrid2 in secs: ");

	//GFLOPS for hybrid version
	db_hybrid = (fop / db_hybrid)/1000000000;
	std::cout << "Dgemm_hybrid GFLOPS: " << db_hybrid << std::endl;
	db_hybrid2 = (fop / db_hybrid2)/1000000000;
	std::cout << "Dgemm_hybrid2 GFLOPS: " << db_hybrid2 << std::endl;
	db_hybrid3 = (fop / db_hybrid3)/1000000000;
	std::cout << "Dgemm_hybrid2 GFLOPS: " << db_hybrid3 << std::endl;

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
	}else if(MD==2){
		test_hybrid_approach(N,Bz);
	}


	return 0;

}
