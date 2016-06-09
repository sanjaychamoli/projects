
#include "Array.h"
#include <iostream>
#include <fstream>

//#include "exp.h"
#include "Timing.h"


#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif

int main(int argc, const char * argv[]) {

    if(argc > 4 || argc < 1) { std :: cout <<" Abort"; abort(); }
    std::string fileA(argv[1]);
    std::string fileB(argv[2]);
    std::string Cname(argv[3]);

	siwir::Timer T;

	Array A(fileA);
	Array B(fileB);
	Array C(A.getSize(0), B.getSize(1));

	int Arow = A.getSize(0);
	int Brow = B.getSize(0);
	int Bcol = B.getSize(1);

	T.reset();


#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "final" );
#endif


   for (int i = 1; i <= Bcol; i++)
   		for (int j = 1; j <= Arow; j++)
   			for (int k = 1; k <= Brow; k++)
   				C(j,i) += A(j,k) * B(k,i);


#ifdef USE_LIKWID
likwid_markerStopRegion( "final" );
#endif

#ifdef USE_LIKWID
likwid_markerClose();
#endif

	std :: cout <<"\nTime:  "<<T.elapsed();


    std::ofstream file;
        file.open(Cname.c_str());
        for (int i = 1; i< A.getSize(0)+1; i++)
            for (int j = 1; j< B.getSize(1) +1; j++){
            file << C(i,j) << std::endl;
        }
    file.close();

	return 0;
}
