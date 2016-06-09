

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

int main(int argc, const char * argv[])
{
//
    if(argc > 4 || argc < 1)
        {
            std :: cout <<" Abort";
            abort();
        }
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
        for (int j = 1; j <= Arow; j++){

            int k = 1;
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            double sum4 = 0;
            double sum5 = 0;
            double sum6 = 0;
            double sum7 = 0;
            double sum8 = 0;
            double unrolling = 0;

            for (; k <= Brow; k = k + 8) {
            	sum1 += A(j,k) * B(k,i);
            	sum2 += A(j,k+1) * B(k+1,i);
            	sum3 += A(j,k+2) * B(k+2,i);
            	sum4 += A(j,k+3) * B(k+3,i);
            	sum5 += A(j,k+4) * B(k+4,i);
            	sum6 += A(j,k+5) * B(k+5,i);
            	sum7 += A(j,k+6) * B(k+6,i);
            	sum8 += A(j,k+7) * B(k+7,i);
            }
            if(Brow%8 != 0)
            {
                for(;k <= Brow ;++k)
                    unrolling += A(j,k)*B(k,i);

            }
           C(j,i) = sum1 +sum2  +sum3 +sum4 +sum5 +sum7 +sum8 +sum6 + unrolling;
         }

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
