#include <iostream>

#include "Array.h"
#include "multiplyfinal.h"
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

	Array A(fileA); Array B(fileB);

#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "final" );
#endif

    {
       siwir::Timer T;
       T.reset();
       Array C(mat::OptimalTransposeMethod1024(A,B));
       std :: cout <<"\nTime:  "<< T.elapsed();
       std::ofstream file;
       file.open(Cname.c_str());
       for (int i = 1; i< A.getSize(0)+1; i++)
    	   for (int j = 1; j< B.getSize(1) +1; j++)
               file << C(i,j) << std::endl;

       file.close();
    }


#ifdef USE_LIKWID
likwid_markerStopRegion( "final" );
#endif

#ifdef USE_LIKWID
likwid_markerClose();
#endif



        return 0;
}
