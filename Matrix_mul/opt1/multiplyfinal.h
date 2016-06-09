#ifndef MULTIPLY_H_
#define MULTIPLY_H_

#include <fstream>
#include "Array.h"

namespace mat {

	Array BlockedTransposeMethod(const Array& A, const Array& B)
	{
		//The product of Matrices is eventually stored in answer Array
        Array answer(A.getSize(0), B.getSize(1));

        //Displays size of A and B
        std :: cout << " \nSize = " << A.getSize(0)<< "x" << A.getSize(1)<< " " << B.getSize(0)<< "x"  << B.getSize(1) << "\n";

        //Blocking two columns of B for getting its transpose.
        //During this blocking, since those particular columns of B are in the cache,
        //we might as well multiply them with rows of A.
        //The following for loop covers even columns of B. In one loop iteration
        //we calculate two rows of the answer.
        if(B.getSize(1)%2 != 0 ) {
            std :: vector<double> temp;
            int j = B.getSize(1);
            int Brow = B.getSize(0);
            int Arow = A.getSize(0);
            for(int i = 1; i <= Brow; i++) {
                temp.push_back(B(i, j));
                answer(1, j) = answer(i, j) + A(1,i) * B(i, j);
            }
            for(int i = 2; i < Arow; i++)
                for(int  k = i ; k< Brow; k++)
                    answer(i,j) = answer(i, j) + A(i,k) * temp[k-1];
        }

        for(int j = 1; j < B.getSize(1); j = j + 2)
        {
                //Temp column vector
                std :: vector<double> temp;
                std :: vector<double> temp2;

                //Extracts columns from B and calculates for first row of A
                double sum1 = 0;
                double sum1a = 0;
                for(int i = 1; i<B.getSize(0) || i==B.getSize(0); i++) {
                        temp.push_back(B(i,j));
                        sum1 = sum1 + A(1,i) * B(i,j);
                        temp2.push_back(B(i,j+1));
                        sum1a = sum1a + A(1,i) * B(i,j+1);
                }
                answer(1,j) = sum1;
                answer(1,j+1) = sum1a;

                //Multiplies the remaining rows of A with temp column vector
                for(int i = 2; i< A.getSize(0) || i== A.getSize(0); i++) {
                    double sum2 = 0.0;
                    double sum21 = 0.0; {

                    int k = 1;
                    int location = A.getSize(1)*(i-1)-1;

                    for(; k < B.getSize(0) || k == B.getSize(0); ++k) {
                    	sum2 = sum2 + A(location + k)*temp[k-1];
                        sum21 = sum21 + A(location + k)*temp2[k-1];
                    }


                    answer(i,j) = sum2;
                    answer(i,j+1) = sum21; }

                    temp.clear();
                    temp2.clear();
                }
        }
      return answer;
    }
}

#endif



