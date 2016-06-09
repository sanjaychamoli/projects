#ifndef MULTIPLY_H_
#define MULTIPLY_H_

#include <fstream>

#include "Array.h"

namespace mat {

Array OptimalTransposeMethod1024 (const Array& A, const Array& B) {
	//The product of Matrices is eventually stored in answer Array
	Array answer(A.getSize(0), B.getSize(1));

	//Displays size of A and B
	std :: cout << " \nSize = " << A.getSize(0)<< "x" << A.getSize(1)<< " " << B.getSize(0)<< "x"  << B.getSize(1) << "\n";

	//The following code is additionally executed when B is odd.
	// Please read algorithm described in the report.
	if(B.getSize(1)%2 != 0 ) {
            std :: vector<double> temp;
            int j = B.getSize(1);
            int Brow = B.getSize(0);
            int Arow = A.getSize(0);
            for(int i = 1; i<Brow || i==Brow; i++) {
                temp.push_back(B(i, j));
                answer(1, j) = answer(i, j) + A(1,i) * B(i, j);
            }
            for(int i = 2; i < Arow; i++)
                for(int  k = i ; k< Brow; k++)
                    answer(i,j) = answer(i, j) + A(i,k) * temp[k-1];
        }

		//Blocking two columns of B for getting its transpose.
		//During this blocking, since those particular columns of B are in the cache,
		//we might as well multiply them with rows of A.
		//The following for loop covers even columns of B. In one loop iteration
		//we calculate two rows of the answer.
		//Loop unrolling and Loop fusion is used.
        for(int j = 1; j < B.getSize(1); j = j + 2) {
        	std :: vector<double> temp;		//Created to store column of B
            std :: vector<double> temp2;	//Will store next column with respect to temp.

            //Extracts columns from B, multiplied with A's first row
            // First row of answer is calculated
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
            int unrollingFactor = 7;
            double sum2 = 0.0;	double sum21 = 0.0;
            double sum3 = 0.0;	double sum31 = 0.0;
            double sum4 = 0.0;	double sum41 = 0.0;
            double sum5 = 0.0;	double sum51 = 0.0;
            double sum6 = 0.0;	double sum61 = 0.0;
            double sum7 = 0.0;	double sum71 = 0.0;
            double sum8 = 0.0;	double sum81 = 0.0; {

            int k = 1;
            int location = A.getSize(1)*(i-1)/*+j*/-1;

            //Loop unrolling
            for(int unroll = B.getSize(0) - B.getSize(0)%unrollingFactor; k < unroll || k == unroll; k = k+unrollingFactor) {
               	sum2 	= sum2 	+ A(location + k)	*temp[k-1];
               	sum21 	= sum21 + A(location + k)	*temp2[k-1];
                sum3 	= sum3 	+ A(location + k+1)	*temp[k];
                sum31 	= sum31 + A(location + k+1)	*temp2[k];
                sum4 	= sum4 	+ A(location + k+2)	*temp[k+1];
                sum41 	= sum41 + A(location + k+2)	*temp2[k+1];
                sum5 	= sum5 	+ A(location + k+3)	*temp[k+2];
                sum51 	= sum51 + A(location + k+3)	*temp2[k+2];
                sum6 = sum6 + A(location + k+4)*temp[k+3];
                sum61 = sum61 + A(location + k+4)*temp2[k+3];
                sum7 = sum7 + A(location + k+5)*temp[k+4];
                sum71 = sum71 + A(location + k+5)*temp2[k+4];
                sum8 = sum8 + A(location + k+6)*temp[k+5];
                sum81 = sum81 + A(location + k+6)*temp2[k+5];
            }

            double unrollSum = 0.0;
            double unrollSum2 = 0.0;
            for(; k < B.getSize(0) || k==B.getSize(0); ++k) {
            	unrollSum 	= unrollSum 	+ A(i,k)*temp[k-1];
                unrollSum2 	= unrollSum2 	+ A(i,k)*temp2[k-1];
                }
                answer(i,j) = sum2 + sum3 + sum4 + sum5 + sum6+ sum7+ sum8 + unrollSum;
                answer(i,j+1) = sum21 + sum31 + sum41 + sum51 + sum61+ sum71+ sum81 + unrollSum2;
            }
            	temp.clear();
                temp2.clear();
            }

        }

      return answer;
}

Array OptimalTranspose2048(const Array& A, const Array& B) {
        //Final answer stored here
        Array answer(A.getSize(0), B.getSize(1));

        //Displays size of A and B
        std :: cout << " \nSize = " << A.getSize(0)<< "x" << A.getSize(1)<< " " << B.getSize(0)<< "x"  << B.getSize(1) << "\n";

        if(B.getSize(1)%2 != 0 ) {
        	std::cout << "odd";
        	std :: vector<double> temp;
            int j = B.getSize(1);
            int Brow = B.getSize(0);
            int Arow = A.getSize(0);
            for(int i = 1; i<Brow || i==Brow; i++) {
                temp.push_back(B(i, j));
                answer(1, j) = answer(i, j) + A(1,i) * B(i, j);
            }
            for(int i = 2; i < Arow; i++)
                for(int  k = i ; k< Brow; k++)
                    answer(i,j) = answer(i, j) + A(i,k) * temp[k-1];
        }

        for(int j = 1; j < B.getSize(1); ++j) {
        	//Temp column vector
        	std :: vector<double> temp;

        	//Extracts columns from B and calculates for first row of A
        	double sum1 = 0;
        	for(int i = 1; i<B.getSize(0) || i==B.getSize(0); i++) {
        		temp.push_back(B(i,j));
                sum1 += sum1 + A(1,i) * B(i,j);
        	}
        	answer(1,j) = sum1;

        	//Multiplies the remaining rows of A with temp column vector
        	for(int i = 1; i< A.getSize(0) || i== A.getSize(0); i++) {
        		int unrollingFactor = 7;
        		double sum2 = 0.0;
                double sum3 = 0.0;
                double sum4 = 0.0;
                double sum5 = 0.0;
                double sum6 = 0.0;
                double sum7 = 0.0;
                double sum8 = 0.0; {

                int k = 1;
                int location = A.getSize(1)*(i-1)-1;
                //Loop unrolling
                for(int unroll = B.getSize(0) - B.getSize(0)%unrollingFactor; k < unroll || k == unroll; k = k+unrollingFactor) {
                	sum2 += sum2 + A(location + k)	*temp[k-1];
                	sum3 += sum3 + A(location + k+1)*temp[k];
                	sum4 += sum4 + A(location + k+2)*temp[k+1];
                	sum5 += sum5 + A(location + k+3)*temp[k+2];
                	sum6 += sum6 + A(location + k+4)*temp[k+3];
                	sum7 += sum7 + A(location + k+5)*temp[k+4];
                	sum8 += sum8 + A(location + k+6)*temp[k+5];
                }

                double unrollSum = 0.0;
                for(; k < B.getSize(0) || k==B.getSize(0); ++k)
                	unrollSum += unrollSum + A(i,k)*temp[k-1];

                	answer(i,j) = sum2 + sum3 + sum4 + sum5 + sum6+ sum7+ sum8 + unrollSum;


                }
                        temp.clear();
                }

        }

      return answer;
}

}

#endif



