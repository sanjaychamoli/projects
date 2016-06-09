/* * multiplyfinal.h
 *
 *  Created on: Oct 28, 2014
 *      Author: sanjay
 */


#ifndef MULTIPLY_H_
#define MULTIPLY_H_

#include <fstream>
#include "Array.h"

namespace mat
{

Array hello(const Array& A, const Array& B/*, std::string& Cname*/)
{
        //Final answer stored here
        Array answer(A.getSize(0), B.getSize(1));

        //Displays size of A and B
        std :: cout << " \nSize = " << A.getSize(0)<< "x" << A.getSize(1)<< " " << B.getSize(0)<< "x"  << B.getSize(1) << "\n";

        // int block  = 16
        //for(int b = 1; b <= block; ++b)
        //
        //Loop for columns of B
    if(B.getSize(1)%2 != 0 )
        {
        std::cout << "odd";
            std :: vector<double> temp;
            int j = B.getSize(1);
            int Brow = B.getSize(0);
            int Arow = A.getSize(0);
            for(int i = 1; i<Brow || i==Brow; i++)
            {
                temp.push_back(B(i, j));
                answer(1, j) = answer(i, j) + A(1,i) * B(i, j);
            }
            for(int i = 2; i < Arow; i++)
                for(int  k = i ; k< Brow; k++)
            {
                    answer(i,j) = answer(i, j) + A(i,k) * temp[k-1];
            }

        }

        for(int j = 1; j < B.getSize(1); j = j + 2)
        {
                //Temp column vector
                std :: vector<double> temp;
                std :: vector<double> temp2;
//              std :: vector<double> temp3;

//              std :: cout << "\n" << j;
                //Extracts columns from B and calculates for first row of A
                double sum1 = 0;
                double sum1a = 0;
//              double sum1b = 0;
                for(int i = 1; i<B.getSize(0) || i==B.getSize(0); i++)
                {

                        temp.push_back(B(i,j));
                        sum1 = sum1 + A(1,i) * B(i,j);
                        temp2.push_back(B(i,j+1));
                        sum1a = sum1a + A(1,i) * B(i,j+1);
//                      temp3.push_back(B(i,j+2));
//                      sum1b = sum1b + A(1,i) * B(i,j+2);
                }
                answer(1,j) = sum1;
                answer(1,j+1) = sum1a;
//              answer(1,j+2) = sum1b;

                //Multiplies the remaining rows of A with temp column vector
                for(int i = 2; i< A.getSize(0) || i== A.getSize(0); i++)
                {
            int unrollingFactor = 7;


            double sum2 = 0.0;
                        double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
                        double sum6 = 0.0;
                        double sum7 = 0.0;
                        double sum8 = 0.0;
//                        double sum9 = 0.0;
//                        double sum10 = 0.0;
//                        double sum11 = 0.0;
//            double sum12 = 0.0;

            double sum21 = 0.0;
                        double sum31 = 0.0;
            double sum41 = 0.0;
            double sum51 = 0.0;
            double sum61 = 0.0;
                        double sum71 = 0.0;
                        double sum81 = 0.0;
//                        double sum91 = 0.0;
//                        double sum101 = 0.0;
//            double sum111 = 0.0;
//            double sum121 = 0.0;

//            std :: vector <double> sum(unrollingFactor, 0.0);
//            std :: vector <double> sum2(unrollingFactor, 0.0);

                                {
                                int k = 1;
                                int location = A.getSize(1)*(i-1)/*+j*/-1;

                                //Loop unrolling
                                for(int unroll = B.getSize(0) - B.getSize(0)%unrollingFactor; k < unroll || k == unroll; k = k+unrollingFactor)
                                {
//                              for(; k < B.getSize(0) || k == B.getSize(0); ++k)
//                              {

                    sum2 = sum2 + A(location + k)*temp[k-1];
                                        sum21 = sum21 + A(location + k)*temp2[k-1];

                                        sum3 = sum3 + A(location + k+1)*temp[k];
                                        sum31 = sum31 + A(location + k+1)*temp2[k];

                                        sum4 = sum4 + A(location + k+2)*temp[k+1];
                    sum41 = sum41 + A(location + k+2)*temp2[k+1];


                    sum5 = sum5 + A(location + k+3)*temp[k+2];
                                        sum51 = sum51 + A(location + k+3)*temp2[k+2];

                                        sum6 = sum6 + A(location + k+4)*temp[k+3];
                                        sum61 = sum61 + A(location + k+4)*temp2[k+3];

                                        sum7 = sum7 + A(location + k+5)*temp[k+4];
                                        sum71 = sum71 + A(location + k+5)*temp2[k+4];

                                        sum8 = sum8 + A(location + k+6)*temp[k+5];
                                        sum81 = sum81 + A(location + k+6)*temp2[k+5];

/*                                        sum9 = sum9 + A(location + k+7)*temp[k+6];
                                        sum91 = sum91 + A(location + k+7)*temp2[k+6];

                                        sum10 = sum10 + A(location + k+8)*temp[k+7];
                                        sum101 = sum101 + A(location + k+8)*temp2[k+7];

                                        sum11 = sum11 + A(location + k+9)*temp[k+8];
                                        sum111 = sum111 + A(location + k+9)*temp2[k+8];

                                        sum12 = sum12 + A(location + k+10)*temp[k+9];
                                        sum121 = sum121 + A(location + k+10)*temp2[k+9];
*/


/*                      for(int m = 1; m<unrollingFactor; ++m)
                            {sum[m] = sum[m] + A(location + k + m)*temp[k-1];
                                sum2[m] = sum2[m] + A(location + k + m) * temp2[k-1];
                               }*/
                }


                                double unrollSum = 0.0;
                                double unrollSum2 = 0.0;
                                for(; k < B.getSize(0) || k==B.getSize(0); ++k)
                                {
                                        unrollSum = unrollSum + A(i,k)*temp[k-1];
                                        unrollSum2 = unrollSum2 + A(i,k)*temp2[k-1];
                }

 /*               for(int count = 0; count<unrollingFactor; ++count)
                {
                    answer(i,j) += sum[count];
                    answer(i, j+1) += sum2[count];
                }*/

                                //Without loop unrolling
//                              for(int k = 1; k < B.getSize(0) || k==B.getSize(0); ++k)
//                              {
//                                      sum2 = sum2 + A(i,k)*temp[k-1];
//                              }
//



                answer(i,j) = sum2 + sum3 + sum4 + sum5 + sum6+ sum7+ sum8 /*+ sum9 + sum10+ sum11+ sum12 */+ unrollSum;
                answer(i,j+1) = sum21 + sum31 + sum41 + sum51 + sum61+ sum71+ sum81 /*+ sum91+ sum101+ sum111+ sum121 */+ unrollSum2;


//                              sum.clear();
                        }
                        temp.clear();
                        temp2.clear();
//                      temp3.clear();
                }

        }

//      std::cout << "\n" << *answer(1024,1024);


//Stores answer
//    std::ofstream file;
//        file.open(Cname.c_str());
//        for (int i = 1; i< A.getSize(0)+1; i++)
//            for (int j = 1; j< B.getSize(1) +1; j++){
//            file << answer(i,j) << std::endl;
//        }
//    file.close();

//      answer.print();
      return answer;
}

}

#endif /* MULTIPLY_H_ */



