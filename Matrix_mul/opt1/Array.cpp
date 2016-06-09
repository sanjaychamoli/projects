#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "Array.h"

//Array constructor which takes file name [in which that array is stored]
Array::Array(const std::string &name){
	std::ifstream infile;
	infile.open(name.c_str());
	infile >> xSize;
	infile >> ySize;
	while(infile) { double temp; infile >> temp; data.push_back(temp); }
}

//Array constructor which takes Number of Rows and Columns as Arguments and sets to zero Array
Array::	Array(const unsigned x, const unsigned y) {
	xSize = x; ySize = y;
	unsigned total = x*y;
	for(unsigned i = 0; i < total; i++)
	data.push_back(0.0);
}

//Array constructor using another Array. Functionality of vector is used.
Array::Array(const Array& s) { *this = s; }

// Print the whole array [This functionality is not used much]
void Array::print() {
   //      -> the line with highest y-value printed first
	if (data.size() == 0)
		std::cout << "No data Available";
	else
		for(unsigned int i = 1; i <= xSize; i++) {
			double maxElementNumber = i * ySize;
			for(unsigned int j = (i-1)*ySize; j < maxElementNumber; j++)
			std::cout << data[j] << " "; std::cout << std::endl;
		}
}

//return size of Row or Column
int Array::getSize( int dimension ) const {
	if(dimension == 0)		// Row
		return xSize;
	else
		if(dimension == 1)	// Column
		return ySize;

	std :: cout << "invalid";
	exit(0);				// Exit if dimension other than 0 or 1 are passed
	return 0;
}

//return total size of the array
int Array::getSize() const {
	if(xSize > 0)
		if(ySize> 0)
			return xSize*ySize;
		return xSize;
	return 0;
}


//Operator*
Array Array:: operator * ( Array& B) {
	if(ySize == B.xSize) {
		Array answer(xSize, B.ySize);
		for(unsigned j = 0; j <ySize; j++) {
			std :: vector < double > temp(B.xSize);
			for(unsigned i = 0; i<B.xSize; i++)
				temp[i] = B.data[i*(B.ySize) + j];

			for(unsigned int i = 0 ; i<temp.size(); i++)
				std::cout << temp[i] << " ";
		}
		return answer;
	}
	else
	{
		exit(0);
		return *this;
	}
}

// Operator=
Array Array:: operator= (const Array& s){
	xSize = s.xSize;
	ySize = s.ySize;
	int number_of_elements = s.getSize();
	data.clear();
	if(xSize>0) {
		data.clear();
		for(int i = 0; i<number_of_elements; i++)
			data.push_back(s.data[i]);
	}
	return *this;
}

