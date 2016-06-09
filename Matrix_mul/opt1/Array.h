
#ifndef ARRAY_H_
#define ARRAY_H_

#include <iostream>
#include <string>
#include <vector>


class Array
{
public:
	Array(const std::string &name);
	Array(const unsigned x, const unsigned j);
    Array(const Array& s);

    inline double& operator () ( int i ,int j );
    inline const double& operator () ( int i ,int j ) const;

    inline double& operator () ( int i );
    inline const double& operator () ( int i ) const;

   Array operator * ( Array& B);
   Array operator = (const Array& s);

   // Print the whole array [This functionality is not used much]
   // Was written for small Array debugging
   void print();

   // Return total size of the array
   int getSize() const;

   // Return xSize for dimension==0 and ySize for dimension==1
   // Other dimension values are not allowed
   int getSize(int dimension ) const;


private:
	unsigned int xSize;
	unsigned int ySize;
	std::vector < double > data;
};

// Operator() 2D
inline double& Array::operator ()(int i,int j) {
	int location = ySize*(i-1)+j-1;
	return data[location];
}

// Operator() 2D const
inline const double& Array::operator ()(int i,int j) const {
	int location = ySize*(i-1)+j-1;
	return data[location];
}

inline double& Array:: operator () ( int i ) { return data[i]; }

inline const double& Array:: operator () ( int i ) const { return data[i]; }

#endif /* ARRAY_H_ */
