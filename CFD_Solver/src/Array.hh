#ifndef ARRAY_HH
#define ARRAY_HH

#include<vector>
#include<iostream>
#include<iomanip>

#include "Types.hh"
#include "Debug.hh"


//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
template<typename T>
class Array
{
private:
    //for storing data
    std::vector<T> data_read;
    int datax_size,datay_size,dataz_size,data_dimension;
	
public:

   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );

   //destructor
   //~Array();

   // Depending on your implementation you might need the following:
   // ~Array();
   // Array(const Array& s);
   // Array& operator= (const Array& s);


   // Access Operators for 1D, 2D and 3D
   inline T & operator () ( int i );
   inline T & operator () ( int i ,int j );
   inline T & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
    inline const T & operator () ( int i ) const;
    inline const T & operator () ( int i ,int j ) const;
    inline const T & operator () ( int i, int j, int k ) const;



   // initialize the whole array with a constant value
   void fill( T value );


   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;


   // Print the whole array ( for debugging purposes )
   void print();

};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
template <typename T>
inline T& Array<T>::operator ()(int i)
{
	CHECK_MSG(0 == data_dimension, "Array is not 1D");
    CHECK_MSG(i <=datax_size, "Invalid access");
	return data_read[i];
}

// Operator() 2D
template <typename T>
inline T& Array<T>::operator ()(int i,int j)
{
	CHECK_MSG(1 == data_dimension, "Array is not 2D");
    CHECK_MSG(i<=datax_size && j<=datay_size, "Invalid access");
	return data_read[i*(datay_size) + j];
}

// Operator() 3D
template <typename T>
inline T& Array<T>::operator ()(int i, int j, int k)
{
	CHECK_MSG(2 == data_dimension, "Array is not 3D");
    CHECK_MSG(i<=datax_size && j<=datay_size && k<=dataz_size, "Invalid access");
	return data_read[i*datay_size*dataz_size + j*dataz_size + k];
}


// Operator() 1D
template <typename T>
inline const T& Array<T>::operator ()(int i) const
{
    CHECK_MSG(0 == data_dimension, "Array is not 1D");
    CHECK_MSG(i <=datax_size, "Invalid access");
    return data_read[i];
}

// Operator() 2D
template <typename T>
inline const T& Array<T>::operator ()(int i,int j) const
{
    CHECK_MSG(1 == data_dimension, "Array is not 2D");
    CHECK_MSG(i<=datax_size && j<=datay_size, "Invalid access");
    return data_read[i*(datay_size) + j];
}

// Operator() 3D
template <typename T>
inline const T& Array<T>::operator ()(int i, int j, int k) const
{
    CHECK_MSG(2 == data_dimension, "Array is not 3D");
    CHECK_MSG(i<=datax_size && j<=datay_size && k<=dataz_size, "Invalid access");
    return data_read[i*datay_size*dataz_size + j*dataz_size + k];
}


#endif //ARRAY_HH

