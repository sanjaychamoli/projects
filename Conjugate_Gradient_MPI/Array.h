/*
 * Template Array class header file
 */

#ifndef ARRAY_H_
#define ARRAY_H_


#include <cstdlib>
#include <fstream>
#include <string>
#include <typeinfo>

#include "Debug.h"
#include "Types.h"

template <class T>
class Array {

public:

  /* Class Constructors for 1D, 2D, and 3D Array */
  Array ( size_t xSize);
  Array ( size_t xSize, size_t ySize);
  Array ( size_t xSize, size_t ySize, size_t zSize);

  /* Copy Constructor for Array */
  Array ( const Array<T> & c);

  /* Access operators for 1D, 2D, and 3D Array */
  inline T &operator()( const size_t i );
  inline T &operator()( const size_t i, const size_t j);
  inline T &operator()( const size_t i, const size_t j, const size_t k);

  /* Access operators for 1D, 2D, and 3D Array */
  inline const T &operator()( const size_t i ) const;
  inline const T &operator()( const size_t i, const size_t j) const;
  inline const T &operator()( const size_t i, const size_t j, const size_t k) const;


  /* initialize the whole Array with a constant value */
  void fill( T value );

  /* vector getter const function */
  const vec_t<T> getConstVector ();

  /* return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
  other dimension values are not allowed */
  size_t getSize( size_t dimension ) const;

  /* Write to file */
  void write2File (std::string filename);

  /* Powerful routine, use with Care !!! */
  inline void deleteRowOfVec ( size_t rowNumber )
  {  v.erase( v.begin() + imax*rowNumber, v.begin() +imax*(rowNumber+1)-1 ); --jmax;}

private:
  vec_t<T>  v;
  size_t imax;
  size_t jmax;
  size_t kmax;
  size_t size;

};

template <class T>
inline T &Array<T>::operator()(const size_t i)
{
  /*Bounds checking error function*/
  CHECK_MSG( jmax==1 && kmax ==1, "Use this operator only for 1D simulation data");
  CHECK_MSG( i < imax && i >=0, "Operator can take values between[inclusive] 0 to imax-1");
  return v[ i ];
}

template <class T>
inline T &Array<T>::operator()(const size_t i, const size_t j)
{
  /*Bounds checking error function*/
  CHECK_MSG( kmax ==1, "Use this operator only for 2D simulation data");
  CHECK_MSG( i < imax && i >=0 && j <jmax && j >= 0,
		  "Operator can take values between[inclusive] 0 to imax-1 and 0 and jmax-1");
  return v[ i + imax*j ];
}

template <class T>
inline T &Array<T>::operator()(const size_t i, const size_t j, const size_t k)
{
  /*Bounds checking error function*/
  CHECK_MSG( i < imax && i >=0 && j <jmax && j >= 0 && k< kmax && k >= 0,
  		"Operator can take values between[inclusive] 0 to imax-1, 0 and jmax-1 and 0 to kmax-1");
  return v[ i + imax*j +k*jmax];
}

template <class T>
inline const T &Array<T>::operator()(const size_t i) const
{
  /*Bounds checking error function*/
  CHECK_MSG( jmax==1 && kmax ==1, "Use this operator only for 1D simulation data");
  CHECK_MSG( i < imax && i >=0, "Operator can take values between[inclusive] 0 to imax-1");
  return v[ i ];
}

template <class T>
inline const T &Array<T>::operator()(const size_t i, const size_t j) const
{
  /*Bounds checking error function*/
  CHECK_MSG( kmax ==1, "Use this operator only for 2D simulation data");
  CHECK_MSG( i < imax && i >=0 && j <jmax && j >= 0,
		  "Operator can take values between[inclusive] 0 to imax-1 and 0 and jmax-1");
  return v[ i + imax*j ];
}

template <class T>
inline const T &Array<T>::operator()(const size_t i, const size_t j, const size_t k) const
{
  /*Bounds checking error function*/
  CHECK_MSG( i < imax && i >=0 && j <jmax && j >= 0 && k< kmax && k >= 0,
    		"Operator can take values between[inclusive] 0 to imax-1, 0 and jmax-1 and 0 to kmax-1");
  return v[ i + imax*j +k*jmax];
}

#endif /* ARRAY_H_ */
