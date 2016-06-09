/*
 * template Array class source file
 */

#include "Array.h"

template <class T>
Array<T>::Array ( size_t xSize )
{
  CHECK_MSG ( xSize > 1, "Size of the simulation data should be more than 1" );
  v.erase( v.begin(), v.end() );
  v.resize( xSize, 0.0 );
  imax = xSize; jmax = 1; kmax = 1;
  size = imax;
  PROGRESS_MESG("An object of the Array class has been created");
}

template <class T>
Array<T>::Array( size_t xSize, size_t ySize )
{
  CHECK_MSG ( xSize > 1 && ySize > 1,
	      "Size of the simulation data should be at least 2x2" );
  v.erase( v.begin(), v.end() );
  v.resize( xSize*ySize, 0.0 );
  imax = xSize; jmax = ySize; kmax = 1;
  size = imax*jmax;
  PROGRESS_MESG("An object of the Array class has been created");
}

template <class T>
Array<T>::Array ( size_t xSize, size_t ySize, size_t zSize )
{
  CHECK_MSG ( xSize > 1 && ySize > 1 && zSize > 1,
		  "Size of the simulation data should be at least 2x2x2" );
  v.erase ( v.begin(), v.end() );
  v.resize( xSize*ySize*zSize, 0.0 );
  imax = xSize; jmax = ySize; kmax = zSize;
  size = imax*jmax*kmax;
  PROGRESS_MESG("An object of the Array class has been created");
}

template <class T>
Array<T>::Array ( const Array<T> & c)
{
  v    = c.getConstVector;
  imax = c.getSize(0);
  jmax = c.getSize(1);
  kmax = c.getSize(2);
  size = imax*jmax*kmax;
}

/*initialize the whole Array with a constant value */
template <class T>
void Array<T>::fill ( T value ) {
	// TODO
	// you might want to use std::fill() here
	std::fill( v.begin(), v.end(), value );
	PROGRESS_MESG("The _ARRAY_ object is populated with presented value");
}

template <class T>
size_t Array<T>::getSize( size_t dimension ) const
{
  if (dimension == 0) return imax; else
  if (dimension == 1) return jmax; else
  if (dimension == 2) return kmax; else
  if (dimension == 3) return size; else
    {
      /*Throw error*/
      PROGRESS_MESG("Invalid request, dimension can be 0, 1, 2, or 3");
      abort();
    }

  return -1;
}

template <class T>
const std::vector<T> Array<T>::getConstVector ()
{
  return v;
}

template <class T>
void Array<T>::write2File ( std::string filename )
{
  std::ofstream file( filename );
  for ( size_t i = 0; i < imax; ++i )
    for ( size_t j = 0; j < jmax; ++j )
      file << i << "\t" << j << "\t" << v[ i + j*imax ]<< "\n";
  file.close();
}
