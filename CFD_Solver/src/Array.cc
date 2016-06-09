#include "Array.hh"


//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

// For 1D Array
template <typename T>
Array<T>::Array( int xSize ):data_read(xSize), datax_size(xSize),  data_dimension(0)
{}

// For 2D Array
template <typename T>
Array<T>::Array( int xSize, int ySize ):data_read(xSize*ySize), datax_size(xSize), datay_size(ySize),  data_dimension(1)
{}

// For 3D Array
template <typename T>
Array<T>::Array( int xSize, int ySize, int zSize ):data_read(xSize*ySize*zSize), datax_size(xSize), datay_size(ySize),
                                                dataz_size(zSize), data_dimension(2)
{}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
template <typename T>
void Array<T>::fill( T value )
{
   data_read.assign(data_read.size(),value);
}


// Print the whole array (for debugging purposes)
template <typename T>
void Array<T>::print()
{
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
   switch(data_dimension)
   {
      case 0: for (int i=0; i<datax_size; i++)
                 std::cout<< std::left<<data_read[i]<<"\t";

              break;

      case 1: for (int j=datay_size-1; j>=0; j--) {
                 for (int i=0; i<datax_size; i++)
                       std::cout<< std::left<<data_read[i*datay_size+j]<<"\t";
                       std::cout<<"\n";}
              break;
      case 2: for (int k=0; k<dataz_size; k++){
                 std::cout<<"\nz = "<< k << "\n";
                 for (int j=datay_size-1; j>=0; j--) {
                    for (int i=0; i<datax_size; i++)
                        std::cout<< std::left<<data_read[i*datay_size*dataz_size + j*dataz_size + k]<<"\t";
                    std::cout<<"\n"; } }
              break;
      default:break;
   }
}

template <typename T>
int Array<T>::getSize( int dimension ) const
{
   CHECK_MSG(dimension <= data_dimension, "Invalid request for dimension");
   switch(dimension)
   {
      case 0: return(datax_size);
      case 1: return(datay_size);
      case 2: return(dataz_size);
      default:return 0;
   }
}

//return total size of the array
template <typename T>
int Array<T>::getSize() const
{
   switch(data_dimension)
   {
      case 0: return(datax_size);
      case 1: return(datax_size*datay_size);
      case 2: return(datax_size*datay_size*dataz_size);
      default:return 0;
   }
}  


template class Array<unsigned char> ;
template class Array<real> ;
