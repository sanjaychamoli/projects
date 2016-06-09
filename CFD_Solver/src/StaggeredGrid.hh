#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include "Types.hh" // for enum directions
#include "Array.hh"
#include "FileReader.hh"

//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
   // Constructors to manually create staggered grid
   StaggeredGrid ( int xSize, int ySize, int zSize, real dx, real dy, real dz);

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid ( const FileReader & configuration );    


   // Getters / Setters for member variables
   Array<real> & p()    { return p_;    } // pressure field
   Array<real> & rhs()  { return rhs_;  } // right handside field
   Array<real> & u()    { return u_;    } // component of velocity in z-direction
   Array<real> & v()    { return v_;    } // component of velocity in z-direction
   Array<real> & w()    { return w_;    } // component of velocity in z-direction
   Array<real> & f()    { return f_;    }
   Array<real> & g()    { return g_;    }
   Array<real> & h()    { return h_;    }

   const Array<real> & p()    const { return p_;    }
   const Array<real> & rhs()  const { return rhs_;  }
   const Array<real> & u()    const { return u_;    }
   const Array<real> & v()    const { return v_;    }
   const Array<real> & w()    const { return w_;    }
   const Array<real> & f()    const { return f_;    }
   const Array<real> & g()    const { return g_;    }
   const Array<real> & h()    const { return h_;    }

   real dx() const { return dx_; }
   real dy() const { return dy_; }
   real dz() const { return dz_; }

   inline int xSize() const {return p_.getSize(0) - 2;}
   inline int ySize() const {return p_.getSize(1) - 2;}
   inline int zSize() const {return p_.getSize(2) - 2;}

   inline bool isFluid(const int x, const int y, const int z) {return isfluid_(x,y,z);}
   inline int getNumFluid();

   //wrapped axis
   inline real u( int x, int y, int z, Direction dir);
   inline real v( int x, int y, int z,Direction dir);
   inline real w( int x, int y, int z,Direction dir);
   inline real p( int x, int y, int z, Direction dir);

protected:
   Array<real> p_;
   Array<real> rhs_;
   Array<real> u_;
   Array<real> v_;
   Array<real> w_;
   Array<real> f_;
   Array<real> g_;
   Array<real> h_;
   Array<unsigned char> isfluid_;

   real xlength,ylength,zlength;
   real dx_;   // distance between two grid points in x direction
   real dy_;   // distance between two grid points in y direction
   real dz_;   // distance between two grid points in z direction

   void setCellToObstacle(int x, int y, int z);
};


inline int StaggeredGrid::getNumFluid()
{
   int numFluid_ = 0;
   for (int i=1; i<isfluid_.getSize(0)-1; i++) // along x-axis
       for (int j=1; j<isfluid_.getSize(1)-1; j++) //along y-axis
           for (int k=1; k<isfluid_.getSize(2)-1; k++) //along z-axis
                numFluid_+=isfluid_(i,j,k);
   return numFluid_;
}


//wrapped axis implmentation.... returning values without storing as per comment

inline real StaggeredGrid::u(int x, int y, int z, Direction dir)
{
    switch(dir)
    {
       case NORTH: if (isfluid_(x,y+1,z))
                      return u_(x,y+1,z);
                   else
                      return -u_(x,y,z);
       case SOUTH: if (isfluid_(x,y-1,z))
                      return u_(x,y-1,z);
                   else
                      return -u_(x,y,z);
       case EAST:  if (isfluid_(x+1,y,z))
                      return u(x+1,y,z,CENTER);
                   else
                      return 0.0;
       case WEST:  if (isfluid_(x-1,y,z))
                      return u_(x-1,y,z);
                   else
                      return 0.0;
       case UP:  if (isfluid_(x,y,z+1))
                      return u_(x,y,z+1);
                   else
                      return 0.0;
       case DOWN:  if (isfluid_(x,y,z-1))
                      return u_(x,y,z-1);
                   else
                      return 0.0;
       default:    if (isfluid_(x+1,y,z)) //covers center
                      return u_(x,y,z);
                   else
                      return 0.0;
    }
}

inline real StaggeredGrid::v(int x, int y, int z, Direction dir)
{
    switch(dir)
     {
        case NORTH: if (isfluid_(x,y+1,z))
                       return v(x,y+1,z,CENTER);
                    else
                       return 0.0;
        case SOUTH: if (isfluid_(x,y-1,z))
                       return v_(x,y-1,z);
                    else
                       return 0.0;
        case EAST:  if (isfluid_(x+1,y,z))
                       return v_(x+1,y,z);
                    else
                       return -v_(x,y,z);
        case WEST:  if (isfluid_(x-1,y,z))
                       return v_(x-1,y,z);
                    else
                       return -v_(x,y,z);
        case UP:  if (isfluid_(x,y,z+1))
                       return v_(x,y,z+1);
                    else
                       return 0.0;
        case DOWN:  if (isfluid_(x,y,z-1))
                       return v_(x,y,z-1);
                    else
                       return 0.0;
        default:    if (isfluid_(x,y+1,z)) //covers center
                       return v_(x,y,z);
                    else
                       return 0.0;
     }
}

inline real StaggeredGrid::w(int x, int y, int z, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (isfluid_(x,y+1,z))
                     return w_(x,y+1,z);
                  else
                     return 0.0;
      case SOUTH: if (isfluid_(x,y-1,z))
                     return w_(x,y-1,z);
                  else
                     return 0.0;
      case EAST:  if (isfluid_(x+1,y,z))
                     return w_(x+1,y,z);
                  else
                     return -w_(x,y,z);
      case WEST:  if (isfluid_(x-1,y,z))
                     return w_(x-1,y,z);
                  else
                     return -w_(x,y,z);
      case UP:  if (isfluid_(x,y,z+1))
                     return w(x,y,z+1,CENTER);
                  else
                     return 0.0;
      case DOWN:  if (isfluid_(x,y,z-1))
                     return w_(x,y,z-1);
                  else
                     return 0.0;
      default:    if (isfluid_(x,y,z+1))  //covers center
                     return w_(x,y,z);
                  else
                     return 0.0;
   }
}

inline real StaggeredGrid::p(int x, int y, int z, Direction dir)
{
    switch(dir)
    {
       case NORTH: if (isfluid_(x,y+1,z))
                      return p_(x,y+1,z);
                   else if (isfluid_(x+1,y+1,z))
                      return 0.5*(p_(x,y,z) + p_(x+1,y+1,z));
                   else if (isfluid_(x-1,y+1,z))
                      return 0.5*(p_(x,y,z) + p_(x-1,y+1,z));
                   else
                      return p_(x,y,z);
       case SOUTH: if (isfluid_(x,y-1,z))
                      return p_(x,y-1,z);
                   else if (isfluid_(x+1,y-1,z))
                      return 0.5*(p_(x,y,z) + p_(x+1,y-1,z));
                   else if (isfluid_(x-1,y-1,z))
                      return 0.5*(p_(x,y,z) + p_(x-1,y-1,z));
                   else
                      return p_(x,y,z);
       case EAST:  if (isfluid_(x+1,y,z))
                      return p_(x+1,y,z);
                   else if (isfluid_(x+1,y+1,z))
                      return 0.5*(p_(x,y,z) + p_(x+1,y+1,z));
                   else if (isfluid_(x+1,y-1,z))
                      return 0.5*(p_(x,y,z) + p_(x+1,y-1,z));
                   else
                      return p_(x,y,z);
       case WEST:  if (isfluid_(x-1,y,z))
                      return p_(x-1,y,z);
                   else if (isfluid_(x-1,y+1,z))
                      return 0.5*(p_(x,y,z) + p_(x-1,y+1,z));
                   else if (isfluid_(x-1,y-1,z))
                      return 0.5*(p_(x,y,z) + p_(x-1,y-1,z));
                   else
                      return p_(x,y,z);
       case UP:  if (isfluid_(x,y,z+1))
                      return p_(x,y,z+1);
                   else if (isfluid_(x+1,y+1,z))
                      return 0.5*(p_(x,y,z) + p_(x+1,y+1,z));
                   else if (isfluid_(x+1,y-1,z))
                      return 0.5*(p_(x,y,z) + p_(x+1,y-1,z));
                   else
                      return p_(x,y,z);
       case DOWN:  if (isfluid_(x,y,z-1))
                      return p_(x,y,z-1);
                   else if (isfluid_(x-1,y+1,z))
                      return 0.5*(p_(x,y,z) + p_(x-1,y+1,z));
                   else if (isfluid_(x-1,y-1,z))
                      return 0.5*(p_(x,y,z) + p_(x-1,y-1,z));
                   else
                      return p_(x,y,z);
       default:    return p_(x,y,z); //covers center
    }
}


#endif //STAGGERED_GRID_HH

