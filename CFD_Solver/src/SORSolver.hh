#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include <cmath>



class SORSolver
{
private:
  //initalizing
  unsigned int imax, jmax, kmax, no_iteration, checkfrequency;
  real relaxfactor, epsilon;

  //setting boundary
  void setboundary(Array<real> & p);

  //calculating residual
  real residual(StaggeredGrid & grid);


public:
   // Constructor to manually create SORSolver
   SORSolver ( real relaxfactor_, real epsilon_,unsigned int no_iteration_,unsigned int checkfrequency_, int imax_, int jmax_, int kmax_);

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );


   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid & grid );
};
#endif //SOR_SOLVER_HH




