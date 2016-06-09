#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "SORSolver.hh"
#include "StaggeredGrid.hh"

class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );


      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid();
      const StaggeredGrid & grid() const;


  private:
      std::string name;
      StaggeredGrid discretization_;
      SORSolver solver_;
      real gamma, Re,tau;
      real gx, gy, gz;
      unsigned int nrofTimesteps, normalizationfrequency, currentstepnumber;
      real currenttime, duration, dt, dtmax, dx, dy, dz;
      unsigned int imax, jmax, kmax, outputinterval;
      real U_init, V_init, W_init, P_init;

      std::string boundary_condition_N, boundary_condition_S, boundary_condition_E, boundary_condition_W, boundary_condition_U, boundary_condition_D;
      real boundary_velocity_N, boundary_velocity_S, boundary_velocity_E, boundary_velocity_W, boundary_velocity_U, boundary_velocity_D;

      Array<real> &u = discretization_.u();
      Array<real> &v = discretization_.v();
      Array<real> &w = discretization_.w();
      Array<real> &f = discretization_.f();
      Array<real> &g = discretization_.g();
      Array<real> &h = discretization_.h();
      Array<real> &p = discretization_.p();
      Array<real> &rhs = discretization_.rhs();

      //computing derivatives
      inline real du2_dx    (int i, int j, int k);
      inline real dv2_dy    (int i, int j, int k);
      inline real dw2_dz    (int i, int j, int k);
      inline real duv_dy    (int i, int j, int k);
      inline real duw_dz    (int i, int j, int k);
      inline real dvu_dx    (int i, int j, int k);
      inline real dvw_dz    (int i, int j, int k);
      inline real dwu_dx    (int i, int j, int k);
      inline real dwv_dy    (int i, int j, int k);
      inline real d2u_dx2   (int i, int j, int k);
      inline real d2u_dy2   (int i, int j, int k);
      inline real d2u_dz2   (int i, int j, int k);
      inline real d2v_dx2   (int i, int j, int k);
      inline real d2v_dy2   (int i, int j, int k);
      inline real d2v_dz2   (int i, int j, int k);
      inline real d2w_dx2   (int i, int j, int k);
      inline real d2w_dy2   (int i, int j, int k);
      inline real d2w_dz2   (int i, int j, int k);

      void computeFGH();
      void composeRHS();
      void updateVelocities();
      void determineNextDT(real const & limit);
      void refreshBoundaries();
      void normalizepressure();

};

inline StaggeredGrid & FluidSimulator::grid()
{
   return discretization_;
}

inline const StaggeredGrid & FluidSimulator::grid() const
{
   return discretization_;
}


/////////////////////////////////////////////////////////////////derivatives function////////////////////////////////////////////////////////////////

inline real FluidSimulator::du2_dx  (int i, int j, int k)
{
return (0.25/dx*(pow(discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,EAST),2) - pow(discretization_.u(i,j,k,WEST)+discretization_.u(i,j,k,CENTER),2)
                    +gamma*( std::abs(discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,EAST)) * (discretization_.u(i,j,k,CENTER)-discretization_.u(i,j,k,EAST))
                            -std::abs(discretization_.u(i,j,k,WEST)+discretization_.u(i,j,k,CENTER)) * (discretization_.u(i,j,k,WEST)-discretization_.u(i,j,k,CENTER)))));
}
inline real FluidSimulator::dv2_dy  (int i, int j, int k)
{
return (0.25/dy*(pow(discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,NORTH),2) - pow(discretization_.v(i,j,k,SOUTH)+discretization_.v(i,j,k,CENTER),2)
                    +gamma*( std::abs(discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,NORTH)) * (discretization_.v(i,j,k,CENTER)-discretization_.v(i,j,k,NORTH))
                            -std::abs(discretization_.v(i,j,k,SOUTH)+discretization_.v(i,j,k,CENTER)) * (discretization_.v(i,j,k,SOUTH)-discretization_.v(i,j,k,CENTER)))));
}
inline real FluidSimulator::dw2_dz  (int i, int j, int k)
{
return (0.25/dy*(pow(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,UP),2) - pow(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,DOWN),2)
                     + gamma*( std::abs(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,UP)) * (discretization_.w(i,j,k,CENTER)-discretization_.w(i,j,k,UP))
                              -std::abs(discretization_.w(i,j,k,DOWN)+discretization_.w(i,j,k,CENTER)) * (discretization_.w(i,j,k,CENTER)-discretization_.w(i,j,k,DOWN)))));
}
inline real FluidSimulator::duv_dy  (int i, int j, int k)
{
return (0.25/dy* ( (discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,EAST))*(discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,NORTH))
                  -(discretization_.v(i,j,k,SOUTH)+discretization_.v(i+1,j,k,SOUTH))*(discretization_.u(i,j,k,SOUTH)+discretization_.u(i,j,k,CENTER))
                  +gamma * (  std::abs(discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,EAST)) * (discretization_.u(i,j,k,CENTER)-discretization_.u(i,j,k,NORTH))
                            - std::abs(discretization_.v(i,j,k,SOUTH)+discretization_.v(i+1,j,k,SOUTH)) * (discretization_.u(i,j,k,SOUTH)-discretization_.u(i,j,k,CENTER)))));
}
inline real FluidSimulator::duw_dz  (int i, int j, int k)
{
return (0.25/dz* ( (discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,EAST))*(discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,UP))
                  -(discretization_.w(i,j,k,DOWN)+discretization_.w(i+1,j,k,DOWN))*(discretization_.u(i,j,k,DOWN)+discretization_.u(i,j,k,CENTER))
                  +gamma * (  std::abs(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,EAST)) * (discretization_.u(i,j,k,CENTER)-discretization_.u(i,j,k,UP))
                            - std::abs(discretization_.w(i,j,k,DOWN)+discretization_.w(i+1,j,k,DOWN)) * (discretization_.u(i,j,k,DOWN)-discretization_.u(i,j,k,CENTER)))));
}
inline real FluidSimulator::dvu_dx  (int i, int j, int k)
{
return (0.25/dx * (  (discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,NORTH))*(discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,EAST))
                   - (discretization_.u(i,j,k,WEST)+discretization_.u(i,j+1,k,WEST))*(discretization_.v(i,j,k,WEST)+discretization_.v(i,j,k,CENTER))
                   + gamma * (  std::abs(discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,NORTH)) * (discretization_.v(i,j,k,CENTER)-discretization_.v(i,j,k,EAST))
                              - std::abs(discretization_.u(i,j,k,WEST)+discretization_.u(i,j+1,k,WEST)) * (discretization_.v(i,j,k,WEST)-discretization_.v(i,j,k,CENTER)))));
}
inline real FluidSimulator::dvw_dz  (int i, int j, int k)
{
return (0.25/dz * (  (discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,NORTH))*(discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,UP))
                   - (discretization_.w(i,j,k,DOWN)+discretization_.w(i,j+1,k,DOWN))*(discretization_.v(i,j,k,DOWN)+discretization_.v(i,j,k,CENTER))
                   + gamma * (  std::abs(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,NORTH)) * (discretization_.v(i,j,k,CENTER)-discretization_.v(i,j,k,UP))
                              - std::abs(discretization_.w(i,j,k,DOWN)+discretization_.w(i,j+1,k,DOWN)) * (discretization_.v(i,j,k,DOWN)-discretization_.v(i,j,k,CENTER)))));
}
inline real FluidSimulator::dwu_dx  (int i, int j, int k)
{
return (0.25/dx * (  (discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,UP))*(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,EAST))
                   - (discretization_.u(i,j,k,WEST)+discretization_.u(i,j,k+1,WEST))*(discretization_.w(i,j,k,WEST)+discretization_.w(i,j,k,CENTER))
                   + gamma * (  std::abs(discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,UP)) * (discretization_.w(i,j,k,CENTER)-discretization_.w(i,j,k,EAST))
                              - std::abs(discretization_.u(i,j,k,WEST)+discretization_.u(i,j,k+1,WEST)) * (discretization_.w(i,j,k,WEST)-discretization_.w(i,j,k,CENTER)))));
}
inline real FluidSimulator::dwv_dy  (int i, int j, int k)
{
return (0.25/dy* ( (discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,UP))*(discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,NORTH))
                  -(discretization_.v(i,j,k,SOUTH)+discretization_.v(i,j,k+1,SOUTH))*(discretization_.w(i,j,k,SOUTH)+discretization_.w(i,j,k,CENTER))
                  +gamma * (  std::abs(discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,UP)) * (discretization_.w(i,j,k,CENTER)-discretization_.w(i,j,k,NORTH))
                            - std::abs(discretization_.v(i,j,k,SOUTH)+discretization_.v(i,j,k+1,SOUTH)) * (discretization_.w(i,j,k,SOUTH)-discretization_.w(i,j,k,CENTER)))));
}
inline real FluidSimulator::d2u_dx2 (int i, int j, int k)
{
   return (( discretization_.u(i,j,k,EAST)-2.0*discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,WEST) )/(dx*dx) );
}
inline real FluidSimulator::d2u_dy2 (int i, int j, int k)
{
   return (( discretization_.u(i,j,k,NORTH)-2.0*discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,SOUTH) )/(dy*dy) );
}
inline real FluidSimulator::d2u_dz2 (int i, int j, int k)
{
   return (( discretization_.u(i,j,k,UP)-2.0*discretization_.u(i,j,k,CENTER)+discretization_.u(i,j,k,DOWN) )/(dz*dz) );
}
inline real FluidSimulator::d2v_dx2 (int i, int j, int k)
{
   return (( discretization_.v(i,j,k,EAST)-2.0*discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,WEST) ) /(dx*dx) );
}
inline real FluidSimulator::d2v_dy2 (int i, int j, int k)
{
   return (( discretization_.v(i,j,k,NORTH)-2.0*discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,SOUTH) )/(dy*dy) );
}
inline real FluidSimulator::d2v_dz2 (int i, int j, int k)
{
   return (( discretization_.v(i,j,k,UP)-2.0*discretization_.v(i,j,k,CENTER)+discretization_.v(i,j,k,DOWN) )/(dz*dz) );
}
inline real FluidSimulator::d2w_dx2 (int i, int j, int k)
{
   return (( discretization_.w(i,j,k,EAST)-2.0*discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,WEST) ) /(dx*dx) );
}
inline real FluidSimulator::d2w_dy2 (int i, int j, int k)
{
   return (( discretization_.w(i,j,k,NORTH)-2.0*discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,SOUTH) )/(dy*dy) );
}
inline real FluidSimulator::d2w_dz2 (int i, int j, int k)
{
   return (( discretization_.w(i,j,k,UP)-2.0*discretization_.w(i,j,k,CENTER)+discretization_.w(i,j,k,DOWN) )/(dz*dz) );
}

#endif
