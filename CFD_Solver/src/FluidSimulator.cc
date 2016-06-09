#include "FluidSimulator.hh"
#include "VTKWriter.hh"
#include <random>
#include <cmath>
#include <algorithm>

FluidSimulator::FluidSimulator( const FileReader & input_data ) : 
                      name(input_data.getStringParameter("name")), discretization_(StaggeredGrid(input_data)),
					  solver_(SORSolver(input_data)), gamma(input_data.getRealParameter("gamma")),
                      Re(input_data.getRealParameter("Re")), tau(input_data.getRealParameter("safetyfactor")),
                      gx(input_data.getRealParameter("GX")), gy(input_data.getRealParameter("GY")), gz(input_data.getRealParameter("GZ")),
                      nrofTimesteps(input_data.getIntParameter("timesteps")),
                      normalizationfrequency(input_data.getIntParameter("normalizationfrequency")),
                      currentstepnumber(1),currenttime(0.0),
                      duration(0.0),dt(input_data.getRealParameter("dt")),dtmax(input_data.getRealParameter("dt")),
                      dx(discretization_.dx()), dy(discretization_.dy()),dz(discretization_.dz()),
                      imax(input_data.getIntParameter("imax")), jmax(input_data.getIntParameter("jmax")), kmax(input_data.getIntParameter("kmax")),
                      outputinterval(1),
                      boundary_condition_N("noslip"), boundary_condition_S("noslip"),
                      boundary_condition_E("noslip"), boundary_condition_W("noslip"),
                      boundary_condition_U("noslip"), boundary_condition_D("noslip"),
                      boundary_velocity_N(0), boundary_velocity_S(0), boundary_velocity_E(0),
                      boundary_velocity_W(0),boundary_velocity_U(0), boundary_velocity_D(0)
{
    U_init = input_data.getRealParameter("U_INIT");
    V_init = input_data.getRealParameter("V_INIT");
    W_init = input_data.getRealParameter("W_INIT");
    P_init = input_data.getRealParameter("P_INIT");

    CHECK_MSG((Re > 0), "Reynolds number (Re) should be other than 0.");
    CHECK_MSG((dt > 0), "dt should me greater than 0.");
    CHECK_MSG((tau > 0 && tau <= 1), "safetyfactor is not in range of 0 to 1.");
    CHECK_MSG(input_data.checkparameter("outputinterval"), "Output Interval not defined exclusively. Default value taken 1");
    outputinterval = input_data.getIntParameter("outputinterval");

    //boundary North
    if (input_data.checkparameter("boundary_condition_N")) {
       boundary_condition_N = input_data.getStringParameter("boundary_condition_N");
       CHECK_MSG ( boundary_condition_N != "outflow" || !(input_data.checkparameter("boundary_velocity_N")),
                   "Program aborted...outflow condition at NOTRH boundary cannot have velocity" );
       CHECK_MSG ( boundary_condition_N != "inflow" || (input_data.checkparameter("boundary_velocity_N")),
                   "Program aborted...inflow condition at NOTRH boundary must have velocity");
       }
    else
        boundary_condition_N = "noslip";

    if (input_data.checkparameter("boundary_velocity_N"))
       boundary_velocity_N = input_data.getRealParameter("boundary_velocity_N");

    //boundary South
    if (input_data.checkparameter("boundary_condition_S")) {
       boundary_condition_S = input_data.getStringParameter("boundary_condition_S");
       CHECK_MSG ( boundary_condition_S != "outflow" || !(input_data.checkparameter("boundary_velocity_S")),
                  "Program aborted...outflow condition at SOUTH boundary cannot have velocity" );
       CHECK_MSG ( boundary_condition_S != "inflow" || (input_data.checkparameter("boundary_velocity_S")),
                  "Program aborted...inflow condition at SOUTH boundary must have velocity" );
       }
    else
       boundary_condition_S = "noslip";

    if (input_data.checkparameter("boundary_velocity_S"))
       boundary_velocity_S = input_data.getRealParameter("boundary_velocity_S");

    //boudary EAST
    if (input_data.checkparameter("boundary_condition_E")) {
       boundary_condition_E = input_data.getStringParameter("boundary_condition_E");
       CHECK_MSG ( boundary_condition_E != "outflow" || !(input_data.checkparameter("boundary_velocity_E")),
                  "Program aborted...outflow condition at EAST boundary cannot have velocity" );
       CHECK_MSG ( boundary_condition_E != "inflow" || (input_data.checkparameter("boundary_velocity_E")),
                   "Program aborted...inflow condition at EAST boundary must have velocity");
       }
    else
       boundary_condition_E = "noslip";

    if (input_data.checkparameter("boundary_velocity_E"))
       boundary_velocity_E = input_data.getRealParameter("boundary_velocity_E");

    //boundary WEST
    if (input_data.checkparameter("boundary_condition_W")) {
       boundary_condition_W = input_data.getStringParameter("boundary_condition_W");
       CHECK_MSG ( boundary_condition_W != "outflow" || !(input_data.checkparameter("boundary_velocity_W")),
                   "Program aborted...outflow condition at WEST boundary cannot have velocity");
       CHECK_MSG ( boundary_condition_W != "inflow" || (input_data.checkparameter("boundary_velocity_W")),
                   "Program aborted...inflow condition at WEST boundary must have velocity" );
       }
    else
       boundary_condition_W = "noslip";

    if (input_data.checkparameter("boundary_velocity_W"))
       boundary_velocity_W = input_data.getRealParameter("boundary_velocity_W");

    //boundary UP
    if (input_data.checkparameter("boundary_condition_U")) {
       boundary_condition_U = input_data.getStringParameter("boundary_condition_U");
       CHECK_MSG ( boundary_condition_U!= "outflow" || !(input_data.checkparameter("boundary_velocity_U")),
                   "Program aborted...outflow condition at UP boundary cannot have velocity" );
       CHECK_MSG ( boundary_condition_U != "inflow" || (input_data.checkparameter("boundary_velocity_U")),
                   "Program aborted...inflow condition at UP boundary must have velocity");
       }
    else
        boundary_condition_U = "noslip";

    if (input_data.checkparameter("boundary_velocity_U"))
       boundary_velocity_U = input_data.getRealParameter("boundary_velocity_U");

    //boundary DOWN
    if (input_data.checkparameter("boundary_condition_D")) {
       boundary_condition_D = input_data.getStringParameter("boundary_condition_D");
       CHECK_MSG ( boundary_condition_D != "outflow" || !(input_data.checkparameter("boundary_velocity_D")),
                   "Program aborted...outflow condition at DOWN boundary cannot have velocity" );
       CHECK_MSG ( boundary_condition_D != "inflow" || (input_data.checkparameter("boundary_velocity_D")),
                   "Program aborted...inflow condition at DOWN boundary must have velocity");
       }
    else
        boundary_condition_D = "noslip";

    if (input_data.checkparameter("boundary_velocity_D"))
       boundary_velocity_D = input_data.getRealParameter("boundary_velocity_D");
}



////////////////////////////////////////////////////////////functions for simulation//////////////////////////////////////////////////////////

//computing FG
void FluidSimulator::computeFGH()
{   
   //for grid points
   for (unsigned int i=1; i<imax; i++)
       for (unsigned int j=1; j<jmax+1; j++)
           for (unsigned int k=1; k<kmax+1; k++)
               if (discretization_.isFluid(i,j,k))
                    f(i,j,k) = u(i,j,k) + dt * ( ((d2u_dx2(i,j,k) + d2u_dy2(i,j,k) + d2u_dz2(i,j,k))/Re)
                                                                         -du2_dx(i,j,k) - duv_dy(i,j,k) - duw_dz(i,j,k) + gx);
   for (unsigned int i=1; i<imax+1; i++)
       for (unsigned int j=1; j<jmax; j++)
           for (unsigned int k=1;k<kmax+1;k++)
                if (discretization_.isFluid(i,j,k))
                    g(i,j,k) = v(i,j,k) + dt * ( ((d2v_dx2(i,j,k) + d2v_dy2(i,j,k) + d2v_dz2(i,j,k))/Re)
                                                                        - dv2_dy(i,j,k) - dvu_dx(i,j,k) - dvw_dz(i,j,k) + gy);

   for (unsigned int i=1; i<imax+1; i++)
       for (unsigned int j=1; j<jmax+1; j++)
           for (unsigned int k=1;k<kmax;k++)
                if (discretization_.isFluid(i,j,k))
                    h(i,j,k) = w(i,j,k) + dt * ( ((d2w_dx2(i,j,k) + d2w_dy2(i,j,k) + d2w_dz2(i,j,k))/Re)
                                                                        - dw2_dz(i,j,k) - dwu_dx(i,j,k) - dwv_dy(i,j,k) + gz);
   //for boundary
   for (unsigned int j=1; j<jmax+1; j++)
       for (unsigned int k=1; k<kmax+1; k++)
       {
           if (discretization_.isFluid(0,j,k))
           f(0,j,k) =u(0,j,k);
           if (discretization_.isFluid(imax,j,k))
           f(imax,j,k) =u(imax,j,k);
       }
   for (unsigned int i=1; i<imax+1; i++)
       for (unsigned int k=0; k<kmax+1; k++)
       {
           if (discretization_.isFluid(i,0,k))
           g(i,0,k) = v(i,0,k);
           if (discretization_.isFluid(i,jmax,k))
           g(i,jmax,k) =v(i,jmax,k);
       }
   for (unsigned int i=0; i<imax+1; i++)
       for (unsigned int j=0; j<jmax+1; j++)
       {
           if (discretization_.isFluid(i,j,0))
           h(i,j,0) = w(i,j,0);
           if (discretization_.isFluid(i,j,kmax))
           h(i,j,kmax) = w(i,j,kmax);
       }
}

void FluidSimulator::composeRHS()
{
    for (unsigned int i=1; i<imax+1; i++)
        for (unsigned int j=1; j<jmax+1; j++)
            for (unsigned int k=1; k<kmax+1; k++)
                if (discretization_.isFluid(i,j,k))
                rhs(i,j,k) = (1.0/dt) * ( (f(i,j,k) - f(i-1,j,k))/dx
                                        + (g(i,j,k) - g(i,j-1,k))/dy
                                        + (h(i,j,k) - h(i,j,k-1))/dz  );

}

void FluidSimulator::updateVelocities() //i am not  pressure part is cultprit... look at him not me :(
{
    for (unsigned int i=1; i<imax; i++)
        for (unsigned int j=1; j<jmax+1; j++)
            for (unsigned int k=1; k<kmax+1; k++)
                if (discretization_.isFluid(i,j,k) && discretization_.isFluid(i+1,j,k))
                    u(i,j,k) = f(i,j,k) - (dt/dx) *  ( p(i+1,j,k) - p(i,j,k) );

    for (unsigned int i=1; i<imax+1; i++)
        for (unsigned int j=1; j<jmax; j++)
            for (unsigned int k=1; k<kmax+1; k++)
                if (discretization_.isFluid(i,j,k) && discretization_.isFluid(i,j+1,k))
                    v(i,j,k) = g(i,j,k) - (dt/dy) *  ( p(i,j+1,k) - p(i,j,k) );

    for (unsigned int i=1; i<imax+1; i++)
        for (unsigned int j=1; j<jmax+1; j++)
            for (unsigned int k=1; k<kmax; k++)
                if (discretization_.isFluid(i,j,k) && discretization_.isFluid(i,j,k+1))
                    w(i,j,k) = h(i,j,k) - (dt/dz) *  ( p(i,j,k+1) - p(i,j,k) );
}

void FluidSimulator::determineNextDT(real const & limit)
{
    real max_u=10e-16;\
    for (int i=0;i<u.getSize(0);i++)
        for(int j=0;j<u.getSize(1);j++)
            for(int k=0;k<u.getSize(2);k++)
                if (u(i,j,k)>max_u)
                    max_u=u(i,j,k);

    real max_v=10e-16;
    for (int i=0;i<v.getSize(0);i++)
        for(int j=0;j<v.getSize(1);j++)
            for(int k=0;k<v.getSize(2);k++)
                if (v(i,j,k)>max_v)
                    max_v=v(i,j,k);

    real stepSize = tau*std::min( (Re/2.0) / (1.0/(dx*dx)+1.0/(dy*dy)),
                                   std::min( dx/std::abs(max_u), dy/std::abs(max_v) ) );
    dt = stepSize<limit ? stepSize : limit;
}

void FluidSimulator::refreshBoundaries()
{
    //boundary NORTH
    if (boundary_condition_N == "noslip")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(i,jmax+1,k) = 2.0 * boundary_velocity_N - u(i,jmax,k);
                v(i,jmax,k) = 0.0;
                w(i,jmax+1,k) = 2.0 * boundary_velocity_N  - w(i,jmax,k); }
    else if (boundary_condition_N == "inflow")
        for (unsigned int i=0; i<=imax; i++)
            for (unsigned int k=0; k<=kmax; k++) {
                u(i,jmax+1,k) = -u(i,jmax,k);
                v(i,jmax,k) = - boundary_velocity_N;
                w(i,jmax+1,k) = -w(i,jmax,k); }
    else if (boundary_condition_N == "outflow")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(i,jmax+1,k) = u(i,jmax,k);
                v(i,jmax,k) = v(i,jmax-1,k);
                w(i,jmax+1,k) = w(i,jmax,k); }

    //boundary SOUTH
    if (boundary_condition_S == "noslip")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(i,0,k) = 2.0 * boundary_velocity_S -u(i,1,k);
                v(i,0,k) = 0.0;
                w(i,0,k) = 2.0 * boundary_velocity_S -w(i,1,k); }
    else if (boundary_condition_S == "inflow")
        for (unsigned int i=0; i<=imax; i++)
            for (unsigned int k=0; k<=kmax; k++) {
                u(i,0,k) = -u(i,1,k);
                v(i,0,k) = boundary_velocity_S;
                w(i,0,k) = -w(i,1,k); }
    else if (boundary_condition_S == "outflow")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(i,0,k) = u(i,1,k);
                v(i,0,k) = v(i,1,k);
                w(i,0,k) = w(i,1,k); }

    //boundary EAST
    if (boundary_condition_E == "noslip")
        for (unsigned int j=1; j<=jmax; j++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(imax,j,k) = 0.0;
                v(imax+1,j,k) = 2.0 * boundary_velocity_E -v(imax,j,k);
                w(imax+1,j,k) = 2.0 * boundary_velocity_E -w(imax,j,k); }
    else if (boundary_condition_E == "inflow")
        for (unsigned int j=0; j<=jmax; j++)
            for (unsigned int k=0; k<=kmax; k++) {
                u(imax,j,k) = -boundary_velocity_E;
                v(imax+1,j,k) = -v(imax,j,k);
                w(imax+1,j,k) = -w(imax,j,k); }
    else if (boundary_condition_E == "outflow")
        for (unsigned int j=1; j<=jmax; j++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(imax,j,k) = u(imax-1,j,k);
                v(imax+1,j,k) = v(imax,j,k);
                w(imax+1,j,k) = w(imax,j,k); }

    //boundary WEST
    if (boundary_condition_W == "noslip")
        for (unsigned int j=1; j<=jmax; j++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(0,j,k) = 0.0;
                v(0,j,k) = 2.0 * boundary_velocity_W -v(1,j,k);
                w(0,j,k) = 2.0 * boundary_velocity_W -w(1,j,k); }
    else if (boundary_condition_W == "inflow")
        for (unsigned int j=0; j<=jmax; j++)
            for (unsigned int k=0; k<=kmax; k++) {
                u(0,j,k) = boundary_velocity_W;
                v(0,j,k) = -v(1,j,k);
                w(0,j,k) = -w(1,j,k); }
    else if (boundary_condition_W == "outflow")
        for (unsigned int j=1; j<=jmax; j++)
            for (unsigned int k=1; k<=kmax; k++) {
                u(0,j,k) = u(1,j,k);
                v(0,j,k) = v(1,j,k);
                w(0,j,k) = w(1,j,k); }

    //boundary UP
    if (boundary_condition_U == "noslip")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int j=1; j<=jmax; j++) {
                u(i,j,kmax+1) = 2.0 * boundary_velocity_U -u(i,j,kmax);
                v(i,j,kmax+1) = 2.0 * boundary_velocity_U -v(i,j,kmax);
                w(i,j,kmax) = 0.0; }
    else if (boundary_condition_U == "inflow")
        for (unsigned int i=0; i<=imax; i++)
            for (unsigned int j=0; j<=jmax; j++) {
                u(i,j,kmax+1) = -u(i,j,kmax);
                v(i,j,kmax+1) = -v(i,j,kmax);
                w(i,j,kmax) = -boundary_velocity_U; }
    else if (boundary_condition_U == "outflow")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int j=1; j<=jmax; j++) {
                u(i,j,kmax+1) = u(i,j,kmax);
                v(i,j,kmax+1) = v(i,j,kmax);
                w(i,j,kmax) = w(i,j,kmax-1); }

    //boundary DOWN
    if (boundary_condition_D == "noslip")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int j=1; j<=jmax; j++) {
                u(i,j,0) = 2.0 * boundary_velocity_D -u(i,j,1);
                v(i,j,0) = 2.0 * boundary_velocity_D -v(i,j,1);
                w(i,j,0) = 0.0; }
    else if (boundary_condition_D == "inflow")
        for (unsigned int i=0; i<=imax; i++)
            for (unsigned int j=0; j<=jmax; j++) {
                u(i,j,0) = -u(i,j,1);
                v(i,j,0) = -v(i,j,1);
                w(i,j,0) = boundary_velocity_D; }
    else if (boundary_condition_D == "outflow")
        for (unsigned int i=1; i<=imax; i++)
            for (unsigned int j=1; j<=jmax; j++) {
                u(i,j,0) = u(i,j,0);
                v(i,j,0) = v(i,j,0);
                w(i,j,0) = w(i,j,0); }
}

void FluidSimulator::normalizepressure()
{
    real mean = 0.0;
    for (unsigned int i=1; i<imax+1; i++)
        for (unsigned int j=1; j<jmax+1; j++)
            for (unsigned int k=1; k<kmax+1; k++)
                 if (discretization_.isFluid(i,j,k))
                 mean += p(i,j,k);
    mean /= (imax*jmax*kmax);
    for (unsigned int i=1; i<imax+1; i++)
        for (unsigned int j=1; j<jmax+1; j++)
            for (unsigned int k=1; k<kmax+1; k++)
                if (discretization_.isFluid(i,j,k))
                p(i,j,k) -= mean;
}

void FluidSimulator::simulateTimeStepCount(unsigned int nrofTimesteps)
{
    p.fill(P_init);
    u.fill(U_init);
    v.fill(V_init);
    w.fill(W_init);
    rhs.fill(0.0);


    VTKWriter vtkWriter(discretization_, name, true, true);

    while (currentstepnumber<=nrofTimesteps)
    {
       determineNextDT( dtmax );

       refreshBoundaries();

       if (currentstepnumber%normalizationfrequency == 0)
       normalizepressure();

       computeFGH();

       composeRHS();

       std::cout<< "\nTime step = "<<currentstepnumber<< "\tTime = "<<currenttime<<"\n";

       solver_.solve(discretization_);

       updateVelocities(); //culprit for error

       currentstepnumber++;
       currenttime +=dt;

       if (currentstepnumber%outputinterval==0)
       vtkWriter.write();
    }
}

void FluidSimulator::simulate(real duration)
{
    VTKWriter vtkWriter(discretization_, name, true, true);
    while (currenttime<=duration)
    {
       determineNextDT( dtmax );

       refreshBoundaries();

       if (currentstepnumber%normalizationfrequency == 0)
          normalizepressure();

       computeFGH();

       composeRHS();

       solver_.solve(discretization_);

       updateVelocities();

       currentstepnumber++;
       currenttime=currenttime + dt;

       if (currentstepnumber%outputinterval==0)
          vtkWriter.write();
    }
}
