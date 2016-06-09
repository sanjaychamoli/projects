#include "SORSolver.hh"



//initialising solver
SORSolver::SORSolver( real relaxfactor_, real epsilon_, unsigned int no_iteration_,unsigned int checkfrequency_, int imax_, int jmax_, int kmax_)
{
 relaxfactor=relaxfactor_;
 epsilon=epsilon_;	
 no_iteration=no_iteration_;
 checkfrequency=checkfrequency_;
 imax=imax_;
 jmax=jmax_;
 kmax=kmax_;
}


SORSolver::SORSolver ( const FileReader & configuration )
{
    imax = configuration.getIntParameter("imax");
    jmax = configuration.getIntParameter("jmax");
    kmax = configuration.getIntParameter("kmax");

    relaxfactor= configuration.getRealParameter("omg");
    CHECK_MSG((relaxfactor <= 1.9 && relaxfactor >= 1.7), "omg should be between 1.7 and 1.9");

    epsilon=configuration.getRealParameter("eps");
    CHECK_MSG((epsilon > 0), "eps should be greater than 0");

    no_iteration=configuration.getIntParameter("itermax");
    CHECK_MSG((no_iteration > 0), "Number of iterations (itermax) should be greater than 0");

    if (!configuration.checkparameter("checkfrequency")) {
        WARN("\nCheck Frequency missing in input...setting it to 1.\n");
       checkfrequency = 1; }
    else
       checkfrequency = configuration.getIntParameter("checkfrequency");
}
 
 
 
// setting boundary
void SORSolver::setboundary(Array<real> & p)
{
    for (unsigned int i=0; i<imax+1;i++)
        for (unsigned int k=0; k<kmax+1;k++)
        {
            p(i,0,k)=p(i,1,k);
            p(i,jmax+1,k)=p(i,jmax,k);
        }

   for(unsigned int j=0;j<jmax+1;j++)
      for (unsigned int k=0; k<kmax+1;k++)
        {
            p(0,j,k)=p(1,j,k);
            p(imax+1,j,k)=p(imax,j,k);
        }
   for(unsigned int i=0;i<imax+1;i++)
      for (unsigned int j=0; j<jmax+1;j++)
        {
            p(i,j,0)=p(i,j,1);
            p(i,j,kmax+1)=p(i,j,kmax);
        }
}


// calculating residual
real SORSolver::residual(StaggeredGrid & grid)
{
    real temp=0.0;
    real resid=0.0;
    real hx2=1.0/(grid.dx()*grid.dx());
    real hy2=1.0/(grid.dy()*grid.dy());
    real hz2=1.0/(grid.dz()*grid.dz());

    for (unsigned int i=1; i<imax+1;i++)
        for (unsigned int j=1; j<jmax+1;j++)
            for (unsigned int k=1; k<kmax+1; k++)
 //               if (grid.isFluid(i,j,k))
                {
                   temp  =   (grid.p(i,j,k,EAST)  + grid.p(i,j,k,WEST)  - 2.0*grid.p(i,j,k,CENTER)) *hx2
                          +  (grid.p(i,j,k,NORTH) + grid.p(i,j,k,SOUTH) - 2.0*grid.p(i,j,k,CENTER)) *hy2
                          +  (grid.p(i,j,k,UP)    + grid.p(i,j,k,DOWN)  - 2.0*grid.p(i,j,k,CENTER)) *hz2
                          -   grid.rhs()(i,j,k) ;
                resid += temp*temp;
                }
    return (sqrt(resid/grid.getNumFluid()));
}


// SOR slover
bool SORSolver::solve( StaggeredGrid & grid )
{
    Array<real> & p = grid.p();
    Array<real> & rhs = grid.rhs();
    real residual_remain = 1e100;

    real hx2 = 1.0/(grid.dx()*grid.dx());
    real hy2 = 1.0/(grid.dy()*grid.dy());
    real hz2=1.0/(grid.dz()*grid.dz());
    real constant = 1/(2.0*(hx2 + hy2 + hz2));

//    setboundary(p);

    for (unsigned int iter =0;iter<no_iteration;iter++)
    {
    for (unsigned int i=1; i<imax+1;i++)
        for (unsigned int j=1; j<jmax+1;j++)
            for (unsigned int k=1; k<kmax+1;k++)
               if (grid.isFluid(i,j,k))
                    p(i,j,k) =      (1.0-relaxfactor)* p(i,j,k)
                                +   relaxfactor * constant * (   (grid.p(i,j,k,EAST)  +   grid.p(i,j,k,WEST)    ) *hx2
                                                               + (grid.p(i,j,k,NORTH) +   grid.p(i,j,k,SOUTH)   ) *hy2
                                                               + (grid.p(i,j,k,UP)    +   grid.p(i,j,k,DOWN)    ) *hz2
                                                               - rhs(i,j,k));
 //  setboundary(p);
	
	if (iter%checkfrequency==0)
	{
        residual_remain = residual(grid);
        //std::cout<<"\nIteration no  "<<iter<< "\tResidual = "<< residual_remain<<"\n";
	}
		
	if (residual_remain<epsilon)
	{
        std::cout<<"!!!!!!!!! Solution converged in "<<iter<<" iterations!!!!!!!!!"<<std::endl;
		return true;
	}
    }
	
   std::cout<<"!!!!!!!!! Solution didn't converge !!!!!!!!!"<<std::endl;
   return false;
   
}
