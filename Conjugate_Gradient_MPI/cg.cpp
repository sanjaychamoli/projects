/*
 *  Conjugated Gradient Method
 *  Group 11
 */

#include <cmath>
#include <fstream>
#include <iostream>

#include <mpi.h>

#include "Array.h"
#include "Array.cpp"
#include "ex03.h"

void updateGhostLayers( Array<real_t> &d, const int &PID, const int &NP );
void crushGhostLayers ( Array<real_t> &u, const int &PID, const int &NP );

int main ( int argc, char *argv[] ) {

  /* Start MPI processes*/
  MPI_Init (& argc, &argv);

  int pID = 0;
  int nP  = 0;

  MPI_Comm_size( MPI_COMM_WORLD, &nP );
  MPI_Comm_rank( MPI_COMM_WORLD, &pID );

  /* Local Variables, to set up our problem*/
  size_t nxGlobal = 1000;
  size_t nyGlobal = 1000;
  size_t maxIter  = 1000;
  real_t eps      = -6;
  real_t l2r      = 0.0;
  real_t beta     = 0.0;
  real_t alpha    = 0.0;

  if ( argc == 5)
  {
    nxGlobal= (size_t) std::stoi( argv[1] );
    nyGlobal= (size_t) std::stoi( argv[2] );
    maxIter = (size_t) std::stoi( argv[3] );
    eps     = (real_t) std::stoi( argv[4] );
  } else if ( pID == ROOT )
      std::cout << "Usage: mpirun -np <N> ./cg <nx> <ny> <c> <eps> \nN is the Number of Processes\n";

  real_t hx = 2.0/nxGlobal;
  real_t hy = 1.0/nyGlobal;

  if (pID == ROOT){
    std::ofstream fileExact("exact.txt");
    for ( size_t i = 0; i < nxGlobal+1; ++i )
      for ( size_t j = 0; j < nyGlobal+1; ++j )
        fileExact << i << "\t" << j << "\t" << U(hx*i, hy*j) << "\n";
    fileExact.close();
  }

  size_t Chunk  = nyGlobal/nP;          /*Amount of Rows  */


  size_t yStart = 1;					/* yStart always from 1 for all processes*/
  size_t yEnd   = (pID == nP-1) ? nyGlobal - Chunk*(nP-1)  : Chunk;
  if (pID == ROOT) yEnd = Chunk-1;
  std::cout << "\t" << yEnd << "\n";
  size_t nx = nxGlobal;
  size_t ny = yEnd+1;

  Array<real_t> r( nx+1, ny+1 ), u( nx+1, ny+1 ), d( nx+1, ny+1 ), z( nx+1, ny+1 );

  /* Imposing boundary conditions on the Upper Chunk alone, all other boundary conditions are zero*/
  if (pID == nP-1)
    for ( size_t i = 0; i < u.getSize(0); ++i )
      u(i, u.getSize(1)-1) = U(hx*i, 1.0);

  real_t delta = 0.0;
  real_t delta1= 0.0;
  for ( size_t i = 1; i < r.getSize(0)-1; ++i ) // Checked
    for ( size_t j = yStart; j < yEnd+1; ++j )
    {
      r(i, j) = f( hx*i, hy*(j+pID*Chunk)) +
                nxGlobal*nxGlobal*( u(i-1, j)   + u(i+1, j)   - 2*u(i, j)) /4.0 - K_2*u(i,j) +
                nyGlobal*nyGlobal*( u(i,   j-1) + u(i,   j+1) - 2*u(i, j)) ;
      delta += r(i ,j)*r(i ,j);
    }

  real_t temp = 0.0;
  /* Gather and add all deltas from processes*/
  MPI_Allreduce ( &delta, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  delta = temp;

  /* If L2 norm of residual is less that 10^eps already, end CG */
  if ( pID == ROOT )
    if ( std::sqrt(  delta/(nxGlobal-1)/(nyGlobal-1) ) < std::pow(10, eps) )
    {
      PROGRESS_MESG("Something Wrong?? ");
      std::cout << std::sqrt(std::fabs( delta/(nxGlobal+1)/(nyGlobal+1) )) << std::endl;
      return 0;
    }

  /* Take the residual as first direction */
  d = r;

  for ( size_t iter = 0; iter < maxIter+1; ++iter)
  {

    updateGhostLayers( d, pID, nP );
    alpha = 0.0;
    temp  = 0.0;

    /* Evaluating inner product (d,z) where z = Ad */
    for ( size_t i = 1; i < nxGlobal; ++i )
      for ( size_t j = yStart; j < yEnd+1; ++j )
      {
	z(i,j) =  nxGlobal*nxGlobal*( -d(i-1, j)   -d(i+1, j)   + 2*d(i, j)) / 4.0 + K_2*d(i,j) +
		  nyGlobal*nyGlobal*( -d(i,   j-1) -d(i,   j+1) + 2*d(i, j));
	temp +=  d(i,j)*z(i,j);
      }

    MPI_Allreduce ( &temp, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    alpha = delta/alpha;

    /* u = u + alpha*d */
    for ( size_t i = 1; i < nxGlobal; ++i )
      for ( size_t j = yStart; j < yEnd+1; ++j )
        u(i,j) += alpha * d(i,j);

    /* r = r - alpha*z */
    for ( size_t i = 1; i < nxGlobal; ++i )
      for ( size_t j = yStart; j < yEnd+1; ++j )
	    r(i, j) -= alpha * z(i, j);

    delta1 = 0.0;
    temp   = 0.0;
    for ( size_t i = 1; i < nxGlobal; ++i )
      for ( size_t j = yStart; j < yEnd+1; ++j )
        temp += r(i ,j)*r(i ,j);

    MPI_Allreduce ( &temp, &delta1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* Let all processes calculate l2r and break iterations if
     * desired residual norm is reached */
    l2r = std::sqrt(std::fabs( delta1/(nxGlobal-1)/(nyGlobal-1)));

    /* Allow only ROOT to print interation and residual norm information */
    if ( pID == ROOT )
    {
      std::cout << "Kome Schon - Iteration number: " <<  iter << "\t";
      std::cout << "L2 norm of residual is : " << l2r << std::endl;
    }

    /* All processes break the loop if l2r reaches desired value */
    if ( l2r <= std::pow(10, eps) ) {   break; }

    beta = delta1/delta;
    for ( size_t i = 1; i < nxGlobal; ++i )
      for ( size_t j = yStart; j < yEnd+1; ++j )
        d(i,j) = r(i ,j) + beta*d(i ,j);

    delta = delta1;

    if (iter == maxIter && pID == ROOT)
      WARN("Maximum iterations reached !! Check if you want to increa");

  }

  crushGhostLayers( u, pID, nP );

  if (pID == ROOT)
    remove("solution.txt");

  /* Write solution to file*/
  for (int PID = 0; PID < nP; ++PID)
  {
    if ( pID == PID )
    {
      std::ofstream fileSolution( "solution.txt", std::ios_base::app | std::ios_base::out);
      for ( size_t i = 0; i < nxGlobal+1; ++i )
	for ( size_t j = PID*Chunk, k=0; k < u.getSize(1)-1; ++j, ++k )
	  fileSolution << i << "\t" << j << "\t" << u(i,k) << "\n";
      fileSolution.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();

  return 0;
}

void updateGhostLayers( Array<real_t> &d, const int &PID, const int &NP)
{
  if (NP == 1)
    return;

  MPI_Status s;
  MPI_Request r[4];

  vec_t<real_t> sendLowerGhostLayer;
  vec_t<real_t> sendUpperGhostLayer;
  vec_t<real_t> recvLowerGhostLayer( d.getSize(0)-2, real_t(0.0)  );
  vec_t<real_t> recvUpperGhostLayer( d.getSize(0)-2, real_t(0.0)  );

  if ( PID == ROOT)
  {
    for ( size_t i = 1; i < d.getSize(0)-1; ++i) {
      //std::cout << " \n \n\n\n\n\n\n\t " << i << " \t" << d.getSize(1)-2 <<std::endl;
      sendLowerGhostLayer.push_back( d(i, d.getSize(1)-2));

    }

    /* TAG for sending data to process 1 is 0 */
    MPI_Isend ( sendLowerGhostLayer.data(), sendLowerGhostLayer.size(), MPI_DOUBLE, 1, PID,         MPI_COMM_WORLD, &r[0]);
    MPI_Irecv ( recvLowerGhostLayer.data(), recvLowerGhostLayer.size(), MPI_DOUBLE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &r[1]);

    MPI_Wait (&r[0], &s);
    MPI_Wait (&r[1], &s);

    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      d(i, d.getSize(1)-1) = recvLowerGhostLayer.at(i-1);

  }
  else if (PID == NP-1)
  {
    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      sendLowerGhostLayer.push_back( d(i, 1) );

    /* TAG for sending data to process 1 is 0 */
    MPI_Isend ( sendLowerGhostLayer.data(), sendLowerGhostLayer.size(), MPI_DOUBLE, PID-1, PID-1,       MPI_COMM_WORLD, &r[0]);
    MPI_Irecv ( recvLowerGhostLayer.data(), recvLowerGhostLayer.size(), MPI_DOUBLE, PID-1, MPI_ANY_TAG, MPI_COMM_WORLD, &r[1]);

    MPI_Wait (&r[0], &s);
    MPI_Wait (&r[1], &s);

    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      d(i, 0) = recvLowerGhostLayer.at(i-1);

  }
  else
  {
    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      sendLowerGhostLayer.push_back( d(i, 1) );

    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      sendUpperGhostLayer.push_back( d(i, d.getSize(1)-2));

    MPI_Isend ( sendLowerGhostLayer.data(), sendLowerGhostLayer.size(), MPI_DOUBLE, PID-1, PID-1,       MPI_COMM_WORLD, &r[0]);
    MPI_Isend ( sendUpperGhostLayer.data(), sendUpperGhostLayer.size(), MPI_DOUBLE, PID+1, PID,         MPI_COMM_WORLD, &r[1]);
    MPI_Irecv ( recvLowerGhostLayer.data(), recvLowerGhostLayer.size(), MPI_DOUBLE, PID-1, MPI_ANY_TAG, MPI_COMM_WORLD, &r[2]);
    MPI_Irecv ( recvUpperGhostLayer.data(), recvUpperGhostLayer.size(), MPI_DOUBLE, PID+1, MPI_ANY_TAG, MPI_COMM_WORLD, &r[3]);

    MPI_Wait (&r[0], &s);
    MPI_Wait (&r[1], &s);
    MPI_Wait (&r[2], &s);
    MPI_Wait (&r[3], &s);

    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      d(i, 0) = recvLowerGhostLayer.at(i-1);

    for ( size_t i = 1; i < d.getSize(0)-1; ++i)
      d(i, d.getSize(1)-1) = recvUpperGhostLayer.at(i-1);

  }

}

/* Delete Ghost Layers of Array Object (run by all processes) */
void crushGhostLayers ( Array<real_t> &u, const int &PID, const int &NP)
{
  if (PID == ROOT)
    u.deleteRowOfVec( u.getSize(1)-1);
  else if (PID == NP-1)
    u.deleteRowOfVec(0);
  else
  {
    u.deleteRowOfVec(0);
    u.deleteRowOfVec( u.getSize(1)-1);
  }
}

