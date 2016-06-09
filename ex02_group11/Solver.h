#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include "Timing.h"


typedef double real_t;
typedef std::vector<real_t> vec_t;

//#define PI  3.141592653589793238462643383279502884
#define PI  M_PI

real_t k2 = 4*PI*PI;


// Function declaration

//Function to write into a file
void writetofile(vec_t &u, const size_t &nx, const size_t&ny, const std::string & filename);

// initialized rhs
void rhs(vec_t &f, const size_t &nx, const size_t &ny) ;

// red black gauss seidel iterations
real_t rbgs(vec_t &u, vec_t &f, vec_t &r, const size_t &nx, const size_t &ny, const size_t &c, int threads) ;

// dirichlet boundary condition
void db(vec_t & u, const size_t &nx, const size_t&ny);



/*********************************************************************************/


// Function defination

void writetofile(vec_t &u, const size_t &nx, const size_t&ny, const std::string & filename) {
        real_t dx=2./nx,dy=1./ny;
        std::ofstream file;
    file.open(filename.c_str());
    file << "# x y u(x,y)\n";
    for (real_t j = 0; j <= ny; file<<std::endl, ++j)
        for (real_t i = 0; i <= nx;++i)
            file << i *dx << " " << j*dy << " " << u[i + (nx + 1)*j] << "\n";
    file.close();


}

void rhs(vec_t &f, const size_t &nx, const size_t &ny) {

    real_t hx = 2.0/real_t(nx); real_t hy = 1.0/real_t(ny);

    // rhs is zero on the boundaries
    for ( size_t i = 0; i < nx+1; ++i)
        for ( size_t j = 0; j < ny+1; ++j)
            f[i + (nx+1)*j] = k2 * std::sin ( 2.0*PI*hx*i ) * std::sinh( 2.0*PI*hy*j );
}


real_t rbgs(vec_t &u, vec_t &f, vec_t &r, const size_t &nx, const size_t &ny, const size_t &c, int threads) {

    std::cout << std::endl << "Threads: " << threads ;
    real_t inv_hx2_  = real_t(nx*nx)/(2.0*2.0);
    real_t inv_hy2_	= real_t(ny*ny);
    real_t h2_		= 1.0 / (2.0* inv_hx2_ + 2.0*inv_hy2_ + k2);
    size_t chunk = ny / threads;

    for (size_t iter = 0; iter < c; ++iter)
    {

            real_t inv_hx2 = inv_hx2_;
            real_t inv_hy2 = inv_hy2_;
            real_t h2 = h2_;
        omp_set_num_threads(threads);
        #pragma omp parallel for firstprivate(inv_hx2_, inv_hy2_, h2_) schedule(static,chunk)
        // Blacks
        for (size_t j = 1; j < ny; ++j)
            for (size_t i = 1+(j+1)%2; i < nx; i += 2)
              u[i + (nx + 1)*j] = h2 * (	inv_hx2 * (u[i-1	+ (nx+1)*j]		+	u[i+1	+	(nx+1)*j] ) +
                                          inv_hy2 * (u[i		+ (nx+1)*(j-1)]	+	u[i		+	(nx+1)*(j+1)]	 ) +
                                                                                  f[i		+	(nx+1)*j] );

        omp_set_num_threads(threads);
        #pragma omp parallel for firstprivate(inv_hx2_, inv_hy2_, h2_) schedule(static,chunk)
        // Reds
        for (size_t j = 1; j < ny; ++j)
          for (size_t i = 1+j%2; i < nx; i += 2)
              u[i + (nx + 1)*j] = h2 * (	inv_hx2 * (u[i-1	+ (nx+1)*j]		+	u[i+1	+ (nx+1)*j]			) +
                                          inv_hy2 * (u[i		+ (nx+1)*(j-1)]	+	u[i		+ (nx+1)*(j + 1)]	) +
                                                                                  f[i		+ (nx+1)*j] );
    }

    real_t res = 0.0;
    omp_set_num_threads(threads);
    #pragma omp parallel for firstprivate(inv_hx2_, inv_hy2_, h2_) reduction(+:res)
    for (size_t j = 1; j < ny; ++j)
        for (size_t i = 1; i < nx; ++i)
        {
            r[i + (nx + 1)*j] =		inv_hx2_ * (u[i-1	+ (nx+1)*j]		+ u[i+1 + (nx+1)*j]		- 2 *	u[i + (nx+1)*j]) +
                                    inv_hy2_ * (u[i		+ (nx+1)*(j-1)] + u[i	+ (nx+1)*(j+1)]	- 2 *	u[i + (nx+1)*j])
                                                                    - k2*u[i + (nx+1)*j] + f[i + (nx + 1)*j];
            res += r[i + (nx + 1)*j] * r[i + (nx + 1)*j];
        }
    return std::sqrt( res / ((nx - 1)*(ny - 1)) );
}

void db(vec_t & u, const size_t &nx, const size_t&ny)
{

    real_t hx = 2.0/real_t(nx);

    // Do nothing for bottom, left and right boundaries
    // Treat top boundary
    for ( size_t i = 1; i < nx; ++i)
        u[i + (nx+1)*ny] =  std::sin ( 2.0*PI*hx*i ) * std::sinh( 2.0*PI );

}



#endif // SOLVER_H
