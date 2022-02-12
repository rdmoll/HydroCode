#ifndef MathOps_h
#define MathOps_h

#include <complex>
#include "Scalar1D.h"
#include "Scalar2D.h"

namespace mathOps
{

void calcDeriv( Scalar1D< std::complex< double > >& f_spec,
                Scalar1D< std::complex< double > >& df_spec );

void calcDerivY( Scalar2D< std::complex< double > >& f_spec,
                 Scalar2D< std::complex< double > >& df_spec,
                 std::size_t Nx, std::size_t Ny, double Ly );

void calcDerivX( Scalar2D< std::complex< double > >& f_spec,
                 Scalar2D< std::complex< double > >& df_spec,
                 std::size_t Nx, std::size_t Ny, double Lx );

void calcDeriv2();

void gradient();

void laplacian();

}

#endif /* mathOps_h */
