#include "MathOps.h"

namespace mathOps
{

void calcDeriv( Scalar1D< std::complex< double > >& f_spec,
                Scalar1D< std::complex< double > >& df_spec )
{
  std::complex< double > fac = std::complex< double >( 0.0, 2.0 * M_PI );
  
  for( size_t i = 0; i < df_spec.size(); ++i )
  {
    df_spec( i ) = fac * std::complex< double >( i, 0.0 ) * f_spec( i );
  }
}

void calcDerivY( Scalar2D< std::complex< double > >& f_spec,
                 Scalar2D< std::complex< double > >& df_spec,
                 std::size_t Nx, std::size_t Ny, double Ly )
{
  std::complex< double > fac = std::complex< double >( 0.0, 2.0 * M_PI / Ly );
  std::size_t nOutX = std::floor( Nx / 2 + 1 );
  std::size_t nOutY = std::floor( Ny / 2 + 1 );
  
  for( std::size_t i = 0 ; i < nOutY ; ++i )
  {
    for( std::size_t j = 0 ; j < nOutX ; ++j )
    {
      df_spec( i, j ) = fac * std::complex< double >( i, 0.0 ) * f_spec( i, j );
    }
  }
  
  for( std::size_t i = nOutY ; i < Ny ; ++i )
  {
    for( std::size_t j = 0 ; j < nOutX ; ++j )
    {
      df_spec( i, j ) = -fac * std::complex< double >( Nx - i, 0.0 ) * f_spec( i, j );
    }
  }
}

void calcDerivX( Scalar2D< std::complex< double > >& f_spec,
                 Scalar2D< std::complex< double > >& df_spec,
                 std::size_t Nx, std::size_t Ny, double Lx )
{
  std::complex< double > fac = std::complex< double >( 0.0, 2.0 * M_PI / Lx );
  std::size_t nOutX = std::floor( Nx / 2 + 1 );
  
  for( std::size_t i = 0 ; i < Ny ; i++ )
  {
    for( std::size_t j = 0 ; j < nOutX ; j++ )
    {
      df_spec( i, j ) = fac * std::complex< double >( j, 0.0 ) * f_spec( i, j );
    }
  }
}

void calcDeriv2()
{
}

void gradient()
{
}

void laplacian()
{
}

} // mathOps
