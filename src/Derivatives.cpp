//
//  Derivatives.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include "Derivatives.h"

namespace mathOps
{

Derivatives::Derivatives()
{
  
}

Derivatives::~Derivatives()
{
  
}

void Derivatives::deriv( std::vector< std::complex< double > > &inOut )
{
  double pi = std::acos(-1.0);
  
  for( size_t i = 0; i < inOut.size(); i++ )
  {
    inOut[ i ] *= 2.0 * pi * i * std::complex< double >( 0.0, 1.0 );
  }
}

void Derivatives::calcDerivX( std::vector< std::vector< std::complex< double > > >& f_spec,
                              std::vector< std::vector< std::complex< double > > >& df_spec,
                              size_t Nx, size_t Ny, double Lx )
{
  double fac = 2.0 * pi / Lx;
  
  size_t nOutX = std::floor( Nx / 2 + 1 );
  size_t nOutY = std::floor( Ny / 2 + 1 );
  
  for( size_t i = 0 ; i < nOutX ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      df_spec[ i ][ j ] = fac * i * std::complex< double >( 0.0, 1.0 ) * f_spec[ i ][ j ];
    }
  }
  
  for( size_t i = nOutX ; i < Nx ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      df_spec[ i ][ j ] = -fac * ( Nx - i ) * std::complex< double >( 0.0, 1.0 ) * f_spec[ i ][ j ];
    }
  }
}

void Derivatives::calcDerivY( std::vector< std::vector< std::complex< double > > >& f_spec,
                              std::vector< std::vector< std::complex< double > > >& df_spec,
                              size_t Nx, size_t Ny, double Ly )
{
  double fac = 2.0 * pi / Ly;
  
  size_t nOutY = std::floor( Ny / 2 + 1 );
  
  for( size_t i = 0 ; i < Nx ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      df_spec[ i ][ j ] = fac * j * std::complex< double >( 0.0, 1.0 ) * f_spec[ i ][ j ];
    }
  }
}

void Derivatives::deriv2( std::vector< std::complex< double > > &inOut )
{
  double pi = std::acos(-1.0);
  
  for( size_t i = 0; i < inOut.size(); i++ )
  {
    inOut[ i ] *= -4.0 * pi * pi * i * i;
  }
}

void Derivatives::gradient()
{
  
}

void Derivatives::laplacian()
{
  
}

} // mathOps
