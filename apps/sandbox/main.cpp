#include <iostream>
#include <complex>

#include "Scalar1D.h"
#include "Scalar2D.h"
#include "Scalar3D.h"
#include "Transforms.h"
#include "MathOps.h"
#include "TestSolver.h"
#include "mpi.h"

int main( int argc, const char * argv[] )
{
  /*
  std::size_t N = 16;
  std::size_t nOut = std::floor( N / 2 + 1 );
  
  Scalar1D< double > in( N );
  Scalar1D< std::complex< double > > out( nOut );
  
  for( std::size_t i = 0; i < N; i++ )
  {
    in( i ) = std::cos( i * ( 2 * M_PI / N ) );
    std::cout << in( i ) << std::endl;
  }
  
  fft::fft_r2c_1d( in, out );
    
  std::cout << std::endl;
  for( std::size_t i = 0; i < nOut; i++ )
  {
    std::cout << out( i ) << std::endl;
  }
   */
  
  /*
  std::size_t N = 32;
  std::size_t nOut = std::floor( N / 2 + 1 );
  
  Scalar1D< double > in( N );
  Scalar1D< std::complex< double > > out( nOut );
  
  out( 1 ) = std::complex< double >( 1.0, 0.0 );
  
  fft::fft_c2r_1d( out, in );
  
  for( std::size_t i = 0; i < N; i++ )
  {
    std::cout << in( i ) / 2 << std::endl;
  }
   */
  
  /*
  std::size_t Nx = 8;
  std::size_t Ny = 8;
  std::size_t nOutX = std::floor( Nx / 2 + 1 );
  
  Scalar2D< double > in( Ny, Nx );
  Scalar2D< std::complex< double > > out( Ny, nOutX );
  
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < Nx; ++j )
    {
      in( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::cos( j * ( 2.0 * M_PI / Nx ) );
    }
  }
  
  fft::fft_r2c_2d_thread( in, out );
  
  std::cout << std::endl;
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < nOutX; ++j )
    {
      std::cout << out( i, j ) << std::endl;
    }
    std::cout << std::endl;
  }
   */
  
  std::size_t Nx = 8;
  std::size_t Ny = 8;
  std::size_t nOutX = std::floor( Nx / 2 + 1 );
  
  Scalar2D< double > in( Ny, Nx );
  Scalar2D< std::complex< double > > out( Ny, nOutX );
  
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < nOutX; ++j )
    {
      std::cout << out( i, j ) << std::endl;
    }
    std::cout << std::endl;
  }
  
  fft::fft_r2c_2d_thread( in, out );
  
  std::cout << std::endl;
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < Nx; ++j )
    {
      in( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::cos( j * ( 2.0 * M_PI / Nx ) );
    }
  }
  
  return 0;
}

