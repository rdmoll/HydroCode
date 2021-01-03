//
//  main.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include <iostream>
#include <netcdf>
#include "TestSuite.h"
#include "FourierTransforms.h"
#include "ioNetCDF.h"

int main( int argc, const char * argv[] )
{
  // do unit testing
  std::cout << "Running unit tests" << std::endl;
  diagnostics::TestSuite unitTests;
  unitTests.testReadParams();
  unitTests.testFourierTransforms1();
  unitTests.testFourierTransforms2();
  unitTests.testDeriv();
  unitTests.testDeriv2();
  std::cout << std::endl;
  
  double pi = std::acos(-1.0);
  
  size_t Nx = 16;
  size_t Ny = 16;
  int nOutX = std::floor( Nx / 2 + 1 );
  int nOutY = std::floor( Ny / 2 + 1 );
  size_t nSteps = 10;
  double fac1 = 0.01;
  double fac2 = 0.0;
  
  std::vector< double > T0_phys( Nx, 0.0 );
  for( int i = 0; i < Nx; i++ )
  {
    T0_phys[ i ] = std::cos( i * ( 2 * pi / Nx ) );
  }
  
  std::vector< double > T1_phys( Nx, 0.0 );
  std::vector< std::complex< double > > T_spec( nOutX, std::complex< double >( 0.0, 0.0 ) );
  
  hydroCode::FourierTransforms fft;
  
  fft.fft_r2c( T0_phys, T_spec );
  
  for( size_t j = 0; j < nSteps; j++ )
  {
    for( size_t i = 1; i < nOutX; i++ )
    {
      T_spec[ i ] /= ( 1.0 + fac1 * 2.0 * pi * i * std::complex< double >( 0.0, 1.0 )
                           + fac2 * 4.0 * pi * pi * i * i );
    }
  }
  
  fft.fft_c2r( T_spec, T1_phys );
  
  for( int i = 0; i < Nx; i++ )
  {
    std::cout << T0_phys[ i ] << " --> " << T1_phys[ i ] / 16.0 << std::endl;
  }
  
  std::vector< std::vector< double > > T_write( Nx, std::vector< double >( Ny, 0.0 ) );
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      T_write[ i ][ j ] = T1_phys[ i ];
    }
  }
  
  // Write netCDF data
  io::ioNetCDF testWriter( "/Users/rmoll/Desktop/test_AdvDiff_xy.nc", "data", Nx, Ny );
  testWriter.write( 0, T_write );
  testWriter.write( 1, T_write );
  testWriter.write( 2, T_write );

  return 0;
}

