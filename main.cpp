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
  unitTests.testFourierTransforms_2D();
  unitTests.testDeriv();
  unitTests.testDeriv2();
  std::cout << std::endl;
  
  // Solve simple equation
  double pi = std::acos(-1.0);
  
  size_t Nx = 64;
  size_t Ny = 64;
  int nOutX = std::floor( Nx / 2 + 1 );
  int nOutY = std::floor( Ny / 2 + 1 );
  size_t nSteps = 40;
  double fac1 = 0.01;
  double fac2 = 0.001;
  
  std::vector< std::vector< double > > T0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > T1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< std::complex< double > > >
    T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  hydroCode::FourierTransforms fft;
  io::ioNetCDF testWriter( "/Users/rmoll/Desktop/test_AdvDiff_xy.nc", "data", Nx, Ny );
  
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      T0_phys[ i ][ j ] = std::cos( i * ( 2 * pi / Nx ) );
    }
  }
  testWriter.write( 0, T0_phys );
  
  // Starting time step loop
  for( size_t t = 0; t < nSteps; t++ )
  {
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, T0_phys, T_spec );
    
    for( int i = 0 ; i < nOutX ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        T_spec[ i ][ j ] /= ( 1.0 - fac1 * 2.0 * pi * i * std::complex< double >( 0.0, 1.0 )
                            + fac2 * 4.0 * pi * pi * i * i );
      }
    }
    
    for( int i = nOutX ; i < Nx ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        T_spec[ i ][ j ] /= ( 1.0 + fac1 * 2.0 * pi * ( Nx - i ) * std::complex< double >( 0.0, 1.0 )
                            + fac2 * 4.0 * pi * pi * ( Nx - i ) * ( Nx - i ) );
      }
    }
    
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, T_spec, T1_phys );
    
    for( int i = 0; i < Nx; i++ )
    {
      for( int j = 0; j < Ny; j++ )
      {
        T0_phys[ i ][ j ] = T1_phys[ i ][ j ] / ( Nx * Ny );
      }
    }
    
    // Write netCDF data
    testWriter.write( t + 1, T0_phys );
  }
  
  return 0;
}

