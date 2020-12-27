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
  
  size_t N = 16;
  int nOut = std::floor( N / 2 + 1 );
  size_t nSteps = 10;
  double fac1 = 0.01;
  double fac2 = 0.0;
  
  std::vector< double > T0_phys( N, 0.0 );
  for( int i = 0; i < N; i++ )
  {
    T0_phys[ i ] = std::cos( i * ( 2 * pi / N ) );
  }
  
  std::vector< double > T1_phys( N, 0.0 );
  std::vector< std::complex< double > > T_spec( nOut, std::complex< double >( 0.0, 0.0 ) );
  
  hydroCode::FourierTransforms fft;
  
  fft.fft_r2c( T0_phys, T_spec );
  
  for( size_t j = 0; j < nSteps; j++ )
  {
    for( size_t i = 1; i < nOut; i++ )
    {
      T_spec[ i ] /= ( 1.0 + fac1 * 2.0 * pi * i * std::complex< double >( 0.0, 1.0 )
                           + fac2 * 4.0 * pi * pi * i * i );
    }
  }
  
  fft.fft_c2r( T_spec, T1_phys );
  
  for( int i = 0; i < N; i++ )
  {
    std::cout << T0_phys[ i ] << " --> " << T1_phys[ i ] / 16.0 << std::endl;
  }
  
  // Write netCDF data
  std::vector< std::vector< double > > testData( 10, std::vector< double >( 10, 0.0 ) );
  io::ioNetCDF testWriter( "/Users/rmoll/Desktop/test_class_xy.nc", "data", 10, 10 );
  testWriter.write( 0, testData );
  testWriter.write( 1, testData );
  testWriter.write( 2, testData );

  return 0;
}

