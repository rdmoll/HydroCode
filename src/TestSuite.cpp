//
//  TestSuite.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include "TestSuite.h"

namespace diagnostics
{

TestSuite::TestSuite()
{
  
}

TestSuite::~TestSuite()
{
  
}

void TestSuite::testReadParams()
{
  std::cout << "testReadParams         : ";
  
  const char* delimiter = "=";
  const char* filename = "/Users/rmoll/Documents/dev/projects/HydroCode/parameters.txt";
  io::ReadParams parameterFile( filename,delimiter );
  
  double value1 = parameterFile.dParam["value1"];
  double value2 = parameterFile.dParam["value2"];
  std::string value3 = parameterFile.strParam["value3"];
  double value4 = parameterFile.dParam["value4"];
  int value5 = parameterFile.iParam["value5"];
  
  if( value1 != 1.5)
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value2 != 2.6)
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value3 != "StringTest")
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value4 != 4.8)
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value5 != 16)
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  std::cout << "PASS" << std::endl;
}

void TestSuite::testFourierTransforms1()
{
  double pi = std::acos(-1.0);
  
  int N = 16;
  int nOut = std::floor( N / 2 + 1 );
  
  std::vector< double > realTruth( N, 0.0 );
  std::vector< std::complex< double > > compTruth( nOut, std::complex< double >( 0.0, 0.0 ) );
  
  hydroCode::FourierTransforms fft;
  
  std::vector< double > in( N, 0.0 );
  std::vector< std::complex< double > > out( nOut, std::complex< double >( 0.0, 0.0 ) );
  
  for( int i = 0; i < N; i++ )
  {
    realTruth[ i ] = std::cos( i * ( 2 * pi / N ) );
    in[ i ] = std::cos( i * ( 2 * pi / N ) );
  }
  
  compTruth[ 1 ] = std::complex< double >( 1.0, 0.0 );
  
  fft.fft_r2c( in, out );
  
  double compMSE = 0.0;
  for( int i = 0 ; i < nOut ; i++ )
  {
    compMSE = compMSE + std::pow( compTruth[ i ].real() - out[ i ].real() / 8.0, 2.0 )
                      + std::pow( compTruth[ i ].imag() - out[ i ].imag() / 8.0, 2.0 );
  }
  compMSE /= 2.0 * nOut;
  
  fft.fft_c2r( out, in );
  
  double realMSE = 0.0;
  for( int i = 0; i < N; i++ )
  {
    realMSE += std::pow( realTruth[ i ] - in[ i ] / 16.0, 2.0 );
  }
  realMSE /= N;
  
  std::cout << "testFourierTransforms1 : complex sigma = " << std::sqrt( compMSE ) << " , " <<
                                         "real sigma = " << std::sqrt( realMSE ) << std::endl;
}

void TestSuite::testFourierTransforms2()
{
  double pi = std::acos(-1.0);
  
  int N = 16;
  int nOut = std::floor( N / 2 + 1 );
  
  std::vector< double > realTruth( N, 0.0 );
  std::vector< std::complex< double > > compTruth( nOut, std::complex< double >( 0.0, 0.0 ) );
  
  hydroCode::FourierTransforms fftNew( N );
  
  for( int i = 0; i < N; i++ )
  {
    realTruth[ i ] = std::cos( i * ( 2 * pi / N ) );
    fftNew.realArray[ i ] = std::cos( i * ( 2 * pi / N ) );
  }
  
  compTruth[ 1 ] = std::complex< double >( 1.0, 0.0 );
  
  fftNew.fft_r2c();
  
  double compMSE = 0.0;
  for( int i = 0 ; i < nOut ; i++ )
  {
    compMSE = compMSE + std::pow( compTruth[ i ].real() - fftNew.compArray[ i ][ 0 ] / 8.0, 2.0 )
                      + std::pow( compTruth[ i ].imag() - fftNew.compArray[ i ][ 1 ] / 8.0, 2.0 );
  }
  compMSE /= 2.0 * nOut;
  
  fftNew.fft_c2r( );
  
  double realMSE = 0.0;
  for( int i = 0; i < N; i++ )
  {
    realMSE += std::pow( realTruth[ i ] - fftNew.realArray[ i ] / 16.0, 2.0 );
  }
  realMSE /= N;
  
  std::cout << "testFourierTransforms2 : complex sigma = " << std::sqrt( compMSE ) << " , " <<
                                         "real sigma = " << std::sqrt( realMSE ) << std::endl;
}

void TestSuite::testDeriv()
{
  double pi = std::acos(-1.0);
  
  int N = 16;
  int nOut = std::floor( N / 2 + 1 );
  
  std::vector< double > realTruth( N, 0.0 );
  std::vector< double > in( N, 0.0 );
  std::vector< std::complex< double > > out( nOut, std::complex< double >( 0.0, 0.0 ) );
  
  hydroCode::FourierTransforms fft;
  mathOps::Derivatives ops;
  
  for( int i = 0; i < N; i++ )
  {
    realTruth[ i ] = -2.0 * pi * std::sin( i * ( 2.0 * pi / N ) );
    in[ i ] = std::cos( i * ( 2.0 * pi / N ) );
  }
  
  fft.fft_r2c( in, out );
  
  ops.deriv( out );
  
  fft.fft_c2r( out, in );
  
  double mse = 0.0;
  for( int i = 0; i < N; i++ )
  {
    mse += std::pow( realTruth[ i ] - in[ i ] / 16.0, 2.0 );
  }
  mse /= N;
  
  std::cout << "testDeriv              : sigma = " << std::sqrt( mse ) << std::endl;
}

void TestSuite::testDeriv2()
{
  double pi = std::acos(-1.0);
  
  int N = 16;
  int nOut = std::floor( N / 2 + 1 );
  
  std::vector< double > realTruth( N, 0.0 );
  std::vector< double > in( N, 0.0 );
  std::vector< std::complex< double > > out( nOut, std::complex< double >( 0.0, 0.0 ) );
  
  hydroCode::FourierTransforms fft;
  mathOps::Derivatives ops;
  
  for( int i = 0; i < N; i++ )
  {
    realTruth[ i ] = -4.0 * pi * pi * std::cos( i * ( 2.0 * pi / N ) );
    in[ i ] = std::cos( i * ( 2.0 * pi / N ) );
  }
  
  fft.fft_r2c( in, out );
  
  ops.deriv2( out );
  
  fft.fft_c2r( out, in );
  
  double mse = 0.0;
  for( int i = 0; i < N; i++ )
  {
    mse += std::pow( realTruth[ i ] - in[ i ] / 16.0, 2.0 );
  }
  mse /= N;
  
  std::cout << "testDeriv2             : sigma = " << std::sqrt( mse ) << std::endl;
}

} // diagnostics
