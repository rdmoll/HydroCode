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
  std::cout << "testReadParams           : ";
  
  const char* delimiter = "=";
  const char* filename = "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/parametersTest.txt";
  io::ReadParams parameterFile( filename, delimiter );
  
  double value1 = parameterFile.dParam[ "value1" ];
  double value2 = parameterFile.dParam[ "value2" ];
  std::string value3 = parameterFile.strParam[ "value3" ];
  double value4 = parameterFile.dParam[ "value4" ];
  int value5 = parameterFile.iParam[ "value5" ];
  
  if( value1 != 1.5 )
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value2 != 2.6 )
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value3 != "StringTest")
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value4 != 4.8 )
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  if( value5 != 16 )
  {
    std::cout << "FAIL" << std::endl;
    return;
  }
  
  std::cout << "PASS" << std::endl;
}

void TestSuite::testFourierTransforms1()
{
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
  
  std::cout << "testFourierTransforms1   : complex sigma = " << std::sqrt( compMSE ) << " , " <<
                                             "real sigma = " << std::sqrt( realMSE ) << std::endl;
}

void TestSuite::testFourierTransforms2()
{
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
  
  std::cout << "testFourierTransforms2   : complex sigma = " << std::sqrt( compMSE ) << " , " <<
                                            "real sigma = " << std::sqrt( realMSE ) << std::endl;
}

void TestSuite::testFourierTransforms_2D()
{
  size_t Nx = 16;
  size_t Ny = 16;
  //int nOutX = std::floor( Nx / 2 + 1 );
  int nOutY = std::floor( Ny / 2 + 1 );
  
  std::vector< std::vector< double > > realTruth( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > realTest( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< std::complex< double > > >
    compTruth( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    compTest( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  hydroCode::FourierTransforms fft;
  
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      realTruth[ i ][ j ] = std::cos( i * ( 2 * pi / Nx ) ) * std::cos( j * ( 2 * pi / Ny ) );
    }
  }
  
  compTruth[ 1 ][ 1 ] = std::complex< double >( 1.0, 0.0 );
  compTruth[ Nx - 1 ][ 1 ] = std::complex< double >( 1.0, 0.0 );
  
  fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, realTruth, compTest );
  
  double compMSE = 0.0;
  for( int i = 0 ; i < Nx ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      compMSE = compMSE + std::pow( compTruth[ i ][ j ].real() - compTest[ i ][ j ].real() / 64, 2.0 )
                        + std::pow( compTruth[ i ][ j ].imag() - compTest[ i ][ j ].imag() / 64, 2.0 );
    }
  }
  compMSE /= 2.0 * Nx * nOutY;
  
  fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, compTest, realTest );
  
  double realMSE = 0.0;
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      realMSE += std::pow( realTruth[ i ][ j ] - realTest[ i ][ j ] / ( Nx * Ny ), 2.0 );
    }
  }
  realMSE /= ( Nx * Ny );
  
  std::cout << "testFourierTransforms_2D : complex sigma = " << std::sqrt( compMSE ) << " , " <<
                                             "real sigma = " << std::sqrt( realMSE ) << std::endl;
}

void TestSuite::testDeriv()
{
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
  
  std::cout << "testDeriv                : sigma = " << std::sqrt( mse ) << std::endl;
}

void TestSuite::testDeriv2()
{
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
  
  std::cout << "testDeriv2               : sigma = " << std::sqrt( mse ) << std::endl;
}

void TestSuite::testReadWriteIO()
{
  testWriteIO();
  testReadIO();
}

void TestSuite::testWriteIO()
{
  size_t Nx = 3;
  size_t Ny = 3;
  size_t Nsteps = 2;
  io::ioNetCDF testIO( "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/test_io.nc", "data", Nx, Ny, 'w' );
  
  std::vector< std::vector< double > > testDataWrite( Nx, std::vector< double >( Ny, 0.0 ) );
  
  for( int ntime = 0; ntime < Nsteps; ntime++ )
  {
    for( int i = 0; i < Nx; i++ )
    {
      for( int j = 0; j < Ny; j++ )
      {
        testDataWrite[ i ][ j ] = 1.0 * ( Ny * i + j + 1.0 );
      }
    }
    
    testIO.write( ntime, testDataWrite );
  }
}

void TestSuite::testReadIO()
{
  bool pass = true;
  size_t Nx = 3;
  size_t Ny = 3;
  size_t Nsteps = 2;
  io::ioNetCDF testIO( "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/test_io.nc", "data", Nx, Ny, 'r' );
  
  std::vector< std::vector< double > > testDataRead( Nx, std::vector< double >( Ny, 0.0 ) );
  
  for( int ntime = 0; ntime < Nsteps; ntime++ )
  {
    testIO.read( ntime, testDataRead );
    
    for( int i = 0; i < Nx; i++ )
    {
      for( int j = 0; j < Ny; j++ )
      {
        if( testDataRead[ i ][ j ] != 1.0 * ( Ny * i + j + 1.0 ) )
        {
          pass = false;
        }
      }
    }
  }
  
  if( pass )
  {
    std::cout << "testReadWriteIO          : PASS" << std::endl;
  }
  else
  {
    std::cout << "testReadParams           : FAIL" << std::endl;
  }
}

void TestSuite::simpleAdvDiff( std::string paramFile )
{
  std::cout << "simpleAdvDiff" << std::endl;
  
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  const char* delimiter = "=";
  io::ReadParams parameterFile( paramFile.c_str(), delimiter );
  
  size_t nSteps = static_cast< size_t >( parameterFile.iParam[ "nSteps" ] );
  double deltaT = parameterFile.dParam[ "deltaT" ];
  size_t Nx = static_cast< size_t >( parameterFile.iParam[ "Nx" ] );
  size_t Ny = static_cast< size_t >( parameterFile.iParam[ "Ny" ] );
  double Lx = parameterFile.dParam[ "Lx" ];
  double Ly = parameterFile.dParam[ "Ly" ];
  double c = parameterFile.dParam[ "c" ];
  double nu = parameterFile.dParam[ "nu" ];
  std::string testWriterFile = parameterFile.strParam[ "testFile" ];
  std::string truthWriterFile = parameterFile.strParam[ "truthFile" ];
  
  int nOutX = std::floor( Nx / 2 + 1 );
  int nOutY = std::floor( Ny / 2 + 1 );
  double advFac = 2.0 * pi * deltaT * c / Lx;
  double diffFac = 4.0 * pi * pi * deltaT * nu / ( Lx * Lx );
  
  std::vector< std::vector< double > > T0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > T1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > truth_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< std::complex< double > > >
    T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  hydroCode::FourierTransforms fft;
  
  io::ioNetCDF testWriter( testWriterFile, "data", Nx, Ny, 'w' );
  io::ioNetCDF truthWriter( truthWriterFile, "data", Nx, Ny, 'w' );
  
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      T0_phys[ i ][ j ] = std::cos( i * ( 2 * pi / Nx ) );
      truth_phys[ i ][ j ] = std::cos( i * ( 2 * pi / Nx ) );
    }
  }
  testWriter.write( 0, T0_phys );
  truthWriter.write( 0, truth_phys );
  
  // Starting time step loop
  for( size_t t = 0; t < nSteps; t++ )
  {
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, T0_phys, T_spec );
    
    for( int i = 0 ; i < nOutX ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        T_spec[ i ][ j ] /= ( 1.0 + advFac * i * std::complex< double >( 0.0, 1.0 )
                            + diffFac * i * i );
      }
    }
    
    for( int i = nOutX ; i < Nx ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        T_spec[ i ][ j ] /= ( 1.0 - advFac * ( Nx - i ) * std::complex< double >( 0.0, 1.0 )
                            + diffFac * ( Nx - i ) * ( Nx - i ) );
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
    
    for( int i = 0; i < Nx; i++ )
    {
      for( int j = 0; j < Ny; j++ )
      {
        truth_phys[ i ][ j ] = std::cos( i * 2.0 * pi / Nx - ( t + 1 ) * advFac );
      }
    }
    
    // Write netCDF data
    testWriter.write( t + 1, T0_phys );
    truthWriter.write( t + 1, truth_phys );
  }
  
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
            << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
            << "Wall clock time passed: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n" << std::endl;
}

void TestSuite::simpleAdvDiffNL( std::string paramFile )
{
  std::cout << "simpleAdvDiffNL" << std::endl;
  
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  const char* delimiter = "=";
  io::ReadParams parameterFile( paramFile.c_str(), delimiter );
  
  size_t nSteps = static_cast< size_t >( parameterFile.iParam[ "nSteps" ] );
  double deltaT = parameterFile.dParam[ "deltaT" ];
  size_t Nx = static_cast< size_t >( parameterFile.iParam[ "Nx" ] );
  size_t Ny = static_cast< size_t >( parameterFile.iParam[ "Ny" ] );
  double Lx = parameterFile.dParam[ "Lx" ];
  double Ly = parameterFile.dParam[ "Ly" ];
  double nu = parameterFile.dParam[ "nu" ];
  std::string testWriterFile = parameterFile.strParam[ "dataFile" ];
  
  int nOutX = std::floor( Nx / 2 + 1 );
  int nOutY = std::floor( Ny / 2 + 1 );
  double advFac = 2.0 * pi / Lx;
  double diffFac = 4.0 * pi * pi * deltaT * nu / ( Lx * Lx );
  
  std::vector< std::vector< double > > T0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > dT0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > NL0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > T1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< std::complex< double > > >
    T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    dT_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    NL_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  hydroCode::FourierTransforms fft;
  
  io::ioNetCDF testWriter( testWriterFile, "data", Nx, Ny, 'w' );
  
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
        dT_spec[ i ][ j ] = advFac * i * std::complex< double >( 0.0, 1.0 ) * T_spec[ i ][ j ];
      }
    }
    
    for( int i = nOutX ; i < Nx ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        dT_spec[ i ][ j ] = -advFac * ( Nx - i ) * std::complex< double >( 0.0, 1.0 ) * T_spec[ i ][ j ];
      }
    }
    
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, dT_spec, dT0_phys );
    
    for( int i = 0; i < Nx; i++ )
    {
      for( int j = 0; j < Ny; j++ )
      {
        dT0_phys[ i ][ j ] /= ( Nx * Ny );
      }
    }
    
    for( int i = 0; i < Nx; i++ )
    {
      for( int j = 0; j < Ny; j++ )
      {
        NL0_phys[ i ][ j ] = T0_phys[ i ][ j ] * dT0_phys[ i ][ j ];
      }
    }
    
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, NL0_phys, NL_spec );
    
    for( int i = 0 ; i < nOutX ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        T_spec[ i ][ j ] = ( T_spec[ i ][ j ] - deltaT * NL_spec[ i ][ j ] )
                           / ( 1.0 + diffFac * i * i );
      }
    }
    
    for( int i = nOutX ; i < Nx ; i++ )
    {
      for( int j = 0 ; j < nOutY ; j++ )
      {
        T_spec[ i ][ j ] = ( T_spec[ i ][ j ] - deltaT * NL_spec[ i ][ j ] )
                           / ( 1.0 + diffFac * ( Nx - i ) * ( Nx - i ) );
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
  
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();
  std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
            << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
            << "Wall clock time passed: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n" << std::endl;
}

} // diagnostics
