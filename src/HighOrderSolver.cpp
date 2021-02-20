//
//  HighOrderSolver.cpp
//  HydroCode
//
//  Created by Ryan Moll on 2/6/21.
//

#include "HighOrderSolver.h"

namespace solvers
{

HighOrderSolver::HighOrderSolver( std::string paramFile )
{
  // Read parameters from file
  const char* delimiter = "=";
  io::ReadParams parameterFile( paramFile.c_str(), delimiter );
  
  nSteps = static_cast< size_t >( parameterFile.iParam[ "nSteps" ] );
  deltaT = parameterFile.dParam[ "deltaT" ];
  Nx = static_cast< size_t >( parameterFile.iParam[ "Nx" ] );
  Ny = static_cast< size_t >( parameterFile.iParam[ "Ny" ] );
  Lx = parameterFile.dParam[ "Lx" ];
  Ly = parameterFile.dParam[ "Ly" ];
  nu = parameterFile.dParam[ "nu" ];
  kappaT = parameterFile.dParam[ "kappaT" ];
  kappaC = parameterFile.dParam[ "kappaC" ];
  testWriterFile = parameterFile.strParam[ "dataFile" ];
  
  nuFac = 4.0 * pi * pi * deltaT * nu;
  kappaFacT = 4.0 * pi * pi * deltaT * kappaT;
  kappaFacC = 4.0 * pi * pi * deltaT * kappaC;
  
  nOutX = std::floor( Nx / 2 + 1 );
  nOutY = std::floor( Ny / 2 + 1 );
  totN = Nx * Ny;
  
  dfdx_spec.resize( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  dfdy_spec.resize( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  dfdx_phys.resize( Nx, std::vector< double >( Ny, 0.0 ) );
  dfdy_phys.resize( Nx, std::vector< double >( Ny, 0.0 ) );
  NL_xPhys.resize( Nx, std::vector< double >( Ny, 0.0 ) );
  NL_yPhys.resize( Nx, std::vector< double >( Ny, 0.0 ) );
  NL_xSpec.resize( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  NL_ySpec.resize( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  wvNumX.resize( Nx );
  wvNumY.resize( nOutY );
  
  for( size_t i = 0 ; i < nOutX ; ++i )
  {
    wvNumX[ i ] = ( i * i ) / ( Lx * Lx );
  }
  
  for( size_t i = nOutX ; i < Nx ; ++i )
  {
    wvNumX[ i ] = ( ( Nx - i ) * ( Nx - i ) ) / ( Lx * Lx );
  }
  
  for( size_t i = 0 ; i < nOutY ; ++i )
  {
    wvNumY[ i ] = ( i * i ) / ( Ly * Ly );
  }
}

HighOrderSolver::~HighOrderSolver()
{
}

void HighOrderSolver::setInitConditions( std::vector< std::vector< double > >& T0_phys,
                                         std::vector< std::vector< double > >& C0_phys,
                                         std::vector< std::vector< double > >& u0_phys,
                                         std::vector< std::vector< double > >& v0_phys )
{
  for( int i = 0; i < Nx; ++i )
  {
    for( int j = 0; j < Ny; ++j )
    {
      T0_phys[ i ][ j ] = std::cos( i * ( 2.0 * pi / Nx ) ) * std::sin( j * ( 2.0 * pi / Ny ) );
      C0_phys[ i ][ j ] = std::cos( i * ( 2.0 * pi / Nx ) ) * std::sin( j * ( 2.0 * pi / Ny ) );
      u0_phys[ i ][ j ] = std::cos( i * ( 2.0 * pi / Nx ) ) * std::cos( j * ( 2.0 * pi / Nx ) );
      v0_phys[ i ][ j ] = 0.0 * std::sin( i * ( 2.0 * pi / Nx ) ) * std::sin( j * ( 4.0 * pi / Ny ) );
    }
  }
}

void HighOrderSolver::calcNonLin( variables::VectorVar& u_phys,
                                  std::vector< std::vector< std::complex< double > > >& f_spec,
                                  std::vector< std::vector< std::complex< double > > >& NL_spec )
{
  ops.calcDerivX( f_spec, dfdx_spec, Nx, Ny, Lx );
  ops.calcDerivY( f_spec, dfdy_spec, Nx, Ny, Ly );
  
  fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, dfdx_spec, dfdx_phys );
  fft.scaleOutput( ( int ) Nx, ( int ) Ny, dfdx_phys );
  
  fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, dfdy_spec, dfdy_phys );
  fft.scaleOutput( ( int ) Nx, ( int ) Ny, dfdy_phys );
  
  for( int i = 0; i < Nx; ++i )
  {
    for( int j = 0; j < Ny; ++j )
    {
      NL_xPhys[ i ][ j ] = u_phys.xTime1[ i ][ j ] * dfdx_phys[ i ][ j ];
      NL_yPhys[ i ][ j ] = u_phys.yTime1[ i ][ j ] * dfdy_phys[ i ][ j ];
    }
  }
  
  fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, NL_xPhys, NL_xSpec );
  fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, NL_yPhys, NL_ySpec );
  
  for( int i = 0 ; i < Nx ; ++i )
  {
    for( int j = 0 ; j < nOutY ; ++j )
    {
      NL_spec[ i ][ j ] = NL_xSpec[ i ][ j ] + NL_ySpec[ i ][ j ];
    }
  }
}

void HighOrderSolver::solve( std::vector< std::vector< std::complex< double > > >& f_spec,
                             std::vector< std::vector< std::complex< double > > >& nl_spec,
                             double diffFac )
{
  for( int i = 0 ; i < nOutX ; ++i )
  {
    for( int j = 0 ; j < nOutY ; ++j )
    {
      f_spec[ i ][ j ] = ( f_spec[ i ][ j ] - deltaT * nl_spec[ i ][ j ] )
                         / ( 1.0 + diffFac * ( wvNumX[ i ] + wvNumY[ j ] ) );
    }
  }
  
  for( int i = nOutX ; i < Nx ; ++i )
  {
    for( int j = 0 ; j < nOutY ; ++j )
    {
      f_spec[ i ][ j ] = ( f_spec[ i ][ j ] - deltaT * nl_spec[ i ][ j ] )
                         / ( 1.0 + diffFac * ( wvNumX[ i ] + wvNumY[ j ] ) );
    }
  }
}

void HighOrderSolver::resetTimePointers( variables::ScalarVar& scalar )
{
  scalar.temp = scalar.time0;
  scalar.time0 = scalar.time1;
  scalar.time1 = scalar.temp;
}

void HighOrderSolver::resetTimePointers( variables::VectorVar& vec )
{
  vec.temp = vec.xTime0;
  vec.xTime0 = vec.xTime1;
  vec.xTime1 = vec.temp;
  
  vec.temp = vec.yTime0;
  vec.yTime0 = vec.yTime1;
  vec.yTime1 = vec.temp;
}

void HighOrderSolver::runSimulation()
{
  std::cout << "Running solver..." << std::endl;
  
  // Start timer
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  // Initialize data writer
  io::ioNetCDF testWriter( testWriterFile, Nx, Ny, 'w' );
  
  // Initialize arrays
  variables::ScalarVar T( Nx, Ny );
  variables::ScalarVar C( Nx, Ny );
  variables::VectorVar u( Nx, Ny );
  
  std::vector< std::vector< std::complex< double > > >
    NL_T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    NL_C_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    NL_ux_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    NL_uy_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  // Set initial conditions
  setInitConditions( T.time1, C.time1, u.xTime1, u.yTime1 );
  testWriter.write_T( 0, T.time1 );
  testWriter.write_C( 0, C.time1 );
  testWriter.write_u( 0, u.xTime1 );
  testWriter.write_v( 0, u.yTime1 );
  
  // Start time step loop
  for( size_t t = 0; t < nSteps; ++t )
  {
    // Transform: phys --> spec
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, T.time1, T.spec );
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, C.time1, C.spec );
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, u.xTime1, u.xSpec );
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, u.yTime1, u.ySpec );
    
    // Calculate Non-linear terms
    calcNonLin( u, T.spec, NL_T_spec );
    calcNonLin( u, C.spec, NL_C_spec );
    calcNonLin( u, u.xSpec, NL_ux_spec );
    calcNonLin( u, u.ySpec, NL_uy_spec );
    
    // Calculate New time step
    solve( T.spec, NL_T_spec, kappaFacT );
    solve( C.spec, NL_C_spec, kappaFacC );
    solve( u.xSpec, NL_ux_spec, nuFac );
    solve( u.ySpec, NL_uy_spec, nuFac );
    
    // Transform: spec --> phys
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, T.spec, T.time0 );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, T.time0 );
    
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, C.spec, C.time0 );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, C.time0 );
    
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, u.xSpec, u.xTime0 );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, u.xTime0 );
    
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, u.ySpec, u.yTime0 );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, u.yTime0 );
    
    // Write to file
    testWriter.write_T( t + 1, T.time0 );
    testWriter.write_C( t + 1, C.time0 );
    testWriter.write_u( t + 1, u.xTime0 );
    testWriter.write_v( t + 1, u.yTime0 );
    
    // Reset pointers
    resetTimePointers( T );
    resetTimePointers( C );
    resetTimePointers( u );
  }
  
  // End timer
  std::clock_t c_end = std::clock();
  auto t_end = std::chrono::high_resolution_clock::now();
  
  // Print completion message and timing information
  std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
            << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
            << "Wall clock time passed: "
            << std::chrono::duration<double, std::milli>(t_end-t_start).count()
            << " ms\n" << std::endl;
}

}
