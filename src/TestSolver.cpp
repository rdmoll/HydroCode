//
//  Solver.cpp
//  HydroCode
//
//  Created by Ryan Moll on 1/22/21.
//

#include "TestSolver.h"

namespace solvers
{

TestSolver::TestSolver( std::string paramFile )
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
  nu = parameterFile.dParam[ "kappa" ];
  testWriterFile = parameterFile.strParam[ "dataFile" ];
  
  nOutX = std::floor( Nx / 2 + 1 );
  nOutY = std::floor( Ny / 2 + 1 );
  totN = Nx * Ny;
}

TestSolver::~TestSolver()
{
}

void TestSolver::setInitConditions( std::vector< std::vector< double > >& T0_phys,
                                    std::vector< std::vector< double > >& u0_phys )
{
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      T0_phys[ i ][ j ] = std::cos( i * ( 2.0 * pi / Nx ) );
      u0_phys[ i ][ j ] = 0.1;
      //u0_phys[ i ][ j ] = 0.1 * std::sin( i * ( 2.0 * pi / Nx ) ) + 0.1;
    }
  }
}

void TestSolver::calcDerivX( std::vector< std::vector< std::complex< double > > >& f_spec,
                             std::vector< std::vector< std::complex< double > > >& df_spec )
{
  for( int i = 0 ; i < nOutX ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      df_spec[ i ][ j ] = 2.0 * pi * i * std::complex< double >( 0.0, 1.0 ) * f_spec[ i ][ j ];
    }
  }
  
  for( int i = nOutX ; i < Nx ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      df_spec[ i ][ j ] = -2.0 * pi * ( Nx - i ) * std::complex< double >( 0.0, 1.0 ) * f_spec[ i ][ j ];
    }
  }
}

void TestSolver::calcNonLin( std::vector< std::vector< double > >& f1_phys,
                             std::vector< std::vector< std::complex< double > > >& f2_spec,
                             std::vector< std::vector< std::complex< double > > >& nl_spec )
{
  std::vector< std::vector< std::complex< double > > >
    df2dx_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< double > > df2dx_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > nl_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  
  calcDerivX( f2_spec, df2dx_spec );
  
  fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, df2dx_spec, df2dx_phys );
  fft.scaleOutput( ( int ) Nx, ( int ) Ny, df2dx_phys );
  
  for( int i = 0; i < Nx; i++ )
  {
    for( int j = 0; j < Ny; j++ )
    {
      nl_phys[ i ][ j ] = f1_phys[ i ][ j ] * df2dx_phys[ i ][ j ];
    }
  }
  
  fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, nl_phys, nl_spec );
}

void TestSolver::solve( std::vector< std::vector< std::complex< double > > >& f_spec,
                        std::vector< std::vector< std::complex< double > > >& nl_spec,
                        double diffusivity )
{
  double diffFac = 4.0 * pi * pi * deltaT * diffusivity / ( Lx * Lx );
  
  for( int i = 0 ; i < nOutX ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      f_spec[ i ][ j ] = ( f_spec[ i ][ j ] - deltaT * nl_spec[ i ][ j ] )
                         / ( 1.0 + diffFac * i * i );
    }
  }
  
  for( int i = nOutX ; i < Nx ; i++ )
  {
    for( int j = 0 ; j < nOutY ; j++ )
    {
      f_spec[ i ][ j ] = ( f_spec[ i ][ j ] - deltaT * nl_spec[ i ][ j ] )
                         / ( 1.0 + diffFac * ( Nx - i ) * ( Nx - i ) );
    }
  }
}

void TestSolver::runSimulation()
{
  std::cout << "Running solver..." << std::endl;
  
  // Start timer
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  // Initialize data writer
  io::ioNetCDF testWriter( testWriterFile, Nx, Ny, 'w' );
  
  // Initialize arrays
  std::vector< std::vector< double > > T0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > u0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > T1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > u1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > temp( Nx, std::vector< double >( Ny, 0.0 ) );
  
  std::vector< std::vector< std::complex< double > > >
    T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    u_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    nl_T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    nl_u_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  // Set initial conditions
  setInitConditions( T0_phys, u0_phys );
  testWriter.write_T( 0, T0_phys );
  testWriter.write_u( 0, u0_phys );
  
  // Start time step loop
  for( size_t t = 0; t < nSteps; t++)
  {
    // Transform: phys --> spec
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, T0_phys, T_spec );
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, u0_phys, u_spec );
    
    // Calculate Non-linear terms
    calcNonLin( u0_phys, T_spec, nl_T_spec );
    calcNonLin( u0_phys, u_spec, nl_u_spec );
    
    // Calculate New time step
    solve( T_spec, nl_T_spec, kappa );
    solve( u_spec, nl_u_spec, nu );
    
    // Transform: spec --> phys
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, T_spec, T1_phys );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, T1_phys );
    
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, u_spec, u1_phys );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, u1_phys );
    
    // Write T1 to file
    testWriter.write_T( t + 1, T1_phys );
    testWriter.write_u( t + 1, u1_phys );
    
    // Reset pointers
    temp = T1_phys;
    T1_phys = T0_phys;
    T0_phys = temp;
    
    temp = u1_phys;
    u1_phys = u0_phys;
    u0_phys = temp;
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
