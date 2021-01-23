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
    }
  }
}

void TestSolver::calcDeriv( std::vector< std::vector< std::complex< double > > >& f_spec,
                            std::vector< std::vector< std::complex< double > > >& df_spec )
{
  
}

void TestSolver::solveU(std::vector< std::vector< std::complex< double > > >& u_spec, double diffFac )
{
  
}

void TestSolver::solveT( std::vector< std::vector< std::complex< double > > >& T_spec,
                         double advFac, double diffFac )
{
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
}

void TestSolver::runSimulation()
{
  std::cout << "Running solver..." << std::endl;
  
  // Start timer
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  // Initialize data writer
  io::ioNetCDF testWriter( testWriterFile, "data", Nx, Ny, 'w' );
  
  // Define constants
  double c = 0.1;
  double advFac = 2.0 * pi * deltaT * c / Lx;
  double diffFac = 4.0 * pi * pi * deltaT * nu / ( Lx * Lx );
  
  // Initialize arrays
  std::vector< std::vector< double > > T0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > u0_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > T1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > u1_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > N_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  std::vector< std::vector< double > > L_phys( Nx, std::vector< double >( Ny, 0.0 ) );
  
  std::vector< std::vector< std::complex< double > > >
    T_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    u_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    dTdx_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    dudx_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    N_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  std::vector< std::vector< std::complex< double > > >
    L_spec( Nx, std::vector< std::complex< double > >( nOutY, std::complex< double >( 0.0, 0.0 ) ) );
  
  // Set initial conditions
  setInitConditions( T0_phys, u0_phys );
  testWriter.write( 0, T0_phys );
  
  // Start time step loop
  for( size_t t = 0; t < nSteps; t++)
  {
    // Transform T0:    phys --> spec
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, T0_phys, T_spec );
    
    // Transform u0:    phys --> spec
    fft.fft_r2c_2d( ( int ) Nx, ( int ) Ny, u0_phys, u_spec );
    
    // Calculate dTdx in spectral space
    //calcDeriv( T_spec, dTdx_spec );
    
    // Calculate dudx in spectral space
    //calcDeriv( u_spec, dudx_spec );
    
    solveT( T_spec, advFac, diffFac );
    
    // Transform du0dx: spec --> phys
    fft.fft_c2r_2d( ( int ) Nx, ( int ) Ny, T_spec, T1_phys );
    fft.scaleOutput( ( int ) Nx, ( int ) Ny, T1_phys );
    
    // Calculate NL = u0*du0dx in physical space
    // Transform NL:    phys --> spec
    // Calculate u1 in spectral space
    // Transform u1:    spec --> phys
    // set: u0 = u1, um1 = u0, um2 = um1, um3 = um2
    
    // Write u1 to file
    testWriter.write( t + 1, T0_phys );
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
