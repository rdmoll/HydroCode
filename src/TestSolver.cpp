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
  kappa = parameterFile.dParam[ "kappa" ];
  testWriterFile = parameterFile.strParam[ "dataFile" ];
  
  nOutX = floor( Nx / 2 + 1 );
  nOutY = floor( Ny / 2 + 1 );
  totN = Nx * Ny;
}

TestSolver::~TestSolver()
{
}

void TestSolver::setInitConditions( Scalar2D< double >& T_phys,
                                    Scalar2D< double >& u_phys,
                                    Scalar2D< double >& v_phys )
{
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < Nx; ++j )
    {
      T_phys( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::sin( j * ( 2.0 * M_PI / Nx ) );
      u_phys( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::cos( j * ( 2.0 * M_PI / Nx ) );
      v_phys( i, j ) = 0.0 * std::sin( i * ( 2.0 * M_PI / Ny ) ) * std::sin( j * ( 4.0 * M_PI / Nx ) );
    }
  }
}

void TestSolver::runSimulation()
{
  std::cout << "Running Simulation" << std::endl;
  std::cout << std::endl;
  
  // Start timer
  std::clock_t c_start = std::clock();
  auto t_start = std::chrono::high_resolution_clock::now();
  
  // Initialize variables
  Scalar2D< double > T_phys( Ny, Nx );
  Scalar2D< double > u_phys( Ny, Nx );
  Scalar2D< double > v_phys( Ny, Nx );
  
  Scalar2D< std::complex< double > > T_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > u_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > v_spec( Ny, nOutX );
  
  // Initialize data writer
  io::ioNetCDF testWriter( testWriterFile, Nx, Ny, 'w' );
  
  // Set initial conditions
  setInitConditions( T_phys, u_phys, v_phys );
  
  // Write initial conditions
  testWriter.write_T( 0, T_phys );
  testWriter.write_u( 0, u_phys );
  testWriter.write_v( 0, v_phys );
  
  // Start time step loop
  for( size_t t = 0; t < nSteps; t++)
  {
    // Transform: phys --> spec
    fft::fft_r2c_2d( T_phys, T_spec );
    fft::fft_r2c_2d( u_phys, u_spec );
    fft::fft_r2c_2d( v_phys, v_spec );
    
    // Transform: spec --> phys
    fft::fft_c2r_2d( T_spec, T_phys );
    fft::scaleOutput( T_phys );
    
    fft::fft_c2r_2d( u_spec, u_phys );
    fft::scaleOutput( u_phys );
    
    fft::fft_c2r_2d( v_spec, v_phys );
    fft::scaleOutput( v_phys );
    
    // Write data to file
    testWriter.write_T( t + 1, T_phys );
    testWriter.write_u( t + 1, u_phys );
    testWriter.write_v( t + 1, v_phys );
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
