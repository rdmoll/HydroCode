#include "TestSolver.h"
#include "MathOps.h"

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
  kappa = parameterFile.dParam[ "kappaT" ];
  testWriterFile = parameterFile.strParam[ "dataFile" ];
  
  nOutX = floor( Nx / 2 + 1 );
  nOutY = floor( Ny / 2 + 1 );
  totN = Nx * Ny;
  
  fourPiSq = 4.0 * M_PI * M_PI;
  
  totWvNum.setSize( Ny, nOutX );
  for( int i = 0 ; i < nOutY ; i++ )
  {
    for( int j = 0 ; j < nOutX ; j++ )
    {
      totWvNum( i, j ) = ( i * i ) / ( Ly * Ly ) + ( j * j ) / ( Lx * Lx );
    }
  }
  
  for( int i = nOutY ; i < Ny ; i++ )
  {
    for( int j = 0 ; j < nOutX ; j++ )
    {
      totWvNum( i, j ) = ( ( Ny - i ) * ( Ny - i ) ) / ( Ly * Ly ) + ( j * j ) / ( Lx * Lx );
    }
  }
}

TestSolver::~TestSolver()
{
}

void TestSolver::setInitConditions( Scalar2D< double >& T_phys,
                                    Scalar2D< double >& u_phys,
                                    Scalar2D< double >& v_phys )
{
  std::srand( time( NULL ) );
  
  Scalar2D< std::complex< double > > init_T_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > init_u_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > init_v_spec( Ny, nOutX );
  
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < Nx; ++j )
    {
      T_phys( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::sin( j * ( 2.0 * M_PI / Nx ) );
      u_phys( i, j ) = -0.5 * std::cos( i * ( 2.0 * M_PI / Ny ) );
      v_phys( i, j ) = 0.05 * std::sin( j * ( 2.0 * M_PI / Nx ) );
    }
  }
  
  fft::fft_r2c_2d( T_phys, init_T_spec );
  fft::fft_r2c_2d( u_phys, init_u_spec );
    
  for( int i = 0; i < Ny; ++i )
  {
    for( int j = 0; j < nOutX; ++j )
    {
      double T_randValR = 2.0 * ( ( double ) rand() / ( RAND_MAX ) ) - 1;
      double T_randValC = 2.0 * ( ( double ) rand() / ( RAND_MAX ) ) - 1;
      init_T_spec( i, j ) += 0.01 * std::complex< double >( T_randValR, T_randValC );
      
      double u_randValR = 2.0 * ( ( double ) rand() / ( RAND_MAX ) ) - 1;
      double u_randValC = 2.0 * ( ( double ) rand() / ( RAND_MAX ) ) - 1;
      init_u_spec( i, j ) += 0.01 * std::complex< double >( u_randValR, u_randValC );
      
      double v_randValR = 2.0 * ( ( double ) rand() / ( RAND_MAX ) ) - 1;
      double v_randValC = 2.0 * ( ( double ) rand() / ( RAND_MAX ) ) - 1;
      init_v_spec( i, j ) = 0.01 * std::complex< double >( v_randValR, v_randValC );
    }
  }
  
  fft::fft_c2r_2d( init_T_spec, T_phys );
  fft::scaleOutput( T_phys );
  fft::fft_c2r_2d( init_u_spec, u_phys );
  fft::scaleOutput( u_phys );
  fft::fft_c2r_2d( init_v_spec, v_phys );
  fft::scaleOutput( v_phys );
}

void TestSolver::calcNL( Scalar2D< double >& f1_phys,
                         Scalar2D< double >& f2_phys,
                         Scalar2D< std::complex< double > >& f3_spec,
                         Scalar2D< std::complex< double > >& nl_spec )
{
  Scalar2D< std::complex< double > > df3dx_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > df3dy_spec( Ny, nOutX );
  Scalar2D< double > df3dx_phys( Ny, Nx );
  Scalar2D< double > df3dy_phys( Ny, Nx );
  Scalar2D< double > nl_x_phys( Ny, Nx );
  Scalar2D< double > nl_y_phys( Ny, Nx );
  Scalar2D< double > nl_phys( Ny, Nx );
  
  mathOps::calcDerivX( f3_spec, df3dx_spec, Nx, Ny, Lx );
  mathOps::calcDerivY( f3_spec, df3dy_spec, Nx, Ny, Ly );
  
  fft::fft_c2r_2d( df3dx_spec, df3dx_phys );
  fft::scaleOutput( df3dx_phys );
  
  fft::fft_c2r_2d( df3dy_spec, df3dy_phys );
  fft::scaleOutput( df3dy_phys );
  
  nl_x_phys = f1_phys * df3dx_phys;
  nl_y_phys = f2_phys * df3dy_phys;
  nl_phys = nl_x_phys + nl_y_phys;
  
  fft::fft_r2c_2d( nl_phys, nl_spec );
}

void TestSolver::solve( Scalar2D< std::complex< double > >& f_spec,
                        Scalar2D< std::complex< double > >& nl_spec,
                        const double timeStep,
                        const double diffusivity )
{
  const std::complex< double > diffFac = std::complex< double >( fourPiSq * deltaT * diffusivity );
  const std::complex< double > timeStepC = std::complex< double >( timeStep, 0.0 );
  const std::complex< double > oneC = std::complex< double >( 1.0, 0 );
  
  f_spec = ( f_spec - ( timeStepC * nl_spec ) ) / ( oneC + ( diffFac * totWvNum ) );
}

void TestSolver::runSimulation()
{
  std::cout << "-- Running Simulation --" << std::endl;
  std::cout << "Writing to file: " << testWriterFile << std::endl;
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
  
  Scalar2D< std::complex< double > > nl_T_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > nl_u_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > nl_v_spec( Ny, nOutX );
  
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
    
    // Calculate nonlinear terms
    calcNL( u_phys, v_phys, T_spec, nl_T_spec );
    calcNL( u_phys, v_phys, u_spec, nl_u_spec );
    calcNL( u_phys, v_phys, v_spec, nl_v_spec );
    
    // Solve next time step
    solve( T_spec, nl_T_spec, deltaT, kappa );
    solve( u_spec, nl_u_spec, deltaT, nu );
    solve( v_spec, nl_v_spec, deltaT, nu );
    
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
    
    // Increment time pointers
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
