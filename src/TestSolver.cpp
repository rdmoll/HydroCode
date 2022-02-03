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
      T_phys( i, j ) = std::cos( i * ( 2.0 * pi / Ny ) ) * std::sin( j * ( 2.0 * pi / Nx ) );
      u_phys( i, j ) = std::cos( i * ( 2.0 * pi / Ny ) ) * std::cos( j * ( 2.0 * pi / Nx ) );
      v_phys( i, j ) = 0.0 * std::sin( i * ( 2.0 * pi / Ny ) ) * std::sin( j * ( 4.0 * pi / Nx ) );
    }
  }
}

void TestSolver::runSimulation()
{
  std::cout << "Run Simulation" << std::endl;
  
  Scalar2D< double > T_phys( Ny, Nx );
  Scalar2D< double > u_phys( Ny, Nx );
  Scalar2D< double > v_phys( Ny, Nx );
  
  Scalar2D< std::complex< double > > T_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > u_spec( Ny, nOutX );
  Scalar2D< std::complex< double > > v_spec( Ny, nOutX );
  
  // Set initial conditions
  setInitConditions( T_phys, u_phys, v_phys );
  
  // Start time step loop
  for( size_t t = 0; t < nSteps; t++)
  {
  }
}

}
