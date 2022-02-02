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
  Nx = static_cast< size_t >( 8 ); //static_cast< size_t >( parameterFile.iParam[ "Nx" ] );
  Ny = static_cast< size_t >( 8 ); //static_cast< size_t >( parameterFile.iParam[ "Ny" ] );
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

void TestSolver::runSimulation()
{
  std::cout << "Run Simulation" << std::endl;
  
  Scalar2D< double > T_phys( Nx, Ny );
  Scalar2D< double > u_phys( Nx, Ny );
  Scalar2D< double > v_phys( Nx, Ny );
  
  Scalar2D< std::complex< double > > T_spec( Nx, nOutY );
  Scalar2D< std::complex< double > > u_spec( Nx, nOutY );
  Scalar2D< std::complex< double > > v_spec( Nx, nOutY );
  
  Scalar2D< double > realTruth( Nx, Ny );
  
  for( int i = 0; i < Ny; i++ )
  {
    for( int j = 0; j < Nx; j++ )
    {
      realTruth( i, j ) = std::cos( j * ( 2 * pi / Nx ) ) * std::cos( i * ( 2 * pi / Ny ) );
    }
  }
  
  for( int i = 0; i < Ny; i++ )
  {
    for( int j = 0; j < Nx; j++ )
    {
      std::cout << realTruth( i, j ) << "  ";
    }
    std::cout << std::endl;
  }
}

}
