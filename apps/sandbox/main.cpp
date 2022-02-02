#include <iostream>
#include <complex>

#include "../../include/Scalar1D.h"
#include "../../include/Transforms.h"
//#include "../../include/MathOps.h"
#include "MathOps.h"
#include "TestSolver.h"

int main( int argc, const char * argv[] )
{
  std::size_t N = 16;
  std::size_t nOut = std::floor( N / 2 + 1 );
  
  Scalar1D< double > realTruth( N );
  Scalar1D< double > in( N );
  Scalar1D< std::complex< double > > out( nOut );
  Scalar1D< std::complex< double > > df_spec( nOut );
  Scalar1D< double > df_phys( N );
  
  for( std::size_t i = 0; i < N; i++ )
  {
    realTruth( i ) = std::cos( i * ( 2 * M_PI / N ) );
    in( i ) = std::cos( i * ( 2 * M_PI / N ) );
    
    std::cout << in( i ) << std::endl;
  }
  
  fft::fft_r2c_1d( in, out );
  
  std::cout << std::endl;
  for( std::size_t i = 0; i < nOut; i++ )
  {
    std::cout << out( i ) << std::endl;
  }
  
  mathOps::calcDeriv( out, df_spec );
  
  fft::fft_c2r_1d( df_spec, df_phys );
  
  std::cout << std::endl;
  for( std::size_t i = 0; i < N; i++ )
  {
    std::cout << df_phys( i ) << std::endl;
  }
  
  Scalar1D< double > test1( 3 );
  Scalar1D< double > test2( 3 );
  Scalar1D< double > test3( 3 );
  
  test1( 0 ) = 1.0;
  test1( 1 ) = 2.0;
  test1( 2 ) = 3.0;
  
  test2( 0 ) = 3.0;
  test2( 1 ) = 4.0;
  test2( 2 ) = 5.0;
  
  test3 = test1 * test2;
  
  std::cout << std::endl;
  
  for( std::size_t i = 0; i < 3; ++i )
  {
    std::cout << test3( i ) << std::endl;
  }
  
  std::cout << std::endl;
  
  std::cout << test3.size() << std::endl;
  
  solvers::TestSolver solver1( "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/advDiffParams.txt" );
  solver1.runSimulation();
  
  return 0;
}

