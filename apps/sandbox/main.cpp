#include <iostream>
#include <complex>

#include "../../include/Scalar1D.h"
#include "../../include/Transforms.h"
#include "../../include/MathOps.h"

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
  
  for( std::size_t i = 0; i < nOut; i++ )
  {
    std::cout << out( i ) << std::endl;
  }
  
  mathOps::calcDeriv( out, df_spec, 1.0 * N );
  
  fft::fft_c2r_1d( df_spec, df_phys );
  
  for( std::size_t i = 0; i < N; i++ )
  {
    std::cout << df_phys( i ) << std::endl;
  }
  
  return 0;
}

