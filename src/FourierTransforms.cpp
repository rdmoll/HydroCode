//
//  FourierTransforms.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/13/20.
//

#include "FourierTransforms.h"

namespace hydroCode
{

FourierTransforms::FourierTransforms()
{
  realArray = NULL;
  compArray = NULL;
  plan_r2c = NULL;
  plan_c2r = NULL;
}

FourierTransforms::FourierTransforms( int N )
{
  int nOut = std::floor( N / 2 + 1 );
  
  realArray = ( double* ) fftw_malloc( sizeof( double ) * N );
  for( size_t i = 0; i < N; i++ )
  {
    realArray[ i ] = 0.0;
  }
  
  compArray = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * nOut );
  for( size_t i = 0; i < nOut; i++ )
  {
    compArray[ i ][ 0 ] = 0.0;
    compArray[ i ][ 1 ] = 0.0;
  }
  
  plan_r2c = fftw_plan_dft_r2c_1d( N, &realArray[0], &compArray[0], FFTW_ESTIMATE );
  plan_c2r = fftw_plan_dft_c2r_1d( N, &compArray[0], &realArray[0], FFTW_ESTIMATE );
}

FourierTransforms::~FourierTransforms()
{
  if( plan_r2c != NULL )
  {
    fftw_destroy_plan( plan_r2c );
  }
  
  if( plan_c2r != NULL )
  {
    fftw_destroy_plan( plan_c2r );
  }
  
  if( realArray != NULL )
  {
    fftw_free( realArray );
  }
  
  if( compArray != NULL )
  {
    fftw_free( compArray );
  }
}

void FourierTransforms::fft_r2c( )
{
  try
  {
    fftw_execute( plan_r2c );
  }
  catch ( std::exception& e )
  {
    std::cout << "Error!" << std::endl;
  }
}

void FourierTransforms::fft_r2c( std::vector< double > &rVec, std::vector< std::complex< double > > &cVec )
{
  int N = static_cast< int >( rVec.size() );
  int nOut = std::floor( N / 2 + 1 );
  
  double* in = ( double* ) fftw_malloc( sizeof( double ) * N );
  fftw_complex* out = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * nOut );
  
  for( int i = 0; i < N; i++ )
  {
    in[ i ] = rVec[ i ];
  }

  fftw_plan plan;
  plan = fftw_plan_dft_r2c_1d( N, &in[0], &out[0], FFTW_ESTIMATE );
  fftw_execute( plan );
  
  for( int i = 0 ; i < nOut ; i++ )
  {
    cVec[ i ] = std::complex< double >( out[ i ][ 0 ], out[ i ][ 1 ] );
  }
  
  fftw_destroy_plan( plan );
  fftw_free( in );
  fftw_free( out );
}

void FourierTransforms::fft_r2c_2d( std::vector< double > &rVec,
                                    std::vector< std::vector< std::complex< double > > > &cVec )
{
}

void FourierTransforms::fft_c2r( )
{
  try
  {
    fftw_execute( plan_c2r );
  }
  catch ( std::exception& e )
  {
    std::cout << "Error!" << std::endl;
  }
}

void FourierTransforms::fft_c2r( std::vector< std::complex< double > > &cVec, std::vector< double > &rVec )
{
  size_t nIn = cVec.size();
  int N = static_cast< int >( 2 * ( nIn - 1 ) );
  
  fftw_complex* in = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * nIn );
  double* out = (double*) fftw_malloc(sizeof(double) * N );
  
  for( int i = 0 ; i < nIn ; i++ )
  {
    in[ i ][ 0 ] = cVec[ i ].real();
    in[ i ][ 1 ] = cVec[ i ].imag();
  }
  
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_1d( N, &in[0], &out[0], FFTW_ESTIMATE );
  fftw_execute( plan );
  
  for( int i = 0 ; i < N ; i++ )
  {
    rVec[ i ] = out[ i ];
  }
  
  fftw_destroy_plan( plan );
  fftw_free( in );
  fftw_free( out );
}

} // hydroCode
