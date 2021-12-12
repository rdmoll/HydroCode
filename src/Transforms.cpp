#include "../include/Transforms.h"

namespace fft
{

void fft_r2c_1d( Scalar1D< double > &rVec, Scalar1D< std::complex< double > > &cVec )
{
  int N = static_cast< int >( rVec.size() );
  int nOut = std::floor( N / 2 + 1 );
  
  double* in = rVec.get();
  fftw_complex* out = reinterpret_cast< fftw_complex* >( cVec.get() );

  fftw_plan plan = fftw_plan_dft_r2c_1d( N, &in[0], &out[0], FFTW_ESTIMATE );
  fftw_execute( plan );
  
  fftw_destroy_plan( plan );
}

void fft_c2r_1d( Scalar1D< std::complex< double > > &cVec, Scalar1D< double > &rVec )
{
  int nIn = static_cast< int >( cVec.size() );
  int N = static_cast< int >( 2 * ( nIn - 1 ) );
  
  fftw_complex* in = reinterpret_cast< fftw_complex* >( cVec.get() );
  double* out = rVec.get();
  
  fftw_plan plan = fftw_plan_dft_c2r_1d( N, &in[0], &out[0], FFTW_ESTIMATE );
  fftw_execute( plan );
  
  fftw_destroy_plan( plan );
}

void fft_r2c_2d( Scalar2D< double > &rVec, Scalar2D< std::complex< double > > &cVec )
{
  int N = static_cast< int >( rVec.rows() );
  int M = static_cast< int >( rVec.cols() );
  
  double* in = rVec.get();
  fftw_complex* out = reinterpret_cast< fftw_complex* >( cVec.get() );
  
  fftw_plan plan = fftw_plan_dft_r2c_2d( N, M, in, out, FFTW_ESTIMATE );
  fftw_execute( plan );
  
  fftw_destroy_plan( plan );
}

void fft_c2r_2d( Scalar2D< std::complex< double > > &cVec, Scalar2D< double > &rVec )
{
  int N = static_cast< int >( rVec.rows() );
  int M = static_cast< int >( rVec.cols() );
  
  fftw_complex* in = reinterpret_cast< fftw_complex* >( cVec.get() );
  double* out = rVec.get();
  
  fftw_plan plan = fftw_plan_dft_c2r_2d( N, M, in, out, FFTW_ESTIMATE );
  fftw_execute( plan );
  
  fftw_destroy_plan( plan );
}

void scaleOutput( Scalar2D< double > &f_phys )
{
  int N = static_cast< int >( f_phys.rows() );
  int M = static_cast< int >( f_phys.cols() );
  
  for( size_t i = 0; i < N; i++ )
  {
    for( size_t j = 0; j < M; j++ )
    {
      f_phys( i, j ) /= ( N * M );
    }
  }
}

} // fft
