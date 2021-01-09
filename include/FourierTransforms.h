//
//  FourierTransforms.h
//  HydroCode
//
//  Created by Ryan Moll on 12/13/20.
//

#ifndef FourierTransforms_h
#define FourierTransforms_h

#include <iostream>
#include <complex>
#include <vector>
#include <exception>
#include "fftw3.h"

namespace hydroCode
{

class FourierTransforms
{
public:
  FourierTransforms();
  FourierTransforms( int N );
  
  ~FourierTransforms();
  
  double* realArray;
  fftw_complex* compArray;
  
  void fft_r2c( std::vector< double > &rVec, std::vector< std::complex< double > > &cVec );
  void fft_r2c( );
  void fft_r2c_2d( int N, int M,
                   std::vector< std::vector< double > > &rVec,
                   std::vector< std::vector< std::complex< double > > > &cVec );
  
  void fft_c2r( std::vector< std::complex< double > > &cVec, std::vector< double > &rVec );
  void fft_c2r( );
  void fft_c2r_2d( int N, int M,
                   std::vector< std::vector< std::complex< double > > > &cVec,
                   std::vector< std::vector< double > > &rVec );
  
protected:
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;
};

} // fluidSim

#endif /* FourierTransforms_h */
