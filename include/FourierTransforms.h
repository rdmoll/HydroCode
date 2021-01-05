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
  
  void fft_r2c( std::vector< double > &inVec, std::vector< std::complex< double > > &outVec );
  void fft_r2c( );
  void fft_r2c_2d( std::vector< double > &rVec, std::vector< std::vector< std::complex< double > > > &cVec );
  
  void fft_c2r( std::vector< std::complex< double > > &inVec, std::vector< double > &outVec );
  void fft_c2r( );
  
protected:
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;
};

} // fluidSim

#endif /* FourierTransforms_h */
