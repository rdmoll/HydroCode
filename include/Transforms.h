#ifndef Transforms_h
#define Transforms_h

#include <iostream>
#include <complex>
#include <vector>
#include <exception>
#include "fftw3.h"
#include "Scalar1D.h"
#include "Scalar2D.h"

namespace fft
{

void fft_r2c_1d( Scalar1D< double > &rVec, Scalar1D< std::complex< double > > &cVec );
void fft_c2r_1d( Scalar1D< std::complex< double > > &cVec, Scalar1D< double > &rVec );
void fft_r2c_2d( Scalar2D< double > &rVec, Scalar2D< std::complex< double > > &cVec );
void fft_c2r_2d( Scalar2D< std::complex< double > > &cVec, Scalar2D< double > &rVec );

void fft_r2c_1d_thread( Scalar1D< double > &rVec, Scalar1D< std::complex< double > > &cVec );
void fft_c2r_1d_thread( Scalar1D< std::complex< double > > &cVec, Scalar1D< double > &rVec );
void fft_r2c_2d_thread( Scalar2D< double > &rVec, Scalar2D< std::complex< double > > &cVec );
void fft_c2r_2d_thread( Scalar2D< std::complex< double > > &cVec, Scalar2D< double > &rVec );

} // fft

#endif /* Transforms_h */
