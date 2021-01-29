//
//  Derivatives.h
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#ifndef Derivatives_h
#define Derivatives_h

#include <iostream>
#include <vector>
#include <complex>

namespace mathOps
{

class Derivatives
{
public:
  Derivatives();
  ~Derivatives();
  
  void deriv( std::vector< std::complex< double > > &inOut );
  void calcDerivX( std::vector< std::vector< std::complex< double > > >& f_spec,
                   std::vector< std::vector< std::complex< double > > >& df_spec,
                   size_t Nx, size_t Ny );
  void calcDerivY( std::vector< std::vector< std::complex< double > > >& f_spec,
                   std::vector< std::vector< std::complex< double > > >& df_spec,
                   size_t Nx, size_t Ny );
  void deriv2( std::vector< std::complex< double > > &inOut );
  void gradient();
  void laplacian();
  
protected:
  double pi = std::acos(-1.0);
  
};

}

#endif /* Derivatives_h */
