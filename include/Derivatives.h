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
  void deriv2( std::vector< std::complex< double > > &inOut );
  void gradient();
  void laplacian();
};

}

#endif /* Derivatives_h */
