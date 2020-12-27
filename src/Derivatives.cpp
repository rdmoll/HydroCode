//
//  Derivatives.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include "Derivatives.h"

namespace mathOps
{

Derivatives::Derivatives()
{
  
}

Derivatives::~Derivatives()
{
  
}

void Derivatives::deriv( std::vector< std::complex< double > > &inOut )
{
  double pi = std::acos(-1.0);
  
  for( size_t i = 0; i < inOut.size(); i++ )
  {
    inOut[ i ] *= 2.0 * pi * i * std::complex< double >( 0.0, 1.0 );
  }
}

void Derivatives::deriv2( std::vector< std::complex< double > > &inOut )
{
  double pi = std::acos(-1.0);
  
  for( size_t i = 0; i < inOut.size(); i++ )
  {
    inOut[ i ] *= -4.0 * pi * pi * i * i;
  }
}

void Derivatives::gradient()
{
  
}

void Derivatives::laplacian()
{
  
}

} // mathOps
