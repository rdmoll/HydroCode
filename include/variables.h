//
//  variables.h
//  HydroCode
//
//  Created by Ryan Moll on 2/4/21.
//

#ifndef variables_h
#define variables_h

#include <iostream>
#include <vector>
#include <complex>

namespace variables
{

class ScalarVar
{
public:
  ScalarVar( size_t Nx, size_t Ny )
  {
    int nOutY = std::floor( Ny / 2 + 1 );
    
    time0.resize( Nx );
    time1.resize( Nx );
    time2.resize( Nx );
    
    spec.resize( Nx );
    
    temp.resize( Nx );
    
    for( size_t i = 0; i < Ny; ++i )
    {
      time0[ i ].resize( Ny, 0.0 );
      time1[ i ].resize( Ny, 0.0 );
      time2[ i ].resize( Ny, 0.0 );
      
      spec[ i ].resize( nOutY, std::complex< double >( 0.0, 0.0 ) );
      
      temp[ i ].resize( Ny, 0.0 );
    }
  }
  ~ScalarVar(){};
  
  std::vector< std::vector< double > > time0;
  std::vector< std::vector< double > > time1;
  std::vector< std::vector< double > > time2;
  std::vector< std::vector< std::complex< double > > > spec;
  
  std::vector< std::vector< double > > temp;
};

class VectorVar
{
public:
  VectorVar( size_t Nx, size_t Ny )
  {
    int nOutY = std::floor( Ny / 2 + 1 );
    
    xTime0.resize( Nx );
    xTime1.resize( Nx );
    xTime2.resize( Nx );
    xSpec.resize( Nx );
    
    yTime0.resize( Nx );
    yTime1.resize( Nx );
    yTime2.resize( Nx );
    ySpec.resize( Nx );
    
    temp.resize( Nx );
    
    for( size_t i = 0; i < Ny; ++i )
    {
      xTime0[ i ].resize( Ny, 0.0 );
      xTime1[ i ].resize( Ny, 0.0 );
      xTime2[ i ].resize( Ny, 0.0 );
      xSpec[ i ].resize( nOutY, std::complex< double >( 0.0, 0.0 ) );
      
      yTime0[ i ].resize( Ny, 0.0 );
      yTime1[ i ].resize( Ny, 0.0 );
      yTime2[ i ].resize( Ny, 0.0 );
      ySpec[ i ].resize( nOutY, std::complex< double >( 0.0, 0.0 ) );
      
      temp[ i ].resize( Ny, 0.0 );
    }
  }
  ~VectorVar(){};
  
  std::vector< std::vector< double > > xTime0;
  std::vector< std::vector< double > > xTime1;
  std::vector< std::vector< double > > xTime2;
  std::vector< std::vector< std::complex< double > > > xSpec;
  
  std::vector< std::vector< double > > yTime0;
  std::vector< std::vector< double > > yTime1;
  std::vector< std::vector< double > > yTime2;
  std::vector< std::vector< std::complex< double > > > ySpec;
  
  std::vector< std::vector< double > > temp;
};

}

#endif /* variables_h */
