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

template< class T >
class ScalarVar
{
public:
  ScalarVar( size_t Nx, size_t Ny )
  {
    time0.resize( Nx );
    time1.resize( Nx );
    time2.resize( Nx );
    temp.resize( Nx );
    for( size_t i = 0; i < Ny; ++i )
    {
      time0[ i ].resize( Ny, 0.0 );
      time1[ i ].resize( Ny, 0.0 );
      time2[ i ].resize( Ny, 0.0 );
      temp[ i ].resize( Ny, 0.0 );
    }
  }
  ~ScalarVar(){};
  
  std::vector< std::vector< T > > time0;
  std::vector< std::vector< T > > time1;
  std::vector< std::vector< T > > time2;
  
protected:
  std::vector< std::vector< T > > temp;
};

template< class T >
class VectorVar
{
public:
  VectorVar( size_t Nx, size_t Ny )
  {
    xTime0.resize( Nx );
    xTime1.resize( Nx );
    xTime2.resize( Nx );
    yTime0.resize( Nx );
    yTime1.resize( Nx );
    yTime2.resize( Nx );
    temp.resize( Nx );
    for( size_t i = 0; i < Ny; ++i )
    {
      xTime0[ i ].resize( Ny, 0.0 );
      xTime1[ i ].resize( Ny, 0.0 );
      xTime2[ i ].resize( Ny, 0.0 );
      yTime0[ i ].resize( Ny, 0.0 );
      yTime1[ i ].resize( Ny, 0.0 );
      yTime2[ i ].resize( Ny, 0.0 );
      temp[ i ].resize( Ny, 0.0 );
    }
  }
  ~VectorVar(){};
  
  std::vector< std::vector< T > > xTime0;
  std::vector< std::vector< T > > xTime1;
  std::vector< std::vector< T > > xTime2;
  
  std::vector< std::vector< T > > yTime0;
  std::vector< std::vector< T > > yTime1;
  std::vector< std::vector< T > > yTime2;
  
protected:
  std::vector< std::vector< T > > temp;
};

}

#endif /* variables_h */
