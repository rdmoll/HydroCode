//
//  Scalar3D.h
//  HydroCode
//
//  Created by Ryan Moll on 2/3/22.
//

#ifndef Scalar3D_h
#define Scalar3D_h

#include <iostream>
#include <memory>

template< class T >
class Scalar3D
{
public:
  Scalar3D(){};
  Scalar3D( std::size_t rowSize, std::size_t colSize, std::size_t sliceSize );
  ~Scalar3D( );
  
  // Container methods
  void set( std::size_t rowIndex, std::size_t colIndex, std::size_t sliceIndex, T value );
  void set( T* arrayPtr );
  void setSize( std::size_t rowSize, std::size_t colSize, std::size_t sliceSize );
  T* get() const;
  std::size_t rows();
  std::size_t cols();
  std::size_t slices();
  
  // Operators
  const T& operator()( const std::size_t rowIndex, const std::size_t colIndex, const std::size_t sliceIndex ) const;
  T& operator()( const std::size_t rowIndex, const std::size_t colIndex, const std::size_t sliceIndex );
  
protected:
  std::unique_ptr< T > _data;
  std::size_t _rowSize;
  std::size_t _colSize;
  std::size_t _sliceSize;
};

#include "Scalar3D.it.h"

#endif /* Scalar3D_h */
