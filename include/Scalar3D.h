#ifndef Scalar3D_h
#define Scalar3D_h

template< class T >
class Scalar3D
{
public:
  Scalar3D(){};
  Scalar3D( std::size_t sliceSize, std::size_t rowSize, std::size_t colSize );
  ~Scalar3D( );
  
  // Container methods
  void set( std::size_t sliceIndex, std::size_t rowIndex, std::size_t colIndex, T value );
  void set( T* arrayPtr );
  void setSize( std::size_t sliceSize, std::size_t rowSize, std::size_t colSize );
  T* get() const;
  std::size_t rows();
  std::size_t cols();
  std::size_t slices();
  
  // Operators
  const T& operator()( const std::size_t sliceIndex, const std::size_t rowIndex, const std::size_t colIndex ) const;
  T& operator()( const std::size_t sliceIndex, const std::size_t rowIndex, const std::size_t colIndex );
  Scalar3D< T >& operator=( Scalar3D< T >& arr );
  
protected:
  std::unique_ptr< T > _data;
  std::size_t _sliceSize;
  std::size_t _rowSize;
  std::size_t _colSize;
};

// Multiplication

template< class T >
Scalar3D< T >& operator*( Scalar3D< T >& arr1, Scalar3D< T >& arr2 );

template< class T >
Scalar3D< T >& operator*( T& val, Scalar3D< T >& arr );

template< class T >
Scalar3D< T >& operator*( Scalar3D< T >& arr, T& val );

// Division

template< class T >
Scalar3D< T >& operator/( Scalar3D< T >& arr1, Scalar3D< T >& arr2 );

template< class T >
Scalar3D< T >& operator/( T& val, Scalar3D< T >& arr );

template< class T >
Scalar3D< T >& operator/( Scalar3D< T >& arr, T& val );

// Addition

template< class T >
Scalar3D< T >& operator+( Scalar3D< T >& arr1, Scalar3D< T >& arr2 );

template< class T >
Scalar3D< T >& operator+( T& val, Scalar3D< T >& arr );

template< class T >
Scalar3D< T >& operator+( Scalar3D< T >& arr, T& val );

// Subtraction

template< class T >
Scalar3D< T >& operator-( Scalar3D< T >& arr1, Scalar3D< T >& arr2 );

template< class T >
Scalar3D< T >& operator-( T& val, Scalar3D< T >& arr );

template< class T >
Scalar3D< T >& operator-( Scalar3D< T >& arr, T& val );

#include "Scalar3D.it.h"

#endif /* Scalar3D_h */
