#ifndef Scalar2D_h
#define Scalar2D_h

template< class T >
class Scalar2D
{
public:
  Scalar2D(){};
  Scalar2D( std::size_t rowSize, std::size_t colSize );
  ~Scalar2D( );
  
  // Container methods
  void set( std::size_t rowIndex, std::size_t colIndex, T value );
  void set( T* arrayPtr );
  void setSize( std::size_t rowSize, std::size_t colSize );
  T* get() const;
  std::size_t rows();
  std::size_t cols();
  
  // Operators
  const T& operator()( const std::size_t rowIndex, const std::size_t colIndex ) const;
  T& operator()( const std::size_t rowIndex, const std::size_t colIndex );
  Scalar2D< T >& operator=( Scalar2D< T >& arr );
  
protected:
  std::unique_ptr< T > _data;
  std::size_t _rowSize;
  std::size_t _colSize;
};

// Multiplication

template< class T >
Scalar2D< T >& operator*( Scalar2D< T >& arr1, Scalar2D< T >& arr2 );

template< class T >
Scalar2D< T >& operator*( T& val, Scalar2D< T >& arr );

template< class T >
Scalar2D< T >& operator*( Scalar2D< T >& arr, T& val );

// Division

template< class T >
Scalar2D< T >& operator/( Scalar2D< T >& arr1, Scalar2D< T >& arr2 );

template< class T >
Scalar2D< T >& operator/( T& val, Scalar2D< T >& arr );

template< class T >
Scalar2D< T >& operator/( Scalar2D< T >& arr, T& val );

// Addition

template< class T >
Scalar2D< T >& operator+( Scalar2D< T >& arr1, Scalar2D< T >& arr2 );

template< class T >
Scalar2D< T >& operator+( T& val, Scalar2D< T >& arr );

template< class T >
Scalar2D< T >& operator+( Scalar2D< T >& arr, T& val );

// Subtraction

template< class T >
Scalar2D< T >& operator-( Scalar2D< T >& arr1, Scalar2D< T >& arr2 );

template< class T >
Scalar2D< T >& operator-( T& val, Scalar2D< T >& arr );

template< class T >
Scalar2D< T >& operator-( Scalar2D< T >& arr, T& val );

#include "Scalar2D.it.h"

#endif /* Scalar2D_h */
