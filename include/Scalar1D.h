#ifndef Scalar1D_h
#define Scalar1D_h

template< class T >
class Scalar1D
{
public:
  Scalar1D( std::size_t arraySize );
  ~Scalar1D( );
  
  // Container methods
  void set( std::size_t index, T value );
  void set( T* arrayPtr );
  T* get() const;
  std::size_t size();
  
  // Operators
  const T& operator()( const std::size_t index ) const;
  T& operator()( const std::size_t index );
  
protected:
  std::unique_ptr< T > _data;
  std::size_t _arraySize;
};

#include "Scalar1D.it.h"

#endif /* Scalar1D_h */
