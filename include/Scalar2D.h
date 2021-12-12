#ifndef Scalar2D_h
#define Scalar2D_h

template< class T >
class Scalar2D
{
public:
  Scalar2D( std::size_t rowSize, std::size_t colSize );
  ~Scalar2D( );
  
  // Container methods
  void set( std::size_t rowIndex, std::size_t colIndex, T value );
  void set( T* arrayPtr );
  T* get() const;
  std::size_t rows();
  std::size_t cols();
  
  // Operators
  const T& operator()( const std::size_t rowIndex, const std::size_t colIndex ) const;
  T& operator()( const std::size_t rowIndex, const std::size_t colIndex );
  
protected:
  std::unique_ptr< T > _data;
  std::size_t _rowSize;
  std::size_t _colSize;
};

#include "Scalar2D.it.h"

#endif /* Scalar2D_h */
