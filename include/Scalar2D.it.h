#include <iostream>
#include <memory>

template< class T >
Scalar2D< T >::Scalar2D( std::size_t rowSize, std::size_t colSize )
{
  _rowSize = rowSize;
  _colSize = colSize;
  
  _data = std::unique_ptr< T >( new T[ _rowSize * _colSize ] );
}

template< class T >
Scalar2D< T >::~Scalar2D()
{
}

template< class T >
void Scalar2D< T >::set( std::size_t rowIndex, std::size_t colIndex, T value )
{
  ( _data.get() )[ _colSize * rowIndex + colIndex ] = value;
}

template< class T >
void Scalar2D< T >::set( T* arrayPtr )
{
  for( size_t i = 0; i < _rowSize; i++ )
  {
    for( size_t j = 0; j < _colSize; j++ )
    {
      ( _data.get() )[ _colSize * i + j ] = arrayPtr[ _colSize * i + j ];
    }
  }
}

template< class T >
T* Scalar2D< T >::get() const
{
  return _data.get();
}

template< class T >
std::size_t Scalar2D< T >::rows()
{
  return _rowSize;
}

template< class T >
std::size_t Scalar2D< T >::cols()
{
  return _colSize;
}

template< class T >
const T& Scalar2D< T >::operator()( const std::size_t rowIndex, const std::size_t colIndex ) const
{
  return ( _data.get() )[ _colSize * rowIndex + colIndex ];
}

template< class T >
T& Scalar2D< T >::operator()( const std::size_t rowIndex, const std::size_t colIndex )
{
  return ( _data.get() )[ _colSize * rowIndex + colIndex ];
}
