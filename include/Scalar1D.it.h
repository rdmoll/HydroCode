#include <iostream>
#include <memory>

template< class T >
Scalar1D< T >::Scalar1D( std::size_t arraySize )
{
  _arraySize = arraySize;
  _data = std::unique_ptr< T >( new T[ arraySize ] );
}

template< class T >
Scalar1D< T >::~Scalar1D()
{
}

template< class T >
void Scalar1D< T >::set( std::size_t index, T value )
{
  ( _data.get() )[ index ] = value;
}

template< class T >
void Scalar1D< T >::set( T* arrayPtr )
{
  for( std::size_t i = 0; i < _arraySize; ++i )
  {
    ( _data.get() )[ i ] = arrayPtr[ i ];
  }
}

template< class T >
T* Scalar1D< T >::get() const
{
  return _data.get();
}

template< class T >
std::size_t Scalar1D< T >::size()
{
  return _arraySize;
}

template< class T >
const T& Scalar1D< T >::operator()( const std::size_t index ) const
{
  return ( _data.get() )[ index ];
}

template< class T >
T& Scalar1D< T >::operator()( const std::size_t index )
{
  return ( _data.get() )[ index ];
}
