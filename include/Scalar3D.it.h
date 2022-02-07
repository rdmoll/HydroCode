#include <iostream>
#include <memory>

template< class T >
Scalar3D< T >::Scalar3D( std::size_t rowSize, std::size_t colSize, std::size_t sliceSize )
{
  _rowSize = rowSize;
  _colSize = colSize;
  _sliceSize = sliceSize;
  
  _data = std::unique_ptr< T >( new T[ _rowSize * _colSize * _sliceSize ] );
}

template< class T >
Scalar3D< T >::~Scalar3D()
{
}

template< class T >
void Scalar3D< T >::set( std::size_t rowIndex, std::size_t colIndex, std::size_t sliceIndex, T value )
{
  ( _data.get() )[ _sliceSize*_colSize*rowIndex + _sliceSize*colIndex + sliceIndex ] = value;
}

template< class T >
void Scalar3D< T >::set( T* arrayPtr )
{
  std::size_t totSize = _rowSize * _colSize * _sliceSize;
  
  for( std::size_t i = 0; i < totSize; ++i )
  {
    ( _data.get() )[ i ] = arrayPtr[ i ];
  }
}

template< class T >
void Scalar3D< T >::setSize( std::size_t rowSize, std::size_t colSize, std::size_t sliceSize )
{
  _rowSize = rowSize;
  _colSize = colSize;
  _sliceSize = sliceSize;
  
  _data = std::unique_ptr< T >( new T[ _rowSize * _colSize * _sliceSize ] );
}

template< class T >
T* Scalar3D< T >::get() const
{
  return _data.get();
}

template< class T >
std::size_t Scalar3D< T >::rows()
{
  return _rowSize;
}

template< class T >
std::size_t Scalar3D< T >::cols()
{
  return _colSize;
}

template< class T >
std::size_t Scalar3D< T >::slices()
{
  return _sliceSize;
}

template< class T >
const T& Scalar3D< T >::operator()( const std::size_t rowIndex, const std::size_t colIndex, const std::size_t sliceIndex ) const
{
  return ( _data.get() )[ _sliceSize*_colSize*rowIndex + _sliceSize*colIndex + sliceIndex ];
}

template< class T >
T& Scalar3D< T >::operator()( const std::size_t rowIndex, const std::size_t colIndex, const std::size_t sliceIndex )
{
  return ( _data.get() )[ _sliceSize*_colSize*rowIndex + _sliceSize*colIndex + sliceIndex ];
}

template< class T >
Scalar3D< T >& Scalar3D< T >::operator=( Scalar3D< T >& arr )
{
  if( arr.rows() == _rowSize && arr.cols() == _colSize && arr.slices() == _sliceSize )
  {
    for( size_t i = 0; i < _rowSize; ++i )
    {
      for( size_t j = 0; j < _colSize; ++j )
      {
        for( size_t k = 0; k < _sliceSize; ++k )
        {
          set( i, j, k, arr( i, j, k ) );
        }
      }
    }
  }
  
  return *this;
}

template< class T >
Scalar3D< T >& operator*( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.rows(), arr1.cols(), arr1.slices() );
  
  for( size_t i = 0; i < arr1.rows(); ++i )
  {
    for( size_t j = 0; j < arr1.cols(); ++j )
    {
      for( size_t k = 0; k < arr1.slices(); ++k )
      {
        arrOut( i, j, k ) = arr1( i, j, k ) * arr2( i, j, k );
      }
    }
  }
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator+( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.rows(), arr1.cols(), arr1.slices() );
  
  for( size_t i = 0; i < arr1.rows(); ++i )
  {
    for( size_t j = 0; j < arr1.cols(); ++j )
    {
      for( size_t k = 0; k < arr1.slices(); ++k )
      {
        arrOut( i, j, k ) = arr1( i, j, k ) + arr2( i, j, k );
      }
    }
  }
  
  return arrOut;
}
