#include <iostream>
#include <memory>

template< class T >
Scalar3D< T >::Scalar3D( std::size_t sliceSize, std::size_t rowSize, std::size_t colSize )
{
  _sliceSize = sliceSize;
  _rowSize = rowSize;
  _colSize = colSize;
  
  _data = std::unique_ptr< T >( new T[ _sliceSize * _rowSize * _colSize ] );
}

template< class T >
Scalar3D< T >::~Scalar3D()
{
}

template< class T >
void Scalar3D< T >::set( std::size_t sliceIndex, std::size_t rowIndex, std::size_t colIndex, T value )
{
  ( _data.get() )[ _rowSize*_colSize*sliceIndex + _colSize*rowIndex + colIndex ] = value;
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
void Scalar3D< T >::setSize( std::size_t sliceSize, std::size_t rowSize, std::size_t colSize )
{
  _rowSize = rowSize;
  _colSize = colSize;
  _sliceSize = sliceSize;
  
  _data = std::unique_ptr< T >( new T[ _sliceSize * _rowSize * _colSize ] );
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
const T& Scalar3D< T >::operator()( const std::size_t sliceIndex, const std::size_t rowIndex, const std::size_t colIndex ) const
{
  return ( _data.get() )[ _rowSize*_colSize*sliceIndex + _colSize*rowIndex + colIndex ];
}

template< class T >
T& Scalar3D< T >::operator()( const std::size_t sliceIndex, const std::size_t rowIndex, const std::size_t colIndex )
{
  return ( _data.get() )[ _rowSize*_colSize*sliceIndex + _colSize*rowIndex + colIndex ];
}

template< class T >
Scalar3D< T >& Scalar3D< T >::operator=( Scalar3D< T >& arr )
{
  if( arr.rows() == _rowSize && arr.cols() == _colSize && arr.slices() == _sliceSize )
  {
    for( size_t i = 0; i < _sliceSize; ++i )
    {
      for( size_t j = 0; j < _rowSize; ++j )
      {
        for( size_t k = 0; k < _colSize; ++k )
        {
          set( i, j, k, arr( i, j, k ) );
        }
      }
    }
  }
  
  return *this;
}

// Multiplication

template< class T >
Scalar3D< T >& operator*( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.slices(), arr1.rows(), arr1.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator*( T& val, Scalar3D< T >& arr )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator*( Scalar3D< T >& arr, T& val )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

// Division

template< class T >
Scalar3D< T >& operator/( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.slices(), arr1.rows(), arr1.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator/( T& val, Scalar3D< T >& arr )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator/( Scalar3D< T >& arr, T& val )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

// Addition

template< class T >
Scalar3D< T >& operator+( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.slices(), arr1.rows(), arr1.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator+( T& val, Scalar3D< T >& arr )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator+( Scalar3D< T >& arr, T& val )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

// Subtraction

template< class T >
Scalar3D< T >& operator-( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.slices(), arr1.rows(), arr1.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator-( T& val, Scalar3D< T >& arr )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}

template< class T >
Scalar3D< T >& operator-( Scalar3D< T >& arr, T& val )
{
  static Scalar3D< T > arrOut( arr.slices(), arr.rows(), arr.cols() );
  
  return arrOut;
}
