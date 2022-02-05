template< class T >
Scalar3D< T >::Scalar3D( std::size_t arraySize )
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
const std::size_t Scalar3D< T >::slices()
{
  return _sliceSize;
}

template< class T >
const T& Scalar3D< T >::operator()( const std::size_t index ) const
{
  return ( _data.get() )[ index ];
}

template< class T >
T& Scalar3D< T >::operator()( const std::size_t index )
{
  return ( _data.get() )[ index ];
}

template< class T >
Scalar3D< T >& Scalar3D< T >::operator=( Scalar3D< T >& arr )
{
  if( arr.size() == _arraySize )
  {
    for( std::size_t i = 0; i < _arraySize; ++i )
    {
      set( i, arr( i ) );
    }
  }
  return *this;
}

template< class T >
Scalar3D< T >& operator*( Scalar3D< T >& arr1, Scalar3D< T >& arr2 )
{
  static Scalar3D< T > arrOut( arr1.size() );
  
  for( std::size_t i = 0; i < arr1.size(); ++i )
  {
    arrOut( i ) = arr1( i ) * arr2( i );
  }
  
  return arrOut;
}
