//
//  ioNetCDF.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/26/20.
//

#include "ioNetCDF.h"

namespace io
{

ioNetCDF::ioNetCDF( std::string fileName, size_t xDimSize, size_t yDimSize, char mode )
{
  _xDimSize = xDimSize;
  _yDimSize = yDimSize;
  
  if( mode == 'r')
  {
    _dataFile.open( fileName, netCDF::NcFile::read );
    
    _data_T = _dataFile.getVar( "data_T" );
    _data_C = _dataFile.getVar( "data_C" );
    _data_u = _dataFile.getVar( "data_u" );
    _data_v = _dataFile.getVar( "data_v" );
  }
  else if( mode == 'w' )
  {
    _dataFile.open( fileName, netCDF::NcFile::replace );
    
    netCDF::NcDim tDim = _dataFile.addDim( "T" );
    netCDF::NcDim xDim = _dataFile.addDim( "X", _xDimSize );
    netCDF::NcDim yDim = _dataFile.addDim( "Y", _yDimSize );
    
    std::vector< netCDF::NcDim > dims;
    dims.push_back( tDim );
    dims.push_back( xDim );
    dims.push_back( yDim );
    
    _data_T = _dataFile.addVar( "data_T", netCDF::ncDouble, dims );
    _data_C = _dataFile.addVar( "data_C", netCDF::ncDouble, dims );
    _data_u = _dataFile.addVar( "data_u", netCDF::ncDouble, dims );
    _data_v = _dataFile.addVar( "data_v", netCDF::ncDouble, dims );
  }
  
  _startp.push_back( 0 );
  _startp.push_back( 0 );
  _startp.push_back( 0 );
  _countp.push_back( 1 );
  _countp.push_back( _xDimSize );
  _countp.push_back( _yDimSize );
}

ioNetCDF::~ioNetCDF()
{
  _dataFile.close();
}

void ioNetCDF::read_T( size_t readIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataIn[ _xDimSize ][ _yDimSize ];
  _startp[ 0 ] = readIndex;
  
  _data_T.getVar( _startp, _countp, dataIn );
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataVector[ i ][ j ] = dataIn[ i ][ j ];
    }
  }
}

void ioNetCDF::read_C( size_t readIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataIn[ _xDimSize ][ _yDimSize ];
  _startp[ 0 ] = readIndex;
  
  _data_C.getVar( _startp, _countp, dataIn );
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataVector[ i ][ j ] = dataIn[ i ][ j ];
    }
  }
}

void ioNetCDF::read_u( size_t readIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataIn[ _xDimSize ][ _yDimSize ];
  _startp[ 0 ] = readIndex;
  
  _data_u.getVar( _startp, _countp, dataIn );
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataVector[ i ][ j ] = dataIn[ i ][ j ];
    }
  }
}

void ioNetCDF::read_v( size_t readIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataIn[ _xDimSize ][ _yDimSize ];
  _startp[ 0 ] = readIndex;
  
  _data_v.getVar( _startp, _countp, dataIn );
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataVector[ i ][ j ] = dataIn[ i ][ j ];
    }
  }
}

void ioNetCDF::write_T( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataArray[ _xDimSize ][ _yDimSize ];
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataArray[ i ][ j ] = dataVector[ i ][ j ];
    }
  }
  
  _startp[ 0 ] = writeIndex;
  _data_T.putVar( _startp, _countp, dataArray );
}

void ioNetCDF::write_C( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataArray[ _xDimSize ][ _yDimSize ];
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataArray[ i ][ j ] = dataVector[ i ][ j ];
    }
  }
  
  _startp[ 0 ] = writeIndex;
  _data_C.putVar( _startp, _countp, dataArray );
}

void ioNetCDF::write_u( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataArray[ _xDimSize ][ _yDimSize ];
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataArray[ i ][ j ] = dataVector[ i ][ j ];
    }
  }
  
  _startp[ 0 ] = writeIndex;
  _data_u.putVar( _startp, _countp, dataArray );
}

void ioNetCDF::write_v( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataArray[ _xDimSize ][ _yDimSize ];
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataArray[ i ][ j ] = dataVector[ i ][ j ];
    }
  }
  
  _startp[ 0 ] = writeIndex;
  _data_v.putVar( _startp, _countp, dataArray );
}

}
