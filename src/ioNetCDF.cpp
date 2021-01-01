//
//  ioNetCDF.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/26/20.
//

#include "ioNetCDF.h"

namespace io
{

ioNetCDF::ioNetCDF( std::string fileName, std::string variable, size_t xDimSize, size_t yDimSize )
{
  _dataFile.open( fileName, netCDF::NcFile::replace );
  
  _xDimSize = xDimSize;
  _yDimSize = yDimSize;
  
  netCDF::NcDim tDim = _dataFile.addDim( "T" );
  netCDF::NcDim xDim = _dataFile.addDim( "X", _xDimSize );
  netCDF::NcDim yDim = _dataFile.addDim( "Y", _yDimSize );
  
  std::vector< netCDF::NcDim > dims;
  dims.push_back( tDim );
  dims.push_back( xDim );
  dims.push_back( yDim );
  
  _data = _dataFile.addVar( variable, netCDF::ncDouble, dims );
  
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

void ioNetCDF::read( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataIn[ _xDimSize ][ _yDimSize ];
  _startp[ 0 ] = writeIndex;
  
  _data.getVar( _startp, _countp, dataIn );
  
  for( size_t i = 0; i < _xDimSize; i++ )
  {
    for( size_t j = 0; j < _yDimSize; j++ )
    {
      dataVector[ i ][ j ] = dataIn[ i ][ j ];
    }
  }
}

void ioNetCDF::write( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
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
  _data.putVar( _startp, _countp, dataArray );
}

}
