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
  dataFile.open( fileName, netCDF::NcFile::replace );
  
  netCDF::NcDim tDim = dataFile.addDim( "T" );
  netCDF::NcDim xDim = dataFile.addDim( "X", xDimSize );
  netCDF::NcDim yDim = dataFile.addDim( "Y", yDimSize );
  
  std::vector< netCDF::NcDim > dims;
  dims.push_back( tDim );
  dims.push_back( xDim );
  dims.push_back( yDim );
  
  data = dataFile.addVar( variable, netCDF::ncDouble, dims );
  
  startp.push_back( 0 );
  startp.push_back( 0 );
  startp.push_back( 0 );
  countp.push_back( 1 );
  countp.push_back( 10 );
  countp.push_back( 10 );
}

ioNetCDF::~ioNetCDF()
{
  dataFile.close();
}

void ioNetCDF::read()
{
  
}

void ioNetCDF::write( size_t writeIndex, std::vector< std::vector< double > >& dataVector )
{
  double dataArray[ dataVector.size() ][ dataVector[ 0 ].size() ];
  startp[ 0 ] = writeIndex;
  data.putVar( startp, countp, dataArray );
}

}
