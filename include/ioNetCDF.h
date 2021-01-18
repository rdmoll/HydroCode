//
//  ioNetCDF.hpp
//  HydroCode
//
//  Created by Ryan Moll on 12/26/20.
//

#ifndef ioNetCDF_h
#define ioNetCDF_h

#include <iostream>
#include <string>
#include <netcdf>

namespace io
{

class ioNetCDF
{
public:
  ioNetCDF( std::string fileName, std::string variable, size_t xDimSize, size_t yDimSize, char mode );
  ~ioNetCDF();
  
  void read( size_t readIndex, std::vector< std::vector< double > >& dataVector );
  void write( size_t writeIndex, std::vector< std::vector< double > >& dataVector );
  
  netCDF::NcFile _dataFile;
  netCDF::NcVar _data;
  size_t _xDimSize, _yDimSize;
  std::vector< size_t > _startp, _countp;
};

}

#endif /* ioNetCDF_h */
