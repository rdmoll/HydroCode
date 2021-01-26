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
  ioNetCDF( std::string fileName, size_t xDimSize, size_t yDimSize, char mode );
  ~ioNetCDF();
  
  void read_T( size_t readIndex, std::vector< std::vector< double > >& dataVector );
  void read_u( size_t readIndex, std::vector< std::vector< double > >& dataVector );
  void write_T( size_t writeIndex, std::vector< std::vector< double > >& dataVector );
  void write_u( size_t writeIndex, std::vector< std::vector< double > >& dataVector );
  
  netCDF::NcFile _dataFile;
  netCDF::NcVar _data_T, _data_u;
  size_t _xDimSize, _yDimSize;
  std::vector< size_t > _startp, _countp;
};

}

#endif /* ioNetCDF_h */
