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
  ioNetCDF( std::string fileName, std::string variable, size_t xDimSize, size_t yDimSize );
  ~ioNetCDF();
  
  void read();
  void write( size_t writeIndex, std::vector< std::vector< double > >& dataVector );
  
  netCDF::NcFile dataFile;
  netCDF::NcVar data;
  std::vector< size_t > startp, countp;
};

}

#endif /* ioNetCDF_h */
