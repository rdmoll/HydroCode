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
#include "Scalar1D.h"
#include "Scalar2D.h"

namespace io
{

class ioNetCDF
{
public:
  ioNetCDF( std::string fileName, size_t xDimSize, size_t yDimSize, char mode );
  ~ioNetCDF();
  
  void read_T( size_t readIndex, Scalar2D< double >& dataArray );
  void read_C( size_t readIndex, Scalar2D< double >& dataArray );
  void read_u( size_t readIndex, Scalar2D< double >& dataArray );
  void read_v( size_t readIndex, Scalar2D< double >& dataArray );
  void write_T( size_t writeIndex, Scalar2D< double >& dataArray );
  void write_C( size_t writeIndex, Scalar2D< double >& dataArray );
  void write_u( size_t writeIndex, Scalar2D< double >& dataArray );
  void write_v( size_t writeIndex, Scalar2D< double >& dataArray );
  
  netCDF::NcFile _dataFile;
  netCDF::NcVar _data_T, _data_C, _data_u, _data_v;
  size_t _xDimSize, _yDimSize;
  std::vector< size_t > _startp, _countp;
};

}

#endif /* ioNetCDF_h */
