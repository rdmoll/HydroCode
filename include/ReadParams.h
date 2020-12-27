//
//  ReadParams.h
//  HydroCode
//
//  Created by Ryan Moll on 12/13/20.
//

#ifndef ReadParams_h
#define ReadParams_h

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

namespace io
{

class ReadParams
{
public:
  ReadParams( const char*,const char* );
  ~ReadParams();
  
  std::map<std::string,double> iParam;
  std::map<std::string,double> dParam;
  std::map<std::string,std::string> strParam;
  
protected:
  int iVal;
  double dVal;
  std::string delim;
  std::string line;
  std::string valStr;
  std::string valName;
  std::ifstream fileHandle;

  void ParseLine( void );
  bool AlphaFound( void );
  void StoreParams( void );
};

} // io

#endif  /* ReadParams_h */
