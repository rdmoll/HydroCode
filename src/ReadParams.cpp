//
//  ReadParams.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include "../include/ReadParams.h"

namespace io
{

ReadParams::ReadParams( const char* paramFile, const char* delimiter)
{
  delim = delimiter;
  fileHandle.open( paramFile );
  if( fileHandle.is_open() )
  {
    while( std::getline(fileHandle,line) )
    {
        ParseLine();
        StoreParams();
    }
    fileHandle.close();
  }
  else
  {
    std::cout << "ERROR: Can't open parameter file" << std::endl;
  }
}

ReadParams::~ReadParams()
{
  
}

void ReadParams::ParseLine( void )
{
  valName = line.substr( 0,line.find(delim) );
  valName.erase( std::remove(valName.begin(),valName.end(),' '), valName.end() );
  valStr = line.substr( line.find(delim) + delim.length() );
}

bool ReadParams::AlphaFound( void )
{
  const char alpha[] = { \
                      'a','b','c','d','f','g','h','i','j','k','l','m','n', \
                      'o','p','q','r','s','t','u','v','w','x','y','z', \
                      'A','B','C','D','F','G','H','I','J','K','L','M','N', \
                      'O','P','Q','R','S','T','U','V','W','X','Y','Z', \
                      };

  for( int i = 0; i < 50; i++)
  {
    if( valStr.find(alpha[i]) < valStr.length() )
    {
      return true;
    }
  }
  
  return false;
}

void ReadParams::StoreParams( void )
{
  if( AlphaFound() )
  {
    valStr.erase( std::remove(valStr.begin(),valStr.end(),' '), valStr.end() );
    strParam[valName] = valStr;
  }
  else
  {
    if( valStr.find('.') < valStr.length() )
    {
        dVal = atof( valStr.c_str() );
        dParam[valName] = dVal;
    }
    else
    {
        iVal = atoi( valStr.c_str() );
        iParam[valName] = iVal;
    }
  }
}

} // io
