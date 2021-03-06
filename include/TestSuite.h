//
//  TestSuite.h
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#ifndef TestSuite_h
#define TestSuite_h

#include <iostream>
#include <chrono>
#include "FourierTransforms.h"
#include "ReadParams.h"
#include "Derivatives.h"
#include "ioNetCDF.h"

namespace diagnostics
{

class TestSuite
{
public:
  TestSuite();
  ~TestSuite();
  
  void testReadParams();
  void testFourierTransforms1();
  void testFourierTransforms2();
  void testFourierTransforms_2D();
  void testDeriv();
  void testDeriv2();
  void testDeriv2D();
  void testReadWriteIO();
  void testWriteIO();
  void testReadIO();
  
  void simpleAdvDiff( std::string paramFile );
  void simpleAdvDiffNL( std::string paramFile );
  
protected:
  const double pi = std::acos( -1.0 );
};

} // diagnostics

#endif /* TestSuite_h */
