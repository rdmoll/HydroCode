//
//  TestSuite.h
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#ifndef TestSuite_h
#define TestSuite_h

#include <iostream>
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
  void testReadWriteIO();
  void testWriteIO();
  void testReadIO();
  
  void simpleAdvDiff( std::string paramFile );
  void simpleAdvDiffNL( std::string paramFile );
};

} // diagnostics

#endif /* TestSuite_h */
