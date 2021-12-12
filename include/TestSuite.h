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
#include "Scalar1D.h"
#include "Scalar2D.h"
#include "Transforms.h"
#include "ReadParams.h"
#include "MathOps.h"
#include "ioNetCDF.h"

namespace diagnostics
{

class TestSuite
{
public:
  TestSuite();
  ~TestSuite();
  
  void testReadParams();
  void testFourierTransform1D();
  void testFourierTransform2D();
  void testDeriv1D();
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
