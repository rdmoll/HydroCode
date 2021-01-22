//
//  main.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include <iostream>
#include <netcdf>
#include "TestSuite.h"
#include "FourierTransforms.h"
#include "ioNetCDF.h"

int main( int argc, const char * argv[] )
{
  // do unit testing
  std::cout << "RUNNING UNIT TESTS" << std::endl;
  diagnostics::TestSuite unitTests;
  unitTests.testReadParams();
  unitTests.testFourierTransforms1();
  unitTests.testFourierTransforms2();
  unitTests.testFourierTransforms_2D();
  unitTests.testDeriv();
  unitTests.testDeriv2();
  unitTests.testReadWriteIO();
  unitTests.simpleAdvDiff( "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/advDiffParams.txt" );
  unitTests.simpleAdvDiffNL( "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/advDiffNLParams.txt" );
  std::cout << std::endl;
  
  return 0;
}

