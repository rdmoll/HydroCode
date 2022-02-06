//
//  main.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include <iostream>
#include <complex>
#include <netcdf>
#include "TestSuite.h"
#include "Transforms.h"
#include "ioNetCDF.h"

int main( int argc, const char * argv[] )
{
  // Do unit testing
  std::cout << "RUNNING UNIT TESTS" << std::endl;
  diagnostics::TestSuite unitTests;
  unitTests.testReadParams();
  unitTests.testFourierTransform1D();
  unitTests.testFourierTransform2D();
  unitTests.testDeriv1D();
  unitTests.testDeriv2D();
  unitTests.testReadWriteIO();
  
  std::cout << std::endl;
  
  unitTests.simpleAdvDiff( "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/advDiffParams.txt" );
  
  std::cout << std::endl;
  
  unitTests.simpleAdvDiffNL( "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/advDiffNLParams.txt" );
  
  return 0;
}

