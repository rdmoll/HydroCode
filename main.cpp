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
  std::cout << std::endl;
  
  //size_t nSteps = 40;
  //double deltaT = 0.005;
  //size_t Nx = 64;
  //size_t Ny = 64;
  //double c = 0.1;
  //double nu = 0.0;
  //unitTests.simpleAdvDiff( nSteps, deltaT, Nx, Ny, c, nu );
  //unitTests.simpleAdvDiffNL( nSteps, deltaT, Nx, Ny, c, nu );
  
  return 0;
}

