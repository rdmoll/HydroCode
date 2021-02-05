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
#include "FourierTransforms.h"
#include "ioNetCDF.h"
#include "TestSolver.h"
#include "variables.h"

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
  unitTests.testDeriv2D();
  
  std::cout << std::endl;
  
  unitTests.simpleAdvDiff( "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/advDiffParams.txt" );
  
  std::cout << std::endl;
  
  unitTests.simpleAdvDiffNL( "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/advDiffNLParams.txt" );
  
  std::cout << std::endl;
  
  std::string testSolverParamFile = "/Users/rmoll/Documents/dev/projects/HydroCode/Testing/testSolverParams.txt";
  solvers::TestSolver simpleAdvDiff( testSolverParamFile );
  simpleAdvDiff.runSimulation();
  
  return 0;
}

