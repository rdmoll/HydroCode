//
//  main.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include <iostream>
#include <complex>
#include <netcdf>
#include "../../include/TestSuite.h"
#include "../../include/Transforms.h"
#include "../../include/ioNetCDF.h"
//#include "../../include/TestSolver.h"
//#include "../../include/HighOrderSolver.h"
//#include "../../include/variables.h"

int main( int argc, const char * argv[] )
{
  // do unit testing
  std::cout << "RUNNING UNIT TESTS" << std::endl;
  diagnostics::TestSuite unitTests;
  unitTests.testReadParams();
  unitTests.testFourierTransform1D();
  unitTests.testFourierTransform2D();
  unitTests.testDeriv1D();
  unitTests.testDeriv2D();
  unitTests.testReadWriteIO();
  
  //std::cout << std::endl;
  
  //unitTests.simpleAdvDiff( "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/advDiffParams.txt" );
  
  //std::cout << std::endl;
  
  //unitTests.simpleAdvDiffNL( "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/advDiffNLParams.txt" );
  
  //std::cout << std::endl;
  
  //std::string testSolverParamFile = "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/testSolverParams.txt";
  //solvers::TestSolver simpleAdvDiff( testSolverParamFile );
  //solvers::HighOrderSolver simpleAdvDiff( testSolverParamFile );
  //simpleAdvDiff.runSimulation();
  
  return 0;
}

