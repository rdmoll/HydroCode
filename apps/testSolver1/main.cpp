//
//  main.cpp
//  HydroCode
//
//  Created by Ryan Moll on 12/21/20.
//

#include <iostream>
#include <complex>
#include <netcdf>
#include "Transforms.h"
#include "ioNetCDF.h"
#include "TestSolver.h"

int main( int argc, const char * argv[] )
{
  std::string testSolverParamFile = "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/testSolverParams.txt";
  solvers::TestSolver testSolver1( testSolverParamFile );
  testSolver1.runSimulation();
  
  return 0;
}

