//
//  TestSolver.h
//  HydroCode
//
//  Created by Ryan Moll on 2/1/22.
//

#ifndef TestSolver_h
#define TestSolver_h

#include <iostream>
#include <string>
#include <math.h>
#include <complex>
#include "ReadParams.h"
#include "Scalar2D.h"
#include "Transforms.h"
#include "ioNetCDF.h"

namespace solvers
{

class TestSolver
{
public:
  TestSolver( std::string paramFile );
  ~TestSolver();
  
  void setInitConditions( Scalar2D< double >& T_phys,
                          Scalar2D< double >& u_phys,
                          Scalar2D< double >& v_phys );
  
  void runSimulation();
  
protected:
  size_t nSteps;
  double deltaT;
  size_t Nx, Ny;
  double Lx, Ly;
  double nu, kappa;
  std::string testWriterFile;
  
  int nOutX;
  int nOutY;
  size_t totN;
  
  const double pi = std::acos( -1.0 );
};

}

#endif /* TestSolver_h */
