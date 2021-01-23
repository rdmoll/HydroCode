//
//  Solver.hpp
//  HydroCode
//
//  Created by Ryan Moll on 1/22/21.
//

#ifndef TestSolver_h
#define TestSolver_h

#include <iostream>
#include <chrono>
#include "FourierTransforms.h"
#include "ReadParams.h"
#include "ioNetCDF.h"

namespace solvers
{

class TestSolver
{
public:
  TestSolver( std::string paramFile );
  ~TestSolver();
  
  void setInitConditions( std::vector< std::vector< double > >& T0_phys,
                          std::vector< std::vector< double > >& u0_phys );
  void calcDeriv( std::vector< std::vector< std::complex< double > > >& f_spec,
                  std::vector< std::vector< std::complex< double > > >& df_spec );
  void solveU( std::vector< std::vector< std::complex< double > > >& u_spec, double diffFac );
  void solveT( std::vector< std::vector< std::complex< double > > >& T_spec, double advFac, double diffFac );
  void runSimulation();
  
protected:
  const double pi = std::acos( -1.0 );
  
  size_t nSteps;
  double deltaT;
  size_t Nx, Ny;
  double Lx, Ly;
  double nu;
  std::string testWriterFile;
  
  int nOutX;
  int nOutY;
  size_t totN;
  
  hydroCode::FourierTransforms fft;
};

}

#endif /* TestSolver_h */
