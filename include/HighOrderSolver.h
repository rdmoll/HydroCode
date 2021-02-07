//
//  HighOrderSolver.h
//  HydroCode
//
//  Created by Ryan Moll on 2/6/21.
//

#ifndef HighOrderSolver_h
#define HighOrderSolver_h

#include <iostream>
#include <chrono>
#include "FourierTransforms.h"
#include "ReadParams.h"
#include "ioNetCDF.h"
#include "Derivatives.h"
#include "variables.h"

namespace solvers
{

class HighOrderSolver
{
public:
  HighOrderSolver( std::string paramFile );
  ~HighOrderSolver();
  
  void setInitConditions( std::vector< std::vector< double > >& T0_phys,
                         std::vector< std::vector< double > >& u0_phys,
                         std::vector< std::vector< double > >& v0_phys );
  void calcNonLin( variables::VectorVar& u,
                   std::vector< std::vector< std::complex< double > > >& f_spec,
                   std::vector< std::vector< std::complex< double > > >& NL_spec );
  void solve( std::vector< std::vector< std::complex< double > > >& f0_spec,
              std::vector< std::vector< std::complex< double > > >& nl_spec,
              double diffusivity );
  void resetTimePointers( variables::ScalarVar& scalar );
  void resetTimePointers( variables::VectorVar& vec );
  void runSimulation();
  
protected:
  const double pi = std::acos( -1.0 );
  
  size_t nSteps;
  double deltaT;
  size_t Nx, Ny;
  double Lx, Ly;
  double nu, kappa;
  std::string testWriterFile;
  
  int nOutX;
  int nOutY;
  size_t totN;
  
  hydroCode::FourierTransforms fft;
  mathOps::Derivatives ops;
};

}

#endif /* HighOrderSolver_h */
