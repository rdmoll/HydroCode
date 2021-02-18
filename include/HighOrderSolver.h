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
                          std::vector< std::vector< double > >& C0_phys,
                          std::vector< std::vector< double > >& u0_phys,
                          std::vector< std::vector< double > >& v0_phys );
  void calcNonLin( variables::VectorVar& u_phys,
                   std::vector< std::vector< std::complex< double > > >& f_spec,
                   std::vector< std::vector< std::complex< double > > >& NL_spec );
  void solve( std::vector< std::vector< std::complex< double > > >& f0_spec,
              std::vector< std::vector< std::complex< double > > >& nl_spec,
              double diffFac );
  void resetTimePointers( variables::ScalarVar& scalar );
  void resetTimePointers( variables::VectorVar& vec );
  void runSimulation();
  
protected:
  const double pi = std::acos( -1.0 );
  
  size_t nSteps;
  double deltaT;
  size_t Nx, Ny;
  double Lx, Ly;
  double nu, kappaT, kappaC, nuFac, kappaFacT, kappaFacC;
  std::string testWriterFile;
  
  int nOutX;
  int nOutY;
  size_t totN;
  
  hydroCode::FourierTransforms fft;
  mathOps::Derivatives ops;
  
  std::vector< std::vector< std::complex< double > > > dfdx_spec;
  std::vector< std::vector< std::complex< double > > > dfdy_spec;
  std::vector< std::vector< double > > dfdx_phys;
  std::vector< std::vector< double > > dfdy_phys;
  std::vector< std::vector< double > > NL_xPhys;
  std::vector< std::vector< double > > NL_yPhys;
  std::vector< std::vector< std::complex< double > > > NL_xSpec;
  std::vector< std::vector< std::complex< double > > > NL_ySpec;
  
  std::vector< double > wvNumX, wvNumY;
};

}

#endif /* HighOrderSolver_h */
