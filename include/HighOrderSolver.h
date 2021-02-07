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
  void calcNonLinDx( std::vector< std::vector< double > >& f1_phys,
                     std::vector< std::vector< std::complex< double > > >& f2_spec,
                     std::vector< std::vector< std::complex< double > > >& nl_spec );
  void calcNonLinDy( std::vector< std::vector< double > >& f1_phys,
                     std::vector< std::vector< std::complex< double > > >& f2_spec,
                     std::vector< std::vector< std::complex< double > > >& nl_spec );
  void solve( std::vector< std::vector< std::complex< double > > >& f0_spec,
              std::vector< std::vector< std::complex< double > > >& nl_spec,
              double diffusivity );
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
