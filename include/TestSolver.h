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
#include "Derivatives.h"

namespace solvers
{

class TestSolver
{
public:
  TestSolver( std::string paramFile );
  ~TestSolver();
  
  void setInitConditions( std::vector< std::vector< double > >& T0_phys,
                          std::vector< std::vector< double > >& u0_phys );
  void calcDerivX( std::vector< std::vector< std::complex< double > > >& f_spec,
                   std::vector< std::vector< std::complex< double > > >& df_spec );
  void calcNonLin( std::vector< std::vector< double > >& f1_phys,
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

#endif /* TestSolver_h */
