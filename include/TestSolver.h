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
  
  void calcNL( Scalar2D< double >& f1_phys,
               Scalar2D< double >& f2_phys,
               Scalar2D< std::complex< double > >& f3_spec,
               Scalar2D< std::complex< double > >& nl_spec );
  
  void solve( Scalar2D< std::complex< double > >& f_spec,
              Scalar2D< std::complex< double > >& nl_spec,
              const double timeStep,
              const double diffusivity );
  
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
  double fourPiSq;
  Scalar2D< std::complex< double > > totWvNum;
};

}

#endif /* TestSolver_h */
