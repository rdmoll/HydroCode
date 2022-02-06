#include <iostream>
#include <complex>

#include "Scalar1D.h"
#include "Scalar2D.h"
#include "Scalar3D.h"
#include "Transforms.h"
#include "MathOps.h"
#include "TestSolver.h"
#include "mpi.h"

int main( int argc, const char * argv[] )
{
  int myRank, p, n;
  
  MPI_Init( NULL, NULL );
  MPI_Comm_size( MPI_COMM_WORLD, &p );
  MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
  
  std::cout << "Total num procs = " << p << " , proc rank = " << myRank << std::endl;
  
  /*
  std::size_t N = 16777216; //16;
  std::size_t nOut = std::floor( N / 2 + 1 );
  
  Scalar1D< double > realTruth( N );
  Scalar1D< double > in( N );
  Scalar1D< std::complex< double > > out( nOut );
  Scalar1D< std::complex< double > > df_spec( nOut );
  Scalar1D< double > df_phys( N );
  
  for( std::size_t i = 0; i < N; i++ )
  {
    realTruth( i ) = std::cos( i * ( 2 * M_PI / N ) );
    in( i ) = std::cos( i * ( 2 * M_PI / N ) );
    
    //std::cout << in( i ) << std::endl;
  }
   */
  
  if( myRank == 0 )
  {
    std::size_t Nx = 4096;
    std::size_t Ny = 4096;
    std::size_t nOutX = std::floor( Nx / 2 + 1 );
    
    Scalar2D< double > realTruth( Ny, Nx );
    Scalar2D< double > in( Ny, Nx );
    Scalar2D< std::complex< double > > out( Ny, nOutX );
    
    for( int i = 0; i < Ny; ++i )
    {
      for( int j = 0; j < Nx; ++j )
      {
        realTruth( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::sin( j * ( 2.0 * M_PI / Nx ) );
        in( i, j ) = std::cos( i * ( 2.0 * M_PI / Ny ) ) * std::cos( j * ( 2.0 * M_PI / Nx ) );
      }
    }
    
    /*
    fft::fft_r2c_1d( in, out );
    
    std::cout << std::endl;
    for( std::size_t i = 0; i < nOut; i++ )
    {
      std::cout << out( i ) << std::endl;
    }
    
    mathOps::calcDeriv( out, df_spec );
    
    fft::fft_c2r_1d( df_spec, df_phys );
    
    std::cout << std::endl;
    for( std::size_t i = 0; i < N; i++ )
    {
      std::cout << df_phys( i ) << std::endl;
    }
    
    Scalar1D< double > test1( 3 );
    Scalar1D< double > test2( 3 );
    Scalar1D< double > test3( 3 );
    
    test1( 0 ) = 1.0;
    test1( 1 ) = 2.0;
    test1( 2 ) = 3.0;
    
    test2( 0 ) = 3.0;
    test2( 1 ) = 4.0;
    test2( 2 ) = 5.0;
    
    test3 = test1 * test2;
    
    std::cout << std::endl;
    
    for( std::size_t i = 0; i < 3; ++i )
    {
      std::cout << test3( i ) << std::endl;
    }
    
    std::cout << std::endl;
    
    std::cout << test3.size() << std::endl;
    
    solvers::TestSolver solver1( "/Users/ryanmoll/Documents/dev/projects/HydroCode/tests/advDiffParams.txt" );
    solver1.runSimulation();
    
    std::cout << std::endl;
    */
    
    {
    // Start timer
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    fft::fft_r2c_2d_thread( in, out );
    
    // End timer
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
    
    // Print completion message and timing information
    std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
              << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n" << std::endl;
    }
    
    {
    // Start timer
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    fft::fft_r2c_2d( in, out );
    
    // End timer
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
    
    // Print completion message and timing information
    std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
              << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed: "
              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n" << std::endl;
    }
  }

  MPI_Finalize();
  
  return 0;
}

