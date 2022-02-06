cmake_minimum_required(VERSION 2.8.12)
project(HydroCode)

#CMake scripts
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#Find FFTW and NetCDF packages
find_package( FFTW )
find_package( NETCDF )
find_package( MPI )

#Bring the headers into the project
list( APPEND includePath "${CMAKE_CURRENT_SOURCE_DIR}/include"
		        ${FFTW_INCLUDE_DIRS} 
                         ${NETCDF_INCLUDE_DIRS}
			${MPI_CXX_INCLUDE_DIRS} )

#Use file(GLOB...) for wildcard additions
file( GLOB SOURCES "src/*.cpp" )

#Use file(GLOB...) for wildcard additions
file( GLOB HEADERS "include/*.h" )

#Link libraries for executables
file( GLOB SANDBOX_MAIN "apps/sandbox/main.cpp" )
add_executable( sandbox ${SANDBOX_MAIN} ${SOURCES} ${HEADERS} )
target_link_libraries( sandbox ${FFTW_DOUBLE_LIB} )
target_link_libraries( sandbox ${FFTW_DOUBLE_THREADS_LIB} )
target_link_libraries( sandbox ${NETCDF_DOUBLE_LIB} )
target_link_libraries( sandbox ${NETCDF_CPP_LIB} )
target_link_libraries( sandbox ${NETCDF_FLOAT_LIB} )
target_link_libraries( sandbox ${MPI_CXX_LIBRARIES} )
target_include_directories( sandbox PUBLIC "$<BUILD_INTERFACE:${includePath}>" )

file( GLOB UNITTESTS_MAIN "apps/unitTests/main.cpp" )
add_executable( unitTests ${UNITTESTS_MAIN} ${SOURCES} ${HEADERS} )
target_link_libraries( unitTests ${FFTW_DOUBLE_LIB} )
target_link_libraries( unitTests ${FFTW_DOUBLE_THREADS_LIB} )
target_link_libraries( unitTests ${NETCDF_DOUBLE_LIB} )
target_link_libraries( unitTests ${NETCDF_CPP_LIB} )
target_link_libraries( unitTests ${NETCDF_FLOAT_LIB} )
target_link_libraries( unitTests ${MPI_CXX_LIBRARIES} )
target_include_directories( unitTests PUBLIC "$<BUILD_INTERFACE:${includePath}>" )

file( GLOB TESTSOLVER1_MAIN "apps/testSolver1/main.cpp" )
add_executable( testSolver1 ${TESTSOLVER1_MAIN} ${SOURCES} ${HEADERS} )
target_link_libraries( testSolver1 ${FFTW_DOUBLE_LIB} )
target_link_libraries( testSolver1 ${FFTW_DOUBLE_THREADS_LIB} )
target_link_libraries( testSolver1 ${NETCDF_DOUBLE_LIB} )
target_link_libraries( testSolver1 ${NETCDF_CPP_LIB} )
target_link_libraries( testSolver1 ${NETCDF_FLOAT_LIB} )
target_link_libraries( testSolver1 ${MPI_CXX_LIBRARIES} )
target_include_directories( testSolver1 PUBLIC "$<BUILD_INTERFACE:${includePath}>" )