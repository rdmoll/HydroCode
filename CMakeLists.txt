cmake_minimum_required(VERSION 2.8.9)
project(HydroCode)

#CMake scripts
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#Find FFTW and NetCDF packages
find_package( FFTW )
find_package( NETCDF )

#Bring the headers into the project
list( APPEND includePath "${CMAKE_CURRENT_SOURCE_DIR}/include"
		        ${FFTW_INCLUDE_DIRS} 
                         ${NETCDF_INCLUDE_DIRS} )

#Use file(GLOB...) for wildcard additions
file( GLOB SOURCES "apps/sandbox/main.cpp" "src/*.cpp" )

#Use file(GLOB...) for wildcard additions
file( GLOB HEADERS "include/*.h" )

#Link libraries
add_executable( hydroTest ${SOURCES} ${HEADERS} )
target_link_libraries( hydroTest ${FFTW_DOUBLE_LIB} )
target_link_libraries( hydroTest ${NETCDF_DOUBLE_LIB} )
target_link_libraries( hydroTest ${NETCDF_CPP_LIB} )
target_link_libraries( hydroTest ${NETCDF_FLOAT_LIB} )
target_include_directories( hydroTest PUBLIC "$<BUILD_INTERFACE:${includePath}>" )