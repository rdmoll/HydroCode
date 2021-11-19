cmake_minimum_required(VERSION 2.8.9)
project(HydroCode)

#message( STATUS ${CMAKE_CURRENT_SOURCE_DIR} )

#cmake scripts
#set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#Include FFTW libraries
set( FFTW_DIR /usr/local )
set( NETCDF_DIR /opt/homebrew )

set( FFTW_LIB_NAME "libfftw3.a" )
find_library( FFTW_LIB NAMES ${FFTW_LIB_NAME} PATHS ${FFTW_DIR}/lib )

set( NETCDF_LIB_NAME1 "libnetcdf.a" )
set( NETCDF_LIB_NAME2 "libnetcdf-cxx4.a" )
set( NETCDF_LIB_NAME3 "libnetcdff.a" )
find_library( NETCDF_LIB1 NAMES ${NETCDF_LIB_NAME1} PATHS ${NETCDF_DIR}/lib )
find_library( NETCDF_LIB2 NAMES ${NETCDF_LIB_NAME2} PATHS ${NETCDF_DIR}/lib )
find_library( NETCDF_LIB3 NAMES ${NETCDF_LIB_NAME3} PATHS ${NETCDF_DIR}/lib )

#Bring the headers into the project
list( APPEND includePath "${CMAKE_CURRENT_SOURCE_DIR}/include" 
                         "${FFTW_DIR}/include"
                         "${NETCDF_DIR}/include" )

#Use file(GLOB...) for wildcard additions
file( GLOB SOURCES "apps/*.cpp" "src/*.cpp" )

add_executable( hydroTest ${SOURCES} )
target_link_libraries( hydroTest ${FFTW_LIB} )
target_link_libraries( hydroTest ${NETCDF_LIB2} )
target_link_libraries( hydroTest ${NETCDF_LIB1} )
target_link_libraries( hydroTest ${NETCDF_LIB3} )
target_include_directories( hydroTest PUBLIC "$<BUILD_INTERFACE:${includePath}>" )