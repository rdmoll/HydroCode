cmake_minimum_required(VERSION 2.8.9)
project(HydroCode)

#message( STATUS ${CMAKE_CURRENT_SOURCE_DIR} )

#cmake scripts
#set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#Include FFTW libraries
set( FFTW_DIR /usr/local )
set( NETCDF_DIR /opt/homebrew/Cellar/netcdf/4.8.1 )

set( FFTW_LIB_NAME "libfftw3.a" )
find_library( FFTW_LIB NAMES ${FFTW_LIB_NAME} PATHS ${FFTW_DIR}/lib )

set( NETCDF_LIB_NAME "libnetcdf.a" )
find_library( NETCDF_LIB NAMES ${NETCDF_LIB_NAME} PATHS ${NETCDF_DIR}/lib )

#Bring the headers into the project
list( APPEND includePath "${CMAKE_CURRENT_SOURCE_DIR}/include" 
                         "${FFTW_DIR}/include"
                         "${NETCDF_DIR}/include" )

#Use file(GLOB...) for wildcard additions
file(GLOB SOURCES "src/*.cpp")

add_executable( hydroTest ${SOURCES} )
target_link_libraries( hydroTest FFTW_LIB )
target_link_libraries( hydroTest NETCDF_LIB )
target_include_directories( hydroTest PUBLIC "$<BUILD_INTERFACE:${includePath}>" )