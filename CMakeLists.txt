cmake_minimum_required(VERSION 2.8.9)
project(HydroCode)

#cmake scripts
set(CMAKE_MODULE_PATH "cmake/findFFTW" "fftw3/findFFTW" )

#Include FFTW libraries
find_package( FFTW )

#Bring the headers into the project
include_directories( include FFTW_INCLUDE_DIRS )

#Manually add the sources using the set command
#set(SOURCES src/[src1].cpp src/[src2].cpp)

#Use file(GLOB...) for wildcard additions
file(GLOB SOURCES "src/*.cpp")

add_executable(hydroTest ${SOURCES})