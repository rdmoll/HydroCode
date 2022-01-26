# Find NetCDF files

if( NOT NETCDF_ROOT AND DEFINED ENV{NETCDFDIR} )
    set( NETCDF_ROOT $ENV{NETCDFDIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT FFTW_ROOT )
    pkg_check_modules( PKG_NETCDF QUIET "netcdf" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

if( ${FFTW_USE_STATIC_LIBS} )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )
endif()

if( NETCDF_ROOT )
    # find libs

    find_library(
	NETCDF_DOUBLE_LIB
	NAMES "netcdf" libnetcdf
	PATHS ${NETCDF_ROOT}
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
	)

    find_library(
	NETCDF_CPP_LIB
	NAMES "netcdf-cxx4" libnetcdf-cxx4
	PATHS ${NETCDF_ROOT}
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
	)

    find_library(
	NETCDF_FLOAT_LIB
	NAMES "netcdff" libnetcdf
	PATHS ${NETCDF_ROOT}
	PATH_SUFFIXES "lib" "lib64"
	NO_DEFAULT_PATH
	)

    #find includes

    find_path(NETCDF_INCLUDE_DIRS
        NAMES "netcdf.h"
        PATHS ${NETCDF_ROOT}
        PATH_SUFFIXES "include"
        NO_DEFAULT_PATH
        )

else()
    # find libs

    find_library(
        NETCDF_DOUBLE_LIB
        NAMES "netcdf"
        PATHS ${PKG_NETCDF_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
        )

    find_library(
        NETCDF_CPP_LIB
        NAMES "netcdf-cxx4"
        PATHS ${PKG_NETCDF_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
        )

    find_library(
        NETCDF_FLOAT_LIB
        NAMES "netcdff"
        PATHS ${PKG_NETCDF_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
        )

    # find includes

    find_path(NETCDF_INCLUDE_DIRS
        NAMES "netcdf.h"
        PATHS ${PKG_NETCDF_LIBRARY_DIRS} ${INCLUDE_INSTALL_DIR}
        )

endif( NETCDF_ROOT )

if (NETCDF_DOUBLE_LIB)
    set(NETCDF_DOUBLE_LIB_FOUND TRUE)
    set(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} ${NETCDF_DOUBLE_LIB})
    add_library(NETCDF::Double INTERFACE IMPORTED)
    set_target_properties(NETCDF::Double
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${NETCDF_DOUBLE_LIB}"
        )
else()
    set(NETCDF_DOUBLE_LIB_FOUND FALSE)
endif()

if (NETCDF_CPP_LIB)
    set(NETCDF_CPP_LIB_FOUND TRUE)
    set(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} ${NETCDF_CPP_LIB})
    add_library(NETCDF::DoubleCpp INTERFACE IMPORTED)
    set_target_properties(NETCDF::DoubleCpp
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${NETCDF_CPP_LIB}"
        )
else()
    set(NETCDF_CPP_LIB_FOUND FALSE)
endif()

if (NETCDF_FLOAT_LIB)
    set(NETCDF_FLOAT_LIB_FOUND TRUE)
    set(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} ${NETCDF_FLOAT_LIB})
    add_library(NETCDF::Float INTERFACE IMPORTED)
    set_target_properties(NETCDF::Float
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${NETCDF_FLOAT_LIB}"
        )
else()
    set(NETCDF_FLOAT_LIB_FOUND FALSE)
endif()

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(NETCDF
    REQUIRED_VARS NETCDF_INCLUDE_DIRS
    HANDLE_COMPONENTS
    )

mark_as_advanced(
    NETCDF_INCLUDE_DIRS
    NETCDF_LIBRARIES
    NETCDF_LIB
    NETCDF_CPP_LIB
    NETCDF_FLOAT_LIB
    )
    
