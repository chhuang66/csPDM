cmake_minimum_required(VERSION 2.6)
project(aeifonethree_semilocal_2)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# enable c++11
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib)

find_package(HDF5)
if(HDF5_FOUND)
    set(HDF5_INCLUDE_DIR "/home/justin/anaconda2/include/")
	include_directories(${HDF5_INCLUDE_DIR})
	add_definitions(-DHAVE_HDF5)
	message("HDF5 include directory: " ${HDF5_INCLUDE_DIR})
	message("HDF5 libraries: " ${HDF5_LIBRARIES})
endif()

find_package(BLAS)
find_package(LAPACK)

add_executable(aEIFONETHREE_semiLocal main.cpp EEG_run.cpp input_mod.cpp csvparser.cc)

install(TARGETS aEIFONETHREE_semiLocal RUNTIME DESTINATION bin)
target_link_libraries (aEIFONETHREE_semiLocal  cst pdm ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
