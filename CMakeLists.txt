cmake_minimum_required(VERSION 3.1)
project(iketama)
#
# armadillo and fgt libraries are explicitly defined
# must install and change file path for new installations
# explicit paths:
# armadillo
# GFORTRAN_LIB
# FGT
# QHULL
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4.1 -fopenmp -I/Users/robertrambo/Downloads/armadillo-5.300.4/include -I/Users/robertrambo/temp/fgt/include")

MESSAGE("** SYSTEM : ${CMAKE_SYSTEM_NAME}")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    MESSAGE("** ${CMAKE_LIBRARY_PATH}")
    # Directory Path will need to change on a different system
    # I have installed GFORTRAN as my compiler but if Intel is installed, then some details below will need to change
    find_path(GFORTRAN_LIB libgfortran.a /sw /sw/lib /sw/lib/gcc4.9/lib /opt /opt/local /opt/local/lib/gcc49)
    MESSAGE("** SETTING  GFORTRAN_LIB=${GFORTRAN_LIB}")

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(CMAKE_VERBOSE_MAKEFILE ON )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

#
# set openMP (does not work for AppleClang g++)
#

set(Boost_USE_MULTITHREADED ON)

find_package(OpenMP)

if(OPENMP_FOUND)
    MESSAGE("** OPENMP FOUND")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(Boost_USE_MULTITHREADED ON)
endif()

SET (QHULL_INCLUDE_DIRS ~/qhull/src ~/qhull/src/libqhull ~/qhull/src/libqhullcpp)
include_directories(${QHULL_INCLUDE_DIRS})
#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)

SET(BOOST_ROOT ~/boost)
find_package(Boost 1.57.0 REQUIRED)

# if you have boost in a non-standard directory, like you home directory, set the following:
#/Users/robertrambo/boost/boost_1_57_0
#/Users/robertrambo/boost/boost_1_57_0/stage/lib

find_library(QHULL_LIB libqhullstatic.a PATHS ~/qhull/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/qhull/lib)
find_library(ARMA_LIB libarmadillo.a HINTS /Users/robertrambo/Downloads/armadillo-5.300.4)
find_library(FGT_LIB libfgt.a HINTS /Users/robertrambo/temp/fgt/build/lib)
# find_library(CPD_LIB libcpd.a PATHS ${CMAKE_CURRENT_SOURCE_DIR}/src/CPD/build)

# MESSAGE("** CPD LIB => ${CPD_LIB}" )
find_library(BOOST_FILESYSTEM libboost_filesystem.a PATHS ~/boost ~/boost/boost_1_57_0/stage)
find_library(BOOST_PO libboost_program_options.a PATHS ~/boost ~/boost/boost_1_57_0/stage)
find_library(BOOST_RANDOM libboost_random.a)
find_library(BOOST_REGEX libboost_regex.a)
find_library(BOOST_SYSTEM libboost_system.a)

set(SOURCE_FILES "")
# order of include may be important
add_subdirectory (src/Base/)

include(src/CMakeLists.txt)
include(src/EulerTour/CMakeLists.txt)
include(src/Annealer/CMakeLists.txt)

MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE("** BOOST FOUND")
endif()

MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIRS}")
MESSAGE("** CURRENT SOURCE DIRECTORY DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

# set(APPEND SOURCE_FILES src/main.cpp)
# create object file then link to make executable
add_executable(iketama ${SOURCE_FILES} src/main.cpp)
SET_TARGET_PROPERTIES(iketama PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES(iketama ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM} base ${QHULL_LIB} ${QHULLCPP_LIB} ${ARMA_LIB} ${FGT_LIB})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE("** with ${CMAKE_CXX_COMPILER_ID}")
    SET(CMAKE_EXE_LINKER_FLAGS "-L${GFORTRAN_LIB} -lgfortran -L/Users/robertrambo/LAPACK/lib -llapack -lblas -ltmglib" )
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("** LINKING WILL BE PERFORMED WITH -lgfortran")
    SET(CMAKE_EXE_LINKER_FLAGS "-L${GFORTRAN_LIB} -lgfortran" )
    # if using icc and gfortran
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

    # if using MS Visual Studio and gfortran
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

endif()

# TARGET_LINK_LIBRARIES(Bubbles ${WETSAXSF90})