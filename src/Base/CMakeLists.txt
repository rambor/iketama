cmake_minimum_required(VERSION 3.1)

PROJECT(base)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
#set(Boost_USE_MULTITHREADED ON)
# find_package(OpenMP)

#if(OPENMP_FOUND)
#    MESSAGE("** OPENMP FOUND")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(Boost_USE_MULTITHREADED ON)
#endif()

SET(BOOST_ROOT ~/boost)
find_package(Boost 1.57.0 REQUIRED)

#find_package(LAPACK PATHS /Users/robertrambo/LAPACK/lib)
find_library(BOOST_FILESYSTEM libboost_filesystem.a PATHS ~/boost ~/boost/boost_1_57_0/stage)
find_library(BOOST_PO libboost_program_options.a PATHS ~/boost ~/boost/boost_1_57_0/stage)
find_library(BOOST_RANDOM libboost_random.a)
find_library(BOOST_REGEX libboost_regex.a)
find_library(BOOST_SYSTEM libboost_system.a)
#find_package(OpenMP)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE("** BOOST FOUND")
endif()

set(SOURCE_FILES
        "${CMAKE_CURRENT_LIST_DIR}/vector3.h"
        "${CMAKE_CURRENT_LIST_DIR}/Phase.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Bead.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Data.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Datum.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Coords.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/PDBModel.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/functions.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/ThreadPool.h"
        "${CMAKE_CURRENT_LIST_DIR}/Component.cpp"
        )


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2  ")

add_library(base STATIC ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(base ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM})

target_include_directories (base PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

MESSAGE("** COMPILE FLAGS ${CMAKE_CXX_FLAGS}")

#list(APPEND SOURCE_FILES
#        "${CMAKE_CURRENT_LIST_DIR}/vector3.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/Phase.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/Bead.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/Data.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/Datum.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/Coords.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/PDBModel.cpp"
#        "${CMAKE_CURRENT_LIST_DIR}/functions.cpp"
#        )