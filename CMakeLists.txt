cmake_minimum_required(VERSION 3.1)
project(iketama)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4.1 -fopenmp -I/Users/robertrambo/Downloads/armadillo-5.300.4/include -I/Users/robertrambo/temp/fgt/include")

set(SOURCE_FILES main.cpp)
add_executable(iketama ${SOURCE_FILES})