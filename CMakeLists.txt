cmake_minimum_required(VERSION 3.16)

project(RiemannIntegrabilityTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(RiemannIntegrabilityTest main.cpp
    src/lstu/rit/Segment.hpp src/lstu/rit/Segment.cpp
    src/lstu/rit/Partition.hpp src/lstu/rit/Partition.cpp
    src/lstu/rit/epsmath.hpp )

target_include_directories(RiemannIntegrabilityTest PUBLIC ./src/lstu/rit)
