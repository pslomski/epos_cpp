cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(epos VERSION 1.0)

configure_file(epos_config.h.in epos_config.h)

include_directories(include epos)
add_subdirectory(epos)
add_subdirectory(nric)
add_subdirectory(step1)
add_subdirectory(test)