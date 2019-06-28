cmake_minimum_required(VERSION 3.14)
project(NewAnalyzer)

set(CMAKE_CXX_STANDARD 14)


find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories(.)



add_executable(NewAnalyzer vDriftAnalyzerV.cpp vDriftAnalyzerV.h)
target_link_libraries(NewAnalyzer PRIVATE ${ROOT_LIBRARIES})
