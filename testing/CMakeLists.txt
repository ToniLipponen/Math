cmake_minimum_required(VERSION 3.20)
project(testing)

set (CMAKE_CXX_STANDARD 14)
add_executable(testing main.cpp)
target_include_directories(testing PUBLIC ..)