cmake_minimum_required(VERSION 3.20)
project(testing)

add_executable(testing main.cpp)
target_include_directories(testing PUBLIC ..)