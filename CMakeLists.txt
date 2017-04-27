cmake_minimum_required(VERSION 3.7)
project(MD)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        control.c
        coord.h
        MD.c
        util.c
        input.dat)

add_executable(MD ${SOURCE_FILES})