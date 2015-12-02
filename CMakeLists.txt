cmake_minimum_required(VERSION 3.3)
project(poset_builder)

set(CMAKE_C_FLAGS "-std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        main.c
        complex_struct.h
        complex_struct.c
        simplex_basic.h
        simplex_basic.c
        file_util.h
        file_util.c
        )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/build/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/build/bin)

add_executable(poset_builder ${SOURCE_FILES})