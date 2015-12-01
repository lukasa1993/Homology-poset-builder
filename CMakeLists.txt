cmake_minimum_required(VERSION 3.3)
project(poset_builder)

set(CMAKE_C_FLAGS "-std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

IF(true)
    message(STATUS "Checking OpenMP")
    find_package(OpenMP)
    IF(OPENMP_FOUND)
        message("Found OpenMP! ^_^")
        # add flags for OpenMP
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE()
        message("Missed OpenMP! x_x")
    ENDIF()
ENDIF()

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