cmake_minimum_required(VERSION 3.16)

project(cynomath LANGUAGES C)

set(CYNOMATH_SOURCES
    BlackScholes.c
    BlackScholesBarrierOption.c
    BlackSholesDoubleBarrier.c
    cynoMathUtility.c
    Integration_fin.c
    integration_GH.c
    integration_GL.c
    matrix.c
    MTRandomNumber.c
    optimizer.c
    pch.c
    randomNumber.c
    rootFinder.c
    sobol.c
    sorting.c
    specialFunctions.c
)

add_library(cynomath SHARED ${CYNOMATH_SOURCES})

target_include_directories(cynomath
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

set_target_properties(cynomath PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME "cynomath"
)
