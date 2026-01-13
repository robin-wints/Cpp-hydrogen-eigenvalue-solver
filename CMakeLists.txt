cmake_minimum_required(VERSION 3.15)
project(HydrogenSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Add executable
add_executable(hydrogen
    main/hydrogen.cpp
    main/functions.cpp
)

# Eigen (header-only)
target_include_directories(hydrogen
    PRIVATE
        ${CMAKE_SOURCE_DIR}/eigen-5.0.0
)

# Include Spectra (header-only)
target_include_directories(hydrogen
    PRIVATE
        ${CMAKE_SOURCE_DIR}/Spectra-master/include
)

