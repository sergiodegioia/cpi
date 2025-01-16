cmake_minimum_required( VERSION 3.31)
project( cpi LANGUAGES CXX VERSION 1.0.0)
set( CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native") 
add_subdirectory( src bin)

