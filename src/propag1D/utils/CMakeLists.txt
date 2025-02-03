cmake_minimum_required( VERSION 3.31)
project( propagation-utils LANGUAGES CXX VERSION 1.0.0)
add_library( utils1D OBJECT cfgldr.cpp)
target_include_directories( utils1D PUBLIC include)

