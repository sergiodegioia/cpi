cmake_minimum_required( VERSION 3.31)

project( reconstruction CXX)
set( CMAKE_CXX_STANDARD 20)
find_package( OpenCV REQUIRED)
add_executable( ghost_image ghost_reconstruction.cpp)
target_link_libraries( ghost_image PRIVATE ${OpenCV_LIBS})

