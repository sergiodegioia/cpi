cmake_minimum_required( VERSION 3.31)

project( propagation CXX)
set( CMAKE_CXX_STANDARD 20)
find_package( Qt6 REQUIRED COMPONENTS Core Widgets)
qt_standard_project_setup()
qt_add_executable( propag1D main.cpp)
list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/module")
find_package( Qwt REQUIRED)
find_package( OpenCV REQUIRED)
add_subdirectory( core1D)
add_subdirectory( utils)
find_package( OpenMP REQUIRED)
#target_link_libraries( propag1D PRIVATE Qt6::Widgets Qwt::Qwt ${OpenCV_LIBS} core utils)
target_link_libraries( propag1D PRIVATE Qt6::Widgets Qwt::Qwt ${OpenCV_LIBS} core1D utils OpenMP::OpenMP_CXX)


