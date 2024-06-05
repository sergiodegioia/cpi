cmake_minimum_required( VERSION 3.20)
project( propagation CXX)
set( CMAKE_CXX_STANDARD 20)
find_package( Qt6 REQUIRED COMPONENTS Core Widgets)
qt_standard_project_setup()
qt_add_executable( propagation main.cpp)
list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake/module")
find_package( Qwt REQUIRED)
find_package( OpenCV REQUIRED)
add_subdirectory( core)
target_link_libraries( propagation PRIVATE Qt6::Widgets Qwt::Qwt ${OpenCV_LIBS} core)


