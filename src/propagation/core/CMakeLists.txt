cmake_minimum_required( VERSION 3.31)
project( propagation-core LANGUAGES CXX VERSION 1.0.0)
add_library( core OBJECT signal.cpp)
target_include_directories( core PUBLIC include)
find_package( OpenCV REQUIRED)
find_package( PkgConfig REQUIRED)
pkg_check_modules( fftw3 REQUIRED IMPORTED_TARGET fftw3)
target_link_libraries( core PRIVATE ${OpenCV_LIBS} PkgConfig::fftw3)

