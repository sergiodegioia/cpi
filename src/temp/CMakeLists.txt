cmake_minimum_required( VERSION 3.31)
set( CMAKE_CXX_STANDARD 20)
project( fft1D LANGUAGES CXX VERSION 1.0.0)
add_executable( fft1D main.cpp)
add_executable( fft2D fft2D.cpp)
find_package( PkgConfig REQUIRED)
pkg_check_modules( fftw3 REQUIRED IMPORTED_TARGET fftw3)
target_link_libraries( fft1D PRIVATE PkgConfig::fftw3)
target_link_libraries( fft2D PRIVATE PkgConfig::fftw3)


