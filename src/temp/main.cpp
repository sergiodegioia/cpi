#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <complex>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fftw3.h>
#include <numbers>
#include <cstring>
#include <string>
#include <random>

using namespace std::complex_literals;
using namespace std::numbers;

int main(){
  double f = 2;
  int N = 1024;
  Eigen::VectorXd t = Eigen::VectorXd::LinSpaced( N, 0, 1);
  Eigen::VectorXcd value = (1i * 2.0 * pi * f * t).array().exp();
  fftw_make_planner_thread_safe();
  //std::cout << "sizeof( fftw_complex) = " << sizeof( fftw_complex) << std::endl;
  //std::cout << "sizeof( double) = " << sizeof( double) << std::endl;
  fftw_complex* in = (fftw_complex*)fftw_malloc( N * sizeof( fftw_complex));
  fftw_complex* out = (fftw_complex*)fftw_malloc( N * sizeof( fftw_complex));
  fftw_plan p_fw = fftw_plan_dft_1d( N, in , out, FFTW_BACKWARD, FFTW_ESTIMATE);
  memcpy( in, value.data(), N * sizeof( fftw_complex));
  fftw_execute( p_fw);
  memcpy( value.data(), out, N * sizeof( fftw_complex));
  std::ofstream file("complex_values.txt");
  if( file.is_open()){
      file << value << std::endl;
  }else{
    std::cout << value << std::endl;
  }
  file.close();
}


