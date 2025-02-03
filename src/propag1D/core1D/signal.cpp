#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <complex>
#include <cmath>
#include <signal.hpp>
#include <eigen3/Eigen/Dense>
#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <fftw3.h>
#include <numbers>
#include <cstring>
#include <string>
#include <random>

using namespace std::complex_literals;
using namespace std::numbers;

std::string matTypeToString(int type) {
    std::string r;

    uchar depth = type & CV_MAT_DEPTH_MASK;
    uchar chans = 1 + (type >> CV_CN_SHIFT);

    switch (depth) {
        case CV_8U:  r = "8U"; break;
        case CV_8S:  r = "8S"; break;
        case CV_16U: r = "16U"; break;
        case CV_16S: r = "16S"; break;
        case CV_32S: r = "32S"; break;
        case CV_32F: r = "32F"; break;
        case CV_64F: r = "64F"; break;
        default:     r = "User"; break;
    }

    r += "C";
    r += std::to_string(chans);

    return r;
}
void Signal::toTextFile( std::filesystem::path filename){
  std::ofstream file( filename);
  if(!file){
    std::cerr << "[core/signal] File {" + filename.string() + "} not created." << std::endl;
  }
  file << value;
  file.close();
}

void Signal::toString(){
  std::cout << value << std::endl;
}

int Signal::illuminate_thermally( double coherence_diameter){
  int N = value.rows();
  Eigen::VectorXd phase( N);
  phase = phase.unaryExpr( [this](double){ return dis( gen);});
  value = Eigen::exp((2i * pi * phase).array());

  double side = 5 * coherence_diameter;
  int knl_size = (int)std::round( N * side / L);
  double sigma = coherence_diameter/2;
  double reslim = side/knl_size;
  Eigen::VectorXcd gauss1D = Eigen::exp(Eigen::VectorXd::LinSpaced( knl_size, -side/2 + reslim, side/2).array().square()/(-2*sigma*sigma)).matrix();

  std::complex< double> *data = value.data();
  cv::Mat sgl = cv::Mat_< std::complex< double>>( 1, N, data);
  std::complex< double> *knl_data = gauss1D.data();
  cv::Mat knl = cv::Mat_<std::complex< double>>( 1, knl_size, knl_data);
  Eigen::VectorXd impulse( 1);
  impulse << 1;
  auto *impulse_data = impulse.data();
  cv::Mat impulse_knl = cv::Mat_<double>( 1, 1, impulse_data);
  //cv::sepFilter2D( sgl, sgl, -1, knl, impulse_knl);
  return knl_size;
}

void Signal::illuminate_uniformly(){
  value.setOnes();
}

void Signal::triple_slit_mask( int w_ratio, int h_ratio, int slits, double w_offset_ratio){
  int N = value.rows();
  int hnw = N/w_ratio;
  if( hnw % 2){
    hnw--;
  }
  int hnh = N/h_ratio;
  if( hnh % 2){
    hnh--;
  }
  int onw = static_cast<int>( hnw * w_offset_ratio);
  Eigen::RowVectorXcd line( N);
  line.setZero();
  int nw = hnw / slits;
  if( nw % 2){
    nw--;
  }
  if( nw < 2){
    nw = 2;
  }
  int offset = nw / 2;
  if( slits % 2){
    line.block( 0, N/2, 1, nw/2) = Eigen::RowVectorXcd::Constant( nw/2, 1);
    offset += nw;
    slits--;
  }
  for( int i = 0; i < slits/2; i++){
    line.block( 0, N/2 + offset, 1, nw) = Eigen::RowVectorXcd::Constant( nw, 1);
    offset =+ 2 * nw;
  }
  line += line.reverse().eval();
  // apply the vertical offset
  if( 0 != onw){
    Eigen::RowVectorXcd temp( N);
    temp.setZero();
    if( onw > 0){
      temp.block( 0, onw, 1, N-onw) = line.block( 0, 0, 1, N-onw);
    }else if( onw < 0){
      onw = -onw;
      temp.block( 0, 0, 1, N-onw) = line.block( 0, onw, 1, N-onw);
    }
    line = temp;
  }
  //
  value = value.array() * line.transpose().array();
}

Signal::Signal( double lambda, double side_length_in_meter, int N):
  lambda{ lambda}, L{ side_length_in_meter}, gen{ std::random_device{}()}, dis{ 0.0, 1.0}
{
  value = Eigen::VectorXcd( N);
  illuminate_uniformly();
}

void Signal::quadratic_phase_lag_shift_fourier( double k, double c){
  int N = value.rows();
  double f0 = 1/L;
  double reslim = f0;
  double f_Nyq = N * f0 / 2.0;
  double dc = f0/2;
  Eigen::VectorXd dom1, dom2;
  if( N%2){ // if N is odd?
    dom1 = Eigen::VectorXd::LinSpaced( (N-1)/2+1, dc, (N-1)/2 * f0 + dc);
    dom2 = Eigen::VectorXd::LinSpaced( (N-1)/2, ((N-1)/2) * f0 + dc, -f0 + dc);
  }else{ // if N is ven
    dom1 = Eigen::VectorXd::LinSpaced( N/2, dc, (N/2 - 1) * f0 + dc);
    dom2 = Eigen::VectorXd::LinSpaced( N/2, -N/2 * f0 + dc, -f0 + dc);
  }
  Eigen::VectorXd dom( dom1.size() + dom2.size());
  dom << dom1, dom2;
  //std::cout << dom << std::endl;
  auto dom_sq = dom.array().square().matrix();
  value = ( 1i * k / 2.0 * c * dom_sq).array().exp() * value.array();
}

void Signal::propagate( double dist){
  if( 0 == dist){
    return;
  }
  int N = value.rows();
  fftw_make_planner_thread_safe();
  fftw_complex* in = (fftw_complex*)fftw_malloc( N * sizeof( fftw_complex));
  fftw_complex* out = (fftw_complex*)fftw_malloc( N * sizeof( fftw_complex));
  fftw_complex* in2 = (fftw_complex*)fftw_malloc( N * sizeof( fftw_complex));
  fftw_complex* out2 = (fftw_complex*)fftw_malloc( N * sizeof( fftw_complex));
  fftw_plan p_fw = fftw_plan_dft_1d( N, in , out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p_bw = fftw_plan_dft_1d( N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
  memcpy( in, value.data(), N * sizeof( fftw_complex));
  fftw_execute( p_fw);
  memcpy( value.data(), out, N * sizeof( fftw_complex));
  double k = 2 * pi / lambda;
  double c = -lambda*lambda * dist;
  quadratic_phase_lag_shift_fourier( k, c);
  memcpy( in2, value.data(), N * sizeof( fftw_complex));
  fftw_execute( p_bw);
  memcpy( value.data(), out2, N * sizeof( fftw_complex));
  fftw_destroy_plan( p_fw);
  fftw_destroy_plan( p_bw);
  fftw_free( in);
  fftw_free( out);
  fftw_free( in2);
  fftw_free( out2);
}

double Signal::max_intensity( Eigen::VectorXcd data){
    return (data.array() * data.conjugate().array()).matrix().real().maxCoeff();
}

double Signal::max_intensity(){
    return (value.array() * value.conjugate().array()).matrix().real().maxCoeff();
}

