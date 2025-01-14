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
#include <random>

using namespace std::complex_literals;
using namespace std::numbers;

class SignalImpl{
  public:
    Eigen::MatrixXd *signal;
};

void Signal::toString(){
//  std::cout << *pSignal->signal << std::endl;
  std::cout << value << std::endl;
}

void Signal::illuminate_thermally( double coherence_diameter){
  int N = value.cols();
  int knl_size = (int)std::round( N * coherence_diameter / L);
  std::cout << "knl_size in pixels = " << knl_size << std::endl;
  //Eigen::MatrixXd phase = Eigen::MatrixXd::Random( N, N);
  Eigen::MatrixXd phase( N, N);
  phase = phase.unaryExpr( [this](double){ return dis( gen);});
  value = Eigen::exp((2i * pi * phase).array());

  double reslim = coherence_diameter/knl_size;
  auto dom_sq = Eigen::VectorXd::LinSpaced( knl_size, -coherence_diameter/2 + reslim, coherence_diameter/2).array().square().matrix();
  auto ones = Eigen::VectorXd::Constant( knl_size, 1);
  Eigen::MatrixXcd disk = (((dom_sq * ones.transpose()) + (ones * dom_sq.transpose())).array() < pow( coherence_diameter/2,2)).cast<double>();
  std::complex< double> *data = value.data();
  cv::Mat sgl = cv::Mat_< std::complex< double>>( N, N, data);
  std::complex< double> *knl_data = disk.data();
  cv::Mat knl = cv::Mat_<std::complex< double>>( knl_size, knl_size, knl_data);
  cv::filter2D( sgl, sgl, -1, knl);
  //std::cout << "cv::filter2D finished!" << std::endl;
  //cv::filter2D( sgl, sgl, -1, knl, cv::Point( -1, -1), 0, cv::BORDER_CONSTANT);
  //Eigen::Map< Eigen::Matrix< std::complex< double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> light( (std::complex< double> *)sgl.data, N, N);
  //detect( light, "light_intensity.tiff");
  //phase_detect( light, "light_phase.tiff");
  //detect( value, "input_intensity.tiff");
  //value = value.array() * light.array();
}

void Signal::illuminate_uniformly(){
  value.setZero();
  int N = value.cols();
  //value = MatrixXi::Random( N, N).cast<double>();
}

void Signal::triple_slit_mask( int w_ratio, int h_ratio, int slits){
  int N = value.cols();
  int hnw = N/w_ratio;
  if( hnw % 2){
    hnw--;
  }
  int hnh = N/h_ratio;
  if( hnh % 2){
    hnh--;
  }
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
  auto ones = Eigen::VectorXd::Constant( hnh, 1);
  value.block( (N-hnh)/2, 0, hnh, N) = value.block( (N-hnh)/2, 0, hnh, N).array() * (ones * line).array();
  value.block( 0, 0, (N-hnh)/2, N).setZero();
  value.block( (N+hnh)/2, 0, (N-hnh)/2, N).setZero();
}

Signal::Signal( double lambda, double side_length_in_meter, int N):
  lambda{ lambda}, L{ side_length_in_meter}, gen{ std::random_device{}()}, dis{ 0.0, 1.0}
{
  value = Eigen::MatrixXcd( N, N);
  illuminate_uniformly();
}
/*
void Signal__Signal( double lambda, double side_length_in_meter, int N, int w_ratio, int h_ratio, int slits):
  lambda{ lambda}, L{ side_length_in_meter}
{
 // pSignal = new SignalImpl();
  value = Eigen::MatrixXcd::Constant( N, N, 0);
  int w = N/w_ratio;
  int h = N/h_ratio;
  int x0 = (N - w)/2;
  int y0 = (N - h)/2;
  int sd_bw = 2; //ratio of distance between two contiguous slits over interslit width
  if( w >= 2 * slits - 1){

    //calculate approximated values for bw and sw
    int bw = w/(slits * sd_bw + 1); // interslit width
    int sw = (w - bw * (slits - 1))/slits; // slit width
 
    //do the best effort to approach the sd_bw condition
    while( (bw+sw)/bw > sd_bw){
      bw++;
      sw--;
    }
    int r = w - bw * (slits - 1) - sw * slits; // remainder
    auto slit = Eigen::MatrixXd::Constant( h, sw, 1);
    int count = 0;
    int offset = r/2;
    for( int x = x0 + offset; count < slits; x += sw + bw){
    std::cout << "x = " << x << std::endl;
      value.block( y0, x, h, sw) = slit;
      count++;
    }
    //pSignal->signal = &signal;
    //std::cout << *pSignal->signal << std::endl;
    std::cout << "sw = " << sw << std::endl;
    std::cout << "bw = " << bw << std::endl;
    std::cout << "r = " << r << std::endl;
    std::cout << "offset = " << offset << std::endl;
  }else{
    //pSignal->signal = &value;
    //std::cout << pSignal->signal << std::endl;
  }
}
*/

void Signal::detect( Eigen::MatrixXcd detecting, std::string filename, int bit_depth){
  Eigen::MatrixXd intensity = (detecting.array() * detecting.conjugate().array()).matrix().real();
  //Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).sqrt().matrix().real();

  store( filename, intensity, bit_depth);
}

double Signal::bucket(){
  Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).matrix().real();
  return intensity.sum();
}

void Signal::picture( std::string filename, int bit_depth){
  Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).matrix().real();
  //Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).sqrt().matrix().real();

  store( filename, intensity, bit_depth);
}

void Signal::phase_detect( Eigen::MatrixXcd detecting, std::string filename, int bit_depth){
  Eigen::MatrixXd phase = detecting.array().arg();

  store( filename, phase, bit_depth);
}

void Signal::phase_picture( std::string filename, int bit_depth){
  Eigen::MatrixXd phase = value.array().arg();

  store( filename, phase, bit_depth);
}

void Signal::bucket( std::string filename){
  std::ofstream file( filename);
  if(!file){
    std::cerr << "[core/signal] File {" + filename + "} not created." << std::endl;
  }
  file << std::setprecision( std::numeric_limits<double>::max_digits10) << bucket();
  file.close();
}

void Signal::store( std::string filename, Eigen::MatrixXd data_to_store, int bit_depth){
  //remap linearly with max intensity to pow( 2, bit_depth)-1
  int maxTiff = 65535;
  if( 64!=bit_depth){
     maxTiff = pow( 2, bit_depth)-1;
  }
  data_to_store /= data_to_store.maxCoeff() / maxTiff;
  int C = data_to_store.cols();
  int R = data_to_store.rows();
  auto typeTiff = CV_64FC1;
  std::vector<int> params;
  params.push_back( cv::IMWRITE_TIFF_COMPRESSION);
  params.push_back( 1);
  if( 8 == bit_depth){
    typeTiff = CV_8UC1;
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> converted_data = data_to_store.cast<uint8_t>();
    uint8_t *data = converted_data.data();
    cv::Mat m( R, C, typeTiff, data);
    cv::imwrite( filename, m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
  }else if( 16 == bit_depth){
    typeTiff = CV_16UC1;
    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> converted_data = data_to_store.cast<uint16_t>();
    uint16_t *data = converted_data.data();
    cv::Mat m( R, C, typeTiff, data);
    cv::imwrite( filename, m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
  }else{
    double *data = data_to_store.data();
    cv::Mat m( R, C, typeTiff, data);
    cv::imwrite( filename, m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
  }
}

void Signal::quadratic_phase_lag_shift( double k, double c){
  int N = value.cols();
  double reslim = L/N;
  auto dom_sq = Eigen::VectorXd::LinSpaced( N, -L/2 + reslim, L/2).array().square().matrix();
  auto ones = Eigen::VectorXd::Constant( N, 1);
  //value = ( 1i * k / 2.0 * c * dom * dom.transpose()).array().exp() * value.array();
  value = ( 1i * k / 2.0 * c * ((dom_sq * ones.transpose()) + (ones * dom_sq.transpose()))).array().exp() * value.array();
}

void Signal::quadratic_phase_lag_shift_fourier( double k, double c){
  int N = value.cols();
  double f0 = 1/L;
  double reslim = f0;
  double f_Nyq = N * f0 / 2.0;
  auto dom1 = Eigen::VectorXd::LinSpaced( N/2+1, 0, N/2 * f0);
  auto dom2 = Eigen::VectorXd::LinSpaced( N/2-1, ((N/2+1)-N) * f0, ((N - 1) - N) * f0);
  Eigen::VectorXd dom( dom1.size() + dom2.size());
  dom << dom1, dom2;
  //std::cout << dom << std::endl;
  auto dom_sq = dom.array().square().matrix();
  auto ones = Eigen::VectorXd::Constant( N, 1);
  //value = ( 1i * k / 2.0 * c * dom * dom.transpose()).array().exp() * value.array();
  value = ( 1i * k / 2.0 * c * ((dom_sq * ones.transpose()) + (ones * dom_sq.transpose()))).array().exp() * value.array();
}

Eigen::MatrixXcd diskmask_it( Eigen::MatrixXcd to_mask, double radius){
  int N = to_mask.cols();
  double reslim = 2 * radius/N;
  auto dom_sq = Eigen::VectorXd::LinSpaced( N, -radius + reslim, radius).array().square().matrix();
  auto ones = Eigen::VectorXd::Constant( N, 1);
  auto masked = (((dom_sq * ones.transpose()) + (ones * dom_sq.transpose())).array() < pow( radius,2)).cast<double>() * to_mask.array();
  //cout << to_mask.real() << std::endl;
  return masked;
}

void Signal::mask( double radius){
  int N = value.cols();
  double reslim = L/N;
  auto dom_sq = Eigen::VectorXd::LinSpaced( N, -L/2 + reslim, L/2).array().square().matrix();
  auto ones = Eigen::VectorXd::Constant( N, 1);
  value = (((dom_sq * ones.transpose()) + (ones * dom_sq.transpose())).array() < pow( radius,2)).cast<double>() * value.array();
  //cout << value.real() << std::endl;
}

void Signal::propagate( double dist){
  if( 0 == dist){
    return;
  }
  int N = value.cols();
  fftw_make_planner_thread_safe();
  //std::cout << "sizeof( fftw_complex) = " << sizeof( fftw_complex) << std::endl;
  //std::cout << "sizeof( double) = " << sizeof( double) << std::endl;
  fftw_complex* in = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_complex* out = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_complex* in2 = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_complex* out2 = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_plan p_fw = fftw_plan_dft_2d( N, N, in , out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p_bw = fftw_plan_dft_2d( N, N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
  memcpy( in, value.data(), N * N * sizeof( fftw_complex));
  fftw_execute( p_fw);
  memcpy( value.data(), out, N * N * sizeof( fftw_complex));
  //Eigen::Map< Eigen::Matrix< complex< double>, N, N, Eigen::RowMajor>> value_ft( reinterpret_cast<complex<double>*>(*out));
  value = value.transpose().eval();
  double k = 2 * pi / lambda;
  double c = -lambda*lambda * dist;
  quadratic_phase_lag_shift_fourier( k, c);
  memcpy( in2, value.data(), N * N * sizeof( fftw_complex));
  fftw_execute( p_bw);
  memcpy( value.data(), out2, N * N * sizeof( fftw_complex));
  value = value.transpose().eval();
  fftw_destroy_plan( p_fw);
  fftw_destroy_plan( p_bw);
  fftw_free( in);
  fftw_free( out);
  fftw_free( in2);
  fftw_free( out2);
}


