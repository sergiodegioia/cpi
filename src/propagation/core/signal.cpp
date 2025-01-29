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

class SignalImpl{
  public:
    Eigen::MatrixXd *signal;
};

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
  //file << std::setprecision( std::numeric_limits<double>::max_digits10) << bucket();
  file << value;
  file.close();
}

void Signal::toString(){
//  std::cout << *pSignal->signal << std::endl;
  std::cout << value << std::endl;
}

int Signal::illuminate_thermally( double coherence_diameter){
  int N = value.cols();
  //Eigen::MatrixXd phase = Eigen::MatrixXd::Random( N, N);
  Eigen::MatrixXd phase( N, N);
  phase = phase.unaryExpr( [this](double){ return dis( gen);});
  value = Eigen::exp((2i * pi * phase).array());

  double side = 5 * coherence_diameter;
  int knl_size = (int)std::round( N * side / L);
  double sigma = coherence_diameter/2;
  double reslim = side/knl_size;
  auto gauss1D = Eigen::exp(Eigen::VectorXd::LinSpaced( knl_size, -side/2 + reslim, side/2).array().square()/(-2*sigma*sigma)).matrix();
  Eigen::MatrixXcd gauss2D = gauss1D * gauss1D.transpose();

  std::complex< double> *data = value.data();
  cv::Mat sgl = cv::Mat_< std::complex< double>>( N, N, data);
  std::complex< double> *knl_data = gauss2D.data();
  cv::Mat knl = cv::Mat_<std::complex< double>>( knl_size, knl_size, knl_data);
  cv::filter2D( sgl, sgl, -1, knl);
  return knl_size;
}

void Signal::illuminate_thermally_flat_knl( double coherence_diameter){
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
  std::cout << "Kernel type: " << matTypeToString( knl.type()) << std::endl;
  std::cout << "Input type: " << matTypeToString( sgl.type()) << std::endl;
  cv::filter2D( sgl, sgl, -1, knl);
}

/*
void Signal::illuminate_thermally_fft( double coherence_diameter){
  int N = value.cols();
  int knl_size = (int)std::round( N * coherence_diameter / L);
  std::cout << "knl_size in pixels = " << knl_size << std::endl;
  //Eigen::MatrixXd phase = Eigen::MatrixXd::Random( N, N);
  Eigen::MatrixXd phase( N, N);
  phase = phase.unaryExpr( [this](double){ return dis( gen);});
  value = Eigen::exp((2i * pi * phase).array());

  fftw_make_planner_thread_safe();
  std::cout << "sizeof( fftw_complex) = " << sizeof( fftw_complex) << std::endl;
  std::cout << "sizeof( double) = " << sizeof( double) << std::endl;
  fftw_complex* in = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_complex* out = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_complex* in2 = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_complex* out2 = (fftw_complex*)fftw_malloc( N * N * sizeof( fftw_complex));
  fftw_plan p_fw = fftw_plan_dft_2d( N, N, in , out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_plan p_bw = fftw_plan_dft_2d( N, N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
  memcpy( in, value.data(), N * N * sizeof( fftw_complex));
  fftw_execute( p_fw);
  memcpy( value.data(), out, N * N * sizeof( fftw_complex));
  value = value.transpose().eval();

  double F = N/L/2;
  auto freq_sq = Eigen::VectorXd::LinSpaced( N, -F + 1/L, F).array().square().matrix();
  auto ones = Eigen::VectorXd::Constant( N, 1);
  Eigen::MatrixXcd knl_ft = Eigen::exp(2*pi*pi*coherence_diameter*coherence_diameter*(freq_sq * ones.transpose() + ones * freq_sq.transpose()).array());

  value = value * knl_ft;
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
*/

void Signal::illuminate_uniformly(){
  value.setOnes();
  //int N = value.cols();
  //value = MatrixXi::Random( N, N).cast<double>();
}

void Signal::triple_slit_mask( int w_ratio, int h_ratio, int slits, double w_offset_ratio){
  int N = value.cols();
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

void Signal::detect( Eigen::MatrixXcd detecting, std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth){
  Eigen::MatrixXd intensity = (detecting.array() * detecting.conjugate().array()).matrix().real();
  Eigen::MatrixXd result = intensity.array().pow( 1./1);
  //Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).sqrt().matrix().real();

  max_value = (result.array() * result.conjugate().array()).matrix().real().maxCoeff();
  double min_value = (result.array() * result.conjugate().array()).matrix().real().minCoeff();
  std::cout << "fft max intensity = " << max_value << std::endl;
  std::cout << "fft min intensity = " << min_value << std::endl;
  store( filename, result, max_value, centered_crop_factor, bit_depth);
}

double Signal::bucket(){
  Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).matrix().real();
  return intensity.sum();
}

void Signal::picture( std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth){
  Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).matrix().real();
  //Eigen::MatrixXd intensity = (value.array() * value.conjugate().array()).sqrt().matrix().real();

  store( filename, intensity, max_value, centered_crop_factor, bit_depth);
}

void Signal::phase_detect( Eigen::MatrixXcd detecting, std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth){
  Eigen::MatrixXd phase = detecting.array().arg();

  store( filename, phase, max_value, bit_depth, centered_crop_factor);
}

void Signal::phase_picture( std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth){
  Eigen::MatrixXd phase = value.array().arg();

  store( filename, phase, max_value, centered_crop_factor, bit_depth);
}

void Signal::bucket( std::filesystem::path filename, double norm_fact){
  std::ofstream file( filename);
  if(!file){
    std::cerr << "[core/signal] File {" + filename.string() + "} not created." << std::endl;
  }
  //file << std::setprecision( std::numeric_limits<double>::max_digits10) << bucket();
  file << static_cast<int>(std::round( bucket()/norm_fact));
  file.close();
}

void Signal::store( std::filesystem::path filename, Eigen::MatrixXd data_to_store, double max_value, double centered_crop_factor, int bit_depth){
  if( centered_crop_factor > 1){
    centered_crop_factor = 1;
  }
  int C = data_to_store.cols();
  int R = data_to_store.rows();
  int i = static_cast<int>(std::round((R - R*centered_crop_factor)/2));
  int j = static_cast<int>(std::round((C - C*centered_crop_factor)/2));
  C = C*centered_crop_factor;
  R = R*centered_crop_factor;
  auto new_data_to_store = data_to_store.block( i, j, R, C);
  //remap linearly with max intensity to pow( 2, bit_depth)-1
  int maxTiff = 65535;
  if( 64!=bit_depth){
     maxTiff = pow( 2, bit_depth)-1;
  }
  //data_to_store /= data_to_store.maxCoeff() / maxTiff;
  //std::cout << "max value of data to store = " << data_to_store.maxCoeff() << "; norm factor = " << max_value << "; max value Tiff = " << maxTiff << std::endl;
  //double min_store = new_data_to_store.minCoeff();
  //double max_store = new_data_to_store.maxCoeff();
  //std::cout << "max_value = " << max_value << std::endl;
  //std::cout << "min_store = " << min_store << std::endl;
  //std::cout << "max_store = " << max_store << std::endl;
  new_data_to_store /= max_value / (maxTiff + .0);
  //C = data_to_store.cols();
  //R = data_to_store.rows();
  auto typeTiff = CV_64FC1;
  std::vector<int> params;
  params.push_back( cv::IMWRITE_TIFF_COMPRESSION);
  params.push_back( 1);
  if( 8 == bit_depth){
    typeTiff = CV_8UC1;
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> converted_data = new_data_to_store.cast<uint8_t>();
    uint8_t *data = converted_data.data();
    cv::Mat m( R, C, typeTiff, data);
    cv::imwrite( filename.string(), m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
  }else if( 16 == bit_depth){
    typeTiff = CV_16UC1;
    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> converted_data = new_data_to_store.cast<uint16_t>();
    uint16_t *data = converted_data.data();
    cv::Mat m( R, C, typeTiff, data);
    cv::imwrite( filename.string(), m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
  }else{
    double *data = new_data_to_store.data();
    cv::Mat m( R, C, typeTiff, data);
    cv::imwrite( filename.string(), m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
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
  double dc = f0/2;
  //auto dom1 = Eigen::VectorXd::LinSpaced( N/2+1, dc, N/2 * f0 + dc);
  //auto dom2 = Eigen::VectorXd::LinSpaced( N/2-1 + (N%2?1:0), ((N/2+1)-N) * f0 + dc, ((N - 1) - N) * f0 + dc);
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

Eigen::MatrixXcd Signal::fft(){
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
  auto fft_data = value.transpose().eval();
  memcpy( in, fft_data.data(), N * N * sizeof( fftw_complex));
  fftw_execute( p_fw);
  memcpy( fft_data.data(), out, N * N * sizeof( fftw_complex));
  fftw_destroy_plan( p_fw);
  fftw_free( in);
  fftw_free( out);
  fft_data /= N;
  /*
  memcpy( in2, fft_data.data(), N * N * sizeof( fftw_complex));
  fftw_execute( p_bw);
  memcpy( fft_data.data(), out2, N * N * sizeof( fftw_complex));
  fft_data = fft_data.transpose().eval();
  fft_data /= N;
  fftw_destroy_plan( p_bw);
  fftw_free( in2);
  fftw_free( out2);
  */
  return fft_data;
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
  memcpy( in, value.transpose().data(), N * N * sizeof( fftw_complex));
  fftw_execute( p_fw);
  memcpy( value.data(), out, N * N * sizeof( fftw_complex));
  value = value.transpose().eval();
  //Eigen::Map< Eigen::Matrix< complex< double>, N, N, Eigen::RowMajor>> value_ft( reinterpret_cast<complex<double>*>(*out));
  double k = 2 * pi / lambda;
  double c = -lambda*lambda * dist;
  quadratic_phase_lag_shift_fourier( k, c);
  memcpy( in2, value.transpose().data(), N * N * sizeof( fftw_complex));
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

double Signal::max_intensity( Eigen::MatrixXcd data){
    return (data.array() * data.conjugate().array()).matrix().real().maxCoeff();
}

double Signal::max_intensity(){
    return (value.array() * value.conjugate().array()).matrix().real().maxCoeff();
}

