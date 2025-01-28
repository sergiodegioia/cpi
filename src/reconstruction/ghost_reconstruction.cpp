#include <iostream>
#include <opencv2/core/mat.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <eigen3/Eigen/Dense>

template< typename T>
void store( Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data_to_store);

int main(){
  std::cout << "ghost reconstruction!!!" << std::endl;
  Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> data_to_store = Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>::Constant( 8,8,255);
  store( data_to_store);
}
template< typename T>
void store( Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data_to_store){
  std::vector<int> params;
  params.push_back( cv::IMWRITE_TIFF_COMPRESSION);
  params.push_back( 1);
  auto typeTiff = CV_8UC1;
  int maxTiff;
  if constexpr (std::is_same< T, uint8_t>::value){
    typeTiff = CV_8UC1;
    maxTiff = 255;
  }else if constexpr (std::is_same< T, uint16_t>::value){
    typeTiff = CV_16UC1;
    maxTiff = 65535;
  }else if constexpr (std::is_same< T, double>::value){
    typeTiff = CV_64FC1;
    maxTiff = 65535;
  }
  //data_to_store /= max_value / (maxTiff + .0);
  uint8_t *data = data_to_store.data();
  int C = data_to_store.cols();
  int R = data_to_store.rows();
  cv::Mat m( R, C, typeTiff, data);
  cv::imwrite( "aTiff.tiff", m.t(), params); // data is in coloum-major order, cv::Mat wants it in row-major
}

