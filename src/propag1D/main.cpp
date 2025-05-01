#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif

#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "QApplication"
#include "qwt_plot.h"
#include "qwt_plot_curve.h"
#include <signal.hpp>
#include <cfgldr.hpp>
#include <numbers>
#include <cmath>
#include <map>
#include <filesystem>

using namespace std::numbers;

int draw( QApplication& a, unsigned N, const double* x_p, const double* t_p);
int draw( QApplication& a, unsigned N, const double* x_p, const double* t_p){
  QwtPlot plot;
  plot.setTitle("signal");
  plot.resize( 400, 600);
  plot.show();
  QwtPlotCurve curve;
  QPolygonF pol;
  curve.setSamples( x_p, t_p, N);
  curve.attach( &plot);
  return a.exec();
}

inline std::string seq( int counter, int frames){
  if( frames == 1){
    return "";
  }
  double tot_figure = std::ceil( std::log10( frames));
  if( tot_figure < 1){
    throw std::runtime_error( "Too few frames: " + std::to_string( frames));
  }
  if( tot_figure == 1){
    return "_" + std::to_string( counter);
  }
  if( tot_figure == 2){
    if( counter < 10){
      return "_0" + std::to_string( counter);
    }
    return std::to_string( counter);
  }
  if( tot_figure == 3){
    if( counter < 10){
      return "_00" + std::to_string( counter);
    }
    if( counter < 100){
      return "_0" + std::to_string( counter);
    }
    return std::to_string( counter);
  }
  if( tot_figure == 4){
    if( counter < 10){
      return "_000" + std::to_string( counter);
    }
    if( counter < 100){
      return "_00" + std::to_string( counter);
    }
    if( counter < 1000){
      return "_0" + std::to_string( counter);
    }
    return std::to_string( counter);
  }
  double pads;
  if( counter == 0){
    pads = tot_figure - 1;
  }else{
    double log_counter = std::log10( counter);
    double figure = std::ceil( log_counter);
    if( figure == log_counter){
      figure++;
    }
    pads = tot_figure - figure;
  }
  std::string prefix = "_";
  for( int i = 0; i < pads; i++){
    prefix += "0";
  }
  return prefix + std::to_string( counter);
}

int main( int argc, char **argv){
  std::cout << "hello, cpi!!" << std::endl;
  QApplication a( argc, argv);
  // it is convenient to load all parameters at startup in order for 
  // a bad formatted config file to abort execution as soon as possible
  ConfigLoader cfg;
  std::cout << "loaded parameters: " << std::endl;
  std::cout << cfg << std::endl;
  constexpr std::string_view CPI = "CPI";
  constexpr std::string_view GI = "GI";
  std::filesystem::path pathname = cfg.get_dirname( "storage.root_dir");
  std::string experiment = cfg.get_experiment();
  int pixel_format = cfg.get_int("tiff.pixel_format");
  //std::string pixel_format = cfg.get_pixel_format();
  int slits = cfg.get_int( "object.slits");
  int w_ratio = cfg.get_int( "object.w_ratio");
  int h_ratio = cfg.get_int( "object.h_ratio");
  double w_offset_ratio = cfg.get_double( "object.vertical_offset_ratio");
  unsigned N = cfg.get_int( "N");
  double side_length_in_meter = cfg.get_double( "side_length_in_meter");
  double lambda = cfg.get_double( "lambda");
  double speckle_diameter = cfg.get_double( "speckle_diameter");
  double fA = cfg.get_double( "lens.focal_length");
  double radius = cfg.get_double("lens.aperture_diameter")/2;
  double object_to_lens = cfg.get_double( "common_arm.object_to_lens");
  double lens_to_detectorA = cfg.get_double( "common_arm.lens_to_detectorA");
  double lens_to_detectorB = cfg.get_double( "common_arm.lens_to_detectorB");
  int frames = cfg.get_int( "frames");

  double intensity_factor = 2.0;
  double max_intens_at_source;
  double max_intens_fft;
  double max_intens_at_lens;
  double max_intens_at_A;
  double max_intens_at_B;

  int i;
    //#pragma omp parallel for default( none) shared( intensity_factor, max_intens_at_source, max_intens_at_lens, max_intens_at_A, max_intens_at_B, lambda, fA, side_length_in_meter, speckle_diameter, N, w_ratio, h_ratio, slits, w_offset_ratio, object_to_lens, radius, lens_to_detectorA, lens_to_detectorB, frames, pathname, pixel_format, std::cout) private( i) schedule( static)
  _Pragma( "omp parallel for default( none) shared( intensity_factor,  max_intens_fft, max_intens_at_source, max_intens_at_lens, max_intens_at_A, max_intens_at_B, lambda, fA, side_length_in_meter, speckle_diameter, N, w_ratio, h_ratio, slits, w_offset_ratio, object_to_lens, radius, lens_to_detectorA, lens_to_detectorB, frames, pathname, pixel_format, std::cout) private( i) schedule( static)")
  for( i = 0; i < frames; i++){
    Signal input( lambda, side_length_in_meter, N);
    input.illuminate_thermally( speckle_diameter);
    input.triple_slit_mask( w_ratio, h_ratio, slits, w_offset_ratio);
    Signal second = input;
    input.propagate( object_to_lens);
    second.propagate( -object_to_lens);
    input.toTextFile( pathname / ("forward" + seq( i, frames) + ".txt"));
    second.toTextFile( pathname / ("backward" + seq( i, frames) + ".txt"));
  }

  //return draw( a, N, x_p, t_p);
}

