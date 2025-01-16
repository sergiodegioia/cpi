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
  int slits = cfg.get_int( "slits");
  int w_ratio = cfg.get_int( "w_ratio");
  int h_ratio = cfg.get_int( "h_ratio");
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
  double max_intens;
  {
    Signal input( lambda, side_length_in_meter, N);
    input.illuminate_thermally( speckle_diameter);
    max_intens = input.max_intensity();
  }

  int i;
    //#pragma omp parallel for default( none) shared( intensity_factor, max_intens, lambda, side_length_in_meter, speckle_diameter, N, w_ratio, h_ratio, slits, object_to_lens, radius, lens_to_detectorA, lens_to_detectorB, frames, std::cout) private( i) schedule( static)
  _Pragma( "omp parallel for default( none) shared( intensity_factor, max_intens, lambda, side_length_in_meter, speckle_diameter, N, w_ratio, h_ratio, slits, object_to_lens, radius, lens_to_detectorA, lens_to_detectorB, frames, std::cout) private( i) schedule( static)")
  for( i = 0; i < frames; i++){
    std::cout << "Thread #" << std::to_string( omp_get_thread_num()) << " is running iteration i=" << std::to_string( i) << std::endl;
    Signal input( lambda, side_length_in_meter, N);
    input.illuminate_thermally( speckle_diameter);
    input.picture("reference" + seq( i, frames) + "_8bit.tiff", intensity_factor * max_intens, 8);
    //input.bucket("ref_bucket" + seq( i, frames) + "_8bit.txt", intensity_factor * max_intens);
    input.triple_slit_mask( w_ratio, h_ratio, slits);
    /*
     //START: comment out for CPI/uncomment for GI
    input.bucket("bucket" + seq( i, frames) + "_8bit.txt", intensity_factor * max_intens);
     //END: comment out for CPI/uncomment for GI
     */
    /*
     //START: uncomment for CPI/comment out for GI
     */
    input.propagate( object_to_lens);
    input.mask( radius);
    Signal secondBeam = input;
    input.propagate( -lens_to_detectorA);
    input.picture("toDetectorA" + seq( i, frames) + ".tiff", intensity_factor * max_intens, 8);
    secondBeam.propagate( - lens_to_detectorB);
    secondBeam.picture("toDetectorB" + seq( i, frames) + ".tiff", intensity_factor * max_intens, 8);
    /*
     //END: uncomment for CPI/comment out for GI
    */
  }
  //return draw( a, N, x_p, t_p);
}

