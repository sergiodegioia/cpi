#include <iostream>
#include <fstream>
#include <sstream>
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

  Signal input( lambda, side_length_in_meter, N);
  input.illuminate_thermally( speckle_diameter);
  input.picture("light.tiff");
  input.triple_slit_mask( w_ratio, h_ratio, slits);
  //input.toString();
  std::cout << "illuminate_thermally finished!" << std::endl;
  input.picture("afterOjbect.tiff");
  double k = 2 * pi / lambda;
  input.propagate( object_to_lens);
  input.picture("propagatedTillLens.tiff");
  input.mask( radius);
  input.picture("afterLens.tiff");
  auto secondBeam = input;
  input.propagate( -lens_to_detectorA);
  input.picture("toDetectorA.tiff");

  
  secondBeam.propagate( - lens_to_detectorB);
  secondBeam.picture("toDetectorB.tiff");
  //return draw( a, N, x_p, t_p);
}

