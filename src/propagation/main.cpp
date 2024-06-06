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

using namespace std;
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
  cout << "hello, cmake!!" << endl;
  QApplication a( argc, argv);
  ConfigLoader cfg;
  std::cout << "loaded parameters: " << cfg << std::endl;
  int slits = cfg.get_int( "slits");
  int w_ratio = cfg.get_int( "w_ratio");
  int h_ratio = cfg.get_int( "h_ratio");
  unsigned N = cfg.get_int( "N");
  double side_length_in_meter = cfg.get_double( "side_length_in_meter");
  double lambda = cfg.get_double( "lambda");
  Signal input( lambda, side_length_in_meter, N);
  double speckle_diameter = cfg.get_double( "speckle_diameter");
  input.illuminate_thermally( speckle_diameter);
  input.picture("input.tiff");
  input.triple_slit_mask( w_ratio, h_ratio, slits);
  //input.toString();
  std::cout << "illuminate_thermally finished!" << std::endl;
  input.picture("tripleslit.tiff");
  double k = 2 * pi / lambda;
  double fA = cfg.get_double( "common_lens.focal_length");
  double c = -1/fA;
  //input.toString();
  //input.quadratic_phase_lag_shift( k, c);
  //input.toString();
  //input.picture("afterLens.tiff");
  double radius = cfg.get_double("common_lens.aperture_diameter")/2;
  input.propagate( cfg.get_double( "common_arm.object_to_lens"));
  input.picture("propagated.tiff");
  input.mask( radius);
  input.picture("mask.tiff");
  input.propagate( - cfg.get_double( "common_arm.lens_to_BS"));
  input.picture("back_propagated.tiff");
  //return draw( a, N, x_p, t_p);
}

