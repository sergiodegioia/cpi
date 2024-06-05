#include <iostream>
#include "QApplication"
#include "qwt_plot.h"
#include "qwt_plot_curve.h"
#include <signal.hpp>
#include <numbers>
#include <cmath>

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
  unsigned N = 2;
  int slits = stoi( string( argv[ 1]));
  int w_ratio = stoi( string( argv[ 2]));
  int h_ratio = 3;
  N = 1116;
  double side_length_in_meter = 1.0e-3;
  double lambda = 632.8e-9;
  Signal input( lambda, side_length_in_meter, N, w_ratio, h_ratio, slits);
  //input.toString();
  input.picture("input.tiff");
  input.illuminate_thermally( 10e-6);
  std::cout << "illuminate_thermally finished!" << std::endl;
  input.picture("convolution.tiff");
  double k = 2 * pi / lambda;
  double fA = 1.0e-3/2;
  double c = -1/fA;
  //input.toString();
  //input.quadratic_phase_lag_shift( k, c);
  //input.toString();
  //input.picture("afterLens.tiff");
  double radius = side_length_in_meter/4;
  input.propagate( 1.0e-3);
  input.picture("propagated.tiff");
  input.mask( radius);
  input.picture("mask.tiff");
  input.propagate( -1.0e-3);
  input.picture("back_propagated.tiff");
  //return draw( a, N, x_p, t_p);
}

