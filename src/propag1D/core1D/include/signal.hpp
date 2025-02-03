#include <string>
#include <eigen3/Eigen/Dense>
#include <random>
#include <filesystem>

class Signal{
  public:
    Signal( double lambda, double side_length_in_meter, int N);
    void toString();
    void toTextFile( std::filesystem::path filename);
    void propagate( double dist);
    void quadratic_phase_lag_shift_fourier( double k, double c);
    int illuminate_thermally( double coherence_diameter);
    void triple_slit_mask( int w_ratio, int h_ratio, int slits, double w_offset_ratio);
    double max_intensity();
    double max_intensity( Eigen::VectorXcd);
  private:
    void illuminate_uniformly();
    double lambda;
    double L;
    Eigen::VectorXcd value;
    //class SignalImpl *pSignal;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
};

