#include <string>
#include <eigen3/Eigen/Dense>
#include <random>

class Signal{
  public:
    Signal( double lambda, double side_length_in_meter, int N);
    void detect( Eigen::MatrixXcd detecting, std::string filename, double max_value, int bit_depth);
    void phase_detect( Eigen::MatrixXcd detecting, std::string filename, double max_value, int bit_depth);
    void picture( std::string filename, double max_value, int bit_depth = 64);
    double bucket();
    void bucket( std::string filename, double norm_fact);
    void phase_picture( std::string filename, double max_value, int bit_depth);
    void toString();
    void quadratic_phase_lag_shift( double k, double c);
    void quadratic_phase_lag_shift_fourier( double k, double c);
    void mask( double radius);
    Eigen::MatrixXcd diskmask_it( Eigen::MatrixXcd to_mask, double radius);
    void propagate( double dist);
    void illuminate_thermally( double coherence_diameter);
    void illuminate_thermally_flat_knl( double coherence_diameter);
    //void illuminate_thermally_fft( double coherence_diameter);
    void triple_slit_mask( int w_ratio, int h_ratio, int slits);
    double max_intensity();
  private:
    void illuminate_uniformly();
    void store( std::string filename, Eigen::MatrixXd data_to_store, double max_value, int bit_depth);
    double lambda;
    double L;
    Eigen::MatrixXcd value;
    Eigen::MatrixXd x, y;
    //class SignalImpl *pSignal;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
};

