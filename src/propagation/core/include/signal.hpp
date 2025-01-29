#include <string>
#include <eigen3/Eigen/Dense>
#include <random>
#include <filesystem>

class Signal{
  public:
    Signal( double lambda, double side_length_in_meter, int N);
    void detect( Eigen::MatrixXcd detecting, std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth);
    void phase_detect( Eigen::MatrixXcd detecting, std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth);
    void picture( std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth = 64);
    double bucket();
    void bucket( std::filesystem::path filename, double norm_fact);
    void phase_picture( std::filesystem::path filename, double max_value, double centered_crop_factor, int bit_depth);
    void toString();
    void toTextFile( std::filesystem::path filename);
    void quadratic_phase_lag_shift( double k, double c);
    void quadratic_phase_lag_shift_fourier( double k, double c);
    void mask( double radius);
    Eigen::MatrixXcd diskmask_it( Eigen::MatrixXcd to_mask, double radius);
    void propagate( double dist);
    Eigen::MatrixXcd  fft();
    int illuminate_thermally( double coherence_diameter);
    void illuminate_thermally_flat_knl( double coherence_diameter);
    //void illuminate_thermally_fft( double coherence_diameter);
    void triple_slit_mask( int w_ratio, int h_ratio, int slits, double w_offset_ratio);
    double max_intensity();
    double max_intensity( Eigen::MatrixXcd);
  private:
    void illuminate_uniformly();
    void store( std::filesystem::path filename, Eigen::MatrixXd data_to_store, double max_value, double centered_crop_factor, int bit_depth);
    double lambda;
    double L;
    Eigen::MatrixXcd value;
    Eigen::MatrixXd x, y;
    //class SignalImpl *pSignal;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
};

