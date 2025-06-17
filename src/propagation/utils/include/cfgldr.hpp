#include <string>
#include <map>
#include <ostream>
#include <filesystem>

class ConfigLoader{
  public:
    ConfigLoader();
    ConfigLoader(const std::string& filename);
    std::string get( std::string param);
    std::filesystem::path get_dirname( std::string param);
    std::filesystem::path get_filename( std::string param);
    std::filesystem::path get_existing_dirname(const std::string& param);
    std::filesystem::path get_existing_file(const std::string& param, const std::string& basedir_param);
    bool get_bool( std::string param);
    int get_int( std::string param);
    double get_double( std::string param);
    std::string get_experiment();
    std::string get_beam_shape();
    std::string get_pixel_format();
    friend auto operator<<( std::ostream& os, ConfigLoader const &cfg) -> std::ostream&;
  private:
    //bool loadconfig();
    bool loadconfig(const std::string& filename = "cpi.cfg");
    std::map<std::string, std::string> kv;
};

class ConfigError: public std::runtime_error
{
  public:
    ConfigError( char const* message) throw();
    ConfigError( const std::string& message) throw();
};

