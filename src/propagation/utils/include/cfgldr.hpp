#include <string>
#include <map>
#include <ostream>
#include <filesystem>

class ConfigLoader{
  public:
    ConfigLoader();
    std::string get( std::string param);
    std::filesystem::path get_pathname( std::string param);
    int get_int( std::string param);
    double get_double( std::string param);
    std::string get_experiment();
    std::string get_pixel_format();
    friend auto operator<<( std::ostream& os, ConfigLoader const &cfg) -> std::ostream&;
  private:
    bool loadconfig();
    std::map<std::string, std::string> kv;
};

class ConfigError: public std::runtime_error
{
  public:
    ConfigError( char const* message) throw();
    ConfigError( const std::string& message) throw();
};

