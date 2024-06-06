#include <string>
#include <map>
#include <ostream>

class ConfigLoader{
  public:
    ConfigLoader();
    std::string get( std::string param);
    int get_int( std::string param);
    double get_double( std::string param);
    friend auto operator<<( std::ostream& os, ConfigLoader const &cfg) -> std::ostream&;
  private:
    bool loadconfig();
    std::map<std::string, std::string> kv;
};


