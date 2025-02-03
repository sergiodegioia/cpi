#include <cfgldr.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numbers>
#include <cmath>
#include <map>
#include <filesystem>

#include <algorithm> 
#include <cctype>
#include <locale>
// trim from start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

inline void trim(std::string &s) {
  rtrim( s);
  ltrim( s);
}


auto operator<<( std::ostream& os, ConfigLoader const &cfg) -> std::ostream&{
  for( const auto& elem: cfg.kv){
    os << "{" << elem.first << "} = {" << elem.second << "}\n";
  }
  return os;
}

ConfigLoader::ConfigLoader(): kv{std::map<std::string, std::string>()}{
  loadconfig();
}

std::string ConfigLoader::get_pixel_format(){
  try{
    std::string pixel_format = get( "tiff.pixel_format");
    if( pixel_format != "8U" && pixel_format != "16U" && pixel_format != "64F"){
      throw ConfigError( "[utils/ConfigLoader] Parameter {pixel_format} is {" + pixel_format + "} and hence is neither 8U nor 16U nor 64F, the only legal values.");
    }
    return pixel_format;
  }catch( ConfigError &cfgerr){
    std::cout << cfgerr.what() << std::endl;
    throw;
  }catch( std::exception &err){
    std::cout << "[utils/ConfigLoader] Parameter {pixel_format} cannot be parsed as a double. Check it in the configuration file for the wrong format and fix it. Previous exception: " << err.what() << std::endl;
    throw;
  }
}

std::string ConfigLoader::get_experiment(){
  try{
    std::string experiment = get( "experiment");
    if( experiment != "CPI" && experiment != "GI" && experiment != "OBJ" && experiment != "NO_OP"){
      throw ConfigError( "[utils/ConfigLoader] Parameter {experiment} is {" + experiment + "} and hence is neither CPI nor GI, the only legal values.");
    }
    return experiment;
  }catch( ConfigError &cfgerr){
    std::cout << cfgerr.what() << std::endl;
    throw;
  }catch( std::exception &err){
    std::cout << "[utils/ConfigLoader] Parameter {experiment} cannot be parsed as a double. Check it in the configuration file for the wrong format and fix it. Previous exception: " << err.what() << std::endl;
    throw;
  }
}

double ConfigLoader::get_double( std::string param){
  try{
    return std::stod( get( param));
  }catch( ConfigError &cfgerr){
    std::cout << cfgerr.what() << std::endl;
    throw;
  }catch( std::exception &err){
    std::cout << "[utils/ConfigLoader] Parameter {" << param << "} cannot be parsed as a double. Check it in the configuration file for the wrong format and fix it. Previous exception: " << err.what() << std::endl;
    throw;
  }
}

std::filesystem::path ConfigLoader::get_pathname( std::string param){
  auto value = get( param);
  try{
    auto path = std::filesystem::path( value);
    if( get_bool( "storage.overwrite")){
      std::filesystem::create_directory( path);
    }else{
      if( !std::filesystem::create_directory( path)){
        throw ConfigError( "[utils/ConfigLoader] Parameter {"  + param +  "} is {" + value + "} which is an existing path. It cannot be used in non-overwrite mode.");
      }
    }
    return path;
  }catch( std::exception &err){
    std::cout << err.what() << std::endl;
    throw;
  }
}

bool ConfigLoader::get_bool(std::string param) {
    try {
        std::string value = get( param);
        std::transform(value.begin(), value.end(), value.begin(), ::tolower); // Convert to lowercase for case-insensitivity
        
        if ( value == "y" || value == "true" || value == "1") {
            return true;
        } else if ( value == "n" || value == "false" || value == "0") {
            return false;
        } else {
            throw ConfigError( "[utils/ConfigLoader] Parameter {" + param + "} is {" + value + "} and hence is neither Y nor N, the only legal values. Check it in the configuration file for the wrong format and fix it. ");
        }
    } catch (ConfigError &cfgerr) {
        std::cout << cfgerr.what() << std::endl;
        throw;
    }
}

int ConfigLoader::get_int( std::string param){
  try{
    return std::stoi( get( param));
  }catch( ConfigError &cfgerr){
    std::cout << cfgerr.what() << std::endl;
    throw;
  }catch( std::exception &err){
    std::cout << "[utils/ConfigLoader] Parameter {" << param << "} cannot be parsed as an int. Check it in the configuration file for the wrong format and fix it. Previous exception: " << err.what() << std::endl;
    throw;
  }
}

std::string ConfigLoader::get( std::string param){
  auto iter = kv.find( param);
  if( iter == kv.end()){
    throw ConfigError( "[utils/ConfigLoader] Parameter {" + param + "} does not exist in the configuration file.");
  }else{
    return iter->second;
  }
}

bool ConfigLoader::loadconfig(){
  std::string filename("cpi.cfg");
  std::ifstream ifs( filename);
  if( !ifs.is_open()){
    perror( ("error while opening file " + filename).c_str());
    return false;
  }
  std::string line;
  while( std::getline( ifs, line)){
    trim( line);
    if( line.empty()){
      continue;
    }
    if( '#' == line[ 0]){
      continue;
    }
    std::istringstream is_line( line);
    std::string key;
    std::getline( is_line, key, '=');
    trim( key);
    if( is_line.eof()){
      std::cout << "[utils/ConfigLoader] While reading config file it was found a line [" << line << "] without an equal sign. This line has been skipped." << std::endl;
      continue;
    }
    std::string value;
    std::getline( is_line, value, '#');
    trim( value);
    if( value.empty()){
      std::cout << "[utils/ConfigLoader] While reading config file it was found parameter [" << key << "] with no value assigned. This parameter has been skipped." << std::endl;
      continue;
    }
    const auto [it, success] = kv.insert({ key, value});
    if( !success){
      std::cout << "[utils/ConfigLoader] While reading config file it was found a duplicate key [" << key << "]. Only the value of the first occurrence of that key has been loaded." << std::endl;
    }
  }
  if( ifs.bad()){
    perror( ("error while reading file " + filename).c_str());
    return false;
  }
  ifs.close();
  return true;
}

ConfigError::ConfigError( const std::string& message) throw()
  : ConfigError( message.c_str())
{
}

ConfigError::ConfigError( char const* message) throw()
  : std::runtime_error( message)
{
}

