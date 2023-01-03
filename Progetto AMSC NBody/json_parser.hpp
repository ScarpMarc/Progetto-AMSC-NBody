#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Constants.h"
// TODO
/*#include "../json/single_include/nlohmann/json.hpp"
using json = nlohmann::json;

// ...

std::ifstream f("example.json");
json data = json::parse(f);*/


class JsonParser
{

    public: 
    // Constructor

    JsonParser(std::string _file_name):
    filename(_file_name)
    {}

    void parse()
    {
        if (filename.empty())
            filename = "settings.json";
        std::string s;
        s = _json_string_from_file(filename);
        _parse_json(s);
    }

    private:

    std::string _json_string_from_file(std::string filepath)
    {
        std::ifstream myfile(filepath);
        std::ostringstream temp;
        while (!myfile.eof())
        {
            std::string _temp;
            getline(myfile, _temp);
            temp << _temp;
        }
        //temp << myfile.rdbuf();
        std::string s = temp.str();
        //std::cout << filepath << std::endl;
        //std::cout << s << std::endl;
        myfile.close();

        return s;
    }

    std::map<std::string, std::string> _parse_json(const std::string& json) {

        std::map<std::string, std::string> values;
        
        // Find the start and end of the object
        size_t object_start = json.find('{');
        size_t object_end = json.rfind('}');
        if (object_start == std::string::npos || object_end == std::string::npos) {
            // Invalid JSON
            return values;
        }

        // Extract the object string
        std::string object_str = json.substr(object_start + 1, object_end - object_start - 1);

        
        // Split the object string into key-value pairs
        size_t key_start = 0, flag = 1;
        while (key_start < object_str.size() && flag != 0) {
            // std::cout << key_start << std::endl;
            // Find the start and end of the key
            size_t key_end = object_str.find(':', key_start);
            if (key_end == std::string::npos) {
            // Invalid JSON
            return values;
            }
            std::string key = object_str.substr(key_start, key_end - key_start);
            //std::cout << key << std::endl;

            // Find the start and end of the value
            size_t value_start = object_str.find('"', key_end + 1);
            if (value_start == std::string::npos) {
            // Invalid JSON
            return values;
            }
            size_t value_end = object_str.find('"', value_start + 1);
            if (value_end == std::string::npos) {
            // Invalid JSON
            return values;
            }
            std::string value = object_str.substr(value_start + 1, value_end - value_start - 1);
            //std::cout << value << std::endl;
            
            if(key.find("max_ticks") != std::string::npos)
            {
                max_ticks = std::stoull(value);
                //std::cout << max_ticks << std::endl;
            }
            if(key.find("ticks_per_second") != std::string::npos)
            {
                ticks_per_second = atof(value.c_str());
                //std::cout << ticks_per_second << std::endl;
            }
            if(key.find("screen_refresh_millis") != std::string::npos)
            {
                screen_refresh_millis = stoul(value);
                //std::cout << screen_refresh_millis << std::endl;
            }
            if(key.find("screenResX") != std::string::npos)
            {
                screenResX = stoul(value);
                //std::cout << screenResX << std::endl;
            }
            if(key.find("screenResY") != std::string::npos)
            {
                screenResY = stoul(value);
                //std::cout << screenResY << std::endl;
            }
            if(key.find("save_status_interval") != std::string::npos)
            {
                save_status_interval = stoul(value);
                //std::cout << save_status_inteval << std::endl;
            }
            if(key.find("total_particles") != std::string::npos)
            {
                total_particles = stoul(value);
                //std::cout << total_particles << std::endl;
            }

            // Store the key-value pair in the map
            values[key] = value;

            // Move to the next key
            key_start = object_str.find(',', value_end) + 1;
            if(key_start == 0){
                flag = 0;
            }
            
        }

        return values;
    }

    //private:
    std::string filename;

};


/*
only of the form
{
  "key1": "value1",
  "key2": "value2",
  ...
}
*/
