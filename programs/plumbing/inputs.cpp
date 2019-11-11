#include<iostream>
#include<fstream>
#include<regex>
#include<type_traits>
#include "inputs.h"

void input::define_essentials(){
    #ifdef NDIM
    #if NDIM >= 4
    add_essential("nt");
    #endif 
    #if NDIM >= 3 
    add_essential("nz");
    #endif 
    #if NDIM >= 2
    add_essential("ny");
    #endif
    #endif
    add_essential("nx");
}

///handle one line of input in parameter file 
void input::handle(const std::string & line){
    std::regex pattern("\\s*([a-zA-Z_-]+[0-9]*)\\s*=\\s*([^\\s]*)\\s*");
    std::smatch results;
    if(!std::regex_match(line, results, pattern)){
        std::cout << "badly formatted line: " + line + "\n";
    };
    std::string variable(results[1]);
    std::string value(results[2]);
    bool is_numeric = (!value.empty() && value.find_first_not_of("0123456789.-") == std::string::npos);
    if (essentials.find(variable)!=essentials.end()) essentials[variable] = true;
    if (is_numeric) {
        if (values.count(variable)==1){
            values[variable] = std::stod(value); 
        } else {
            values.insert(std::pair<std::string, double>(variable, std::stod(value)));
        }
    } else {
        if (names.count(variable)==1){
            names[variable] = value; 
        } else {
            names.insert(std::pair<std::string, std::string>(variable, value));
        }
    }
}

input::input(const std::string & fname) {
    define_essentials(); 
    std::ifstream inputfile;
    inputfile.open(fname);
    if (inputfile.is_open()){
        std::string line;
        getline(inputfile, line, '\n');
        while(!inputfile.eof()){
            handle(line);
            getline(inputfile, line, '\n');
        }
    } else {
        std::cout << "input file couldn't be opened. Checking default params...\n";
    }
    inputfile.close();
    check_essentials();
}

void input::add_essential(const std::string & var) {
    essentials.insert(std::pair<std::string, bool>(var, false));
}

template<typename T>
void input::add_essential(const std::string & var, const T & default_value) {
    bool paramok = true;
    switch (typeid(T))
    {
    case typeid(std::string):
        names[var] = default_value;
        break;

    case typeid(int):
        values[var] = default_value;
        break;

    case typeid(float):
        values[var] = default_value;
        break;

    case typeid(double):
        values[var] = default_value;
        break;

    default:
        std::cout << "type of " + var + " not recognized (try int, float, double or string)";
        paramok = false;
        break;
    }
    essentials.insert(std::pair<std::string, bool>(var, paramok));
}

void input::check_essentials(){
    bool fail = false;
    for (auto i = essentials.begin(); i != essentials.end(); ++i){
        if (!(*i).second){
            std::cout << "required parameter " + (*i).first + " not found\n"; 
        }
        fail = true;
    }
    if (fail){
        std::cout << "exiting...";
        exit(1);
    }
}

input::returntype input::get(const std::string & variable){
    return returntype(variable, this);
}

