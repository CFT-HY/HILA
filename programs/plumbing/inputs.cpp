#include "inputs.h"
#include<typeinfo>
#include<iostream>
#include<fstream>
#include<regex>

void inputs::define_essentials(){
    #ifdef USE_MPI
    //call add_essential(variable) here to make variable necessary for an mpi run 
    #endif
    #ifdef CUDA
    //same thing here, for variables that need to be in params file for cuda runs etc.
    #endif
}

///handle one line of input in parameter file 
void inputs::handle(std::string & line){
    std::regex pattern("\\s*([a-zA-Z_-]+[0-9]*)\\s*=\\s*([^\\s]*)\\s*");
    std::smatch results;
    if(!std::regex_match(line, results, pattern)) return;
    std::string variable(results[1]);
    std::string right(results[2]);
    bool is_numeric = (!right.empty() && right.find_first_not_of("0123456789.-") == std::string::npos);
    if (essentials.find(variable)!=essentials.end()) essentials[variable] = true;
    if (is_numeric) {
        if (values.find(variable)!=values.end()){
            values[variable] = std::stod(right); 
        } else {
            values.insert(std::pair<std::string, double>(variable, std::stod(right)));
        }
    } else {
        if (names.find(variable)!=names.end()){
            names[variable] = right; 
        } else {
            names.insert(std::pair<std::string, std::string>(variable, right));
        }
    }
}

inputs::inputs(std::string & fname) {
    construct_essentials(); 
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
