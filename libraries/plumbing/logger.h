#ifndef LOGGER_H_
#define LOGGER_H_

#include <iostream>
#include <fstream>

/// Define the logger class here.
class logger_class {
  /// The output stream logs go to
  std::ostream output;
  /// An output file stream, for opening and closing files
  std::ofstream output_file;
  /// Default logging level is 1, only print only at main function level
  int logging_level = 1; 
  /// Current logging level
  int current_level = 0;

  public:
    /// Constructor without parameters, set the logger to std::out
    logger_class() : output(NULL) {
      output.rdbuf( std::cout.rdbuf() );
    }

    /// Construct with filename
    logger_class(std::string filename) : output(NULL) {
      set_output(filename);
    }

    /// Construct with stream
    logger_class(std::ostream stream) : output(NULL) {
      output.rdbuf( stream.rdbuf() );
    }

    /// Close the file in the destructor
    ~logger_class() {
      if(output_file.is_open()){
        output_file.close();
      }
    }

    /// Set logging level
    void set_verbosity(int level){
      logging_level = level;
    }

    /// Increase logging level. Should be called when entering a specialized
    /// area of code.
    void increase_level(){
      current_level += 1;
    }

    /// Decrease logging level. Should be called when exiting a specialized
    /// area of code.
    void decrease_level(){
      current_level += 1;
    }

    /// Set the output stream
    void set_output(std::ostream stream){
      output.rdbuf( stream.rdbuf() );
      if(output_file.is_open()){
        output_file.close();
      }
    }

    /// Open stream to a file. Interprets strings as filenames
    void set_output(std::string filename){
      if(output_file.is_open()){
        output_file.close();
      }
      output_file.open(filename, std::ios::out);
      output.rdbuf( output_file.rdbuf() );
    }

    /// Log function
    void log(std::string text){
      output << text;
    }

    template<typename T>
    friend logger_class& operator<<(logger_class& logger, T rhs);
};

/// Streaming operator for the logger class. Only node 0 prints.
/// This should also accept formatting correctly.
template<typename T>
logger_class& operator<<(logger_class& logger, T rhs)
{
    if(hila::myrank() == 0){
      logger.output << rhs;
    }
    return logger;
}

#endif
