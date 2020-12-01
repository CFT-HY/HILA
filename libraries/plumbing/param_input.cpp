#include <sstream>
#include<iostream>
#include<fstream>
#include<regex>
#include<type_traits>
#include "defs.h"
#include "param_input.h"
#include "globals.h"
#include <errno.h>
#include <iomanip>

#define COMMENT_CHAR '#'

//////////////////////////////////////////////////////////////////////////////
/// Parameter file input system
/// Check param_input.h for user instructions
//////////////////////////////////////////////////////////////////////////////


input::returntype input::get(const std::string & variable){
    return returntype(variable, this);
}

void input::init(const std::string &fname) {
  bool got_error;
  if (hila::my_rank == 0) {
    filename = fname;
    inputfile.open(fname);
    if (inputfile.is_open()){
      is_initialized = true;
      hila::output << "----- Reading file '" << filename << "' -----\n";

    } else {
      hila::output << "Error: input file '" << fname << "' could not be opened\n";
      got_error = true;
    }
  }
  broadcast(got_error);
  if (got_error) terminate(0);
}

void input::close() {
  if (is_initialized) {
    inputfile.close();
    is_initialized = false;
    output0 << "----------------------------------------\n";
  }
  // automatic cleaning of other vars
}

// read one line skipping comments and initial whitespace
bool input::get_line() {
  if (hila::my_rank == 0) {
    do {
      inputfile >> std::ws;  // remove initial whitespace
      if (!std::getline(inputfile, linebuffer)) return false;
    } while (linebuffer.at(0) == COMMENT_CHAR);
    size_t i = linebuffer.find(COMMENT_CHAR);
    if (i != std::string::npos) linebuffer.resize(i);
    lb_start = 0;
    print_linebuf();
  }
  return true;  // sync this at another spot
}

// print the read-in line with a bit special formatting
void input::print_linebuf() {
  if (hila::my_rank == 0) {
    std::string out;
    size_t i;
    for (i=0; i<linebuffer.size() && !std::isspace(linebuffer[i]); i++) 
      out.push_back(linebuffer[i]);
    out.push_back(' ');
    while (out.size() < 16) out.push_back(' ');
    while (i<linebuffer.size() && std::isspace(linebuffer[i])) i++;
    if (i < linebuffer.size()) out.append( linebuffer.substr(i) );
    hila::output << out << '\n';
  }
}

// find tokens separated by whitespace or , - 
// inside of " .. " is one token too.
// returns true if the next token (on the same line) is found

bool input::peek_token( std::string & tok ) {
  while (lb_start < linebuffer.size() && std::isspace(linebuffer[lb_start])) lb_start++;
  if (lb_start == linebuffer.size()) return false;
  size_t i;
  bool in_quotes = false;
  for (i=lb_start; i<linebuffer.size() && 
       ((!std::isspace(linebuffer[i]) && linebuffer[i] != ',') || in_quotes);
       i++) {
    if (linebuffer[i] == '"') {
      in_quotes = !in_quotes;
    }
  }
  if (i == lb_start && linebuffer[i] == ',') i++;  // this happens for only ,

  if (in_quotes) {
    hila::output << "Error: unbalanced quotes\n";
    exit(0);   // unclean exit
  }
  tok = linebuffer.substr(lb_start,i-lb_start);
  return true;
}

// consumes the token too

bool input::get_token( std::string & tok) {
  if (peek_token(tok)) {
    lb_start += tok.size();
    return true;
  }
  return false;
}

// match_token returns true and consumes the token if it matches the argument

bool input::match_token(const std::string & tok) {
  std::string s;
  if (peek_token(s) && s == tok) {
    lb_start += tok.size();
    return true;
  }
  return false;
}

// require the (typically beginning of line) label for parameters

bool input::handle_label(const std::string &label) {
  if (hila::my_rank==0) {
    // check the linebuffer for stuff
    std::string l;
    if (!get_token(l)) {
      get_line();
      get_token(l);
    }

    if (l != label) {
      hila::output << "Error: expecting item '" << label << "'\n";
      return false;
    }
  }
  return true;
}

// is the input string int/double/string and return it

bool input::is_value(const std::string & str, int &val) {
  char *lp;
  val = std::strtol(str.c_str(),&lp,10);
  if (str.c_str() == lp || *lp) return false;
  return true;
}

bool input::is_value(const std::string & str, double &val) {
  char *lp;
  val = std::strtod(str.c_str(),&lp);
  if (str.c_str() == lp || *lp) return false;
  return true;
}

// a trivial function, useful for template
bool input::is_value(const std::string & str, std::string & val) {
  val = remove_quotes(str);
  return true;
}

// require "label <int>"

int input::get_int(const std::string &label) {
  int val = 0;
  bool no_error = handle_label(label);

  if (hila::my_rank==0 && no_error) {
    std::string tok;
    if (!(get_token(tok) && is_value(tok,val))) {
      hila::output << "Error: expecting integer after '" << label << "'\n";
      no_error = false;
    }
  }
  broadcast(no_error);
  if (!no_error) terminate(0);

  broadcast(val);
  return (int)val;
}

// "label <double>"

double input::get_double(const std::string &label) {
  double val = 0;
  bool no_error = handle_label(label); // removes whitespace

  if (hila::my_rank==0 && no_error) {
    std::string tok;
    if (!(get_token(tok) && is_value(tok,val))) {
      hila::output << "Error: expecting double after '" << label << "'\n";
      no_error = false;
    }
  }
  broadcast(no_error);
  if (!no_error) terminate(0);

  broadcast(val);
  return val;
}

std::string input::remove_quotes(const std::string &val) {
  size_t i,j;
  std::string res;
  res = val;
  for (j=i=0; i<val.size(); i++) if (val[i] != '"') res[j++] = val[i];
  res.resize(j);
  return res;
}

// find "label <string>"

std::string input::get_string(const std::string &label) {
  std::string val;
  bool no_error = handle_label(label);

  if (hila::my_rank==0 && no_error) {
    no_error = get_token(val);
    if (no_error) {
      val = remove_quotes(val);  
    } else {
      hila::output << "Error: expecting a string value after '" << label << "'\n";
    }
  }
  broadcast(no_error);
  if (!no_error) terminate(0);

  broadcast(val);
  return val;
}

//  expects "label   <item>"  -line, where <item> matches one of the std::strings in items.
//  returns the index of the item. If not found, errors out
 
int input::get_item(const std::string &label, const std::vector<std::string> &items) {
  bool no_error = handle_label(label);
  int i=0;

  if (hila::my_rank == 0 && no_error) {
    std::string s;
    no_error = get_token(s);
    if (no_error) {
      for (i=0; i<items.size(); i++) {
        if (items[i] == "%i") { 
          if (is_value(s,item_int_val)) break;
        } else if (items[i] == "%f") { 
          if (is_value(s,item_double_val)) break;
        } else if (items[i] == "%s") {
          // %s always matches
          item_string_val = remove_quotes(s);
          break;
        } else if (items[i] == s) {
          // got the item
          break;
        }
      }
    } 
    if (!no_error || i == items.size()) {
      hila::output << "Input '" << label << "' must be one of: ";
      for (int j=0; j<items.size(); j++) {
        if (items[j] == "%f") hila::output << "(double) ";
        else if (items[j] == "%i") hila::output << "(integer) ";
        else hila::output << '\'' << items[j] << "' ";
      }
      hila::output << '\n';
      no_error = false;
    }
  }

  broadcast(no_error);
  if (!no_error) terminate(0);

  broadcast(i);
  return i;
}

// utility function to get a comma-separated list of ints/doubles

template <typename T>
void input::get_type_list(const std::string & label, std::vector<T> & res, 
                          const char * name) {
  bool no_error = handle_label(label);
  const char * rest = linebuffer.c_str() + lb_start;
  while (std::isspace(*rest)) rest++;
  res.clear();

  if (hila::my_rank==0 && no_error) {
    std::string s;
    T val;
    do {
      if (get_token(s) && is_value(s,val)) {
        res.push_back(val);
      } else {
        no_error = false;
      }
    } while (no_error && match_token(","));

    if (!no_error) {
      hila::output << "Error: expecting comma-separated list of " << name <<"s after '" << label << "'\n";
    }
  }
  broadcast(no_error);
  if (!no_error) terminate(0);

  // broadcast the vector too
  broadcast(res);
}

std::vector<int> input::get_int_list(const std::string & label) {
  std::vector<int> res;
  get_type_list(label,res,"integer");
  return res;
}

std::vector<double> input::get_double_list(const std::string & label) {
  std::vector<double> res;
  get_type_list(label,res,"double");
  return res;
}

std::vector<std::string> input::get_string_list(const std::string & label) {
  std::vector<std::string> res;
  get_type_list(label,res,"strings");
  return res;
}
