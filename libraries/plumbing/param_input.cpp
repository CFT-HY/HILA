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

static std::string empty_key("");

input::returntype input::get(const std::string & key){
    return returntype(key, this);
}

input::returntype input::get() {
  return returntype(empty_key,this);
}

void input::open(const std::string &fname) {
  bool got_error;
  if (hila::myrank() == 0) {
    if (is_initialized) {
      hila::output << "Error: file '" << fname << "' cannot be opened because '" << filename 
                   << "' is open in this input variable\n";
      got_error = true;
    } else {
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
  }
  broadcast(got_error);
  if (got_error) hila::terminate(0);
}

void input::close() {
  if (is_initialized) {
    inputfile.close();
    is_initialized = false;
    print_dashed_line();
  }
  // automatic cleaning of other vars
}

// read one line skipping comments and initial whitespace
bool input::get_line() {
  if (hila::myrank() == 0) {
    do {
      inputfile >> std::ws;  // remove initial whitespace
      if (!std::getline(inputfile, linebuffer)) {
        linebuffer.clear();
        lb_start = 0;
        return false;
      }
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
  if (hila::myrank() == 0) {
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

// remove leading whitespace, incl. lines
bool input::remove_whitespace() {
  if (hila::myrank() == 0) {
    while (lb_start < linebuffer.size() && std::isspace(linebuffer[lb_start])) lb_start++;
    if (lb_start == linebuffer.size()) return get_line();
  }
  return true;    // do not broadcast yet, should be used only in node 0
}

// find tokens separated by whitespace or , - 
// inside of " .. " is one token too.
// returns true if the next token (on the same line) is found

bool input::peek_token( std::string & tok ) {
  if (!remove_whitespace()) return false;
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

// require the (typically beginning of line) key for parameters

bool input::handle_key(const std::string &key) {
  if (hila::myrank()==0) {
    // check the linebuffer for stuff
    remove_whitespace();
    if (key.size() == 0) return true;
    std::string l;
    if (!get_token(l) || l != key) {
      hila::output << "Error: expecting key '" << key << "', got '" << l << "'\n";
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
  bool no_error = handle_key(label);

  if (hila::myrank()==0 && no_error) {
    std::string tok;
    if (!(get_token(tok) && is_value(tok,val))) {
      hila::output << "Error: expecting integer after '" << label << "'\n";
      no_error = false;
    }
  }

  // helper to do it in a single broadcast
  struct bc_helper {
    int val;
    bool noerr;
  } bcdata = {val, no_error};

  broadcast(bcdata);
  if (!bcdata.noerr) hila::terminate(0);
  return bcdata.val;
}

int input::get_int() { return get_int(empty_key); }


// "label <double>"

double input::get_double(const std::string &label) {
  double val = 0;
  bool no_error = handle_key(label); // removes whitespace

  if (hila::myrank()==0 && no_error) {
    std::string tok;
    if (!(get_token(tok) && is_value(tok,val))) {
      hila::output << "Error: expecting double after '" << label << "'\n";
      no_error = false;
    }
  }

  // helper to do it in a single broadcast
  struct bc_helper {
    double val;
    bool noerr;
  } bcdata = {val, no_error};

  broadcast(bcdata);
  if (!bcdata.noerr) hila::terminate(0);
  return bcdata.val;
}

double input::get_double() { return get_double(empty_key); }

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
  bool no_error = handle_key(label);

  if (hila::myrank()==0 && no_error) {
    no_error = get_token(val);
    if (no_error) {
      val = remove_quotes(val);  
    } else {
      hila::output << "Error: expecting a string value after '" << label << "'\n";
    }
  }
  broadcast(no_error);
  if (!no_error) hila::terminate(0);

  broadcast(val);
  return val;
}

std::string input::get_string() { return get_string(empty_key); }

//  expects "label   <item>"  -line, where <item> matches one of the std::strings in items.
//  returns the index of the item. If not found, errors out
 
int input::get_item(const std::string &label, const std::vector<std::string> &items, double *dval) {
  bool no_error = handle_key(label);
  int i=0;
  double d;

  if (hila::myrank() == 0 && no_error) {
    std::string s;
    no_error = get_token(s);
    if (no_error) {
      for (i=0; i<items.size() && items[i] != s; i++) ;
      if (i == items.size() && dval != nullptr) {
        no_error = is_value(s,d);
        i = -1;
      }
     }
    if (!no_error || i == items.size()) {
      hila::output << "Input '" << label << "' must be one of: ";
      for (int j=0; j<items.size(); j++) {
        hila::output << '\'' << items[j] << "' ";
      }
      if (dval != nullptr) {
        hila::output << " (double value)";
      }
      hila::output << '\n';
      no_error = false;
    }
  }

  struct bcast_helper {
    double d;
    int i;
    bool noerror;
  } bcdata = {d, i, no_error};     // helper struct to do everything in a single broadcast

  broadcast(bcdata);
  if (!bcdata.noerror) hila::terminate(0);
  if (bcdata.i == -1) *dval = bcdata.d;
  return bcdata.i;
}

int input::get_item(const std::string &label, const std::vector<std::string> &items, double & dval) {
  return get_item(label, items, &dval);
}


// utility function to get a comma-separated list of ints/doubles

template <typename T>
void input::get_type_vector(const std::string & label, std::vector<T> & res, 
                            const char * name) {
  bool no_error = handle_key(label);
  const char * rest = linebuffer.c_str() + lb_start;
  while (std::isspace(*rest)) rest++;
  res.clear();

  if (hila::myrank()==0 && no_error) {
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
  if (!no_error) hila::terminate(0);

  // broadcast the vector too
  broadcast(res);
}

std::vector<int> input::get_int_vector(const std::string & label) {
  std::vector<int> res;
  get_type_vector(label,res,"integer");
  return res;
}

std::vector<double> input::get_double_vector(const std::string & label) {
  std::vector<double> res;
  get_type_vector(label,res,"double");
  return res;
}

std::vector<std::string> input::get_string_vector(const std::string & label) {
  std::vector<std::string> res;
  get_type_vector(label,res,"strings");
  return res;
}
