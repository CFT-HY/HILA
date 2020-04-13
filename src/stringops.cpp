#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include "stringops.h"

/// this routine changes the input to alphanumeric + _, for naming purposes
std::string clean_name(const std::string & s) {
  
  std::string r = s;
  size_t j=0;
  for (size_t i=0; i<s.length(); i++,j++) {
    char c = s[i];
    if (std::isalnum(c) || c == '_') r[j] = c;
    else {
      switch(c) {
        case '+': r[j] = 'P'; break;
        case '-': r[j] = 'M'; break;
        case '*': r[j] = 'X'; break;
        case '/': r[j] = 'D'; break;
        case '=': r[j] = 'E'; break;
        default: r[j] = '_' ;
      }
    }
  }
  return r;
}


std::string remove_initial_whitespace(const std::string & line) {
  // clear whitespace at the beginning of string
  size_t j = 0;
  for (char p : line) {
    if (!(p == ' ' || p == '\t')) break;
    j++;
  }
  if (j > 0) return line.substr(j,std::string::npos);
  return line;
}

std::string remove_all_whitespace(const std::string & line) {
  std::string out = line; // init string 
  int j=0;
  for (char p : line) {
    if (!std::isspace(p)) out[j++] = p;
  }
  out.resize(j);
  return out;
}

/// True if string contains word (note: this word is c++ alphanumeric word, ie. split as in )

int find_word(const std::string & in, const std::string & pattern, int pos) {
  int i = in.find(pattern,pos);
  if (i == std::string::npos) return std::string::npos;  // not found

  if ( i>0 && std::isalnum(in[i-1]) ) 
    return std::string::npos;   // is at the end of a longer word
  if ( i<in.length()-pattern.length()-1 && std::isalnum(in[i+pattern.length()]) ) 
    return std::string::npos;

  return i;
}



// returns true if line contains the word list at the beginning
bool contains_word_list(const std::string & line, const std::vector<std::string> & list) {
  const char *p = line.c_str();
  for (const std::string & r : list) {
    while (isspace(*p)) p++;
    const char *rp = r.c_str();
    for (int i=0; i<r.length(); i++) {
      if (*rp != *p) return false;
      ++rp;
      ++p;
    }
  }
  return true;
}

// remove extra whitespace chars
std::string remove_extra_whitespace(const std::string & line) {
  std::string out = line; // init string 
  int j=0;
  bool previous_char_space = true;  // guarantees leading space removed
  for (char p : line) {

    if (!std::isspace(p)) {
      out[j++] = p;
      previous_char_space = false;
    } else {
      if (!previous_char_space) out[j++] = ' ';   // substitute w. real ' '
      previous_char_space = true;
    }
  }
  if (j>0 && out[j-1] == ' ') j--;  // remove trailing space
  out.resize(j);
  return out;
}   

std::string indent_string(const std::string & s) {
  
  std::string indentstr = "  ";
  std::string res = "";
  std::string line;
  
  int lev = 0;
  size_t current = 0, i = 0;
  while (i < s.length()) {
    i = s.find('\n',current);
    if (i > s.length()) {
      // no more \n, print rest
      line = s.substr(current,s.length()-current);
    } else {
      line = s.substr(current,i-current+1);
    }
    line = remove_initial_whitespace(line);
    // do the actual indent
    for (char p : line) if (p == '}') lev--;
    for (int j=0; j<lev; j++) res += indentstr;
    for (char p : line) if (p == '{') lev++;

    res += line;
    current = i+1;
    
  }
  return res;
}


std::string comment_string(const std::string & s) {

  const std::string comment("//--  ");
  std::string res = comment + s;
  
  size_t i = 0;
  while (i < res.length()) {
    i = res.find('\n',i);
    if (i < res.length()) {
      i++;
      res.insert(i,comment);
    }
  }
  return res;
}


/// From parity + dir -expression remove X, i.e.
/// X + dir -> dir
/// X - dir -> -dir
std::string remove_X(const std::string & s, bool * was_there) {
  std::string r = remove_extra_whitespace(s);
  if (r.size() == 0 || r[0] != 'X') {
    if (was_there != nullptr) *was_there = false;
    return r;
  }
  if (was_there != nullptr) *was_there = true;
  int i = 1;
  if (i < r.size() && std::isspace(r[i])) i++;
  if (i < r.size() && s[i] == '+') {
    i++;
    if (i < r.size() && std::isspace(r[i])) i++;
  }
  return r.substr(i,std::string::npos );
}

/// Types ofen seem to have "class name" -names, harmful
std::string remove_class_from_type(const std::string & s) {
  size_t i = s.find("class ",0);
  if (i < s.size() && ( i == 0 || !std::isalnum(s[i-1]))) {
    return s.substr(0,i) + s.substr(i+6,std::string::npos);
  } else    
    return s;
}
  


// Check if the cmdline has -I<include> or -D<define> -
// arguments and move these after -- if that exists on the command line.
// Clang's optionparser expects these "generic compiler and linker"
// args to be after --
// return value new argc
int rearrange_cmdline(int argc, const char **argv, const char **av) {

  bool found_ddash = false;
  av[argc+1] = nullptr;      // I read somewhere that in c++ argv[argc] = 0
  static char s[3] = "--";   // needs to be static because ptrs
  int ddashloc = 0;

  for (int i=0; i<argc; i++) {
    av[i] = argv[i];
    if (strcmp(av[i],s) == 0) {
      found_ddash = true;
      ddashloc = i;
    }
  }
  if (!found_ddash) {
    // add ddash, does not hurt in any case
    av[argc] = s;
    ddashloc = argc;
    argc++;
  }

  // now find -I and -D -options and move them after --
  for (int i=0; i<ddashloc; ) {
    if (i < ddashloc-1 && (strcmp(av[i],"-D") == 0 || strcmp(av[i],"-I") == 0)) {
      // type -D define
      const char * a1 = av[i];
      const char * a2 = av[i+1];
      for (int j=i+2; j<argc; j++) av[j-2] = av[j];
      av[argc-2] = a1;
      av[argc-1] = a2;
      ddashloc -= 2;
    } else if (strncmp(av[i],"-D",2) == 0 || strncmp(av[i],"-I",2) == 0) {
      // type -Ddefine
      const char * a1 = av[i];
      for (int j=i+1; j<argc; j++) av[j-1] = av[j];
      av[argc-1] = a1;
      ddashloc--;
    } else {
      i++;
    }
  }

  return argc;
}
