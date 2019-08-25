#include <string>
#include <cstring>
#include "stringops.h"



std::string clean_name(const std::string & s) {
  
  std::string r = s;
  size_t j=0;
  for (size_t i=0; i<s.length(); i++) {
    char c = s[i];
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
        c == '_' ||(i > 0 && c >= '0' && c <= '9'))
      {
        r[j++] = c;
      }
  }
  r.resize(j);
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
  


// Check if the cmdline has -I<include> or -D<define> -
// arguments and move these after -- if that exists on the command line.
// Clang's optionparser expects these "generic compiler and linker"
// args to be after --
// return value new argc
int rearrange_cmdline(int argc, const char **argv, const char **av) {

  bool found_ddash = false;
  av[argc+1] = nullptr;  // I read somewhere that in c++ argv[argc] = 0
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
  for (int i=0; i<ddashloc; i++) {
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
    }      
  }

  return argc;
}
