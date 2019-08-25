#include <string>
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
  
