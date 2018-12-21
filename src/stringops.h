// -*- mode: c++ -*-
#ifndef STRINGOPS_H
#define STRINGOPS_H

#include <string>

/// Convert the string so that it is a valid variable name in c++
std::string clean_name(const std::string & s);

/// Remove whitespace at the beginning of the string (meant as a single line
std::string remove_whitespace(const std::string & line);

/// Indent the (multiline) string according to the levels of { in it
std::string indent_string(const std::string & s);

/// Convert string to a valid c++ comment.  Indentation is not changed
std::string comment_string(const std::string & s);

#endif
