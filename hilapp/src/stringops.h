// -*- mode: c++ -*-
#ifndef STRINGOPS_H
#define STRINGOPS_H

#include <string>

/// Return char * to git sha value, if defined
std::string git_sha_value();

/// Convert the string so that it is a valid variable name in c++
std::string clean_name(const std::string &s);

/// Remove whitespace at the beginning of the string (meant as a single line
std::string remove_initial_whitespace(const std::string &line);

/// Remove extra whitespace from a string, string will contain words separated by space
std::string remove_extra_whitespace(const std::string &line);

/// Remove all whitespace, including \n
std::string remove_all_whitespace(const std::string &line);

/// if in contains the pattern returns its first location.
/// pattern must be delimited by non-alphanumeric chars (c++ symbol name rules)
std::string::size_type find_word(const std::string &in, const std::string &pattern,
                                 int pos = 0, bool reverse = false);

/// Check whether line contains the list of space-separated words in list in this order.
/// Extra whitespaces skipped if remainder != nullptr and return value == true, the string
/// pointed by remainder contains the rest of the string.
bool contains_word_list(const std::string &line, const std::string &list,
                        std::string *remainder = nullptr);

/// Indent the (multiline) string according to the levels of { in it
std::string indent_string(const std::string &s);

/// Convert string to a valid c++ comment.  Indentation is not changed
std::string comment_string(const std::string &s);

/// remove initial X from X+dir string.  Optional is_X gives true if X was there
std::string remove_X(const std::string &s, bool *is_X = nullptr);

/// Remove "class " keyword from type which sometimes pops there
std::string remove_class_from_type(const std::string &s);

/// Rearranges command line args according to the wishes of clang tool
int rearrange_cmdline(int argc, const char **argv, const char **av);

#endif
