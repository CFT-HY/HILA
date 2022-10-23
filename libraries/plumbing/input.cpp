#include <sstream>
#include <iostream>
#include <fstream>
#include <regex>
#include <type_traits>
#include "defs.h"
#include "input.h"
#include <errno.h>
#include <iomanip>

#define COMMENT_CHAR '#'

//////////////////////////////////////////////////////////////////////////////
/// Parameter file input system
/// Check input.h for user instructions
//////////////////////////////////////////////////////////////////////////////

namespace hila {

static std::string empty_key("");

bool input::open(const std::string &file_name, bool use_cmdline, bool exit_on_error) {
    extern const char *input_file;

    std::string fname;

    if (use_cmdline && input_file != nullptr) {

        fname = input_file;
        input_file = nullptr; // set this null, so that it is not used again
    } else {
        fname = file_name;
    }

    bool got_error = false;
    if (hila::myrank() == 0) {
        if (is_initialized) {
            if (speaking)
                hila::out << "Error: file '" << fname
                             << "' cannot be opened because '" << filename
                             << "' is open in this input variable\n";

            got_error = true;
        } else {
            filename = fname;
            if (fname == "-") {
                use_cin = true;
                if (speaking)
                    print_dashed_line("Reading from standard input");
            } else {
                use_cin = false;
                inputfile.open(fname);
                if (inputfile.is_open()) {
                    is_initialized = true;

                    if (speaking)
                        print_dashed_line("Reading file " + filename);

                } else {

                    if (speaking)
                        hila::out << "Error: input file '" << fname
                                     << "' could not be opened\n";

                    got_error = true;
                }
            }
        }
    }
    hila::broadcast(got_error);
    if (got_error && exit_on_error) {
        hila::finishrun();
    }

    return !got_error;
}

void input::close() {
    if (is_initialized && !use_cin) {
        inputfile.close();
        is_initialized = false;
    }
    if (speaking)
        print_dashed_line();

    // automatic cleaning of other vars
}

// read one line skipping comments and initial whitespace
bool input::get_line() {
    if (hila::myrank() == 0) {
        do {
            std::istream *inptr;

            if (use_cin)
                inptr = &std::cin;
            else
                inptr = &inputfile;

            *inptr >> std::ws; // remove initial whitespace
            if (!std::getline(*inptr, linebuffer)) {
                linebuffer.clear();
                lb_start = 0;
                return false;
            }
        } while (linebuffer.at(0) == COMMENT_CHAR);
        size_t i = linebuffer.find(COMMENT_CHAR);
        if (i != std::string::npos)
            linebuffer.resize(i);
        lb_start = 0;

        is_line_printed = false; // not yet printed
    }
    return true; // sync this at another spot
}

// print the read-in line with a bit special formatting
void input::print_linebuf(int end_of_key) {
    if (hila::myrank() == 0 && speaking) {

        if (is_line_printed)
            return;
        is_line_printed = true;

        int i = 0;
        while (std::isspace(linebuffer[i]))
            i++;
        for (; i < linebuffer.size() && i < end_of_key; i++) {
            hila::out << linebuffer[i];
        }
        hila::out << ' ';
        if (end_of_key > 0) {
            for (int j = i; j < 20; j++)
                hila::out << ' ';
        }

        while (i < linebuffer.size() && std::isspace(linebuffer[i]))
            i++;
        if (i < linebuffer.size()) {
            hila::out << linebuffer.substr(i);
        }
        hila::out << std::endl; // use endl to show output here
    }
}

// remove leading whitespace, incl. lines
bool input::remove_whitespace() {
    if (hila::myrank() == 0) {
        while (lb_start < linebuffer.size() && std::isspace(linebuffer[lb_start]))
            lb_start++;
        if (lb_start == linebuffer.size())
            return get_line();
    }
    return true; // do not broadcast yet, should be used only in node 0
}

// returns true if line contains the word list at the beginning of line.  list
// contains the word separated by whitespace.  If match found, advances the
// lb_start to new position.  end_of_key is the index where key match on line ends

bool input::contains_word_list(const std::string &list, int &end_of_key) {
    const char *p = linebuffer.c_str() + lb_start;
    const char *q = list.c_str();
    while (std::isspace(*p))
        p++;
    while (std::isspace(*q))
        q++;

    while (*p && *q) {
        // compare non-space chars
        while (*p && *q && *p == *q && !std::isspace(*q)) {
            p++;
            q++;
        }
        if (std::isspace(*q) && std::isspace(*p)) {
            // matching spaces, skip
            while (std::isspace(*p))
                p++;
            while (std::isspace(*q))
                q++;
        }

        if (*p != *q)
            break;
    }
    // if line contained the words in list, *q = 0.
    while (std::isspace(*q))
        q++;
    if (*q != 0)
        return false;

    end_of_key = p - linebuffer.c_str();

    while (std::isspace(*p))
        p++;
    lb_start = p - linebuffer.c_str();
    return true;
}

// find tokens separated by whitespace or ,
// inside of " .. " is one token too.
// returns true if the next token (on the same line) is found

bool input::peek_token(std::string &tok) {
    if (!remove_whitespace())
        return false;
    size_t i;
    bool in_quotes = false;
    for (i = lb_start;
         i < linebuffer.size() &&
         ((!std::isspace(linebuffer[i]) && linebuffer[i] != ',') || in_quotes);
         i++) {
        if (linebuffer[i] == '"') {
            in_quotes = !in_quotes;
        }
    }
    if (i == lb_start && linebuffer[i] == ',')
        i++; // this happens for only ,

    if (in_quotes) {
        hila::out << "Error: unbalanced quotes\n";
        exit(1); // unclean exit
    }
    tok = linebuffer.substr(lb_start, i - lb_start);
    return true;
}

// consumes the token too

bool input::get_token(std::string &tok) {
    if (peek_token(tok)) {
        lb_start += tok.size();
        return true;
    }
    return false;
}

// match_token returns true and consumes the token if it matches the argument

bool input::match_token(const std::string &tok) {
    std::string s;
    if (peek_token(s) && s == tok) {
        lb_start += tok.size();
        return true;
    }
    return false;
}

// require the (typically beginning of line) key for parameters

bool input::handle_key(const std::string &key) {
    if (hila::myrank() == 0) {
        // check the linebuffer for stuff
        remove_whitespace();

        int end_of_key = 0;
        if (key.size() > 0 && !contains_word_list(key, end_of_key)) {
            if (speaking) {
                print_linebuf(0);
                hila::out << "Error: expecting key '" << key << "'\n";
            }
            return false;
        }

        print_linebuf(end_of_key);
    }
    return true;
}

// is the input string int/double/string and return it


// a trivial function, useful for template
bool input::is_value(const std::string &str, std::string &val) {
    val = remove_quotes(str);
    return true;
}

std::string input::remove_quotes(const std::string &val) {
    size_t i, j;
    std::string res;
    res = val;
    for (j = i = 0; i < val.size(); i++)
        if (val[i] != '"')
            res[j++] = val[i];
    res.resize(j);
    return res;
}

//  expects "label   <item>"  -line, where <item> matches one of the std::strings in
//  items. returns the index of the item. If not found, errors out

int input::get_item(const std::string &label, const std::vector<std::string> &items,
                    bool bcast) {

    bool no_error = handle_key(label);
    int item = -1;
    double d;
    std::string s;

    if (hila::myrank() == 0) {
        if (no_error && peek_token(s)) {

            for (int i = 0; i < items.size() && item < 0; i++) {
                double dval;
                long lval;
                int end_of_key;

                // clang-format off
                if ((items[i] == "%s") ||
                    (items[i] == "%f" && is_value(s,dval)) ||
                    (items[i] == "%i" && is_value(s,lval)) ||
                    contains_word_list(items[i],end_of_key)) {
                    item = i;
                }
                // clang-format on
            }
        }

        if (item < 0) {
            // error, nothing was found
            no_error = false;
            if (speaking) {
                hila::out << "Input '" << label << "' must be one of: ";
                for (int i = 0; i < items.size(); i++) {
                    if (items[i] == "%s")
                        hila::out << "<string> ";
                    else if (items[i] == "%f")
                        hila::out << "<float/double> ";
                    else if (items[i] == "%i")
                        hila::out << "<int/long> ";
                    else
                        hila::out << '\'' << items[i] << "' ";
                }
                hila::out << '\n';
            }
        }
    }

    if (bcast) {
        hila::broadcast2(item, no_error);

        // with broadcast exit on error
        if (!no_error)
            hila::finishrun();
    }

    return item;
}

} // namespace hila
