#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include "stringops.h"

/////////////////////////////////////////////////////////////////////////////////////////
/// String manipulation operations, not dependent on the AST analysis routines
/////////////////////////////////////////////////////////////////////////////////////////

// If file has been compiled with -DGIT_SHA_VALUE=<val>, return the value

std::string git_sha_value() {
#if defined(GIT_SHA_VALUE)
#define xstr(s) makestr(s)
#define makestr(s) #s
    return xstr(GIT_SHA_VALUE);
#else
    return "not available";
#endif
}

/// this routine changes the input to alphanumeric + _, for naming purposes
std::string clean_name(const std::string &s) {

    std::string r = s;
    size_t j = 0;
    for (size_t i = 0; i < s.length(); i++, j++) {
        char c = s[i];
        if (std::isalnum(c) || c == '_')
            r[j] = c;
        else {
            switch (c) {
            case '+':
                r[j] = 'P';
                break;
            case '-':
                r[j] = 'M';
                break;
            case '*':
                r[j] = 'X';
                break;
            case '/':
                r[j] = 'D';
                break;
            case '=':
                r[j] = 'E';
                break;
            default:
                r[j] = '_';
            }
        }
    }
    return r;
}

std::string remove_initial_whitespace(const std::string &line) {
    // clear whitespace at the beginning of string
    size_t j = 0;
    for (char p : line) {
        if (!(p == ' ' || p == '\t'))
            break;
        j++;
    }
    if (j > 0)
        return line.substr(j, std::string::npos);
    return line;
}

std::string remove_all_whitespace(const std::string &line) {
    std::string out = line; // init string
    int j = 0;
    for (char p : line) {
        if (!std::isspace(p))
            out[j++] = p;
    }
    out.resize(j);
    return out;
}

/// True if string contains word (note: this word is c++ alphanumeric word, ie. split as
/// in )

std::string::size_type find_word(const std::string &in, const std::string &pattern, int pos,
                                 bool reverse) {
    int i;
    if (!reverse)
        i = in.find(pattern, pos);
    else
        i = in.rfind(pattern, pos);

    if (i == std::string::npos)
        return std::string::npos; // not found

    if (i > 0 && (std::isalnum(in[i - 1]) || in[i - 1] == '_'))
        return std::string::npos; // is at the end of a longer word
    if (i < in.length() - pattern.length() - 1) {
        char c = in[i + pattern.length()];
        if (std::isalnum(c) || c == '_')
            return std::string::npos; // word continues
    }
    return i;
}

// returns true if line contains the word list at the beginning of line.  list
// contains the word separated by whitespace.  Remainder, if non-nullptr, will contain
// the rest of the line if return is true.
bool contains_word_list(const std::string &line, const std::string &list, std::string *remainder) {
    const char *p = line.c_str();
    const char *q = list.c_str();
    while (*p && *q) {
        while (std::isspace(*p))
            p++;
        while (std::isspace(*q))
            q++;

        if (*p != *q)
            break;

        while (*p && *q && *p == *q) {
            p++;
            q++;
        }
    }
    // if line contained the words in list, *q = 0.
    while (std::isspace(*q))
        q++;
    if (*q != 0)
        return false;

    while (std::isspace(*p))
        p++;
    if (remainder != nullptr)
        *remainder = p;
    return true;
}

// remove extra whitespace chars
std::string remove_extra_whitespace(const std::string &line) {
    std::string out = line; // init string
    int j = 0;
    bool previous_char_space = true; // guarantees leading space removed
    for (char p : line) {

        if (!std::isspace(p)) {
            out[j++] = p;
            previous_char_space = false;
        } else {
            if (!previous_char_space)
                out[j++] = ' '; // substitute w. real ' '
            previous_char_space = true;
        }
    }
    if (j > 0 && out[j - 1] == ' ')
        j--; // remove trailing space
    out.resize(j);
    return out;
}

std::string indent_string(const std::string &s) {

    std::string indentstr = "  ";
    std::string res = "";
    std::string line;

    int lev = 0;
    size_t current = 0, i = 0;
    while (i < s.length()) {
        i = s.find('\n', current);
        if (i > s.length()) {
            // no more \n, print rest
            line = s.substr(current, s.length() - current);
        } else {
            line = s.substr(current, i - current + 1);
        }
        line = remove_initial_whitespace(line);
        // do the actual indent
        for (char p : line)
            if (p == '}')
                lev--;
        for (int j = 0; j < lev; j++)
            res += indentstr;
        for (char p : line)
            if (p == '{')
                lev++;

        res += line;
        current = i + 1;
    }
    return res;
}

std::string comment_string(const std::string &s) {

    const std::string comment("//--  ");
    std::string res = comment + s;

    size_t i = 0;
    while (i < res.length()) {
        i = res.find('\n', i);
        if (i < res.length()) {
            i++;
            res.insert(i, comment);
        }
    }
    return res;
}

/// From parity + dir -expression remove X, i.e.
/// X + dir -> dir
/// X - dir -> -dir
std::string remove_X(const std::string &s, bool *was_there) {
    std::string r = remove_extra_whitespace(s);
    if (r.size() == 0 || r[0] != 'X') {
        if (was_there != nullptr)
            *was_there = false;
        return r;
    }
    if (was_there != nullptr)
        *was_there = true;
    int i = 1;
    if (i < r.size() && std::isspace(r[i]))
        i++;
    if (i < r.size() && s[i] == '+') {
        i++;
        if (i < r.size() && std::isspace(r[i]))
            i++;
    }
    return r.substr(i, std::string::npos);
}

/// Types ofen seem to have "class name" -names, harmful
std::string remove_class_from_type(const std::string &s) {
    size_t i = s.find("class ", 0);
    if (i < s.size() && (i == 0 || !std::isalnum(s[i - 1]))) {
        return s.substr(0, i) + s.substr(i + 6, std::string::npos);
    } else
        return s;
}


/////////////////////////////////////////////////////////////////////////////////////
/// Get includes using system gcc - only if static compile!
/////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_COMPILER_INCLUDES

#include <fstream>

int get_includes_from_gcc(std::vector<const char *> &av) {

    std::string compiler = "g++";

    // scan args for opt
    bool got_compiler = false;
    char opt[] = "--" GET_INCLUDES_WITH "=";
    int optlen = strlen(opt);
    for (auto &r : av) {
        if (strncmp(r, opt, optlen) == 0) {
            if (strlen(r) <= optlen) {
                std::cerr << "ERROR: --get-includes-with=<compiler>: no compiler name given\n";
                exit(1);
            }
            compiler = r + optlen;
            got_compiler = true;
            break;
        }
    }

    if (compiler != "none") {
        char pipebuf[20000];
        FILE *pipe;
        if (!got_compiler) {
            // this pipe cmd gives bash completions for g++-<tab>
            pipe = popen("echo -n 'compgen -c g++-' | bash", "r");

            int version = 0;
            while (pipe && !feof(pipe) && fgets(pipebuf, sizeof(pipebuf), pipe) != NULL) {

                int v = 0, i = 4;
                while (pipebuf[i] >= '0' && pipebuf[i] <= '9') {
                    i++;
                    v *= 10;
                    v += pipebuf[i] - '0';
                }
                if (pipebuf[i] == '\n')
                    pipebuf[i] = 0;
                if (pipebuf[i] == 0 && v > version) {
                    compiler = pipebuf;
                    version = v;
                }
            }

            pclose(pipe);
        }

        // std::cerr << "FOUND COMPILER " << compiler << '\n';

        // The following commmand makes the compiler to produce list of include dirs
        std::string pipecmd =
            "echo | " + compiler + " -c -xc++ --std=c++17 -Wp,-v - 2>&1 | grep '^ '";

        // std::cerr << pipecmd << '\n';
        pipe = popen(pipecmd.c_str(), "r");

        // this needs to be static to store the strings!
        static std::vector<std::string> includedirs;

        if (pipe) {
            while (!feof(pipe) && fgets(pipebuf, sizeof(pipebuf), pipe) != NULL) {

                includedirs.push_back("-I");
                includedirs.back() += pipebuf + 1;
                includedirs.back().back() = 0;
            }
            for (auto &r : includedirs) {
                av.push_back(r.c_str());
                // std::cerr << r << '\n';
            }
        }
        pclose(pipe);
    }

    // for (auto &r : av) {
    //     std::cerr << r << '\n';
    // }

    return av.size();
}

#endif // USE_COMPILER_INCLUDES


// Check if the cmdline has -I<include> or -D<define> -
// arguments and move these after -- if that exists on the command line.
// Clang's optionparser expects these "generic compiler and linker"
// args to be after --
// return value new argc
int rearrange_cmdline(int argc, const char **argv, std::vector<const char *> &avvect) {

    bool found_ddash = false;
    // av[argc + 1] = nullptr;  // I read somewhere that in c++ argv[argc] = 0
    static char ddash[3] = "--"; // needs to be static because ptrs
    int ddashloc = 0;

    avvect.clear();

    for (int i = 0; i < argc; i++) {

        avvect.push_back(argv[i]);
        if (strcmp(avvect[i], ddash) == 0) {
            found_ddash = true;
            ddashloc = i;
        }
    }
    if (!found_ddash) {
        // add ddash, does not hurt in any case
        avvect.push_back(ddash);
        ddashloc = argc;
        argc++;
    }

    // now find -I and -D -options and move them after --
    for (int i = 0; i < ddashloc;) {
        if (i < ddashloc - 1 && (strcmp(avvect[i], "-D") == 0 || strcmp(avvect[i], "-I") == 0)) {
            // type -D define
            const char *a1 = avvect[i];
            const char *a2 = avvect[i + 1];
            for (int j = i + 2; j < argc; j++)
                avvect[j - 2] = avvect[j];
            avvect[argc - 2] = a1;
            avvect[argc - 1] = a2;
            ddashloc -= 2;
        } else if (strncmp(avvect[i], "-D", 2) == 0 || strncmp(avvect[i], "-I", 2) == 0) {
            // type -Ddefine
            const char *a1 = avvect[i];
            for (int j = i + 1; j < argc; j++)
                avvect[j - 1] = avvect[j];
            avvect[argc - 1] = a1;
            ddashloc--;
        } else {
            i++;
        }
    }

#ifdef USE_COMPILER_INCLUDES
    argc = get_includes_from_gcc(avvect);
#endif

    avvect.resize(avvect.size() + 3);
    avvect[argc++] = "-std=c++17"; // use c++17 std
    avvect[argc++] = "-DHILAPP";   // add global defn
    avvect[argc] = nullptr;

    return argc;
}
