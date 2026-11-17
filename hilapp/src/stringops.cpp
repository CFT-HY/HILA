#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include "stringops.h"

/////////////////////////////////////////////////////////////////////////////////////////
/// String manipulation operations, not dependent on the AST analysis routines
/////////////////////////////////////////////////////////////////////////////////////////

// global to store the cpp standard
const char *cpp_standard = "c++17";

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
/// Locate the system headers for the clang parser:
///  - clang's builtin headers (stddef.h, stdarg.h, limits.h, intrinsics ...) from
///    a clang resource dir:  -resource-dir=<dir>, headers in <dir>/include.
///    These belong to the clang version hilapp is built with.
///  - (static build) the gcc installation whose C++ library (libstdc++) is used:
///    --gcc-install-dir=<dir>.  Clang then adds the libstdc++ and libc include
///    dirs itself, as system include dirs.
/// gcc's own builtin header dir (<gcc-install-dir>/include) is not used: those
/// headers are gcc specific.
/////////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <fstream>
#include "llvm/Config/llvm-config.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/Path.h"

namespace {

/// args to be added to av.  These are pushed to av only after all have been added,
/// because av stores the c_str() pointers
std::vector<std::string> added_args;

void add_arg(const std::string &arg) {
    added_args.push_back(arg);
}

/// is there an arg starting with prefix (user given on the command line)
bool has_arg_prefix(const std::vector<const char *> &av, const char *prefix) {
    for (const char *a : av)
        if (strncmp(a, prefix, strlen(prefix)) == 0)
            return true;
    return false;
}

/// run a shell command, return the first line of its stdout (empty on failure)
std::string first_output_line(const std::string &cmd) {
    std::string res;
    FILE *pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr)
        return res;
    char buf[4096];
    if (fgets(buf, sizeof(buf), pipe) != nullptr)
        res = buf;
    pclose(pipe);
    while (!res.empty() && (res.back() == '\n' || res.back() == '\r' || res.back() == ' '))
        res.pop_back();
    return res;
}

/// Check the clang resource dir candidate: return 1 if <dir>/include has the builtin
/// headers of this clang version, 0 if there are no headers, -1 if they belong to
/// another clang version (found_version).  hilapp/Makefile and get_hilapp.sh write
/// the clang major version to <dir>/include/.clang_version.
int check_resource_dir(const std::string &dir, int &found_version) {
    if (dir.empty() || !llvm::sys::fs::exists(dir + "/include/stddef.h"))
        return 0;
    std::ifstream vf(dir + "/include/.clang_version");
    if (vf >> found_version && found_version != LLVM_VERSION_MAJOR)
        return -1;
    return 1;
}

/// Add -resource-dir=<dir>.  Candidates:
///  1. <hilapp-dir>/include, where the hilapp binary is <hilapp-dir>/bin/hilapp or
///     <hilapp-dir>/build/hilapp.  Installed by hilapp/Makefile and get_hilapp.sh.
///  2. the resource dir of the clang hilapp was built with (HILAPP_CLANG_RESOURCE_DIR)
/// If neither is found, clang uses its default (relative to the binary).
void add_resource_dir(const char *argv0, std::vector<const char *> &av) {
    if (has_arg_prefix(av, "-resource-dir"))
        return;

    std::string exe = llvm::sys::fs::getMainExecutable(argv0, (void *)&add_resource_dir);
    std::string hilapp_dir =
        llvm::sys::path::parent_path(llvm::sys::path::parent_path(exe)).str();

    int version = 0;
    int status = check_resource_dir(hilapp_dir, version);
    if (status == 1) {
        add_arg("-resource-dir=" + hilapp_dir);
        return;
    }

#ifdef HILAPP_CLANG_RESOURCE_DIR
    int v2;
    if (check_resource_dir(HILAPP_CLANG_RESOURCE_DIR, v2) == 1) {
        add_arg(std::string("-resource-dir=") + HILAPP_CLANG_RESOURCE_DIR);
        return;
    }
#endif

    // Not found.  A non-static hilapp may still find the headers of the clang installation
    // (clang default), thus "not found" is a warning only for the static hilapp
    const char *hint = "  Reinstall hilapp with hilapp/get_hilapp.sh, or build it with "
                       "'make -C hilapp static' or 'make -C hilapp install'\n";
    if (status == -1) {
        std::cerr << "hilapp WARNING: clang builtin headers in " << hilapp_dir
                  << "/include are for clang " << version << ", but hilapp uses clang "
                  << LLVM_VERSION_MAJOR << ", ignoring them.\n"
                  << hint;
    }
#ifdef USE_COMPILER_INCLUDES
    else {
        std::cerr << "hilapp WARNING: clang builtin headers (stddef.h etc.) not found in "
                  << hilapp_dir << "/include.\n"
                  << hint;
    }
#endif
}

} // namespace

#ifdef USE_COMPILER_INCLUDES

/// Static build: hilapp uses the C++ standard library (libstdc++) of the gcc installation
/// of the compiler given with --get-includes-with=<compiler>, by default the newest g++-N
/// or g++.  The installation dir is from '<compiler> -print-libgcc-file-name', and it is
/// passed to clang with --gcc-install-dir.  With 'none' clang searches for gcc itself.
/// Returns the compiler name.
std::string add_gcc_installation(std::vector<const char *> &av) {

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

    if (compiler == "none")
        return compiler;

    if (!got_compiler) {
        // this pipe cmd gives bash completions for g++-<tab>
        // the code below selects the largest version number
        FILE *pipe = popen("echo -n 'compgen -c g++-' | bash", "r");
        char pipebuf[1000];

        int version = 0;
        while (pipe && !feof(pipe) && fgets(pipebuf, sizeof(pipebuf), pipe) != NULL) {

            int v = 0, i = 4;
            while (pipebuf[i] >= '0' && pipebuf[i] <= '9') {
                v *= 10;
                v += pipebuf[i] - '0';
                i++;
            }
            if (pipebuf[i] == '\n')
                pipebuf[i] = 0;
            if (pipebuf[i] == 0 && v > version) {
                compiler = pipebuf;
                version = v;
            }
        }

        if (pipe)
            pclose(pipe);
    }

    // gcc installation given explicitly to clang (after --)
    if (has_arg_prefix(av, "--gcc-install-dir") || has_arg_prefix(av, "--gcc-toolchain"))
        return compiler;

    // libgcc.a is in the gcc installation dir <prefix>/lib/gcc/<triple>/<version>
    std::string libgcc = first_output_line(compiler + " -print-libgcc-file-name 2>/dev/null");
    std::string dir = llvm::sys::path::parent_path(libgcc).str();

    if (dir.empty() || !llvm::sys::fs::exists(dir + "/crtbegin.o")) {
        std::cerr << "hilapp WARNING: could not find the gcc installation of compiler \""
                  << compiler << "\" (with '" << compiler
                  << " -print-libgcc-file-name'), letting clang search for it\n";
        return compiler;
    }

#if LLVM_VERSION_MAJOR >= 16
    add_arg("--gcc-install-dir=" + dir);
#else
    // older clang: give the installation prefix, clang chooses the version under it
    add_arg("--gcc-toolchain=" + dir + "/../../../..");
#endif

    return compiler;
}

#endif // USE_COMPILER_INCLUDES

// Check if the cmdline has -I<include> or -D<define> -
// arguments and move these after -- if that exists on the command line.
// Clang's optionparser expects these "generic compiler and linker"
// args to be after --
// return value new argc
std::string handle_cmdline_args(int argc, const char **argv, std::vector<const char *> &avvect) {

    bool found_ddash = false;
    // av[argc + 1] = nullptr;  // I read somewhere that in c++ argv[argc] = 0
    static char ddash[3] = "--"; // needs to be static because ptrs
    int ddashloc = 0;
    const char *cppstd = "-std=c++17";

    avvect.clear();

    for (int i = 0; i < argc; i++) {

        avvect.push_back(argv[i]);
        if (strcmp(avvect[i], ddash) == 0) {
            found_ddash = true;
            ddashloc = i;
        } else if (strncmp(avvect[i], "-std", 4) == 0) {
            if (strcmp(avvect[i], "-std=c++17") == 0 || strcmp(avvect[i], "-std=c++20") == 0 ||
                strcmp(avvect[i], "-std=c++23") == 0 || strcmp(avvect[i], "-std=c++26") == 0) {
                cppstd = avvect[i];
                cpp_standard = cppstd + 5;
            } else {
                std::cerr << "ERROR: unknown c++ standard \'" << avvect[i] << "\'\n";
                std::cerr << "allowed values -std=c++17, c++20, c++23, c++26\n";
                exit(1);
            }
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

    std::string compiler;
#ifdef USE_COMPILER_INCLUDES
    compiler = add_gcc_installation(avvect);
#endif
    add_resource_dir(argv[0], avvect);
    for (auto &a : added_args)
        avvect.push_back(a.c_str());

    avvect.push_back(cppstd);
    avvect.push_back("-DHILAPP");
    avvect.push_back(nullptr); // last element (avvect[argc]) nullptr

    return compiler;
}
