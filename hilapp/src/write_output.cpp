#include <fstream>
#include <ctime>
#include "hilapp.h"
#include "stringops.h"

bool write_output_file(const std::string &filename, const std::string &buf) {

    auto tv = std::time(nullptr);
    std::string header = "// File generated by " + program_name + " at " +
                         std::asctime(std::localtime(&tv));
    header += "// Git version: " + git_sha_value() + '\n';

    header += "// cmd: ";
    std::string line;
    for (int i = 0; i < cmdline::argc; i++) {
        line += cmdline::argv[i];
        line += ' ';
        if (line.size() > 65) {
            header += line + '\n';
            line = "//        ";
        }
    }
    header += line + "\n\n";

    if (filename == "-") {
        llvm::outs() << header << buf;
    } else {
        std::string name;
        if (filename == "") {
            // use now main_file_name, strip dirs
            size_t s = global.main_file_name.rfind('/', std::string::npos);
            if (s < global.main_file_name.size()) {
                name = global.main_file_name.substr(s + 1, std::string::npos);
            } else {
                name = global.main_file_name;
            }

            // remove last max 3 letter extension TODO: is this safe?
            s = name.rfind('.', std::string::npos);
            if (s < name.size() && s >= name.size() - 4 && s > 0) {
                name = name.substr(0, s);
            }

            name += "." + default_output_suffix;
        } else {
            name = filename;
        }

        std::ofstream f(name, std::ios_base::out | std::ios_base::trunc);
        if (f.is_open()) {
            f << header << buf;
        }
        f.close();

        if (f.bad() || f.fail()) {
            llvm::errs() << program_name << ": output error for file " << name << '\n';
            exit(1);
        }
    }
    return true;
}
