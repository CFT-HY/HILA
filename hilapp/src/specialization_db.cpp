//------------------------------------------------------------------------------
//
// hilapp tools to convert "lattice loops" into
// hardware-dependent "kernels".
//
// database for generated specializations
//
//
//------------------------------------------------------------------------------
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>

// c-type includes for stat call
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "clang/AST/AST.h"

#include "hilapp.h"
#include "specialization_db.h"
#include "stringops.h"

static std::vector<spec> specializations = {};

static bool db_loaded = false;

void load_spec_db() {
    std::ifstream db(specialization_db_filename);
    if (db.is_open()) {
        std::string line;
        while (std::getline(db, line)) {
            // skip commented lines (start with #)

            if (line.length() > 0 && remove_initial_whitespace(line)[0] != '#') {
                spec s;
                std::istringstream iss(line);
                iss >> s.decl;
                iss >> s.file;
                iss >> s.timestamp;

                if (iss.fail()) {
                    llvm::errs() << program_name << ": I/O error in "
                                 << specialization_db_filename << '\n';
                    if (iss.eof())
                        llvm::errs()
                            << "Too few arguments on a line \"" << line << "\"\n";
                    llvm::errs() << "Suggest removing " << specialization_db_filename
                                 << " and recompile\n";
                    exit(1);
                }
                // do not include own specializations in database, will be added again

                if (s.file != global.main_file_name)
                    specializations.push_back(s);
            }
        }
        db.close();
    } else {
        // db does not exist, nothing to do
        specializations.clear();
    }
}

spec *search_spec_db(std::string &decl) {
    for (auto &d : specializations) {
        if (d.decl == decl)
            return (&d);
    }
    return nullptr;
}

// checks whether this decl exists in database.  If so, returns true and filename in
// "here". if not, returns false and adds the item in db
bool in_specialization_db(const std::string &decl_in, std::string &here) {

    std::string decl = remove_all_whitespace(decl_in);

    if (!db_loaded)
        load_spec_db();
    db_loaded = true;

    spec *s = search_spec_db(decl);
    if (s == nullptr) {
        // add new item
        spec sp = {decl, global.main_file_name, std::time(nullptr)};
        specializations.push_back(sp);
    } else if (s->file != global.main_file_name) {
        // somebody else has specialization already -- however, it is possible that
        // it is waiting for compilation and it may not include it any more.

        struct stat status;
        time_t mtime;
#if defined(__APPLE__)
        mtime = status.st_mtimespec;
#else
        mtime = status.st_mtim.tv_sec;
#endif
        if (stat(s->file.c_str(), &status) == 0 && mtime > s->timestamp) {
            // If the mod time of the file is later than timestamp, future compilation
            // expected: steal the specialization
            s->timestamp = std::time(nullptr);
            s->file = global.main_file_name;
        } else {
            // "found" exit branch
            here = s->file;
            return true;
        }
    } else {
        // now main_file == s->file -- should not happen
        llvm::errs() << program_name << ": specialization db handling error\n";
        exit(1);
    }

    here = global.main_file_name;
    return false;
}

void write_specialization_db() {

    // have to load db in any case, because need to strip possible specializatons from
    // this file

    if (!db_loaded)
        load_spec_db();

    std::ofstream db(specialization_db_filename,
                     std::ios_base::out | std::ios_base::trunc);
    if (db.is_open()) {
        db << "# Automatically generated specialization database by " << program_name
           << '\n';
        db << "# Do not edit unless you know what you are doing\n";
        db << "# Removing this file should cause full recompilation\n#\n";
        db << "# Legend: specialization  compilation-unit  timestamp\n";
        db << "# ---------------------------------------------------------\n";
        for (auto &d : specializations) {
            db << d.decl << " " << d.file << " " << d.timestamp << '\n';
        }
        db << "# end\n";
        db.close();
    }
    if (db.bad() || db.fail()) {
        llvm::errs() << program_name << ": specialization db write error\n";
        exit(1);
    }

    specializations.clear();
    db_loaded = false;
}
