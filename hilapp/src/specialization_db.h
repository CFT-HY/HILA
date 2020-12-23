#ifndef SPECIALIZATION_DB_H
#define SPECIALIZATION_DB_H

#include <string>

/// File info for stored in the specialization database
struct spec {
  std::string decl;
  std::string file;
  std::time_t timestamp;
};

void load_spec_db();
spec * search_spec_db( std::string & decl );
bool in_specialization_db( const std::string & decl_in, std::string & here );
void write_specialization_db();

#endif
