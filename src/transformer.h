// -*- mode: c++ -*-
#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include <string>
#include <vector>
#include <list>

#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
//#include "llvm/Support/raw_ostream.h"


// set namespaces globally
using namespace clang;
//using namespace clang::driver;
using namespace clang::tooling;

// constant names for program
const std::string program_name("Transformer");
const std::string specialization_db_filename("specialization_db.txt");
const std::string default_output_suffix("cpt");

bool write_output_file( const std::string & name, const std::string & buf ) ;

struct codetype {
  bool kernelize;
  bool CUDA;
  bool openacc;
};

enum class parity { none, even, odd, all, x };

struct loop_parity_struct {
  const Expr * expr;
  parity value;
  std::string text;
};

struct global_state {
  std::string main_file_name = "";
  bool assert_loop_parity = false;
  std::string full_loop_text = "";
  bool in_func_template = false;
  bool in_class_template = false;
  TemplateParameterList *function_tpl = nullptr;
//  std::vector<const TemplateParameterList *> class_templ_params = {};
//  std::vector<const TemplateArgumentList *> class_templ_args = {};
  FunctionDecl * currentFunctionDecl = nullptr;
  struct location_struct {
    SourceLocation function;
    SourceLocation top;
    SourceLocation bot;
    SourceLocation loop;
  } location;
};


// struct field_info;

struct field_ref {
  Expr * fullExpr;
  Expr * nameExpr;
  Expr * parityExpr;
  Expr * dirExpr;
  struct field_info * info;
  // unsigned nameInd, parityInd;
  int direction;
  bool is_written, is_read;
};

struct dir_ptr {
  Expr * e;
  unsigned count;
};
  
struct field_info {
  std::string type_template;             // This will be the <T> part of field<T>
  std::string old_name;                  // "name" of field variable, can be an expression
  std::string new_name;                  // replacement field name
  std::string loop_ref_name;             // var which refers to payload, loop_ref_name v = new_name->fs.payload
  bool is_written, is_read;              // is the field read from or written to in this loop
  std::vector<dir_ptr> dir_list;         // nb directions TODO: more general gather ptr
  std::vector<field_ref *> ref_list;     // where the var is referred at

  field_info() {
    type_template = old_name = new_name = loop_ref_name = "";
    is_written = is_read = false;
    dir_list = {};
    ref_list = {};
  }

  ~field_info() {
    dir_list.clear();
    ref_list.clear();
    type_template.clear();
    old_name.clear();
    new_name.clear();
    loop_ref_name.clear();
  }
};

struct var_ref {
  DeclRefExpr *ref;
  // unsigned ind;
  std::string assignop;
  bool is_assigned;
};

enum class reduction { NONE, SUM, PRODUCT };

struct var_info {
  std::vector<var_ref> refs;
  VarDecl * decl;
  struct var_decl * var_declp;
  std::string type;
  std::string name;
  bool is_loop_local;
  reduction reduction_type;
  bool is_assigned;
};


struct var_decl {
  VarDecl *decl;
  std::string name;
  std::string type;
  int scope;
};



#endif
