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
#include "srcbuf.h"

// set namespaces globally
using namespace clang;
//using namespace clang::driver;
using namespace clang::tooling;

// constant names for program
const std::string program_name("Transformer");
const std::string specialization_db_filename("specialization_db.txt");
const std::string default_output_suffix("cpt");
enum class reduction { NONE, SUM, PRODUCT };
enum class parity { none, even, odd, all, x };

// variables describing the type of code to be generated
struct codetype {
  bool kernelize;
  bool CUDA;
  bool openacc;
  bool flag_loop_function;
};


// collection of variables holding the state of parsing - definition in transformer.cpp
namespace state {
  extern unsigned skip_children; //= 0;
  extern unsigned scope_level; // = 0;
  extern int skip_next; // = 0;
  extern bool in_loop_body; // = false;
  extern bool accept_field_parity; // = false;
  extern bool loop_found; // = false;
  extern bool dump_ast_next; // = false;
  extern bool compile_errors_occurred; // = false;
  extern bool check_loop; // = false;
  extern bool no_device_code; // = false;
};

extern llvm::cl::OptionCategory TransformerCat;

namespace cmdline {
  // command line options
  extern llvm::cl::opt<bool> dump_ast;
  extern llvm::cl::opt<bool> no_include;
  extern llvm::cl::opt<std::string> dummy_def;
  extern llvm::cl::opt<std::string> dummy_incl;
  extern llvm::cl::opt<bool> function_spec_no_inline;
  extern llvm::cl::opt<bool> method_spec_no_inline; 
  extern llvm::cl::opt<bool> funcinfo;
  extern llvm::cl::opt<bool> no_output;
  extern llvm::cl::opt<bool> syntax_only;
  extern llvm::cl::opt<std::string> output_filename;
  extern llvm::cl::opt<bool> kernel;
  extern llvm::cl::opt<bool> vanilla;
  extern llvm::cl::opt<bool> CUDA;
  extern llvm::cl::opt<bool> openacc;
  extern llvm::cl::opt<bool> func_attribute;
};

struct loop_parity_struct {
  const Expr * expr;
  parity value;
  std::string text;
};

//class storing global variables
struct global_state {
  std::string main_file_name = "";
  bool assert_loop_parity = false;
  std::string full_loop_text = "";
  bool in_func_template = false;
  // bool in_class_template = false;
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


struct special_function_call {
  Expr * fullExpr;
  std::string full_expr;
  std::string replace_expression;
  bool add_loop_var;
  int scope;
};

bool write_output_file( const std::string & name, const std::string & buf ) ;
reduction get_reduction_type(bool, std::string &, var_info &);
void set_fid_modified(const FileID FID);
bool search_fid(const FileID FID);
srcBuf * get_file_buffer(Rewriter & R, const FileID fid);

// take global CI just in case
extern CompilerInstance *myCompilerInstance;
extern global_state global;
extern loop_parity_struct loop_parity;
extern codetype target;

/// global variable declarations - definitions on transformer.cpp

extern ClassTemplateDecl * field_decl;   // Ptr to field primary def in AST
extern ClassTemplateDecl * field_storage_type_decl;   // Ptr to field primary def in AST
extern const std::string field_element_type;
extern const std::string field_type;

// TODO: THESE SHOULD PROBABLY BE CHANGED INTO vectors,
// but they contain pointers to list elements.  pointers to vector elems are not good!
extern std::list<field_ref> field_ref_list;
extern std::list<field_info> field_info_list;
extern std::list<var_info> var_info_list;
extern std::list<var_decl> var_decl_list;
extern std::list<special_function_call> special_function_call_list;
extern std::vector<Expr *> remove_expr_list;
extern std::vector<FunctionDecl *> loop_functions;

#endif
