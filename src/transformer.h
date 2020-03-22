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


/// The following section contains command line options and functions
/// for implementing a backend

// variables describing the type of code to be generated
struct codetype {
  bool CUDA=false;
  bool VECTORIZE=false;
  int vector_size=1;
  bool openacc=false;
};

extern codetype target;  // make this global var

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
  extern llvm::cl::opt<bool> AVX512;
  extern llvm::cl::opt<bool> AVX;
  extern llvm::cl::opt<bool> SSE;
  extern llvm::cl::opt<bool> openacc;
  extern llvm::cl::opt<bool> func_attribute;
  extern llvm::cl::opt<int> VECTORIZE;
};

namespace state {
  extern bool loop_found; // = false;
  extern bool compile_errors_occurred; // = false;
};

extern llvm::cl::OptionCategory TransformerCat;


// Stores the parity of the current loop: Expr, value (if known), Expr as string

struct loop_parity_struct {
  const Expr * expr;
  parity value;
  std::string text;
};

//class storing global state variables used in parsing

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


// field_ref contains info about references to field vars within loops

struct field_ref {
  Expr * fullExpr;              // full expression a[X+d]
  Expr * nameExpr;              // name "a"
  Expr * parityExpr;            // expr within [], here "X+d" or "X+XUP+YUP"
  // Expr * dirExpr;               // expr of the directon -- non-null only for nn-dirs! NOT VERY USEFUL
  std::string direxpr_s;        // original dir expr: "d" or "XUP+YUP" etc.
  struct field_info * info;     // ptr to field info struct
  // unsigned nameInd, parityInd;
  int  sequence;                // sequence of the full stmt where ref appears
  bool is_written, is_read;
  bool is_direction;            // true if ref contains nn OR offset direction - used as a general flag
  bool is_constant_direction;   // true if dir is const. XUP etc.
  bool is_offset;               // true if dir is for offset instead of simple direction
  unsigned constant_value;

  field_ref() {
    fullExpr = nameExpr = parityExpr = nullptr;
    direxpr_s.clear();
    info = nullptr;
    is_written = is_read = is_offset = is_direction = is_constant_direction = false;
    sequence = 0;
    
  }

  ~field_ref() {
    direxpr_s.clear();
  }
};


// dir_ptr is a "subfield" of field_info, storing direction/offset of field ref
// There may be several equivalent field[dir] -references within the loop, 
// these are described the same dir_ptr struct

struct dir_ptr {
  Expr *parityExpr;         // pointer to full parity+dir expression (1st of equivalent ones)
  // Expr *dirExpr;            // non-null only for non-const. nn  NOT USEFUL
  std::string direxpr_s;     // direction expression string (1st of equivalent ones)
  std::vector<field_ref *> ref_list;  // pointers references equivalent to this field[dir]
  unsigned count;           // how many genuine direction refs?  if count==0 this is offset
  bool is_offset;           // is this dir offset?
  bool is_constant_direction;  // if constant nn
  unsigned constant_value;  // value of it
  std::string name_with_dir;         // new name for this field[X+dir] -variable

  dir_ptr() {
    ref_list = {};
    parityExpr = nullptr;
    count = 0;
    is_offset = is_constant_direction = false;
    name_with_dir.clear();
    direxpr_s.clear();
  }

  ~dir_ptr() {
    ref_list.clear();
    name_with_dir.clear();
    direxpr_s.clear();
  }
};


// main struct for storing info about each field variable inside loops
// one field_info for each loop variable
  
struct field_info {
  std::string type_template;             // This will be the <T> part of field<T>
  std::string old_name;                  // "name" of field variable, can be an expression
  std::string new_name;                  // replacement field name
  std::string loop_ref_name;             // var which refers to payload, loop_ref_name v = new_name->fs.payload
  std::vector<dir_ptr> dir_list;         // nb directions TODO: more general gather ptr
  std::vector<field_ref *> ref_list;     // where the var is referred at
  bool is_written;                       // is the field written to in this loop
  bool is_read_atX;                      // local read, i.e. field[X]
  bool is_read_nb;                       // read using nn-neighbours 
  bool is_read_offset;                  // read with an offset (non-nn) index
  int  first_assign_seq;                 // the sequence of the first assignment

  field_info() {
    type_template = old_name = new_name = loop_ref_name = "";
    is_written = is_read_nb = is_read_atX = is_read_offset = false;
    first_assign_seq = 0;
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
  std::string new_name;
  bool is_loop_local;
  reduction reduction_type;
  std::string reduction_name;
  bool is_assigned;
};

// Stores onformation for a single reference to an array
// These are similar to variable references, but often
// need to be handled differently
struct array_ref {
  ArraySubscriptExpr *ref;
  std::string new_name;
  std::string type;
};


// store necessary information for vector reductions
struct vector_reduction_ref {
  CXXOperatorCallExpr *ref;
  std::string type;
  std::string index_name;
  std::string vector_name;
  std::string new_vector_name;
  reduction reduction_type;
  std::string reduction_name;
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
  std::string name;
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
extern ClassTemplateDecl * field_storage_decl;   // Ptr to field primary def in AST
extern const std::string field_storage_type;
extern const std::string field_type;

// global lists used in modifying the field loops
// but they contain pointers to list elements.  pointers to vector elems are not good!
extern std::list<field_ref> field_ref_list;
extern std::list<field_info> field_info_list;
extern std::list<var_info> var_info_list;
extern std::list<var_decl> var_decl_list;
extern std::list<array_ref> array_ref_list;
extern std::list<vector_reduction_ref> vector_reduction_ref_list;
extern std::list<special_function_call> special_function_call_list;
extern std::vector<Expr *> remove_expr_list;

#endif
