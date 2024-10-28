#ifndef HILAPP_H
#define HILAPP_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

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
// using namespace clang::driver;
using namespace clang::tooling;

// constant names for program
const std::string program_name("hilapp");
const std::string specialization_db_filename("specialization_db.txt");
const std::string default_output_suffix("cpt");
const std::string out_only_keyword("out_only");
const std::string const_function_keyword("const_function");

const std::string name_prefix("_HILA_");
const std::string var_name_prefix(name_prefix + "var_");
const std::string field_name_prefix(name_prefix + "field_");
const std::string type_name_prefix(name_prefix + "type_");
const std::string kernel_name_prefix(name_prefix + "kernel_");


enum class reduction { NONE, SUM, PRODUCT }; // TBD:  MIN MAX MINLOC MAXLOC
enum class Parity { none, even, odd, all };

/// The following section contains command line options and functions
/// for implementing a backend

// variables describing the type of code to be generated
struct codetype {
    bool cuda = false;
    bool hip = false;
    bool vectorize = false;
    int vector_size = 1;
    bool openacc = false;
    bool openmp = false;
    bool kernelize = false;
    bool GPU = false;
};

extern codetype target; // make this global var

namespace cmdline {
// command line options
extern llvm::cl::opt<bool> dump_ast;
extern llvm::cl::opt<bool> no_include;
extern llvm::cl::opt<std::string> dummy_def;
extern llvm::cl::opt<std::string> dummy_incl;
extern llvm::cl::opt<bool> funcinfo;
extern llvm::cl::opt<bool> no_output;
extern llvm::cl::opt<bool> syntax_only;
extern llvm::cl::opt<bool> check_initialization;
extern llvm::cl::opt<std::string> output_filename;
extern llvm::cl::opt<bool> vanilla;
extern llvm::cl::opt<bool> CUDA;
extern llvm::cl::opt<bool> HIP;
extern llvm::cl::opt<bool> AVX512;
extern llvm::cl::opt<bool> AVX;
extern llvm::cl::opt<bool> SSE;
extern llvm::cl::opt<bool> openacc;
extern llvm::cl::opt<bool> c_openmp;
// extern llvm::cl::opt<bool> func_attribute;
extern llvm::cl::opt<int> vectorize;
extern llvm::cl::opt<bool> no_interleaved_comm;
// extern llvm::cl::opt<bool> no_mpi;
extern llvm::cl::opt<int> verbosity;
extern llvm::cl::opt<int> avx_info;
extern llvm::cl::opt<bool> comment_pragmas;
extern llvm::cl::opt<bool> insert_includes;
extern llvm::cl::opt<bool> slow_gpu_reduce;

extern llvm::cl::opt<bool> allow_func_globals;

// save also the original argc and argv
extern int argc;
extern const char **argv;
}; // namespace cmdline

namespace state {
extern bool compile_errors_occurred; // = false;
};

extern llvm::cl::OptionCategory HilappCategory;

/**
 * @brief class storing global state variables used in parsing
 *
 */
struct global_state {
    std::string main_file_name = "";
    bool assert_loop_parity = false;
    std::string full_loop_text = "";
    bool in_func_template = false;
    // bool in_class_template = false;
    TemplateParameterList *function_tpl = nullptr;
    //  std::vector<const TemplateParameterList *> class_templ_params = {};
    //  std::vector<const TemplateArgumentList *> class_templ_args = {};
    FunctionDecl *currentFunctionDecl = nullptr;
    int namespace_level = 0;
    SourceRange namespace_range;

    /// Describes a code location
    struct location_struct {
        SourceLocation function;
        SourceLocation top;
        SourceLocation bot;
        SourceLocation loop;
        SourceLocation kernels;
    } location;
};

/// field_ref contains info about references to field vars within loops
struct field_ref {
    Expr *fullExpr;          // full expression a[X+d]
    Expr *nameExpr;          // name "a"
    Expr *parityExpr;        // expr within [], here "X+d" or "X+e_x+e_y"
    std::string direxpr_s;   // original dir expr: "d" or "e_x+e_y" etc.
    struct field_info *info; // ptr to field info struct
    int sequence;            // sequence of the full stmt where ref appears
    bool is_written, is_read;
    bool is_direction;          // true if ref contains nn OR offset Direction - used as a
                                // general flag
    bool is_constant_direction; // true if dir is const. e_x etc.
    bool is_offset;             // true if dir is for offset instead of simple Direction
    bool is_loop_local_dir;     // if dir depends on loop local var - not known outside loop
    unsigned constant_value;

    field_ref() {
        fullExpr = nameExpr = parityExpr = nullptr;
        direxpr_s.clear();
        info = nullptr;
        is_written = is_read = is_offset = is_direction = is_constant_direction =
            is_loop_local_dir = false;
        sequence = 0;
    }

    ~field_ref() {
        direxpr_s.clear();
    }
};

/// dir_ptr is a "subfield" of field_info, storing Direction/offset of field ref
/// There may be several equivalent field[dir] -references within the loop,
/// these are described the same dir_ptr struct
struct dir_ptr {
    Expr *parityExpr; // pointer to full parity+dir expression (1st of equivalent ones)
    // Expr *dirExpr;            // non-null only for non-const. nn  NOT USEFUL
    std::string direxpr_s;             // Direction expression string (1st of equivalent ones)
    std::vector<field_ref *> ref_list; // pointers references equivalent to this Field[dir]
    unsigned count;             // how many genuine Direction refs?  if count==0 this is offset
    bool is_offset;             // is this dir offset?
    bool is_loop_local_dir;     // is direction expr loop local?
    bool is_constant_direction; // if constant nn
    unsigned constant_value;    // value of it
    std::string name_with_dir;  // new name for this Field[X+dir] -variable

    dir_ptr() {
        ref_list = {};
        parityExpr = nullptr;
        count = 0;
        is_offset = is_constant_direction = is_loop_local_dir = false;
        name_with_dir.clear();
        direxpr_s.clear();
    }

    ~dir_ptr() {
        ref_list.clear();
        name_with_dir.clear();
        direxpr_s.clear();
    }
};

/// field_type_info contains information of the template argument of
/// Field<type> -expression.  Type is vectorizable if:
/// a) just float, double, int  or
/// b) is templated type, with float/double in template and  implements
///    the method using base_type = hila::arithmetic_type<T>;
enum class number_type { INT, UNSIGNED, LONG, UNSIGNED_LONG, FLOAT, DOUBLE, LONG_DOUBLE, UNKNOWN };

// Collection of legal number types
struct legal_types {
    std::vector<std::string> vector = {"int",     "unsigned", "long",
                                       "int64_t", "uint64_t", "unsigned long",
                                       "float",   "double",   "long double"};

    bool check_if_legal(std::string is_legal_type) {
        for (const auto &vector_element : this->vector) {
            if (is_legal_type == vector_element)
                return true;
        }
        return false;
    }

    std::string as_string() {
        const std::string delimiter = ", ";
        std::string vector_as_str;
        for (const auto &vector_element : this->vector)
            vector_as_str += vector_element + delimiter;
        return vector_as_str.erase(vector_as_str.size() - delimiter.size(), delimiter.size());
    }

    void add_type(std::string type) {
        auto temporary(this->vector);
        for (auto &element : temporary) {
            this->vector.push_back(type + "<" + element + ">");
        }
    }
};

/// Stores information about how a field is vectorized
struct vectorization_info {
    bool is_vectorizable;
    int vector_size;
    number_type basetype;
    std::string basetype_str;    // name of the base var
    std::string vectortype;      // basetype -> vectorized, e.g. "Vec4d" etc.
    std::string vectorized_type; // var.type now as vectorized
};

/// main struct for storing info about each field variable inside loops
/// one field_info for each loop variable
struct field_info {
    std::string type_template;         // This will be the <T> part of Field<T>
    std::string element_type;          // type of the element of Field
    std::string old_name;              // "name" of Field variable, can be an expression
    std::string new_name;              // replacement Field name
    std::string loop_ref_name;         // var which refers to payload, loop_ref_name v =
                                       // new_name->fs.payload
    std::vector<dir_ptr> dir_list;     // nb directions TODO: more general gather ptr
    std::vector<field_ref *> ref_list; // where the var is referred at
    Expr *nameExpr;                    // first of the name exprs to this Field
    vectorization_info vecinfo;        // info of the type in Field<type>

    bool is_written;        // is the Field written to in this loop
    bool is_read_atX;       // local read, i.e. Field[X]
    bool is_read_nb;        // read using nn-neighbours
    bool is_read_offset;    // read with an offset (non-nn) index
    bool is_loop_local_dir; // is read with loop local direction?
    int first_assign_seq;   // the sequence of the first assignment

    field_info() {
        type_template = old_name = new_name = loop_ref_name = "";
        is_written = is_read_nb = is_read_atX = is_read_offset = is_loop_local_dir = false;
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

////////////////////////////////////////////////////////////////////////////////////////
/// Stores information about a single reference
/// to a variable
struct var_ref {
    DeclRefExpr *ref;
    // unsigned ind;
    std::string assignop;
    bool is_assigned;
};

/// This struct keeps track of all variables appearing in loops
/// variable can be external or defined inside loop (loop_local)
/// is_site_dependent means that the variable value can be site dependent,
/// which has implications for vectorization
struct var_info {
    std::vector<var_ref> refs;  // references of this var in loop
    std::string type;           // type as string
    std::string name;           // variable name
    std::string new_name;       // name to be used in loop
    VarDecl *decl;              // declaration of this var
    std::string reduction_name; // name of reduction variable
    std::vector<var_info *>
        dependent_vars;             // vector of var_infos which may affect is_site_dependent
    reduction reduction_type;       // what type of reduction
    vectorization_info vecinfo;     // info about vectorization
    bool is_loop_local;             // true if defined inside loop
    bool is_assigned;               // is the var assigned to
    bool is_site_dependent;         // is the value of variable site dependent
    bool is_special_reduction_type; // is variable defined with Reduction<T>
    bool is_raw;                    // is it raw access var?
    bool is_lambda_var_ref;         // is it a var ref to a lambda function

    var_info() {
        is_loop_local = is_assigned = is_site_dependent = is_special_reduction_type = is_raw =
            is_lambda_var_ref = false;
        decl = nullptr;
        reduction_type = reduction::NONE;
    }
};

///////////////////////////////////////////////////////////////////////////////////////////
/// Store reference to constant expression (field of a struct, array element)

struct loop_const_expr_ref {
    std::vector<Expr *> refs; // references of this expression in loop
    std::string type;         // type as string
    std::string expression;   // (1st) expression as string verbatim
    std::string exprstring;   // compressed expression as a string, for comparison
    std::string new_name;     // name to be used in loop
    reduction reduction_type; // is the expr used in reduction?  allow it
};


///////////////////////////////////////////////////////////////////////////////////////////
//  Store "array-type" references inside site loops
struct bracket_ref_t {
    Expr *E;                 // full bracket expression "a[i]"
    Expr *BASE;              // "a"  - can be DeclRefExpr, MemberExpr or CXXThisExpr
    Expr *DRE;               // Variable where this is found.  == BASE in simple cases
    std::vector<Expr *> Idx; // "i", one element per index - note outermost as first element
    Stmt *assign_stmt;       // if it is reduction ref, full assign stmt here
};

/// Stores information for a single reference to an loop-extern array or related
/// variable These are similar to variable references, but often need to be handled
/// differently
struct array_ref {
    std::vector<bracket_ref_t> refs; // vector of array references
    VarDecl *vd;                     // declaration of referred to variable
    std::string name;
    std::string new_name;
    std::string element_type;       // type of element
    std::string wrapper_type;       // if need to wrap this var
    size_t size;                    // size of array (if known)
    std::vector<size_t> dimensions; // multidimensional array size
    std::string size_expr;          // if only size expression known
    std::string data_ptr;           // == name for array, name.data() for others
    using reftype = enum { REPLACE, ARRAY, STD_VECTOR, STD_ARRAY, REDUCTION };
    reduction reduction_type;
    reftype type;

    array_ref() {
        refs = {};
        size = 0;
        dimensions.clear();
    }
};

/// store necessary information for vector reductions
// struct vector_reduction_ref {
//     CXXOperatorCallExpr *ref;
//     std::string type;
//     std::string index_name;
//     std::string vector_name;
//     std::string new_vector_name;
//     reduction reduction_type;
//     std::string reduction_name;
// };

/// store information about variable declarations
struct var_decl {
    VarDecl *decl;
    std::string name;
    std::string type;
    int scope;
};

/// store information about calls to special Hila functions
struct special_function_call {
    Expr *fullExpr;
    std::string full_expr;
    std::string name;
    std::string replace_expression;
    std::string args_string;
    Expr *argsExpr;
    SourceRange replace_range;
    bool add_loop_var;
    int scope;
};

/// Stores the parity of the current loop: Expr, value (if known), Expr as string
struct loop_info_struct {
    const Expr *parity_expr;
    std::string parity_text; // parity text in source
    std::string parity_str;  // what string to use
    Parity parity_value;

    bool has_pragma_novector;
    bool has_pragma_access;
    bool has_pragma_safe;
    bool has_pragma_omp_parallel_region;
    const char *pragma_access_args;
    const char *pragma_safe_args;
    bool has_site_dependent_cond_or_index;    // if, for, while w. site dep. cond?
    bool contains_random;                     // does it contain rng (also in loop functions)?
    bool has_conditional;                     // if, for, while, switch, ternary in loop
    std::vector<var_info *> conditional_vars; // may depend on variables
    Expr *condExpr;

    SourceRange range;

    inline void clear_except_external() { // do not remove parity values, may be set in loop init
        has_site_dependent_cond_or_index = contains_random = has_conditional = false;
        conditional_vars.clear();
        condExpr = nullptr;
    }
};

/// Stores information about loop function arguments
struct argument_info {
    Expr *E;
    ParmVarDecl *PV;
    std::vector<var_info *> dependent_vars;
    bool is_modifiable;
    bool is_site_dependent;
    bool is_const;
    bool is_out_only;
    bool is_const_function;

    argument_info() {
        is_modifiable = is_site_dependent = is_out_only = is_const = is_const_function = false;
        E = nullptr;
        PV = nullptr;
        dependent_vars = {};
    }
};

/// Stores information about loop function calls
struct call_info_struct {
    FunctionDecl *funcdecl;
    CXXConstructorDecl *ctordecl;
    CallExpr *call;
    CXXConstructExpr *constructor;
    std::vector<argument_info> arguments;
    Expr *condExpr;
    number_type vector_type_only; // This is number_type::UNKNOWN if the type is not restricted
    argument_info object;
    bool decl_only;
    bool is_operator;
    bool is_method;
    bool is_vectorizable;
    bool is_site_dependent;
    bool has_site_dependent_conditional;
    bool contains_random;
    bool is_defaulted;         // We assume that Clang "default" classification functions
                               // do not need to be handled (incl. default methods)
    bool is_loop_local_lambda; // if it calls to a lambda function declared inside the loop

    call_info_struct() {
        call = nullptr;
        constructor = nullptr;
        funcdecl = nullptr;
        ctordecl = nullptr;
        arguments.clear();
        is_method = is_operator = is_site_dependent = contains_random = false;
        has_site_dependent_conditional = is_defaulted = is_loop_local_lambda = false;
        decl_only = false;
        is_vectorizable = true;
        vector_type_only = number_type::UNKNOWN;
    }
};

/// store information about SiteSelect and SiteValueSelect operations
struct selection_info {
    CXXMemberCallExpr *MCE; // select-expr in loop (a.select())
    Expr *ref;              // var expression 'a'
    Expr *assign_expr;      // assignment if value select
    std::string val_type;
    std::string new_name; //
    selection_info
        *previous_selection; // pointer to first ref to the same variable (otherwise nullptr)
    std::string valname;     // names of tmp variables - used for GPU code
    std::string maskname;
    std::string sitename;
};


/// accumulate reductions (except vector reductions).  This is done just before
/// code generation in codegen
struct reduction_expr {
    std::vector<Expr *> refs;   // refrences to this variable / expression
    std::string name;           // variable name or expression to be assigned
    std::string reduction_name; // name generated for use in loop to collect reduction
    std::string loop_name;      // name used in loop for avx or gpu targets
    std::string type;           // type as string
    var_info *variable; // if this reduction is a variable in var_info_list, point to it, otherwise
                        // nullptr
    reduction reduction_type;   // NONE or SUM
    bool is_special_reduction;  // true if special reduction var
    vectorization_info vecinfo; // is type avx vectorizable?
};


enum class pragma_hila {
    SKIP,
    AST_DUMP,
    LOOP_FUNCTION,
    NOVECTOR,
    NOT_VECTORIZABLE,
    CONTAINS_RNG,
    ACCESS,
    SAFE,
    IN_OMP_PARALLEL_REGION
};

/// Pragma handling things
bool has_pragma_hila(const SourceManager &SM, SourceLocation l0, pragma_hila pragma,
                     SourceLocation &pragmaloc, const char **arg = nullptr);

/// Some sourceloc utilities

SourceLocation getNextLoc(const SourceManager &SM, SourceLocation sl, bool forward = true);
char getChar(const SourceManager &SM, SourceLocation sl);
SourceLocation findChar(const SourceManager &SM, SourceLocation sloc, char ct);
SourceLocation skipParens(const SourceManager &SM, SourceLocation sl, const char partype = '(');
SourceLocation skipString(const SourceManager &SM, SourceLocation sl);
std::string getNextWord(const SourceManager &SM, SourceLocation sl, SourceLocation *end = nullptr);
std::string getRangeText(const SourceManager &SM, SourceLocation begin, SourceLocation end);

bool write_output_file(const std::string &name, const std::string &buf);
reduction get_reduction_type(bool, const std::string &, var_info &);
void set_fid_modified(const FileID FID);
bool search_fid(const FileID FID);
srcBuf *get_file_buffer(Rewriter &R, const FileID fid);

/// reset the status of vectorizable types
void reset_vectorizable_types();

/// and clear loop function info
void clear_loop_functions_in_compilation_unit();

/// check if macro "name" is defined (non-function macros only)
bool is_macro_defined(const char *name, std::string *arg = nullptr);

/// Utility to check if typename is native number
number_type get_number_type(const std::string &s);

/// take global CI just in case
extern CompilerInstance *myCompilerInstance;
extern global_state global;
extern loop_info_struct loop_info;
extern codetype target;

/// global variable declarations - definitions on hilapp.cpp

extern ClassTemplateDecl *field_decl; // Ptr to field primary def in AST
extern const std::string field_storage_type;
extern const std::string field_type;

/// global lists used in modifying the site loops
/// but they contain pointers to list elements.  pointers to vector elems are not good!
extern std::list<field_ref> field_ref_list;
extern std::list<field_info> field_info_list;
extern std::list<array_ref> array_ref_list;
extern std::list<loop_const_expr_ref> loop_const_expr_ref_list;
extern std::list<special_function_call> special_function_call_list;
extern std::list<selection_info> selection_info_list;
extern std::vector<reduction_expr> reduction_list;


// struct func_prototype_struct {
//     FunctionDecl *decl;   // this is prototype decl
//     FunctionDecl *def;    // this should be the definition
//     int flag;             // use this for convenient information
// };

// extern std::vector<func_prototype_struct> prototype_vector;


#endif
