#ifndef MYASTVISITOR_H
#define MYASTVISITOR_H

#include <string>
#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"

#include "srcbuf.h" //srcbuf class interface 
#include "hilapp.h" //global vars needed 
#include "generalvisitor.h"  // Definitions for the general visitor case

class TopLevelVisitor;

extern TopLevelVisitor * globalTopLevelVisitor;

//////////////////////////////////////////////
/// toplevelvisitor.h : overloaded ASTVisitor for 
/// generating code from AST
///
/// Used in:
/// - myastvisitor.cpp
/// - loop_function.cpp
/// - codegen.cpp and its derivatives
/// 
//////////////////////////////////////////////



class TopLevelVisitor : public GeneralVisitor, public RecursiveASTVisitor<TopLevelVisitor> {

private:
  srcBuf *writeBuf;
  srcBuf *toplevelBuf;

  //flags used during AST parsing 
  struct {
    unsigned skip_children;     // if > 0 skip children of this ast node
    unsigned scope_level;       // level of variable scoping: {}
    int  ast_depth;             // depth of ast nodes within loop body.  ast_depth = 0 at top level
    int  stmt_sequence;         // sequence number of full statements in loops.  Full stmts separated by ;
    bool in_loop_body;          // true if in site loop
    bool accept_field_parity;   // if parity of loop not resolved yet
    bool loop_function_next;
  } parsing_state;

public:
  using GeneralVisitor::GeneralVisitor;

  void reset_parsing_state() {
    parsing_state.skip_children = 0;
    parsing_state.scope_level = 0;
    parsing_state.ast_depth = 1;
    parsing_state.stmt_sequence = 0;
    parsing_state.in_loop_body = false;
    parsing_state.accept_field_parity = false;
    parsing_state.loop_function_next = false;
  }

  bool shouldVisitTemplateInstantiations() const { return true; }
  // is false by default, but still goes?

  /// TraverseStmt is called recursively for each level in the AST
  /// We can keep track of the level here
  bool TraverseStmt(Stmt *S);

  /// TraverseDecl is called recursively for each declaration in the AST
  /// We can keep track of the level here
  bool TraverseDecl(Decl *S);
  
  /// VisitStmt is called for each statement in AST.  Thus, when traversing the
  /// AST or part of it we always start from here
  bool VisitStmt(Stmt *s);

  bool VisitVarDecl(VarDecl *var);

  bool VisitDecl( Decl * D);
  bool VisitType( Type * T);
  
  /// Visit function declarations
  bool VisitFunctionDecl(FunctionDecl *f);

  /// typealiases are used to determine if class is vectorizable
  bool VisitTypeAliasDecl(TypeAliasDecl * ta);

  /// True if the decl is preceded by "#pragma hila <string>" where s is the string
  bool has_pragma(Decl *d, const char *s);
  bool has_pragma(Stmt *S, const char *s);
  bool has_pragma(const SourceLocation sl, const char *s);

  const char * get_pragma_value(SourceLocation sl);

  /// true if function contains parity loop
  bool does_function_contain_loop( FunctionDecl *f );

  /// check if there's field reference in the Expr.
  bool does_expr_contain_field(Expr *E);
  
  /// same for function templates
  // bool VisitFunctionTemplateDecl(FunctionTemplateDecl *tf);

  void specialize_function_or_method( FunctionDecl *f );

  int get_param_substitution_list( CXXRecordDecl * r,
                                   std::vector<std::string> & par,
                                   std::vector<std::string> & arg,
                                   std::vector<const TemplateArgument *> & typeargs );

  void make_mapping_lists( const TemplateParameterList * tpl, 
                           const TemplateArgumentList & tal,
                           std::vector<std::string> & par,
                           std::vector<std::string> & arg,
                           std::vector<const TemplateArgument *> & typeargs,
                           std::string *al );

  void check_spec_insertion_point(std::vector<const TemplateArgument *> & typeargs,
                                  SourceLocation ip, 
                                  FunctionDecl *f);
 
  // bool VisitCXXRecordDecl( CXXRecordDecl * D);

  /// and a hook for getting templated class template params
  bool VisitClassTemplateDecl(ClassTemplateDecl *D);

  /// handle the templated class specializations
  // int handle_class_specializations(ClassTemplateDecl *D);
  
  /// special handler for Field<>
  int handle_field_specializations(ClassTemplateDecl *D);
  
  bool is_array_expr(Expr *E); 
  
  void check_allowed_assignment(Stmt * s);
      
  bool check_field_ref_list();

  void check_var_info_list();
  
  bool handle_field_X_expr(Expr *e, bool is_assign, bool is_compound, bool is_X, bool is_func_arg = false);
  
  var_info * handle_var_ref(DeclRefExpr *E, bool is_assign, const std::string & op, Stmt * assign_stmt = nullptr);

  int handle_array_var_ref(ArraySubscriptExpr *E, bool is_assign, std::string & op);

  var_info * new_var_info(VarDecl *decl);

  /// Check that the addressof-operators and reference vars are OK
  void check_addrofops_and_refs(Stmt * S);

  // void handle_function_call_in_loop(Stmt * s, bool is_assignment, bool is_compund);
  void handle_function_call_in_loop(Stmt * s);

  call_info_struct handle_loop_function_args(FunctionDecl *D, CallExpr *Call, 
                                             bool sitedep, bool main_level = true);

  bool handle_call_argument( Expr *E, const ParmVarDecl * pv, bool sitedep,
                             std::vector<var_info *> * out_variables, 
                             std::vector<var_info *> * dep_variables,
                             argument_info & ai, bool main_level = true );

  void handle_member_call_in_loop(Stmt * s);

  void handle_constructor_in_loop(Stmt * s);

  bool attach_dependent_vars( std::vector<var_info *> & variables, bool sitedep,
                              std::vector<var_info *> & dep_variables );

  bool loop_function_check(Decl *fd);

  bool handle_loop_function_if_needed(FunctionDecl *fd);
  
  void process_loop_functions();

  void visit_loop_functions( std::vector<call_info_struct> & calls );

  bool handle_special_loop_function(CallExpr *Call);

  void clear_loop_function_calls();

  // check if stmt is lf[par] = ... -type
  bool is_field_parity_assignment( Stmt *s );

  /// Does ; follow the statement?

  SourceRange getRangeWithSemicolon(Stmt * S, bool flag_error = true);
  
  void requireGloballyDefined(Expr * e);

  /// Entry point for the full site loop
  bool handle_full_loop_stmt(Stmt *ls, bool field_parity_ok );

  /// Function for each stmt within loop body
  bool handle_loop_body_stmt(Stmt * s);

  void remove_vars_out_of_scope(unsigned level);
  
  // add handle to get rewriter too - for source control
  Rewriter &getRewriter() { return TheRewriter; }

  // shofthand for obtaining file buffer within this class
  srcBuf * get_file_srcBuf( SourceLocation sl ) {
    return get_file_buffer( TheRewriter, TheRewriter.getSourceMgr().getFileID(sl) ); 
  }

  /// Code generation headers start here
  /// Starting point for new code
  void generate_code(Stmt *S);
  void handle_field_plus_offsets(std::stringstream &code, srcBuf & loopbuf, std::string & par );

  std::string backend_generate_code(Stmt *S, bool semicolon_at_end, srcBuf & loopBuf, bool generate_wait);
  void backend_handle_loop_function(FunctionDecl *fd);
  void backend_handle_loop_constructor(CXXConstructorDecl *fd);

  bool check_loop_vectorizable(Stmt *S, int & vector_size, std::string & diag);

  /// Generate a header for starting communication and marking fields changed
  std::string generate_code_cpu(Stmt *S, bool semicolon_at_end, srcBuf &sb, bool generate_wait);
  std::string generate_code_cuda(Stmt *S, bool semicolon_at_end, srcBuf &sb, bool generate_wait);
  void generate_openacc_loop_header(std::stringstream & code);
  //   std::string generate_code_openacc(Stmt *S, bool semicolon_at_end, srcBuf &sb);
  std::string generate_code_avx(Stmt *S, bool semicolon_at_end, srcBuf &sb, bool generate_wait);

  /// Handle functions called in a loop
  void handle_loop_function_cuda(FunctionDecl *fd);
  void handle_loop_function_openacc(FunctionDecl *fd);
  void handle_loop_function_avx(FunctionDecl *fd);

  /// Handle functions called in a loop
  void handle_loop_constructor_cuda(CXXConstructorDecl *fd);
  void handle_loop_constructor_openacc(CXXConstructorDecl *fd);
  void handle_loop_constructor_avx(CXXConstructorDecl *fd);


  /// inspect if the type name is vectorizable
  /// returns the vectorized type in vectorized_type, if it is
  bool is_vectorizable_type(const std::string & type_name, vectorization_info & vi);
  bool is_vectorizable_type(const QualType & QT, vectorization_info & vi);

  /// Check if the field type is vectorizable and how
  vectorization_info inspect_field_type(Expr *fE);


  /// Generate a candidate for a kernel name
  std::string make_kernel_name();

  /// Change field references within loops
  void replace_field_refs_and_funcs(srcBuf &sb);
  
  /// utility used in finding pragmas on the previous line
  bool is_preceded_by_pragma( SourceLocation l, std::string & args, SourceLocation & ploc );

  void set_writeBuf(const FileID fid);

  SourceRange get_func_decl_range(FunctionDecl *f);

  void ast_dump(const Stmt *s);
  void ast_dump(const Decl *d);
  void ast_dump_header(const char *s, const SourceRange sr);

};


/// An AST Visitor for checking constraints for a field
/// reference expression. Walks the tree to check each
/// variable reference
class FieldRefChecker : public GeneralVisitor, public RecursiveASTVisitor<FieldRefChecker> {
public:
  using GeneralVisitor::GeneralVisitor;

  bool TraverseStmt(Stmt *s);
  bool VisitDeclRefExpr(DeclRefExpr * e);
};

/// An AST Visitor for checking constraints for assigments
/// in lattice loops
class LoopAssignChecker : public GeneralVisitor, public RecursiveASTVisitor<LoopAssignChecker> {
public:
  using GeneralVisitor::GeneralVisitor;

  bool TraverseStmt(Stmt *s);
  bool VisitDeclRefExpr(DeclRefExpr * e);
};





#endif
