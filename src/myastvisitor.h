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
#include "transformer.h" //global vars needed 
#include "generalvisitor.h"  // Definitions for the general visitor case

//////////////////////////////////////////////
/// myastvisitor.h : overloaded ASTVisitor for 
/// generating code from AST
///
/// Used in:
/// - myastvisitor.cpp
/// - loop_function.cpp
/// - codegen.cpp and its derivatives
/// 
//////////////////////////////////////////////



class MyASTVisitor : public GeneralVisitor, public RecursiveASTVisitor<MyASTVisitor> {

private:
  srcBuf *writeBuf;
  srcBuf *toplevelBuf;

  //flags used during AST parsing 
  struct {
    unsigned skip_children;     // if > 0 skip children of this ast node
    unsigned scope_level;       // level of variable scoping: {}
    int  ast_depth;             // depth of ast nodes within loop body.  ast_depth = 0 at top level
    int  stmt_sequence;         // sequence number of full statements in loops.  Full stmts separated by ;
    bool in_loop_body;          // true if in field loop
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

  /// True if the decl is preceded by "#pragma transformer <string>" where s is the string
  bool has_pragma(Decl *d, const char *s);
  bool has_pragma(Stmt *S, const char *s);
  bool has_pragma(const SourceLocation sl, const char *s);

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
  
  /// special handler for field<>
  int handle_field_specializations(ClassTemplateDecl *D);
  
  bool is_array_expr(Expr *E); 
  
  /// Checks if expr points to a variable defined in the same loop
  var_decl * is_loop_local_var_ref(Expr *E);

  bool is_assignment_expr(Stmt * s, std::string * opcodestr, bool & is_compound);
  
  bool is_loop_extern_var_ref(Expr *E);

  void check_allowed_assignment(Stmt * s);
  
  parity get_parity_val(const Expr *pExpr);
      
  bool check_field_ref_list();

  void check_var_info_list();
  
  bool handle_field_parity_X_expr(Expr *e, bool is_assign, bool is_compound, bool is_X, bool is_func_arg = false);
  
  void handle_var_ref(DeclRefExpr *E, bool is_assign, std::string & op, Stmt * assign_stmt = nullptr);
  void handle_array_var_ref(ArraySubscriptExpr *E, bool is_assign, std::string & op);

  var_info * new_var_info(VarDecl *decl);

  /// check the dependency chain of variables in assignments
  bool check_rhs_of_assignment(Stmt *s, std::vector<var_info *> * dependent = nullptr);

  /// this checks if the statement s is site-dependent inside site loops
  /// if return is false, vi (if non-null) will contain a list of variables
  /// which may turn out to be dependent on site later.  Check after loop complete!
  bool is_site_dependent(Expr *e, std::vector<var_info *> * vi = nullptr);

  /// Check that the addressof-operators and reference vars are OK
  void check_addrofops_and_refs(Stmt * S);


  // void handle_function_call_in_loop(Stmt * s, bool is_assignment, bool is_compund);
  void handle_function_call_in_loop(Stmt * s);
  void handle_loop_function_args(FunctionDecl *D, CallExpr *Call);

  void handle_member_call_in_loop(Stmt * s);

  void handle_constructor_in_loop(Stmt * s);

  bool check_argument( Expr * E, bool is_lvalue, bool output_only,
       argument_info & ai, std::vector<var_info *> & lvalue_refs);


  bool loop_function_check(Decl *fd);

  bool handle_loop_function_if_needed(FunctionDecl *fd);
  
  void handle_loop_function(FunctionDecl *fd);

  bool handle_special_loop_function(CallExpr *Call);

  // check if stmt is lf[par] = ... -type
  bool is_field_parity_assignment( Stmt *s );

  /// Does ; follow the statement?

  SourceRange getRangeWithSemicolon(Stmt * S, bool flag_error = true);
  
  void requireGloballyDefined(Expr * e);

  /// Entry point for the full field loop
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

  bool check_loop_vectorizable(Stmt *S, int & vector_size, std::string & diag);

  /// Generate a header for starting communication and marking fields changed
  std::string generate_code_cpu(Stmt *S, bool semicolon_at_end, srcBuf &sb, bool generate_wait);
  std::string generate_code_cuda(Stmt *S, bool semicolon_at_end, srcBuf &sb);
  void generate_openacc_loop_header(std::stringstream & code);
  //   std::string generate_code_openacc(Stmt *S, bool semicolon_at_end, srcBuf &sb);
  std::string generate_code_avx(Stmt *S, bool semicolon_at_end, srcBuf &sb, bool generate_wait);

  /// Handle functions called in a loop
  void handle_loop_function_cuda(FunctionDecl *fd);
  void handle_loop_function_openacc(FunctionDecl *fd);
  void handle_loop_function_avx(FunctionDecl *fd);


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
  
  /// utility used in inserting stuff after new line in buffer
  SourceLocation getSourceLocationAtEndOfLine( SourceLocation l );
  /// another utility (cannot trust r.getEnd())
  SourceLocation getSourceLocationAtEndOfRange( SourceRange r );

  // get next char and loc, while skipping comments
  SourceLocation getNextLoc(SourceLocation sl,bool forward = true);

  char getChar(SourceLocation sl);

  // get next word or symbol, if it is not a legal name symbol
  std::string getNextWord(SourceLocation sl);
  std::string getPreviousWord(SourceLocation sl);

  /// jump over following () expr
  SourceLocation skipParens( SourceLocation sl);

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
