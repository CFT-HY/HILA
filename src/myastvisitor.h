// -*- mode: c++ -*-
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

//////////////////////////////////////////////
/// myastvisitor.h : overloaded ASTVisitor for 
/// generating code from AST
///
/// Implemented by:
/// - myastvisitor.cpp
/// - codegen.cpp 
//////////////////////////////////////////////


class GeneralVisitor {
protected:

  Rewriter &TheRewriter;
  ASTContext *Context;

  //flags used during AST parsing 
  struct {
    unsigned skip_children;
    unsigned scope_level; 
    bool in_loop_body;
    bool accept_field_parity;
    bool dump_ast_next;
    bool check_loop;
    bool loop_function_next;
  } parsing_state;
  
public:

  GeneralVisitor(Rewriter &R) : TheRewriter(R) {
    parsing_state.skip_children = 0;
    parsing_state.scope_level = 0;
    parsing_state.in_loop_body = false;
    parsing_state.accept_field_parity = false;
    parsing_state.dump_ast_next = false;
    parsing_state.check_loop = false;
    parsing_state.loop_function_next = false;
  }

  GeneralVisitor(Rewriter &R, ASTContext *C) : TheRewriter(R) { 
    parsing_state.skip_children = 0;
    parsing_state.scope_level = 0;
    parsing_state.in_loop_body = false;
    parsing_state.accept_field_parity = false;
    parsing_state.dump_ast_next = false;
    parsing_state.check_loop = false;
    parsing_state.loop_function_next = false;
    Context=C; 
  }

  template <unsigned N>
  void reportDiag(DiagnosticsEngine::Level lev, const SourceLocation & SL,
                  const char (&msg)[N],
                  const char *s1 = nullptr,
                  const char *s2 = nullptr,
                  const char *s3 = nullptr );


  std::string get_stmt_str(const Stmt *s) {
    return TheRewriter.getRewrittenText(s->getSourceRange());
  }

  std::string get_expr_type(const Expr *e) {
    // This is somehow needed for printing type without "class" id
    // TODO: perhaps reanalyse and make more robust?
    PrintingPolicy pp(Context->getLangOpts());
    return e->getType().getUnqualifiedType().getAsString(pp);
  }
};


class MyASTVisitor : public GeneralVisitor, public RecursiveASTVisitor<MyASTVisitor> {

private:
  srcBuf *writeBuf;
  srcBuf *toplevelBuf;

public:
  using GeneralVisitor::GeneralVisitor;

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
  
  /// Visit function declarations
  bool VisitFunctionDecl(FunctionDecl *f);

  /// True if the decl is preceded by "#pragma transformer <string>" where s is the string
  bool has_pragma(Decl *d,const char *s);

  /// true if function contains parity loop
  bool does_function_contain_loop( FunctionDecl *f );
  
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

  // void VisitTypeAliasTemplateDecl(TypeAliasTemplateDecl *D);
  
  // bool VisitClassTemplateSpecalializationDeclImpl(ClassTemplateSpecializationDecl *D);
  // bool VisitClassTemplateSpecalializationDecl(ClassTemplateSpecializationDecl *D);
  
  bool is_field_storage_expr(Expr *E);
  bool is_field_expr(Expr *E);
  bool is_field_decl(ValueDecl *D);

  /// allowed index types: parity, parity_plus_direction, parity_plus_offset
  bool is_parity_index_type(Expr *E);

  // catches field[parity-type] expressions, incl. _plus -versions
  bool is_field_parity_expr(Expr *E);

  bool is_array_expr(Expr *E); 
  
  /// this tries to "fingerprint" expressions and see if they're duplicate
  bool is_duplicate_expr(const Expr * a, const Expr * b);
  
  /// Checks if expr points to a variable defined in the same loop
  var_decl * is_loop_local_var_ref(Expr *E);

  bool is_assignment_expr(Stmt * s, std::string * opcodestr, bool & is_compound);
  
  bool is_function_call_stmt(Stmt * s);

  bool is_constructor_stmt(Stmt * s);

  bool is_loop_extern_var_ref(Expr *E);

  void check_allowed_assignment(Stmt * s);
  
  parity get_parity_val(const Expr *pExpr);
    
  void require_parity_X(Expr * pExpr);
  
  bool check_field_ref_list();

  void check_var_info_list();
  
  bool handle_field_parity_expr(Expr *e, bool is_assign, bool is_compound);
  
  void handle_var_ref(DeclRefExpr *E, bool is_assign, std::string & op);
  void handle_array_var_ref(ArraySubscriptExpr *E);

  void handle_function_call_in_loop(Stmt * s, bool is_assignment, bool is_compund);
  void handle_function_call_in_loop(Stmt * s);
  
  void handle_constructor_in_loop(Stmt * s);

  bool loop_function_check(Decl *fd);

  bool handle_loop_function_if_needed(FunctionDecl *fd);
  
  void handle_loop_function(FunctionDecl *fd);

  bool handle_special_loop_function(CallExpr *Call);

  // check if stmt is lf[par] = ... -type
  bool is_field_parity_assignment( Stmt *s );

  /// Does ; follow the statement?
  bool isStmtWithSemi(Stmt * S);  
  SourceRange getRangeWithSemi(Stmt * S, bool flag_error = true);
  
  void requireGloballyDefined(Expr * e);

  /// Entry point for the full field loop
  bool handle_full_loop_stmt(Stmt *ls, bool field_parity_ok );

  /// Function for each stmt within loop body
  bool handle_loop_body_stmt(Stmt * s);

  void remove_vars_out_of_scope(unsigned level);
  
  // add handle to get rewriter too - for source control
  Rewriter &getRewriter() { return TheRewriter; }


  /// Code generation headers start here
  /// Starting point for new code
  void generate_code(Stmt *S, codetype & target);
  std::string backend_generate_code(Stmt *S, bool semi_at_end, srcBuf & loopBuf);
  void backend_handle_loop_function(FunctionDecl *fd);

  /// Generate a header for starting communication and marking fields changed
  std::string generate_code_cpu(Stmt *S, bool semi_at_end, srcBuf &sb);
  std::string generate_code_cuda(Stmt *S, bool semi_at_end, srcBuf &sb);
  std::string generate_code_openacc(Stmt *S, bool semi_at_end, srcBuf &sb);
  std::string generate_code_avx(Stmt *S, bool semi_at_end, srcBuf &sb);

  /// Handle functions called in a loop
  void handle_loop_function_cuda(FunctionDecl *fd);
  void handle_loop_function_openacc(FunctionDecl *fd);
  void handle_loop_function_avx(FunctionDecl *fd);

  /// Generate a candidate for a kernel name
  std::string make_kernel_name();

  /// Change field references within loops
  void replace_field_refs_and_funcs(srcBuf &sb);

  /// shortcut for "pragma"-like transformer_control("cmd")-functin
  // bool handle_control_stmt(Stmt *s);
  bool control_command(VarDecl *var);

  
  /// utility used in inserting stuff after new line in buffer
  SourceLocation getSourceLocationAtEndOfLine( SourceLocation l );
  /// another utility (cannot trust r.getEnd())
  SourceLocation getSourceLocationAtEndOfRange( SourceRange r );

  /// utility used in finding pragmas on the previous line
  SourceRange getSourceRangeAtPreviousLine( SourceLocation l );

  void set_writeBuf(const FileID fid);

  SourceRange get_templatefunc_decl_range(FunctionTemplateDecl *tf, FunctionDecl *f); 
  SourceRange get_func_decl_range(FunctionDecl *f);

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
