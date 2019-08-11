// -*- mode: c++ -*-
#ifndef MYASTVISITOR_H
#define MYASTVISITOR_H

#include <string>
#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
//#include "llvm/Support/raw_ostream.h"

#include "stringops.h"
#include "srcbuf.h"
#include "transformer.h"
// #include "codegen.h"

class MyASTVisitor : public RecursiveASTVisitor<MyASTVisitor> {

private:  
  Rewriter &TheRewriter;
  ASTContext *Context;
  srcBuf Buf;

public:
  MyASTVisitor(Rewriter &R) : TheRewriter(R) {}
  MyASTVisitor(Rewriter &R, ASTContext *C) : TheRewriter(R) { Context=C; }

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

  /// same for function templates
  bool VisitFunctionTemplateDecl(FunctionTemplateDecl *tf);

  /// same for function templates
  bool VisitCXXMethodDecl(CXXMethodDecl *tf);

  /// and a hook for getting templated class template params
  bool VisitClassTemplateDecl(ClassTemplateDecl *D);

  int handle_field_specializations(ClassTemplateDecl *D);

  // void VisitTypeAliasTemplateDecl(TypeAliasTemplateDecl *D);
  
  bool VisitClassTemplateSpecalializationDeclImpl(ClassTemplateSpecializationDecl *D);
  bool VisitClassTemplateSpecalializationDecl(ClassTemplateSpecializationDecl *D);
  
  bool is_field_element_expr(Expr *E);
  bool is_field_expr(Expr *E);
  bool is_field_decl(ValueDecl *D);
  
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
    PrintingPolicy pp(Context->getLangOpts());
    return e->getType().getUnqualifiedType().getAsString(pp);
  }
  
  void replace_expr(Expr *e, const std::string &s) {
    TheRewriter.ReplaceText(e->getSourceRange(),s);
  }  

  bool is_duplicate_expr(const Expr * a, const Expr * b);
  
  // catches both parity and parity_plus_direction 
  bool is_field_parity_expr(Expr *e);

  /// Checks if expr points to a variable defined in the same loop
  var_decl * is_loop_local_var_ref(Expr *E);

  bool is_assignment_expr(Stmt * s, std::string * opcodestr);
  
  bool is_loop_extern_var_ref(Expr *E);
  
  parity get_parity_val(const Expr *pExpr);
    
  void require_parity_X(Expr * pExpr);
  
  bool check_field_ref_list();

  void check_var_info_list();
  
  bool handle_field_parity_expr(Expr *e, bool is_assign);
  
  void handle_var_ref(DeclRefExpr *E, bool is_assign, std::string & op);

  // check if stmt is lf[par] = ... -type
  bool is_field_parity_assignment( Stmt *s );

  /// Does ; follow the statement?
  bool isStmtWithSemi(Stmt * S);  
  SourceRange getRangeWithSemi(Stmt * S, bool flag_error = true);
  
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

  std::string generate_kernel(Stmt *S, bool semi_at_end);
  std::string generate_in_place(Stmt *S, bool semi_at_end);

  /// Generate a candidate for a kernel name
  std::string make_kernel_name();

  /// Change field references within loops
  void replace_field_refs();

  /// shortcut for "pragma"-like transformer_control("cmd")-functin
  bool handle_control_stmt(Stmt *s);

  /// utility used in inserting stuff after new line in buffer
  SourceLocation getSourceLocationAfterNewLine( SourceLocation l );
  
};

#endif
