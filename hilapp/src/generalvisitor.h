#ifndef GENERALVISITOR_H
#define GENERALVISITOR_H

#include <string>
#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"

#include "srcbuf.h" //srcbuf class interface 
#include "hilapp.h" //global vars needed  

/////////////////////////////////////////////////////////////////////////////////////
/// This implements a "general visitor", a placeholder to many utility functions 
/// which can be used to visit AST 
///
/// Example:  
///
/// class exampleVisitor : public GeneralVisitor, public RecursiveASTVisitor<exampleVisitor> {
/// 
///    // if you need a special constructor, do it like this:
///   public:
///    exampleVisitor(Rewriter & R, ASTContext *C, ... other vars... ) : GeneralVisitor(R,C) {
///       ... custom constructor stmts ...
///    }       
///    ... whatever definitions ...
///  
///    // if constructor is not needed, use the GeneralVisitor constructor:
///    using GeneralVisitor::GeneralVisitor;
///
/// };
/// 
/// exampleVisitor ev(TheRewriter,Context, .... );
///    // here TheRewriter is of type Rewriter &, Context of type ASTContext *  
///    // these are defined in MyASTVisitor, so this works directly there
/// ev.TraverseStmt( .. some Stmt * variable ..);   // or TraverseDecl or some other Traverse
///
/////////////////////////////////////////////////////////////////////////////////////
class GeneralVisitor {
protected:

  /// Reference to the Clang rewriter object
  Rewriter &TheRewriter;
  /// Store the context pointer
  ASTContext *Context;
  /// Store a printing policy. It is required quite often
  PrintingPolicy PP;
  
public:
  /// Construct with rewriter and context
  GeneralVisitor(Rewriter &R, ASTContext *C) : TheRewriter(R), PP(C->getLangOpts()) { 
    Context=C;
  }

  Rewriter & getRewriter() { return TheRewriter; }
  ASTContext * getASTContext() { return Context; }

  /// Report diagnostic info
  template <unsigned N>
  void reportDiag(DiagnosticsEngine::Level lev, const SourceLocation & SL,
                  const char (&msg)[N],
                  const char *s1 = nullptr,
                  const char *s2 = nullptr,
                  const char *s3 = nullptr) {
    // we'll do reporting only when output is on, avoid double reports
    auto & DE = Context->getDiagnostics();    
    auto ID = DE.getCustomDiagID(lev, msg );
    auto DB = DE.Report(SL, ID);
    if (s1 != nullptr) DB.AddString(s1);
    if (s2 != nullptr) DB.AddString(s2);
    if (s3 != nullptr) DB.AddString(s3);
  }

  /// Get the written expression of a C++ statement
  std::string get_stmt_str(const Stmt *s) {
    return TheRewriter.getRewrittenText(s->getSourceRange());
  }

  /// Get the type of an expression (i.e. int, double...) as a string
  std::string get_expr_type(const Expr *e) {
    // This is somehow needed for printing type without "class" id
    // TODO: perhaps reanalyse and make more robust?
    return e->getType().getUnqualifiedType().getAsString(PP);
  }


  /// a list of utility inspection functions
  /// getCanonicalType takes away typedefs, getUnqualifiedType() qualifiers, pp just in case
  /// the type string needs to begin with the string
  bool is_field_storage_expr(Expr *E) {
    return( E && 
      E->getType().getCanonicalType().getUnqualifiedType().getAsString(PP).find(field_storage_type) == 0);
  }

  /// Check the expression of a field expression
  bool is_field_expr(Expr *E) {
    return( E && 
      E->getType().getCanonicalType().getUnqualifiedType().getAsString(PP).find(field_type) == 0);
  }

  /// Check if declaration declate a field
  bool is_field_decl(ValueDecl *D) {
    return( D && 
      D->getType().getCanonicalType().getUnqualifiedType().getAsString(PP).find(field_type) == 0);
  }

  /// try to figure out whether expressions are duplicates
  bool is_duplicate_expr(const Expr * a, const Expr * b) {
    // Use the Profile function in clang, which "fingerprints"
    // statements
    llvm::FoldingSetNodeID IDa, IDb;
    a->Profile(IDa, *Context, true);
    b->Profile(IDb, *Context, true);
    return ( IDa == IDb );
  }

  /// Just check that the expression is of type parity
  bool is_parity_index_type(Expr *E) {
    return (get_expr_type(E) == "parity");
  }

  /// Checks if E is of type Field[parity]  parity=EVEN,ODD,ALL
  bool is_field_parity_expr(Expr *E) {
    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC &&
        strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
        is_field_expr(OC->getArg(0))) {

      return is_parity_index_type(OC->getArg(1));

    } else {
      // DON'T DO TEMPLATES NOW!  ONLY SPECIALIZATIONS
      #if 0
      // This is for templated expressions
      // for some reason, expr a[X] "getBase() gives X, getIdx() a...
      if (ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(E)) {
        Expr * lhs = ASE->getLHS()->IgnoreParens();
        
        if (is_field_expr(ASE->getLHS()->IgnoreParens())) {
          // llvm::errs() << " FP: and field\n";
          return is_parity_index_type(ASE->getRHS());
        }
      }
      #endif
    }
    return false;   
  }

  /// true if X, X+dir, X+offset -type expression
  bool is_X_type(Expr *E) {
    std::string s = get_expr_type(E);
    return (s == "X_index_type");
  }


  /// true if X, X+dir, X+offset -type expression
  bool is_X_index_type(Expr *E) {
    std::string s = get_expr_type(E);
    return (s == "X_index_type" || s == "X_plus_direction" || s == "X_plus_offset");
  }

  /// true if X+dir - type
  bool is_X_and_dir_type(Expr *E) {
    std::string s = get_expr_type(E);
    return (s == "X_plus_direction" || s == "X_plus_offset");
  }


  /// Checks if E is parity of a field (for example f[X]).
  /// Catches both parity and X_plus_direction 
  bool is_field_with_X_expr(Expr *E) {
    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC &&
        strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
        is_field_expr(OC->getArg(0))) {

      return is_X_index_type(OC->getArg(1));

    }
    return false;   
  }

  /// Checks if E is parity plus direction of a field (for example f[X+dir]).
  bool is_field_with_X_and_dir(Expr *E) {
    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC &&
        strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
        is_field_expr(OC->getArg(0))) {
      
      return is_X_and_dir_type(OC->getArg(1));

    }
    return false;   
  }

  /// is the stmt pointing now to a function call
  bool is_function_call_stmt(Stmt * s) {
    if (auto *Call = dyn_cast<CallExpr>(s)){
      // llvm::errs() << "Function call found: " << get_stmt_str(s) << '\n';
      return true;
    }
    return false;
  }

  /// is the stmt pointing now to a member call
  bool is_member_call_stmt(Stmt * s) {
    if (auto *Call = dyn_cast<CXXMemberCallExpr>(s)){
      // llvm::errs() << "Member call found: " << get_stmt_str(s) << '\n';
      return true;
    }
    return false;
  }

  /// is the stmt pointing now to a constructor
  bool is_constructor_stmt(Stmt * s) {
    if (auto *Call = dyn_cast<CXXConstructExpr>(s)){
      // llvm::errs() << "Constructor found: " << get_stmt_str(s) << '\n';
      return true;
    }
    return false;
  }

  /// is the stmt pointing now to a constructor
  bool is_user_cast_stmt(Stmt * s) {
    if (auto *ce = dyn_cast<ImplicitCastExpr>(s)) {
    // if (CastExpr *ce = dyn_cast<CastExpr>(s)) {
      if (ce->getCastKind() == CK_UserDefinedConversion) {
      // llvm::errs() << "Constructor found: " << get_stmt_str(s) << '\n';
        return true;
      }
    }
    return false;
  }


  /// Does the statement end with a semicolon
  bool isStmtWithSemicolon(Stmt * S) {
    SourceLocation l = Lexer::findLocationAfterToken(S->getEndLoc(),
                                                    tok::semi,
                                                    TheRewriter.getSourceMgr(),
                                                    Context->getLangOpts(),
                                                    false);
    if (l.isValid()) {
      //    llvm::errs() << "; found " << get_stmt_str(S) << '\n';
      return true;
    }
    return false;
  }

};


#endif // ifdef GENERALVISITOR_H