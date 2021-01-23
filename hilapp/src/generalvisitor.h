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
///    // these are defined in TopLevelVisitor, so this works directly there
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

  /// check if stmt contains random number generator
  bool contains_random(Stmt *s);


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
  bool is_duplicate_expr(const Expr * a, const Expr * b);

  /// Just check that the expression is of type parity
  bool is_parity_index_type(Expr *E) {
    return (get_expr_type(E) == "parity");
  }

  /// Checks if E is of type Field[parity]  parity=EVEN,ODD,ALL
  bool is_field_parity_expr(Expr *E);

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


  bool is_field_with_X_expr(Expr *E);

  /// Checks if E is parity plus direction of a field (for example f[X+dir]).
  bool is_field_with_X_and_dir(Expr *E);

  bool is_assignment_expr(Stmt * s, std::string * opcodestr, bool &iscompound);

  bool is_site_dependent(Expr * e, std::vector<var_info *> * dependent_var);

  bool is_rhs_site_dependent(Stmt *s, std::vector<var_info *> * vi);

  /// is the stmt pointing now to a function call
  bool is_function_call_stmt(Stmt * s);
  
  /// is the stmt pointing now to a member call
  bool is_member_call_stmt(Stmt * s);

  /// is the stmt pointing now to a constructor
  bool is_constructor_stmt(Stmt * s);

  /// is the stmt pointing now to a constructor
  bool is_user_cast_stmt(Stmt * s);

  /// Does the statement end with a semicolon
  bool isStmtWithSemicolon(Stmt * S);
  
  parity get_parity_val(const Expr *pExpr);

  /// utility used in inserting stuff after new line in buffer
  // SourceLocation getSourceLocationAtEndOfLine( SourceLocation l );
  /// another utility (cannot trust r.getEnd())
  SourceLocation getSourceLocationAtEndOfRange( SourceRange r );

  // get next char and loc, while skipping comments
  SourceLocation getNextLoc(SourceLocation sl,bool forward = true);

  char getChar(SourceLocation sl);
  SourceLocation findChar(SourceLocation sl, char ch);

  // get next word or symbol, if it is not a legal name symbol
  std::string getNextWord(SourceLocation sl, SourceLocation *end = nullptr);
  std::string getPreviousWord(SourceLocation sl, SourceLocation *start = nullptr);

  /// jump over following () expr
  SourceLocation skipParens( SourceLocation sl);


};


#endif // ifdef GENERALVISITOR_H