#include <sstream>
#include <iostream>
#include <string>

#include "myastvisitor.h"
#include "transformer.h"


//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if the function body contains a field loop
/// Logic: 
///   - Find if it contains X (X_index_type), which appears inside loops
///   - Find if it has field[parity] -stmt.  This can appear without X in
///      statemets like  f[EVEN] = 1;  etc.
///
//////////////////////////////////////////////////////////////////////////////

class containsFieldLoopChecker : public GeneralVisitor, public RecursiveASTVisitor<containsFieldLoopChecker> {

public:
  using GeneralVisitor::GeneralVisitor;   // use general visitor constructor

  bool found_X;
  bool found_field_parity;


  // bool VisitStmt(Stmt *s) { llvm::errs() << "In stmt\n"; return true; }

  bool VisitDeclRefExpr(DeclRefExpr * e) {
    /// if we visit X
    if (is_X_type(e)) {
      found_X = true;
      // llvm::errs() << "FOUND X index!\n";
      return false;  // we do not need to go further, do we?
    }
    return true;
  }

  bool VisitCXXOperatorCallExpr(CXXOperatorCallExpr * OC) {
    if (is_field_parity_expr(OC)) {
      // llvm::errs() << "FOUND FIELD WITH PARITY! \n";
      found_field_parity = true;
      return false;
    }
    return true;
  }
};

///////////////////////////////////////////////////////////////////////////////////
/// And loop checker interface here
///////////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::does_function_contain_loop(FunctionDecl * f) {

  if (f->hasBody()) {
    containsFieldLoopChecker flc(TheRewriter,Context);
    flc.TraverseStmt(f->getBody());
    return (flc.found_X || flc.found_field_parity);
  }
  return false;
}
