#include <sstream>
#include <iostream>
#include <string>

#include "myastvisitor.h"
#include "hilapp.h"


//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if a statement (inside site loop)
/// is dependent on site (X).
/// Logic: find if it contains X (X_index_type), which takes care of field 
/// and methods of X (check!), or random number generator.
//////////////////////////////////////////////////////////////////////////////

// this vector keeps track of FunctionDecls already checked in this search 
// - prevents infinite loop in recursion
static std::vector<FunctionDecl *> fdecls;

class containsRandomChecker : public GeneralVisitor, public RecursiveASTVisitor<containsRandomChecker> {

public:
  bool found_random;

  containsRandomChecker(Rewriter &R, ASTContext *C) : GeneralVisitor(R,C) {
    found_random = false;
  }

  // check random number calls 
  bool VisitStmt(Stmt *S) {
    if (CallExpr * CE = dyn_cast<CallExpr>(S)) {
      if (FunctionDecl * FD = CE->getDirectCallee()) {
        std::string name = FD->getNameInfo().getAsString();
        if (name == "hila_random") {
          found_random = true;
          return false;
        } else if (FD->hasBody()) {
          // some other function, go in it hierarchially
          for (auto d : fdecls) if (d == FD) {
            // was checked, return and continue
            return true;
          }
          // llvm::errs() << "RNG CHECK FUNCTION NAME: " << name << "\n";
          fdecls.push_back(FD);

          if (contains_random(FD->getBody())) {
            // llvm::errs() << "RNG CHECK FUNCTION NAME: " << name << " HAS RNG INSIDE\n";
            found_random = true;
            return false;
          }
        }
      }
    }
    return true;
  }

  // need to do new "starter pack" here, because the MyASTVisitor version not callable here
  bool contains_random(Stmt *s) {
    containsRandomChecker chkagain(TheRewriter,Context);
    chkagain.TraverseStmt(s);
    return (chkagain.found_random);
  }

};

////////////////////////////////////////////////////////////////////////////////////
/// Check site dependence.  If return status is false, the optional vector 
/// dependent_var will contains pointers to variables whose status may change
/// later.  These should be checked if needed.  The vector is unmodified if 
/// return is true.
////////////////////////////////////////////////////////////////////////////////////


bool MyASTVisitor::contains_random(Stmt * s) {

  containsRandomChecker checker(TheRewriter,Context);
  fdecls.clear();

  checker.TraverseStmt(s);

  fdecls.clear();
  return (checker.found_random);

}
