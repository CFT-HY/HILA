#include <sstream>
#include <iostream>
#include <string>

#include "toplevelvisitor.h"
#include "hilapp.h"

//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if a statement (inside site loop)
/// contains a random number generator call.  Goes inside functions as far as it can
//////////////////////////////////////////////////////////////////////////////

// this vector keeps track of FunctionDecls already checked in this search
// - prevents infinite loop in recursion
static std::vector<FunctionDecl *> fdecls;

class containsNovectorChecker : public GeneralVisitor,
                                public RecursiveASTVisitor<containsNovectorChecker> {

  public:
    bool found_novector;

    template <typename visitor_type>
    containsNovectorChecker(visitor_type &v) : GeneralVisitor(v) {
        found_novector = false;
    }

    // check function call exprs
    bool VisitStmt(Stmt *S) {
        if (CallExpr *CE = dyn_cast<CallExpr>(S)) {
            if (FunctionDecl *FD = CE->getDirectCallee()) {

                // llvm::errs() << "checking func " << FD->getNameAsString();
                if (has_pragma(FD, pragma_hila::NOVECTOR)) {
                    found_novector = true;
                    // llvm::errs() << " - has pragma novector\n";
                    return false; // stop here, found
                }

                if (FD->hasBody()) {
                    // some other function, go in it hierarchially
                    for (auto d : fdecls) {
                        if (d == FD) {
                            return true; // was checked, return and continue
                        }
                    }
                    // llvm::errs() << "NOVECTOR CHECK FUNCTION NAME: " << name << "\n";
                    fdecls.push_back(FD);

                    if (contains_novector(FD->getBody())) {
                        // llvm::errs() << "novector CHECK FUNCTION NAME: " << name << " HAS
                        // novector INSIDE\n";
                        found_novector = true;
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // need to do new "starter pack" here, because the TopLevelVisitor version not
    // callable here
    bool contains_novector(Stmt *s) {
         containsNovectorChecker chkagain(*this);
         chkagain.TraverseStmt(s);
         return (chkagain.found_novector);
    }
};

////////////////////////////////////////////////////////////////////////////////////
/// Check site dependence.  If return status is false, the optional vector
/// dependent_var will contains pointers to variables whose status may change
/// later.  These should be checked if needed.  The vector is unmodified if
/// return is true.
////////////////////////////////////////////////////////////////////////////////////

bool GeneralVisitor::contains_novector(Stmt *s) {

    containsNovectorChecker checker(*this);
    fdecls.clear();

    checker.TraverseStmt(s);

    fdecls.clear();
    return (checker.found_novector);
}
