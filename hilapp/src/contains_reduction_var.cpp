#include <sstream>
#include <iostream>
#include <string>

#include "toplevelvisitor.h"
#include "hilapp.h"

//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if a statement (inside site loop)
/// contains reference to a special reduction variable, of type
/// Reduction<> or ReductionVector<>
//////////////////////////////////////////////////////////////////////////////

class containsReductionChecker : public GeneralVisitor,
                                 public RecursiveASTVisitor<containsReductionChecker> {

  public:
    bool found_reduction_var;


    template <typename visitortype>
    containsReductionChecker(visitortype &v) : GeneralVisitor(v) {

        found_reduction_var = false;
    }

    // bool VisitStmt(Stmt *s) { llvm::errs() << "In stmt\n"; return true; }

    bool VisitExpr(Expr *e) {
        std::string typ = e->getType().getUnqualifiedType().getNonReferenceType().getAsString();
         llvm::errs() << " REDUCTION CANDIDATE TYPE " << typ << '\n';
        if (typ.find("Reduction<") == 0 || typ.find("ReductionVector<") == 0) {
            found_reduction_var = true;
            return false; // can stop now
        }
        return true;
    }
};

////////////////////////////////////////////////////////////////////////////////////
/// Check loop local variables in expr. If variables found, returns the optional
/// var_info pointers in loop_var.
////////////////////////////////////////////////////////////////////////////////////

bool GeneralVisitor::contains_special_reduction_var(Expr *e) {

    containsReductionChecker checker(*this);

    checker.TraverseStmt(e);
    return checker.found_reduction_var;
}
