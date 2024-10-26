#include <sstream>
#include <iostream>
#include <string>

#include "toplevelvisitor.h"
#include "hilapp.h"

//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if a statement (inside site loop)
/// is dependent on site (X).
/// Logic: find if it contains X (X_index_type), which takes care of field
/// and methods of X (check!), or random number generator.
///
/// No need to descent into functions, it is sufficient to scan the arguments
/// and check random
//////////////////////////////////////////////////////////////////////////////

class isSiteDependentChecker : public GeneralVisitor,
                               public RecursiveASTVisitor<isSiteDependentChecker> {

  public:
    bool found_X;
    bool found_var_depends_on_site;
    bool found_X_method;
    bool found_dependent_var;
    bool found_local_var;
    bool found_assign;
    bool check_local_var;
    std::vector<var_info *> *depends_on_var;

    template <typename visitortype>
    isSiteDependentChecker(visitortype &v, std::vector<var_info *> *dep_var, bool check_local)
        : GeneralVisitor(v) {

        depends_on_var = dep_var; // we do not clear the vector, because it may contain
                                  // earlier dependencies

        check_local_var = check_local; // Exit if loop local var is found

        found_X = found_var_depends_on_site = found_X_method = found_dependent_var =
            found_local_var = found_assign = false;
    }

    // bool VisitStmt(Stmt *s) { llvm::errs() << "In stmt\n"; return true; }

    bool VisitDeclRefExpr(DeclRefExpr *e) {
        /// if we visit X
        if (is_X_type(e)) {
            found_X = true;
            // llvm::errs() << "FOUND X index!\n";
            return false; // we do not need to go further, do we?

        } else {
            if (isa<VarDecl>(e->getDecl())) {
                VarDecl *decl = dyn_cast<VarDecl>(e->getDecl());
                // it's a variable and found the decl, check if we have it here
                bool found = false;
                for (var_info &v : var_info_list) {
                    if (v.decl == decl && !v.is_raw) {
                        found = true;
                        if (check_local_var && v.is_loop_local) {
                            found_local_var = true;
                            return false; // stop
                        }
                        if (v.is_site_dependent) {
                            found_var_depends_on_site = true;
                            return false; // again, the inspection can be stopped
                        }

                        // now it is possible that v later becomes site dependent
                        found_dependent_var = true;
                        if (depends_on_var != nullptr)
                            depends_on_var->push_back(&v);

                        // continue, may find more
                    }
                }
            }
        }
        return true;
    }

    bool VisitCXXMemberCallExpr(CXXMemberCallExpr *MCall) {
        // if we have an X-method, i.e. X.something() which is site dependent
        // it appears that the X DeclRefExpr finds the X-methods before this, but no
        // matter.
        if (is_X_type(MCall->getImplicitObjectArgument())) {
            found_X_method = true;
            return false;
        }
        return true;
    }

    // get assignments also
    bool VisitCXXOperatorCallExpr(CXXOperatorCallExpr *OpCall) {
        if (check_local_var) {
            if (OpCall->isAssignmentOp() ||
                OpCall->getOperator() == OverloadedOperatorKind::OO_PlusPlus ||
                OpCall->getOperator() == OverloadedOperatorKind::OO_MinusMinus) {

                found_assign = true;
                return false; // stop
            }
        }
        return true;
    }

    bool VisitBinaryOperator(BinaryOperator *O) {
        if (check_local_var && O->isAssignmentOp()) {
            found_assign = true;
            return false;
        }
        return true;
    }

    bool VisitUnaryOperator(UnaryOperator *UO) {
        if (check_local_var && (UO->isIncrementOp() || UO->isDecrementOp())) {
            found_assign = true;
            return false;
        }
        return true;
    }
};

////////////////////////////////////////////////////////////////////////////////////
/// Check site dependence.  If return status is false, the optional vector
/// dependent_var will contains pointers to variables whose status may change
/// later.  These should be checked if needed.  The vector is unmodified if
/// return is true.
////////////////////////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_site_dependent(Expr *e, std::vector<var_info *> *dependent_var) {

    isSiteDependentChecker checker(*this, dependent_var, false);

    checker.TraverseStmt(e);
    if (checker.found_X || checker.found_var_depends_on_site || checker.found_X_method) {
        return true;
    }

    // if nothing else found, check if there's rng generator
    return contains_random(e);
}


////////////////////////////////////////////////////////////////////////////////////
/// Check whether the expression is constant within the site loop.  This implies
/// it is not site dependent and does not depend on loop local variable.
/// This can be determined on a single pass.
/// Expr also cannot change the value of any variable.
////////////////////////////////////////////////////////////////////////////////////


bool GeneralVisitor::is_loop_constant(Expr *e) {

    // Init checker with loop local variable trap - we do not care about dependent_var list
    isSiteDependentChecker checker(*this, nullptr, true);

    checker.TraverseStmt(e);
    if (checker.found_X || checker.found_var_depends_on_site || checker.found_X_method ||
        checker.found_local_var || checker.found_assign) {
        return false;
    }

    // if nothing else found, check if there's rng generator - if so, not loop_constant
    return !contains_random(e);
}

////////////////////////////////////////////////////////////////////////////////////
/// Check whether the expression contains X or X method, i.e. ref to field var
////////////////////////////////////////////////////////////////////////////////////

bool GeneralVisitor::contains_field_ref(Expr *e) {

    // Init checker with loop local variable trap - we do not care about dependent_var list
    isSiteDependentChecker checker(*this, nullptr, true);

    checker.TraverseStmt(e);
    return (checker.found_X || checker.found_X_method);

}

