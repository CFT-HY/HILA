#include <sstream>
#include <iostream>
#include <string>

#include "toplevelvisitor.h"
#include "hilapp.h"

//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if a statement (inside site loop)
/// contains loop local variables, i.e. variables defined within the
/// loop body.  If not and expression is not site dependent, the expression is constant 
/// within the loop
//////////////////////////////////////////////////////////////////////////////

class isLoopLocalChecker : public GeneralVisitor,
                           public RecursiveASTVisitor<isLoopLocalChecker> {

  public:

    std::vector<var_info *> *loop_vars = nullptr;
    bool found_local_var;



    template <typename visitortype>
    isLoopLocalChecker(visitortype &v, std::vector<var_info *> *loop_var_list)
        : GeneralVisitor(v) {

        loop_vars = loop_var_list; // we do not clear the vector, because it may contain
                                   // earlier dependencies
        found_local_var = false;

    }

    // bool VisitStmt(Stmt *s) { llvm::errs() << "In stmt\n"; return true; }

    bool VisitDeclRefExpr(DeclRefExpr *e) {
        if (isa<VarDecl>(e->getDecl())) {
            VarDecl *decl = dyn_cast<VarDecl>(e->getDecl());
                // it's a variable and found the decl, check if we have it
                for (var_info &v : var_info_list)
                    if (v.decl == decl && v.is_loop_local) {

                        // found loop local variable - insert in list
                        found_local_var = true;
                        if (loop_vars != nullptr) {
                            for (var_info *rp : *loop_vars) {
                                // if it is already in list, return
                                if (rp->decl == v.decl) return true;
                            }
                            // now it was not in list, add to it
                            loop_vars->push_back(&v);
                            return true;
                        } else {
                            return false;  // might as well stop if not interested in list
                        }
                    }
        }
        return true;
    }
};

////////////////////////////////////////////////////////////////////////////////////
/// Check loop local variables in expr. If variables found, returns the optional
/// var_info pointers in loop_var.  
////////////////////////////////////////////////////////////////////////////////////

bool GeneralVisitor::contains_loop_local_var(Expr *e, std::vector<var_info *> *loop_var) {

    isLoopLocalChecker checker(*this, loop_var);

    checker.TraverseStmt(e);
    if (checker.found_local_var) return true;
    else return false;
}
