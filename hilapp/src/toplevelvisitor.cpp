#include "stringops.h"
#include "hilapp.h"
#include "toplevelvisitor.h"
#include "specialization_db.h"
#include <sstream>
#include <iostream>
#include <string>

/////////
/// Implementation of most toplevelvisitor methods
/////////

// define the ptr to the visitor here
TopLevelVisitor *g_TopLevelVisitor;

// function used for development
std::string print_TemplatedKind(const enum FunctionDecl::TemplatedKind kind) {
    switch (kind) {
    case FunctionDecl::TemplatedKind::TK_NonTemplate:
        return "TK_NonTemplate";
    case FunctionDecl::TemplatedKind::TK_FunctionTemplate:
        return "TK_FunctionTemplate";
    case FunctionDecl::TemplatedKind::TK_MemberSpecialization:
        return "TK_MemberSpecialization";
    case FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization:
        return "TK_FunctionTemplateSpecialization";
    case FunctionDecl::TemplatedKind::TK_DependentFunctionTemplateSpecialization:
        return "TK_DependentFunctionTemplateSpecialization";
    default:
        return "unknown";
    }
}

/// Check the validity a variable reference in a loop
bool FieldRefChecker::VisitDeclRefExpr(DeclRefExpr *e) {
    // It must be declared already. Get the declaration and check
    // the variable list. (If it's not in the list, it's not local)
    // llvm::errs() << "LPC variable reference: " <<  get_stmt_str(e) << "\n" ;
    for (var_info &vi : var_info_list)
        if (vi.is_loop_local && vi.decl == dyn_cast<VarDecl>(e->getDecl())) {
            // It is local! set status
            found_loop_local_var = true;
            vip = &vi;
            break;
        }

    return true;
}

// Walk the tree recursively
bool FieldRefChecker::TraverseStmt(Stmt *s) {
    RecursiveASTVisitor<FieldRefChecker>::TraverseStmt(s);
    return true;
}

// Walk the tree recursively
bool LoopAssignChecker::TraverseStmt(Stmt *s) {
    RecursiveASTVisitor<LoopAssignChecker>::TraverseStmt(s);
    return true;
}

/// Check the validity a variable reference in a loop
bool LoopAssignChecker::VisitDeclRefExpr(DeclRefExpr *e) {
    std::string type = e->getType().getAsString();
    type = remove_extra_whitespace(type);
    if (type.rfind("element<", 0) != std::string::npos) {
        reportDiag(DiagnosticsEngine::Level::Error, e->getSourceRange().getBegin(),
                   "Cannot assign a Field element to a non-element type");
    }
    return true;
}

/// Check if an assignment is allowed -- IS THIS NOW SUPERFLUOUS?
// void TopLevelVisitor::check_allowed_assignment(Stmt *s) {
//     if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
//         if (OP->getNumArgs() == 2) {
//             // Walk the right hand side to check for element types. None are allowed.
//             std::string type = OP->getArg(0)->getType().getAsString();
//             type = remove_extra_whitespace(type);
//             if (type.rfind("element<", 0) == std::string::npos) {

//                 LoopAssignChecker lac(*this);
//                 lac.TraverseStmt(OP->getArg(1));
//             } else {
//                 // llvm::errs() << " ** Element type : " << type << '\n';
//                 // llvm::errs() << " ** Canonical type without keywords: " <<
//                 // OP->getArg(0)->getType().getCanonicalType().getAsString(PP) <<
//                 '\n';
//             }
//         }
//     }
// }

/// -- Handler utility functions --

//////////////////////////////////////////////////////////////////////////////
/// Go through one field reference within parity loop and store relevant info
/// is_assign: assignment, is_compound: compound assign, is_X: argument is X,
/// is_func_arg: expression is a lvalue-argument (non-const. reference) to function
//////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::handle_field_X_expr(Expr *e, bool &is_assign, bool is_also_read,
                                          bool is_X, bool is_func_arg) {

    e = e->IgnoreParens();
    field_ref lfe;

    // we know here that Expr is of field-parity type
    if (CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(e)) {
        lfe.fullExpr = OC;
        // take name
        lfe.nameExpr = OC->getArg(0)->IgnoreImplicit();
        lfe.parityExpr = OC->getArg(1)->IgnoreImplicit();
    } else if (ArraySubscriptExpr *ASE = dyn_cast<ArraySubscriptExpr>(e)) {
        // In template definition TODO: should be removed?

        lfe.fullExpr = ASE;
        lfe.nameExpr = ASE->getLHS();
        lfe.parityExpr = ASE->getRHS();
        // llvm::errs() << lfe.fullExpr << " " << lfe.nameExpr << " " << lfe.parityExpr
        // <<
        // "\n";
    } else {
        llvm::errs() << "Should not happen! Error in Field parity\n";
        exit(1);
    }

    if (is_assign && lfe.nameExpr->getType().isConstQualified()) {
        // llvm::errs() << " ******** CANNOT ASSIGN TO A CONST QUAL VAR, NAME " <<
        // get_stmt_str(lfe.nameExpr) << '\n';
        is_assign = false;
    }

    // Check if the expression is already handled
    for (field_ref r : field_ref_list)
        if (r.fullExpr == lfe.fullExpr) {
            return (true);
        }

    // lfe.nameInd    = writeBuf->markExpr(lfe.nameExpr);
    // lfe.parityInd  = writeBuf->markExpr(lfe.parityExpr);

    lfe.is_written = is_assign;
    lfe.is_read = (is_also_read || !is_assign);
    lfe.sequence = parsing_state.stmt_sequence;

    if (is_assign && !lfe.nameExpr->isLValue()) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   lfe.nameExpr->getSourceRange().getBegin(),
                   "Cannot assign to non-lvalue Field expression");
    }

    std::string parity_expr_type = get_expr_type(lfe.parityExpr);

    if (parity_expr_type == "Parity") {
        if (is_X) {
            llvm::errs() << "Internal error in handle_loop_parity\n";
            exit(1);
        }
        if (parsing_state.accept_field_parity) {
            // 1st parity statement on a single line lattice loop
            loop_info.parity_expr = lfe.parityExpr;
            loop_info.parity_value = get_parity_val(loop_info.parity_expr);
            loop_info.parity_text = get_stmt_str(loop_info.parity_expr);
        } else {
            reportDiag(DiagnosticsEngine::Level::Error,
                       lfe.parityExpr->getSourceRange().getBegin(),
                       "Field[Parity] not allowed here, use Field[X] -type instead");
        }
    }

    // next ref must have wildcard parity
    parsing_state.accept_field_parity = false;

    if (parity_expr_type == "X_plus_direction" || parity_expr_type == "X_plus_offset") {

        if (is_assign && !is_func_arg) {
            reportDiag(
                DiagnosticsEngine::Level::Error,
                lfe.parityExpr->getSourceRange().getBegin(),
                "Cannot assign to Field expression with [X + dir] -type argument.");
        }
        if (is_assign && is_func_arg) {
            reportDiag(
                DiagnosticsEngine::Level::Error,
                lfe.parityExpr->getSourceRange().getBegin(),
                "Cannot use a non-const. reference to Field expression with [X + "
                "dir] -type argument.");
        }

        // Now need to split the expr to parity and dir-bits
        // Because of offsets this is pretty complicated to do in AST.
        // We now know that the expr is of type
        // [X+Direction]  or  [X+CoordinateVector] -- just
        // use the textual form of the expression!

        bool has_X;
        lfe.direxpr_s = remove_X(get_stmt_str(lfe.parityExpr), &has_X);

        if (!has_X) {
            reportDiag(DiagnosticsEngine::Level::Fatal,
                       lfe.parityExpr->getSourceRange().getBegin(),
                       "Internal error: index should have been X");
            exit(1);
        }

        // llvm::errs() << "Direxpr " << lfe.direxpr_s << '\n';

        lfe.is_direction = true;

        if (parity_expr_type == "X_plus_offset") {

            // It's an offset, no checking here to be done
            lfe.is_offset = true;

            FieldRefChecker frc(*this);
            frc.TraverseStmt(lfe.parityExpr);
            if (frc.isLoopLocal()) {
                reportDiag(DiagnosticsEngine::Level::Error,
                           lfe.parityExpr->getSourceRange().getBegin(),
                           "Non-nearest neighbour reference cannot depend on variable "
                           "'%0' defined inside site loop",
                           frc.getLocalVarInfo()->name.c_str());
            }

        } else {

            // Now make a check if the reference is just constant (e_x etc.)
            // Need to descend quite deeply into the expr chain
            Expr *e = lfe.parityExpr->IgnoreParens()->IgnoreImplicit();
            CXXOperatorCallExpr *Op = dyn_cast<CXXOperatorCallExpr>(e);
            if (!Op) {
                if (CXXConstructExpr *Ce = dyn_cast<CXXConstructExpr>(e)) {
                    // llvm::errs() << " ---- got Ce, args " << Ce->getNumArgs() <<
                    // '\n';
                    if (Ce->getNumArgs() == 1) {
                        e = Ce->getArg(0)->IgnoreImplicit();
                        Op = dyn_cast<CXXOperatorCallExpr>(e);
                    }
                }
            }

            if (!Op) {
                reportDiag(
                    DiagnosticsEngine::Level::Fatal,
                    lfe.parityExpr->getSourceRange().getBegin(),
                    "Internal error: could not parse X + Direction/offset -statement");
                exit(1);
            }

            Expr *dirE = Op->getArg(1)->IgnoreImplicit();
            if (dirE->isIntegerConstantExpr(*Context)) {
                llvm::APSInt result;

                // Got constant -- interface changes in clang 12 or 13(!)
#if defined(__clang_major__) && (__clang_major__ <= 11)
                dirE->isIntegerConstantExpr(result, *Context);
#else
                auto res = dirE->getIntegerConstantExpr(*Context);
                result = res.getValue();
#endif

                // Op must be + or - --get the sign
                const char *ops = getOperatorSpelling(Op->getOperator());

                int offset = 0;
                if (strcmp(ops, "+") == 0)
                    offset = 0;
                else if (strcmp(ops, "-") == 0)
                    offset = 50;
                else {
                    llvm::errs() << "This cannot happen, direction op " << ops << '\n';
                    exit(1);
                }

                lfe.is_constant_direction = true;
                // constant_value is used to uniquely label directions.
                // Give negative dirs an offset of 50.
                lfe.constant_value = result.getExtValue() + offset;


                // llvm::errs() << " GOT DIR CONST, value " << lfe.constant_value << "
                // expr " << lfe.direxpr_s << '\n';
            } else {
                lfe.is_constant_direction = false;
                // llvm::errs() << "GOT DIR NOT-CONST " << lfe.direxpr_s << '\n';

                // If the Direction is a variable, add it to the list
                // DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(lfe.dirExpr);
                // static std::string assignop;
                // if(DRE && isa<VarDecl>(DRE->getDecl())) {
                //   handle_var_ref(DRE, false, assignop);
                // }

                // traverse the dir-expression to find var-references etc.
                is_assign = false;
                TraverseStmt(lfe.parityExpr);

                // do it again with fieldrefchecker to see if the dir depends on
                // internal var
                FieldRefChecker frc(*this);
                frc.TraverseStmt(lfe.parityExpr);
                if (frc.isLoopLocal()) {
                    lfe.is_loop_local_dir = true;
                }
            }
        }
    } // end of "Direction"-branch

    // llvm::errs() << "field expr " << get_stmt_str(lfe.nameExpr)
    //              << " Parity " << get_stmt_str(lfe.parityExpr)
    //              << "\n";

    // Check that there are no local variable references up the AST of the name
    FieldRefChecker frc(*this);
    frc.TraverseStmt(lfe.nameExpr);
    if (frc.isLoopLocal()) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   lfe.nameExpr->getSourceRange().getBegin(),
                   "Field reference cannot depend on loop-local variable '%0'",
                   frc.getLocalVarInfo()->name.c_str());
    }

    if (contains_random(lfe.fullExpr)) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   lfe.fullExpr->getSourceRange().getBegin(),
                   "Field reference cannot contain random number generator");
    }

    field_ref_list.push_back(lfe);

    return (true);
}

//////////////////////////////////////////////////////////////////

///  Utility to find the reduction type

reduction get_reduction_type(bool is_assign, const std::string &assignop,
                             var_info &vi) {
    if (is_assign && (!vi.is_loop_local)) {
        if (assignop == "+=")
            return reduction::SUM;
        if (assignop == "*=")
            return reduction::PRODUCT;
    }
    return reduction::NONE;
}

///////////////////////////////////////////////////////////////////
/// Find the the base of a compound variable expression
///////////////////////////////////////////////////////////////////

DeclRefExpr *TopLevelVisitor::find_base_variable(Expr *E) {
    Expr *RE = E;

    while (!dyn_cast<DeclRefExpr>(RE)) {
        // RE may be a compound expression. We want the base variable.
        if (ArraySubscriptExpr *ASE = dyn_cast<ArraySubscriptExpr>(RE)) {
            RE = ASE->getBase()->IgnoreImplicit();
        } else if (MemberExpr *ME = dyn_cast<MemberExpr>(RE)) {
            RE = ME->getBase()->IgnoreImplicit();
        } else if (CXXOperatorCallExpr *OCE = dyn_cast<CXXOperatorCallExpr>(RE)) {
            if (strcmp(getOperatorSpelling(OCE->getOperator()), "[]") == 0) {
                RE = OCE->getArg(0)->IgnoreImplicit();
            } else {
                // It's not a variable
                return nullptr;
            }
        } else {
            // It's not a variable
            return nullptr;
        }
    }
    return dyn_cast<DeclRefExpr>(RE);
}

///////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::is_variable_loop_local(VarDecl *decl) {
    for (var_decl &d : var_decl_list) {
        if (d.scope >= 0 && decl == d.decl) {
            return true;
        }
    }
    return false;
}

/// handle an array subscript expression.  Operations depend
/// on whether the array is defined loop extern or intern:
///  - if intern, nothing special needs to be done
///  - if extern, check index:
///     - if site dependent, mark loop_info.has_site_dep_cond_or_index
///     - if contains loop local var ref or is site dependent:
///       whole array is input to loop: need to know array size.
///       If size not knowable, flag error
///     - If index is not site dependent or loop local:
///          -  it is sufficient to read in this array element only,
///             and remove var references to variables in index
///
int TopLevelVisitor::handle_bracket_var_ref(bracket_ref_t &ref, array_ref::reftype type,
                                            bool &is_assign, std::string &assignop) {

    if (ref.DRE == nullptr) {

        reportDiag(DiagnosticsEngine::Level::Warning,
                   ref.E->getSourceRange().getBegin(),
                   "array brackets '[]' applied to a non-variable, and hilapp does not "
                   "(yet) have logic to analyse this."
                   " Fingers crossed ...");

        return 0;
    }

    // if this has #pragma direct_access don't do anything
    if (loop_info.has_pragma_access &&
        find_word(loop_info.pragma_access_args, get_stmt_str(ref.DRE)) !=
            std::string::npos) {

        // no need to handle the expression here, bubble up and handle the raw ptr on
        // next visit
        return 0;
    }

    VarDecl *vd;
    DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(ref.DRE);
    if (DRE == nullptr) {
        llvm::errs() << "hilapp error: array refs to member variables not yet "
                        "implemented - probably does not work\n";
        llvm::errs() << "Expression: " << get_stmt_str(ref.E) << '\n';
        return 0;
    }
    vd = dyn_cast<VarDecl>(DRE->getDecl());

    // Check if it's local
    if (is_variable_loop_local(vd)) {
        if (type == array_ref::ARRAY) {

            // if an array is declared within the loop, nothing special needs to
            // be done.  Let us anyway put it through the std handler in order to
            // flag it.

            handle_var_ref(DRE, is_assign, assignop);

            // and traverse whatever is in the index in normal fashion
            is_assign = false; // we're not assigning to the index
            TraverseStmt(ref.Idx);
            parsing_state.skip_children = 1; // it's handled now
            return 1;

        } else {

            reportDiag(DiagnosticsEngine::Level::Error,
                       ref.E->getSourceRange().getBegin(),
                       "Cannot define this variable inside site loop");
            parsing_state.skip_children = 1;
            return 1;
        }
    }

    // Now array is declared outside the loop

    if (is_assign && type != array_ref::REDUCTION) {
        reportDiag(DiagnosticsEngine::Level::Error, ref.E->getSourceRange().getBegin(),
                   "Cannot assign to an array, std::vector or std::array here.  Use "
                   "VectorReduction if reduction is needed.");
        return 1;
    }

    // If index is site dep or contains local variables, we need
    // the whole array in the loop

    bool site_dep = is_site_dependent(ref.Idx, &loop_info.conditional_vars);

    // if it is assignment = reduction, don't vectorize
    if (site_dep || is_assign)
        loop_info.has_site_dependent_cond_or_index = true;

    reduction reduction_type;
    if (is_assign) {
        if (assignop == "+=")
            reduction_type = reduction::SUM;
        else
            reduction_type = reduction::PRODUCT;

    } else {
        reduction_type = reduction::NONE;
    }

    // has this array been referred to already?
    // if so, mark and return
    for (array_ref &ar : array_ref_list) {
        if (vd == ar.vd && ar.type != array_ref::REPLACE) {

            if ((ar.type == array_ref::REDUCTION) ^ (type == array_ref::REDUCTION)) {
                reportDiag(
                    DiagnosticsEngine::Level::Error, ref.E->getSourceRange().getBegin(),
                    "VectorReduction cannot be used on RHS and LHS simultaneously.");
                parsing_state.skip_children = 1;
                return 1;
            }

            if (ar.type == array_ref::REDUCTION &&
                ar.reduction_type != reduction_type) {
                reportDiag(DiagnosticsEngine::Level::Error,
                           ref.E->getSourceRange().getBegin(),
                           "Cannot use '+=' and '*=' reduction to the same variable "
                           "simultaneously.");
                parsing_state.skip_children = 1;
                return 1;
            }

            ar.refs.push_back(ref);
            is_assign = false;     // reset assign for index traversal
            TraverseStmt(ref.Idx); // need to traverse index normally
            parsing_state.skip_children = 1;
            return 1;
        }
    }

    // now it is a new array ref
    array_ref ar;
    ar.refs.push_back(ref);
    ar.vd = vd;
    ar.name = vd->getNameAsString();

    // get type of the element of the array
    ar.element_type = ref.E->getType().getCanonicalType().getAsString(PP);

    // and save the type of reduction
    ar.reduction_type = reduction_type;

    if (!site_dep && !contains_loop_local_var(ref.Idx, nullptr) && !is_assign) {

        // now it is a fixed index - move whole arr ref outside the loop
        // Note: multiple refrences are not checked, thus, same element can be
        // referred more than once.  TODO? (small optimization)

        // use the array_ref_list to note this
        ar.type = array_ref::REPLACE;
        array_ref_list.push_back(ar);

        parsing_state.skip_children = 1; // no need to look inside the replacement
                                         // variable refs inside should not be recorded
        return 1;
    }

    // Now there is site dep/loop local stuff in index.  Whole array has to be taken
    // "in".

    ar.type = type;

    if (type == array_ref::ARRAY) {
        const ConstantArrayType *cat = Context->getAsConstantArrayType(vd->getType());
        if (cat) {
            llvm::errs() << " %%% Found constant array type expr, size "
                         << cat->getSize() << "\n";

            ar.size = cat->getSize().getZExtValue(); // get (extended) size value
            ar.size_expr = std::to_string(ar.size);
            ar.data_ptr = ar.name;

        } else {
            // Do not accept arrays with variable size!  The size expression does not
            // nacessarily evaluate correctly due to problems with C(++) arrays
            // (TODO: use pragma?)
            //
            // const VariableArrayType *vat =
            // Context->getAsVariableArrayType(vd->getType()); if (vat) {
            //     llvm::errs() << " %%% Found variable array type expr, size "
            //                  << get_stmt_str(vat->getSizeExpr()) << "\n";

            //     ar.size = 0;
            //     ar.size_expr = vat->getSizeExpr();

            // } else {

            // Now different array type - flag as error
            reportDiag(DiagnosticsEngine::Level::Error,
                       ref.E->getSourceRange().getBegin(),
                       "Array size is unknown - recommend using Vector<>, "
                       "std::array<> or std::vector<> instead");

            parsing_state.skip_children = 1;
            return 1;
        }
    } else {

        std::string typestr = vd->getType().getCanonicalType().getAsString(PP);

        if (type == array_ref::STD_ARRAY) {

            // find the last arg in template

            int i = typestr.rfind('>');
            int j = typestr.rfind(',', i);
            ar.size = std::stoi(typestr.substr(j + 1, i - j - 1));
            ar.size_expr = std::to_string(ar.size);
            ar.data_ptr = ar.name + ".data()";

        } else {

            // Now it must be std::vector<> or VectorReduction.  Size is dynamic

            ar.size = 0;
            ar.size_expr = ar.name + ".size()";
            ar.data_ptr = ar.name + ".data()";

            if (type == array_ref::REDUCTION)
                ar.reduction_type = reduction_type;
        }
    }

    array_ref_list.push_back(ar);

    // traverse whatever is in the index in normal fashion
    is_assign = false; // not assign to index
    TraverseStmt(ref.Idx);

    parsing_state.skip_children = 1;
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int TopLevelVisitor::handle_array_var_ref(ArraySubscriptExpr *ASE, bool &is_assign,
                                          std::string &assignop) {

    bracket_ref_t br = {ASE, find_base_variable(ASE), ASE->getIdx()};
    return handle_bracket_var_ref(br, array_ref::ARRAY, is_assign, assignop);
}

///////////////////////////////////////////////////////////////////////////////
/// is_vector_reference() true if expression is of type var[], where
/// var is std::vector<> or std::array<>
///////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::is_vector_reference(Stmt *s) {
    Expr *E = dyn_cast<Expr>(s);
    if (E) {
        E = E->IgnoreParens();
        CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);
        if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0) {

            // Arg(1) is the "root"
            std::string type =
                OC->getArg(0)->getType().getCanonicalType().getAsString(PP);

            if (type.find("std::vector<") == 0 || type.find("std::array<") == 0 ||
                type.find("VectorReduction") == 0)
                return true;
        }
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
/// handle_vector_reference() processes references like v[index], where
/// v is std::vector or std::array.  Called only if is_vector_reference() is true
///////////////////////////////////////////////////////////////////////////////
bool TopLevelVisitor::handle_vector_reference(Stmt *s, bool &is_assign,
                                              std::string &assignop,
                                              Stmt *assign_stmt) {

    bracket_ref_t br;

    br.E = dyn_cast<Expr>(s);
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(br.E);

    br.DRE = find_base_variable(br.E);
    br.Idx = OC->getArg(1)->IgnoreImplicit();

    if (is_assign)
        br.assign_stmt = assign_stmt;

    std::string type = OC->getArg(0)->getType().getCanonicalType().getAsString(PP);

    array_ref::reftype rt;
    if (type.find("std::vector<") == 0)
        rt = array_ref::STD_VECTOR;
    else if (type.find("std::array<") == 0)
        rt = array_ref::STD_ARRAY;
    else {
        // now we know it is ArrayReduction.  Treat is as reduction only if it is
        // assignment
        if (is_assign)
            rt = array_ref::REDUCTION;
        else
            rt = array_ref::STD_ARRAY;
    }

    handle_bracket_var_ref(br, rt, is_assign, assignop);

    return true;
}

///////////////////////////////////////////////////////////////////////////////
///  Handle constant expressions referred to in loops
///  Const reference is left as is, EXCEPT if:
///   - target is kernelized
///   - const is defined outside loop body
///  If this is true, const ref is substituted with the const value
///  This may cause problems if types are different: e.g. const is enum value,
///  and it is substituted with an int literal.  Try to help with type casting!
///////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::handle_constant_ref(Expr *E) {

    APValue val;
    if (!E->isCXX11ConstantExpr(*Context, &val, nullptr))
        return false; // nothing

    // no need to do anything if not kernelized
    if (!target.kernelize)
        return true;

    E = E->IgnoreImplicit();
    // If it is not a declrefexpr continue to next node in ast
    DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(E);
    if (DRE == nullptr)
        return true;

    SourceLocation sl = DRE->getDecl()->getSourceRange().getBegin();
    if (sl.isValid() && get_FileId(sl) == get_FileId(global.location.loop)) {

        // if const is defined on file scope earlier or
        // const defined within the loop
        if (sl > global.location.top ||
            (sl < global.location.loop && sl > E->getBeginLoc()))
            // leave as is
            parsing_state.skip_children = 1;
        return true;

        // what is the type of the const?
        QualType ty = DRE->getType().getCanonicalType();
        std::string typestr = ty.getAsString();

        // replace int const expr by the value
        srcBuf *buf = get_file_srcBuf(DRE->getBeginLoc());
        if (val.isInt()) {
            std::string repl = std::to_string(val.getInt().getExtValue());
            if (typestr != "int")
                repl = "(" + typestr + ")" + repl;
            buf->replace(DRE->getSourceRange(), repl);

        } else if (val.isFloat()) {
            buf->replace(DRE->getSourceRange(),
                         std::to_string(val.getFloat().convertToDouble()));
        } else {
            // Not int or float, treat as var?
            return false;
        }
    }
    parsing_state.skip_children = 1;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
/// handle_constant_expr_ref handles (non-variable) expressions which are site
/// loop constants.  This includes struct/class members (and array refs?)
///////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::handle_loop_const_expr_ref(Expr *E) {

    // First, get the string rep of the expression
    std::string expstr = remove_all_whitespace(get_stmt_str(E));

    // Did we already have it?
    for (loop_const_expr_ref &cer : loop_const_expr_ref_list) {
        if (cer.exprstring == expstr) {
            cer.refs.push_back(E);
            return;
        }
    }

    // New ref
    loop_const_expr_ref eref;

    eref.refs.push_back(E);
    eref.exprstring = expstr;

    // Get the type of the expr
    clang::QualType typ =
        E->getType().getUnqualifiedType().getCanonicalType().getNonReferenceType();
    typ.removeLocalConst();
    eref.type = typ.getAsString(PP);

    loop_const_expr_ref_list.push_back(eref);
}

///////////////////////////////////////////////////////////////////////////////
/// handle_full_loop_stmt() is the starting point for the analysis of all
/// "parity" -loops
///////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::handle_full_loop_stmt(Stmt *ls, bool field_parity_ok) {
    // init edit buffer
    // Buf.create( &TheRewriter, ls );

    field_ref_list.clear();
    special_function_call_list.clear();
    var_info_list.clear();
    var_decl_list.clear();
    array_ref_list.clear();
    loop_const_expr_ref_list.clear();
    loop_function_calls.clear();

    global.location.loop = ls->getSourceRange().getBegin();
    loop_info.clear_except_external();
    loop_info.range = ls->getSourceRange();
    parsing_state.accept_field_parity = field_parity_ok;

    // the following is for taking the parity from next elem
    parsing_state.scope_level = 0;
    parsing_state.in_loop_body = true;
    parsing_state.ast_depth = 0; // renormalize to the beginning of loop
    parsing_state.stmt_sequence = 0;

    // code analysis starts here
    TraverseStmt(ls);

    parsing_state.in_loop_body = false;
    parsing_state.ast_depth = 0;

    // check and analyze the field expressions
    check_var_info_list();
    check_addrofops_and_refs(ls); // scan through the full loop again
    check_field_ref_list();
    process_loop_functions(); // revisit functions when vars are fully resolved

    // check here also if conditionals are site dependent through var dependence
    // because var_info_list was checked above, once is enough
    if (loop_info.has_site_dependent_cond_or_index == false) {
        for (auto *n : loop_info.conditional_vars)
            if (n->is_site_dependent)
                loop_info.has_site_dependent_cond_or_index = true;
    }

    // if (loop_info.has_site_dependent_conditional) llvm::errs() << "Cond is site
    // dep!\n";

    // and now generate the appropriate code

    generate_code(ls);

    // Buf.clear();

    // Emit the original command as a commented line
    writeBuf->insert(get_real_range(ls->getSourceRange()).getBegin(),
                     comment_string(global.full_loop_text) + "\n", true, true);

    global.full_loop_text = "";

    // don't go again through the arguments
    parsing_state.skip_children = 1;

    return true;
}

////////////////////////////////////////////////////////////////////////////////
///  act on statements within the parity loops.  This is called
///  from VisitStmt() if the status state::in_loop_body is true
////////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::handle_loop_body_stmt(Stmt *s) {

    // This keeps track of the assignment to field
    // must remember the set value across calls
    static bool is_assignment = false;
    static bool is_compound = false;
    static Stmt *assign_stmt = nullptr;
    static std::string assignop;
    static bool is_field_assign = false;

    // depth = 1 is the "top level" statement, should give fully formed
    // c++ statements separated by ';'.  These act as sequencing points
    // This is used to obtain assignment and read ordering
    if (parsing_state.ast_depth == 1)
        parsing_state.stmt_sequence++;

    // Need to recognize assignments lf[X] =  or lf[X] += etc.
    // And also assignments to other vars: t += norm2(lf[X]) etc.
    Expr *assignee;
    if (is_assignment_expr(s, &assignop, is_compound, &assignee)) {

        // check_allowed_assignment(s);
        assign_stmt = s;
        is_assignment = true;

        // Need to mark/handle the assignment method if necessary
        if (is_constructor_stmt(s)) {
            handle_constructor_in_loop(s);
        } else if (is_function_call_stmt(s)) {
            handle_function_call_in_loop(s);
        }

        // next visit here will be to the assigned to variable
        // Check type, if it is f[X] or f[parity] - then we do not (necessarily)
        // need to read in the variable
        is_field_assign =
            (is_field_parity_expr(assignee) || is_field_with_X_expr(assignee));

        return true;
    }

    // if (isa<MemberExpr>(s)) {
    //     if (is_assignment) is_member_expr = true;
    // }

    if (is_constructor_stmt(s)) {
        handle_constructor_in_loop(s);
        // return true;
    }

    // Check for std::vector or std::array references
    if (is_vector_reference(s)) {
        handle_vector_reference(s, is_assignment, assignop, assign_stmt);
        parsing_state.skip_children = 1;
        is_assignment = false;
        return true;
    }

    // Check c++ methods  -- HMM: it seems above function call stmt catches these first
    if (0 && is_member_call_stmt(s)) {
        handle_member_call_in_loop(s);
        // let this ripple trough, for now ...
        // return true;
    }

    // Check for function calls parameters. We need to determine if the
    // function can assign to the a field parameter (is not const).
    if (is_function_call_stmt(s)) {
        handle_function_call_in_loop(s);
        // let this ripple trough, for - expr f[X] is a function call and is trapped
        // below too
        // return true;
    }

    if (is_user_cast_stmt(s)) {
        // llvm::errs() << "GOT USER CAST " << get_stmt_str(s) << '\n';
    }

    // catch then expressions
    if (Expr *E = dyn_cast<Expr>(s)) {

        // Avoid treating constexprs as variables
        if (handle_constant_ref(E)) {
            return true;
        }

        // if (is_field_element_expr(E)) {
        // run this expr type up until we find field variable refs
        if (is_field_with_X_expr(E)) {
            // It is Field[X] reference
            // get the expression for field name
            handle_field_X_expr(E, is_assignment, is_compound || !is_field_assign,
                                true);
            is_assignment = false; // next will not be assignment
            // (unless it is a[] = b[] = c[], which is OK)

            parsing_state.skip_children = 1;
            return true;
        }

        if (is_field_parity_expr(E)) {
            // Now we know it is a field parity reference
            // get the expression for field name

            handle_field_X_expr(E, is_assignment, is_compound || !is_field_assign,
                                false);
            is_assignment = false; // next will not be assignment
            // (unless it is a[] = b[] = c[], which is OK)

            parsing_state.skip_children = 1;
            return true;
        }

        if (is_field_expr(E)) {
            // field without [X], bad usually (TODO: allow  scalar func(field)-type?)
            reportDiag(DiagnosticsEngine::Level::Error, E->getSourceRange().getBegin(),
                       "Field expressions without [X] not allowed within site loop");
            parsing_state.skip_children = 1; // once is enough
            return true;
        }

        // check in general array-type access ops - may prevent vectorization
        // if index is site dependent
        if (is_site_dependent_access_op(E)) {
            loop_info.has_site_dependent_cond_or_index = true;
        }

        // if (UnaryOperator * UO = dyn_cast<UnaryOperator>(E)) {
        //   if (UO->getOpcode() == UnaryOperatorKind::UO_AddrOf &&
        //       does_expr_contain_field( UO->getSubExpr() ) ) {
        //     reportDiag(DiagnosticsEngine::Level::Error,
        //                E->getSourceRange().getBegin(),
        //                "Taking address of '%0' is not allowed, suggest using
        //                references. " "If a pointer is necessary, copy first: 'auto v
        //                = %1; auto *p = &v;'", get_stmt_str(UO->getSubExpr()).c_str(),
        //                get_stmt_str(UO->getSubExpr()).c_str() );

        //     parsing_state.skip_children = 1;  // once is enough
        //     return true;
        //   }
        // }

        if (DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(E)) {
            if (isa<VarDecl>(DRE->getDecl())) {
                // now it should be var ref non-field

                // check if it is raw ptr access, mark if so
                if (loop_info.has_pragma_access &&
                    find_word(loop_info.pragma_access_args,
                              DRE->getDecl()->getNameAsString()) != std::string::npos) {

                    handle_var_ref(DRE, is_assignment, assignop, assign_stmt, true);

                } else {

                    handle_var_ref(DRE, is_assignment, assignop, assign_stmt);
                }

                is_assignment = false;

                // llvm::errs() << "Variable ref: "
                //              << TheRewriter.getRewrittenText(E->getSourceRange()) <<
                //              '\n';

                parsing_state.skip_children = 1;
                return true;
            }
            // TODO: function ref?
        }

        if (MemberExpr *ME = dyn_cast<MemberExpr>(E)) {

            // nothing to do if not loop const (PROBABLY?)
            // Check also the type: if it is not trivial, don't know what to do here
            if (!is_assignment && is_loop_constant(E) &&
                ME->getType().isTrivialType(*Context)) {

                handle_loop_const_expr_ref(E);

                parsing_state.skip_children = 1;
                return true;
            }
        }

        if (isa<ArraySubscriptExpr>(E)) {
            // llvm::errs() << "  It's array expr "
            //              << TheRewriter.getRewrittenText(E->getSourceRange()) <<
            //              "\n";
            auto a = dyn_cast<ArraySubscriptExpr>(E);

            // At this point this should be an allowed expression?
            // llvm::errs() << " ARRAY EXPRS " << get_stmt_str(E) << '\n';
            int is_handled = handle_array_var_ref(a, is_assignment, assignop);

            is_assignment = false;
            // We don't want to handle the array variable or the index separately
            parsing_state.skip_children = is_handled;
            return true;
        }

        if (0) {

            // not field type non-const expr
            // llvm::errs() << "Non-const other Expr: " << get_stmt_str(E) << '\n';
            // loop-local variable refs inside? If so, we cannot evaluate this as
            // "whole"

            // check_local_loop_var_refs = 1;

            // TODO: find really uniq variable references
            // var_ref_list.push_back( handle_var_ref(E) );

            parsing_state.skip_children = 1;
            return true;
        }
    } // Expr checking branch - now others...

    // This reached only if s is not Expr

    // start {...} -block or other compound
    if (isa<CompoundStmt>(s) || isa<ForStmt>(s) || isa<IfStmt>(s) ||
        isa<WhileStmt>(s) || isa<DoStmt>(s) || isa<SwitchStmt>(s) ||
        isa<ConditionalOperator>(s)) {

        static bool passthrough = false;
        // traverse each stmt - use passthrough trick if needed
        if (passthrough) {
            passthrough = false;
            return true;
        }

        parsing_state.scope_level++;
        passthrough = true; // next visit will be to the same node, skip

        // Reset ast_depth, so that depth == 0 again for the block.
        if (isa<CompoundStmt>(s))
            parsing_state.ast_depth = -1;

        TraverseStmt(s);

        // check also the conditionals - are these site dependent?
        if (!loop_info.has_site_dependent_cond_or_index) {
            Expr *condexpr = nullptr;
            if (IfStmt *IS = dyn_cast<IfStmt>(s))
                condexpr = IS->getCond();
            else if (ForStmt *FS = dyn_cast<ForStmt>(s))
                condexpr = FS->getCond();
            else if (WhileStmt *WS = dyn_cast<WhileStmt>(s))
                condexpr = WS->getCond();
            else if (DoStmt *DS = dyn_cast<DoStmt>(s))
                condexpr = DS->getCond();
            else if (SwitchStmt *SS = dyn_cast<SwitchStmt>(s))
                condexpr = SS->getCond();
            else if (ConditionalOperator *CO = dyn_cast<ConditionalOperator>(s))
                condexpr = CO->getCond();

            if (condexpr != nullptr) {
                loop_info.has_site_dependent_cond_or_index =
                    is_site_dependent(condexpr, &loop_info.conditional_vars);
                if (loop_info.has_site_dependent_cond_or_index)
                    loop_info.condExpr = condexpr;

                // Flag general cond expression
                loop_info.has_conditional = true;
            }
        }

        parsing_state.ast_depth = 0;

        parsing_state.scope_level--;
        remove_vars_out_of_scope(parsing_state.scope_level);
        parsing_state.skip_children = 1;
        return true;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////
///  List Field<> specializations
///  This does not currently do anything necessary, specializations are
///  now handled as Field<> are used in loops
///
////////////////////////////////////////////////////////////////////////////

int TopLevelVisitor::handle_field_specializations(ClassTemplateDecl *D) {
    // save global, perhaps needed (perhaps not)
    field_decl = D;

    if (cmdline::verbosity >= 2)
        llvm::errs() << "Field<type> specializations in this compilation unit:\n";

    int count = 0;
    for (auto spec = D->spec_begin(); spec != D->spec_end(); spec++) {
        count++;
        auto &args = spec->getTemplateArgs();

        if (args.size() != 1) {
            llvm::errs() << " *** Fatal: More than one type arg for Field<>\n";
            exit(1);
        }
        if (TemplateArgument::ArgKind::Type != args.get(0).getKind()) {
            reportDiag(DiagnosticsEngine::Level::Error, D->getSourceRange().getBegin(),
                       "Expecting type argument in \'Field\' template");
            return (0);
        }

        // Get typename without class, struct... qualifiers
        std::string typestr = args.get(0).getAsType().getAsString(PP);

        if (cmdline::verbosity >= 2) {
            llvm::errs() << "  Field < " << typestr << " >";
            if (spec->isExplicitSpecialization())
                llvm::errs() << " explicit specialization\n";
            else
                llvm::errs() << '\n';
        }
    }
    return (count);

} // end of "field"

////////////////////////////////////////////////////////////////////////////////////
/// These are the main traverse methods
/// By overriding these methods in TopLevelVisitor we can control which nodes are
/// visited. These are control points for the depth of the traversal;
///  skip_children,  ast_depth
////////////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::TraverseStmt(Stmt *S) {

    // if state::skip_children > 0 we'll skip all until return to level up
    if (parsing_state.skip_children > 0)
        parsing_state.skip_children++;

    // go via the original routine...
    if (!parsing_state.skip_children) {
        parsing_state.ast_depth++;
        RecursiveASTVisitor<TopLevelVisitor>::TraverseStmt(S);
        if (parsing_state.ast_depth > 0)
            parsing_state.ast_depth--;
    }

    if (parsing_state.skip_children > 0)
        parsing_state.skip_children--;

    return true;
}

bool TopLevelVisitor::TraverseDecl(Decl *D) {

    // if state::skip_children > 0 we'll skip all until return to level up
    if (parsing_state.skip_children > 0)
        parsing_state.skip_children++;

    // go via the original routine...
    if (!parsing_state.skip_children) {
        parsing_state.ast_depth++;
        RecursiveASTVisitor<TopLevelVisitor>::TraverseDecl(D);
        if (parsing_state.ast_depth > 0)
            parsing_state.ast_depth--;
    }

    if (parsing_state.skip_children > 0)
        parsing_state.skip_children--;

    return true;
}

//  Obsolete when X is new type
// void TopLevelVisitor::require_parity_X(Expr * pExpr) {
//   // Now parity has to be X (or the same as before?)
//   if (get_parity_val(pExpr) != Parity::x) {
//     reportDiag(DiagnosticsEngine::Level::Error,
//                pExpr->getSourceRange().getBegin(),
//                "Use wildcard parity \"X\" or \"Parity::x\"" );
//   }
// }

//////////////////////////////////////////////////////////////////////////////
/// Process the Field<> -references appearing in this loop, and
/// construct the field_info_list
//////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::check_field_ref_list() {

    bool no_errors = true;

    global.assert_loop_parity = false;

    field_info_list.clear();

    for (field_ref &p : field_ref_list) {

        std::string name = get_stmt_str(p.nameExpr);

        field_info *fip = nullptr;

        // search for duplicates: if found, lfip is non-null

        for (field_info &li : field_info_list) {
            if (is_duplicate_expr(li.nameExpr, p.nameExpr)) {
                //  if (name.compare(li.old_name) == 0) {
                fip = &li;
                break;
            }
        }

        if (fip == nullptr) {
            field_info lfv;
            lfv.old_name = name;
            lfv.type_template = get_expr_type(p.nameExpr);
            if (lfv.type_template.find("Field", 0) != 0) {
                reportDiag(DiagnosticsEngine::Level::Error,
                           p.nameExpr->getSourceRange().getBegin(),
                           "Confused: type of Field expression?");
                no_errors = false;
            }
            lfv.type_template.erase(0, 5); // Remove "Field"  from Field<T>

            // get also the fully canonical Field<T>  type.
            lfv.element_type = p.nameExpr->getType()
                                   .getUnqualifiedType()
                                   .getCanonicalType()
                                   .getAsString(PP);
            int a = lfv.element_type.find('<') + 1;
            int b = lfv.element_type.rfind('>') - a;
            lfv.element_type = lfv.element_type.substr(a, b);

            lfv.nameExpr = p.nameExpr; // store the first nameExpr to this field

            field_info_list.push_back(lfv);
            fip = &field_info_list.back();
        }
        // now fip points to the right info element
        // copy that to lf reference
        p.info = fip;

        if (p.is_written && !fip->is_written) {
            // first write to this field var
            fip->first_assign_seq = p.sequence;
            fip->is_written = true;
        }

        // note special treatment of file[X] -read: it is real read only if it
        // comes before or at the same time than assign
        if (p.is_read) {
            if (p.is_direction) {
                fip->is_read_nb = true;
            } else if (!fip->is_written || fip->first_assign_seq >= p.sequence) {
                fip->is_read_atX = true;
            }
        }

        if (p.is_offset)
            fip->is_read_offset = true;

        // save expr record
        fip->ref_list.push_back(&p);

        if (p.is_direction) {

            if (p.is_loop_local_dir)
                fip->is_loop_local_dir = true;

            // does this dir with this field name exist before?
            // Use is_duplicate_expr() to resolve the ptr, it has (some) intelligence to
            // find equivalent constant expressions
            // TODO: use better method?
            bool found = false;
            for (dir_ptr &dp : fip->dir_list) {
                if (p.is_constant_direction) {
                    found = (dp.is_constant_direction &&
                             dp.constant_value == p.constant_value);
                } else {
                    found = is_duplicate_expr(dp.parityExpr, p.parityExpr);
                }

                if (found) {
                    dp.count += (p.is_offset == false); // only nn in count
                    dp.ref_list.push_back(&p);
                    break;
                }
            }

            if (!found) {
                dir_ptr dp;
                dp.parityExpr = p.parityExpr;
                dp.count = (p.is_offset == false);
                dp.is_offset = p.is_offset;
                dp.is_constant_direction = p.is_constant_direction;
                dp.constant_value = p.constant_value;
                dp.is_loop_local_dir = p.is_loop_local_dir;
                dp.direxpr_s = p.direxpr_s; // copy the string expr of Direction

                dp.ref_list.push_back(&p);

                fip->dir_list.push_back(dp);
            }
        } // Direction
    }     // p-loop

    for (field_info &l : field_info_list) {

        // Check if the field can be vectorized

        l.vecinfo = inspect_field_type(l.nameExpr);

        // check for f[ALL] = f[X+dir] -type use, which is undefined

        if (l.is_written && l.dir_list.size() > 0) {

            // There may be error, find culprits
            bool found_error = false;
            for (field_ref *p : l.ref_list) {
                if (p->is_direction && !p->is_written && !p->is_offset) {
                    if (loop_info.parity_value == Parity::all) {

                        reportDiag(DiagnosticsEngine::Level::Error,
                                   p->parityExpr->getSourceRange().getBegin(),
                                   "Simultaneous access '%0' and assignment '%1' not "
                                   "allowed with parity ALL",
                                   get_stmt_str(p->fullExpr).c_str(),
                                   l.old_name.c_str());
                        no_errors = false;
                        found_error = true;

                    } else if (loop_info.parity_value == Parity::none) {
                        reportDiag(DiagnosticsEngine::Level::Remark,
                                   p->parityExpr->getSourceRange().getBegin(),
                                   "Simultaneous access '%0' and assignment to '%1' is "
                                   "allowed "
                                   "only when parity %2 is EVEN or ODD.  Inserting "
                                   "assertion to ensure that",
                                   get_stmt_str(p->fullExpr).c_str(),
                                   l.old_name.c_str(), loop_info.parity_text.c_str());
                        found_error = true;
                    }
                }
            }

            if (found_error) {
                for (field_ref *p : l.ref_list) {
                    if (p->is_written) {
                        reportDiag(DiagnosticsEngine::Level::Note,
                                   p->fullExpr->getSourceRange().getBegin(),
                                   "Location of assignment");
                    }
                }
            }
        }
    }
    return no_errors;
}

/////////////////////////////////////////////////////////////////////////////
/// Check now that the references to variables are as required
/////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::check_var_info_list() {

    for (var_info &vi : var_info_list) {
        if (!vi.is_loop_local && !vi.is_raw) {
            if (vi.reduction_type != reduction::NONE) {
                if (vi.refs.size() > 1) {
                    // reduction only once
                    int i = 0;
                    for (auto &vr : vi.refs) {
                        if (vr.assignop == "+=" || vr.assignop == "*=") {
                            reportDiag(DiagnosticsEngine::Level::Error,
                                       vr.ref->getSourceRange().getBegin(),
                                       "Reduction variable \'%0\' used more than once "
                                       "within one site loop",
                                       vi.name.c_str());
                            break;
                        }
                        i++;
                    }
                    int j = 0;
                    for (auto &vr : vi.refs) {
                        if (j != i)
                            reportDiag(DiagnosticsEngine::Level::Remark,
                                       vr.ref->getSourceRange().getBegin(),
                                       "Other reference to \'%0\'", vi.name.c_str());
                        j++;
                    }
                }

            } else if (vi.is_special_reduction_type) {
                // Use Reduction<> -type vars only as reductions
                for (auto &vr : vi.refs) {
                    reportDiag(DiagnosticsEngine::Level::Error,
                               vr.ref->getSourceRange().getBegin(),
                               "Reduction variables are restricted only for reductions "
                               "(on the lhs of \'+=\' or \'*=\')");
                }

            } else if (vi.is_assigned) {
                // now not reduction
                for (auto &vr : vi.refs) {
                    if (vr.is_assigned)
                        reportDiag(
                            DiagnosticsEngine::Level::Error,
                            vr.ref->getSourceRange().getBegin(),
                            "Cannot assign to variable defined outside site loop "
                            "(unless reduction \'+=\' or \'*=\')");
                }
            }

            // Check if product reduction is done for legal variables
            if (vi.reduction_type == reduction::PRODUCT) {
                legal_types default_legal_types;
                std::string var_type =
                    vi.decl->getType().getCanonicalType().getAsString();
                const bool allowed_reduction_type =
                    default_legal_types.check_if_legal(var_type);
                if (!allowed_reduction_type) {
                    for (auto &vr : vi.refs) {
                        reportDiag(DiagnosticsEngine::Level::Error,
                                   vr.ref->getSourceRange().getBegin(),
                                   "\nProduct reduction variable of type \'%0\' not "
                                   "allowed. \nMust be of type: \'%1\'",
                                   var_type.c_str(),
                                   default_legal_types.as_string().c_str());
                    }
                }
            }
        }
    }

    // iterate through var_info_list until no more is_site_dependent -relations found
    // this should not leave any corner cases behind

    int found;
    do {
        found = 0;
        for (var_info &vi : var_info_list) {
            if (vi.is_site_dependent == false) {
                for (var_info *d : vi.dependent_vars)
                    if (d->is_site_dependent) {
                        vi.is_site_dependent = true;
                        found++;
                        break; // go to next var
                    }
            }
        }
    } while (found > 0);

    // and also get the vectorized type for them, to be prepared...

    if (target.vectorize) {
        for (var_info &vi : var_info_list)
            if (!vi.is_raw) {
                vi.vecinfo.is_vectorizable = is_vectorizable_type(vi.type, vi.vecinfo);
            }
    }
}

////////////////////////////////////////////////////////////////////////////////////////
/// flag_error = true by default in toplevelvisitor.h

SourceRange TopLevelVisitor::getRangeWithSemicolon(Stmt *S, bool flag_error) {
    return getRangeWithSemicolon(S->getSourceRange(), flag_error);
}

SourceRange TopLevelVisitor::getRangeWithSemicolon(SourceRange SR, bool flag_error) {
    SourceRange range(SR.getBegin(),
                      Lexer::findLocationAfterToken(SR.getEnd(), tok::semi,
                                                    TheRewriter.getSourceMgr(),
                                                    Context->getLangOpts(), false));
    if (!range.isValid()) {
        if (flag_error) {
            reportDiag(DiagnosticsEngine::Level::Fatal, SR.getEnd(),
                       "Expecting ';' after expression");
        }
        // put a valid value in any case
        range = SR;
    }

    // llvm::errs() << "Range w semi: " << TheRewriter.getRewrittenText(range) << '\n';
    return range;
}

bool TopLevelVisitor::hasSemicolonAfter(SourceRange sr) {

    SourceLocation s = getSourceLocationAtEndOfRange(sr);
    do {
        s = s.getLocWithOffset(1);
    } while (std::isspace(getChar(s)));
    return getChar(s) == ';';
}


/////////////////////////////////////////////////////////////////////////////
/// Variable decls inside site loops, but also some outside
/////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitVarDecl(VarDecl *var) {

    if (var->getName().str() == "X") {
        static bool second_def = false;
        if (second_def) {
            reportDiag(DiagnosticsEngine::Level::Warning,
                       var->getSourceRange().getBegin(),
                       "Declaring variable 'X' may shadow the site index X");
        }
        second_def = true;
    }


    // and then loop body statements
    if (parsing_state.in_loop_body) {
        // for now care only loop body variable declarations

        if (!var->hasLocalStorage()) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       var->getSourceRange().getBegin(),
                       "Static or external variable declarations not allowed within "
                       "site loops");
            return true;
        }

        if (var->isStaticLocal()) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       var->getSourceRange().getBegin(),
                       "Cannot declare static variables inside site loops");
            return true;
        }

        if (is_field_decl(var)) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       var->getSourceRange().getBegin(),
                       "Cannot declare Field<> variables within site loops");
            parsing_state.skip_children = 1;
            return true;
        }

        // Now it should be automatic local variable decl

        add_var_to_decl_list(var, parsing_state.scope_level);
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::ast_dump_header(const char *s, const SourceRange sr_in,
                                      bool is_function) {
    SourceRange sr = sr_in;
    unsigned linenumber = srcMgr.getSpellingLineNumber(sr.getBegin());
    std::string name = srcMgr.getFilename(sr.getBegin()).str();

    // check if it is macro
    if (sr.getBegin().isMacroID()) {
        CharSourceRange CSR =
            TheRewriter.getSourceMgr().getImmediateExpansionRange(sr.getBegin());
        sr = CSR.getAsRange();
    }

    if (!is_function) {
        std::string source = TheRewriter.getRewrittenText(sr);
        auto n = source.find('\n');

        if (n == std::string::npos) {
            llvm::errs() << "**** AST dump of " << s << " \'" << source << "\' on line "
                         << linenumber << " in file " << name << '\n';
        } else {
            llvm::errs() << "**** AST dump of " << s << " starting with \'"
                         << source.substr(0, n) << "\' on line " << linenumber
                         << " in file " << name << '\n';
        }
    } else {
        llvm::errs() << "**** AST dump of declaration of function '" << s
                     << "' on line " << linenumber << " in file " << name << '\n';
    }
}

void TopLevelVisitor::ast_dump(const Stmt *S) {
    ast_dump_header("statement", S->getSourceRange(), false);
    S->dumpColor();
    llvm::errs() << "*****************************\n";
}

void TopLevelVisitor::ast_dump(const Decl *D) {
    ast_dump_header("declaration", D->getSourceRange(), false);
    D->dumpColor();
    llvm::errs() << "*****************************\n";
}

void TopLevelVisitor::ast_dump(const FunctionDecl *D) {
    ast_dump_header(D->getQualifiedNameAsString().c_str(), D->getSourceRange(), true);
    D->dumpColor();
    llvm::errs() << "*****************************\n";
}

///////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::remove_vars_out_of_scope(unsigned level) {
    while (var_decl_list.size() > 0 && var_decl_list.back().scope > level)
        var_decl_list.pop_back();
}

///////////////////////////////////////////////////////////////////////////////
/// VisitStmt is called for each statement in AST.  Thus, when traversing the
/// AST or part of it we start here, and branch off depending on the statements
/// and parsing_state.lags
///////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitStmt(Stmt *s) {

    if (parsing_state.ast_depth <= 1 && has_pragma(s, pragma_hila::AST_DUMP)) {
        ast_dump(s);
    }

    // Entry point when inside Field[par] = .... body
    if (parsing_state.in_loop_body) {
        return handle_loop_body_stmt(s);
    }

    // loop of type "onsites(p)"
    // Defined as a macro, needs special macro handling
    if (isa<ForStmt>(s)) {

        ForStmt *f = cast<ForStmt>(s);
        SourceLocation startloc = f->getSourceRange().getBegin();

        if (startloc.isMacroID()) {
            Preprocessor &pp = myCompilerInstance->getPreprocessor();
            static std::string loop_call("onsites");
            if (pp.getImmediateMacroName(startloc) == loop_call) {
                // Now we know it is onsites-macro

                CharSourceRange CSR =
                    TheRewriter.getSourceMgr().getImmediateExpansionRange(startloc);
                std::string macro = TheRewriter.getRewrittenText(CSR.getAsRange());
                bool internal_error = true;

                // llvm::errs() << "MACRO STRING " << macro << '\n';

                loop_info.has_pragma_novector = has_pragma(s, pragma_hila::NOVECTOR);
                loop_info.has_pragma_access =
                    has_pragma(s, pragma_hila::ACCESS, &loop_info.pragma_access_args);
                loop_info.has_pragma_omp_parallel_region =
                    has_pragma(s, pragma_hila::IN_OMP_PARALLEL_REGION);

                DeclStmt *init = dyn_cast<DeclStmt>(f->getInit());
                if (init && init->isSingleDecl()) {
                    VarDecl *vd = dyn_cast<VarDecl>(init->getSingleDecl());
                    if (vd) {
                        const Expr *ie = vd->getInit();
                        if (ie) {
                            loop_info.parity_expr = ie;
                            loop_info.parity_value =
                                get_parity_val(loop_info.parity_expr);
                            loop_info.parity_text = remove_initial_whitespace(
                                macro.substr(loop_call.length(), std::string::npos));

                            global.full_loop_text =
                                macro + " " + get_stmt_str(f->getBody());

                            // Delete "onsites()" -text

                            // TheRewriter.RemoveText(CSR);
                            writeBuf->remove(CSR);

                            handle_full_loop_stmt(f->getBody(), false);
                            internal_error = false;
                        }
                    }
                }
                if (internal_error) {
                    reportDiag(DiagnosticsEngine::Level::Error,
                               f->getSourceRange().getBegin(),
                               "\'onsites\'-macro: not a Parity type argument");
                    return true;
                }
            }
        }
        return true;
    }

    //  Starting point for fundamental operation
    //  Field[par] = ....  version with Field<class>
    //  Arg(0)  is the LHS of assignment

    CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s);
    bool found = false;
    if (OP && OP->isAssignmentOp() &&
        is_field_parity_expr(OP->getArg(0)->IgnoreImplicit())) {
        found = true;
    } else {
        // check also Field<double> or some other non-class var
        BinaryOperator *BO = dyn_cast<BinaryOperator>(s);
        if (BO && BO->isAssignmentOp() &&
            is_field_parity_expr(BO->getLHS()->IgnoreImplicit()))
            found = true;
    }

    if (found) {

        loop_info.has_pragma_novector = has_pragma(s, pragma_hila::NOVECTOR);
        loop_info.has_pragma_access =
            has_pragma(s, pragma_hila::ACCESS, &loop_info.pragma_access_args);

        SourceRange full_range = getRangeWithSemicolon(s, false);
        global.full_loop_text = TheRewriter.getRewrittenText(full_range);

        handle_full_loop_stmt(s, true);
        return true;
    }

    // And, for correct level for pragma handling - turns to 0 for stmts inside
    if (isa<CompoundStmt>(s))
        parsing_state.ast_depth = -1;

    Expr *E = dyn_cast<Expr>(s);

    // new stuff: if there is field[coordinate], modify these to appropriate
    // functions

    if (is_field_with_coordinate_stmt(s))
        return true;

    //  Finally, if we get to a Field[Parity] -expression without a loop or assignment
    //  flag error

    if (E && is_field_parity_expr(E)) {
        reportDiag(
            DiagnosticsEngine::Level::Error, E->getSourceRange().getBegin(),
            "Field[Parity] -expression is allowed only in LHS of Field assignment "
            "statements (Field[par] = ...)");
        parsing_state.skip_children = 1;
        return true;

    } else if (E && is_field_with_X_expr(E)) {
        reportDiag(DiagnosticsEngine::Level::Error, E->getSourceRange().getBegin(),
                   "Field[X] -expressions allowed only in site loops");
        parsing_state.skip_children = 1;
        return true;
    }

    //   else if (E && is_X_index_type(E)) {
    //   reportDiag(DiagnosticsEngine::Level::Error,
    //               E->getSourceRange().getBegin(),
    //               "Use of \"X\" is allowed only in site loops");
    //   parsing_state.skip_children = 1;
    // }

    // and add special handling for special function calls here

    if (CallExpr *CE = dyn_cast<CallExpr>(s)) {
        if (FunctionDecl *FD = CE->getDirectCallee()) {

            // THIS DOES NOT WORK BECAUSE REWRITE-BUFFERS ARE NOT CHANGED ANY MORE
            // is it memalloc(size) -call -> substitute with
            // memalloc( size, filename, linenumber)
            // don't know if this is really useful

            // if (FD->getNameAsString() == "memalloc" && CE->getNumArgs() == 1) {
            //   SourceLocation sl = CE->getRParenLoc();

            //   // generate new args
            //   std::string name(srcMgr.getFilename(sl));
            //   std::size_t i = name.rfind('/');
            //   if (i != std::string::npos) name = name.substr(i);

            //   std::string args(", \"");
            //   args.append(name).append( "\", " ).append(
            //       std::to_string( srcMgr.getSpellingLineNumber(sl)));

            // if (!writeBuf->is_edited(sl)) writeBuf->insert(sl,args);
            // }
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////
/// Handle Field[Coordinate] -expressions
////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::is_field_with_coordinate_stmt(Stmt *s) {

    static bool was_previously_assigned = false;

    // Check first if this is field[c] = ... -stmt

    // assigment through operator=()

    CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s);
    if (OP && OP->isAssignmentOp() &&
        is_field_with_coordinate(OP->getArg(0)->IgnoreImplicit())) {
        const char *sp = getOperatorSpelling(OP->getOperator());
        if (sp[0] != '=') {
            // it's a compound assignment, not allowed
            reportDiag(
                DiagnosticsEngine::Level::Error, OP->getOperatorLoc(),
                "Only direct assignment '=' allowed for Field[CoordinateVector]");
            return false;
        }

        field_with_coordinate_assign(OP->getArg(0)->IgnoreImplicit(),
                                     OP->getArg(1)->IgnoreImplicit(),
                                     OP->getOperatorLoc());
        was_previously_assigned = true;

        return true;
    }

    // check also Field<double> or some other non-class assign
    BinaryOperator *BO = dyn_cast<BinaryOperator>(s);
    if (BO && BO->isAssignmentOp() &&
        is_field_with_coordinate(BO->getLHS()->IgnoreImplicit())) {
        if (BO->isCompoundAssignmentOp()) {
            reportDiag(
                DiagnosticsEngine::Level::Error, BO->getOperatorLoc(),
                "Only direct assignment '=' allowed for Field[CoordinateVector]");
            return false;
        }
        field_with_coordinate_assign(BO->getLHS()->IgnoreImplicit(),
                                     BO->getRHS()->IgnoreImplicit(),
                                     BO->getOperatorLoc());
        was_previously_assigned = true;

        return true;
    }

    // Handle now just field[coord] -expressions, these are not assignments.
    // If assginment was previously seen, the expr has been handled already

    Expr *E = dyn_cast<Expr>(s);
    if (E && is_field_with_coordinate(E)) {
        if (!was_previously_assigned) {
            field_with_coordinate_read(E);
        }
        was_previously_assigned = false;

        return true;
    }

    return false;
}

void TopLevelVisitor::field_with_coordinate_assign(Expr *lhs, Expr *rhs,
                                                   SourceLocation oploc) {

    // lhs is field[par] -expr - we know here that arg is of the
    // right type, so this succeeds
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(lhs);

    if (writeBuf->get(oploc, 1) == ",")
        return; // this has already been changed

    // now change = -> ,
    writeBuf->replace(SourceRange(oploc, oploc), ",");
    // Remove ]
    writeBuf->remove(SourceRange(OC->getRParenLoc(), OC->getRParenLoc()));
    // find [ and insert method call
    SourceLocation sl = OC->getArg(1)->getBeginLoc();
    while (sl.isValid() && getChar(sl) != '[')
        sl = sl.getLocWithOffset(-1);
    writeBuf->replace(SourceRange(sl, sl), ".set_element_at(");
    // and insert closing )
    sl = getSourceLocationAtEndOfRange(rhs->getSourceRange());
    sl = sl.getLocWithOffset(1);
    writeBuf->insert(sl, ")", true);
}

void TopLevelVisitor::field_with_coordinate_read(Expr *E) {

    // expr is field[par] -expr - we know here that arg is of the
    // right type, so this succeeds
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    // Replace ] with )
    writeBuf->replace(SourceRange(OC->getRParenLoc(), OC->getRParenLoc()), ")");
    // find [ and insert method call
    SourceLocation sl = OC->getArg(1)->getBeginLoc();
    while (sl.isValid() && getChar(sl) != '[')
        sl = sl.getLocWithOffset(-1);
    writeBuf->replace(SourceRange(sl, sl), ".get_element(");
}

////////////////////////////////////////////////////////////////////////
/// This is visited for every function declaration and specialization
/// (either function specialization or template class method specialization)
/// If needed, specializations are "rewritten" open
////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitFunctionDecl(FunctionDecl *f) {

    // operate (usually) only non-templated functions and template specializations
    if (has_pragma(f, pragma_hila::AST_DUMP))
        ast_dump(f);

    // For the pragmas, we need to check the "previousdecl" -- prototype, if there is a
    // pragma there Automatic now in has_pragma

    if (has_pragma(f, pragma_hila::LOOP_FUNCTION)) {
        // This function can be called from a loop,
        // mark as noticed -- note that we do not have call arguments
        loop_function_check(f);
    }

    // llvm::errs() << "Function " << f->getNameInfo().getName() << "\n";

    if (f->isThisDeclarationADefinition() && f->hasBody()) {

        // llvm::errs() << "LOOKING AT FUNC " << f->getQualifiedNameAsString() << " on
        // line "
        //              << srcMgr.getSpellingLineNumber(f->getBeginLoc()) << " file "
        //              << srcMgr.getFilename( f->getBeginLoc() ) << '\n';

        global.currentFunctionDecl = f;

        Stmt *FuncBody = f->getBody();

        // Type name as string
        std::string TypeStr = f->getReturnType().getAsString();

        // Function name
        std::string FuncName = f->getNameInfo().getName().getAsString();

        // llvm::errs() << " - Function "<< FuncName << "\n";

        // if (does_function_contain_field_access(f)) {
        //   loop_callable = false;
        // }

        switch (f->getTemplatedKind()) {
        case FunctionDecl::TemplatedKind::TK_NonTemplate:
            // Normal, non-templated class method -- nothing here

            if (f->isCXXClassMember()) {
                CXXMethodDecl *method = dyn_cast<CXXMethodDecl>(f);
                CXXRecordDecl *parent = method->getParent();
                if (parent->isTemplated()) {
                    if (does_function_contain_field_access(f)) {
                        // Skip children here. Loops in the template may be
                        // incorrect before they are specialized
                        parsing_state.skip_children = 1;
                    }
                }
            }
            break;

        case FunctionDecl::TemplatedKind::TK_FunctionTemplate:
            // not descent inside templates
            parsing_state.skip_children = 1;
            break;

            // do all specializations here: either direct function template,
            // specialized through class template, or specialized
            // because it's a base class function specialized by derived class

        case FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization:
        case FunctionDecl::TemplatedKind::TK_MemberSpecialization:
        case FunctionDecl::TemplatedKind::TK_DependentFunctionTemplateSpecialization:

            if (does_function_contain_field_access(f)) {
                specialize_function_or_method(f);
            } else {
                parsing_state.skip_children = 1; // no reason to look at it further
            }
            break;

        default:
            // do nothing
            break;
        }

        SourceLocation ST = f->getSourceRange().getBegin();
        global.location.function = ST;

        if (cmdline::funcinfo) {
            // Add comment before
            std::stringstream SSBefore;
            SSBefore << "// hilapp info:\n"
                     << "//   begin function " << FuncName << " returning " << TypeStr
                     << '\n'
                     << "//   of template type "
                     << print_TemplatedKind(f->getTemplatedKind()) << '\n';
            writeBuf->insert(ST, SSBefore.str(), true, true);
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////
/// And do the visit for constructors too - used here just for flag
////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitCXXConstructorDecl(CXXConstructorDecl *c) {

    if (has_pragma(c, pragma_hila::LOOP_FUNCTION)) {
        // This function can be called from a loop,
        // mark as noticed -- note that we do not have call arguments
        loop_function_check(c);
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////
/// This does the heavy lifting of specializing function templates and
/// methods defined within template classes.  This is needed if there are
/// site loops within the functions
////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::specialize_function_or_method(FunctionDecl *f) {
    // This handles all functions and methods. Parent is non-null for methods,
    // and then is_static gives the static flag

    bool no_inline;
    bool is_static = false;
    CXXRecordDecl *parent = nullptr;

    /* Check if the function is a class method */
    if (f->isCXXClassMember()) {
        // method is defined inside template class.  Could be a chain of classes!
        CXXMethodDecl *method = dyn_cast<CXXMethodDecl>(f);
        parent = method->getParent();
        is_static = method->isStatic();
        no_inline = cmdline::method_spec_no_inline;
    } else {
        no_inline = cmdline::function_spec_no_inline;
    }

    srcBuf *writeBuf_saved = writeBuf;
    srcBuf funcBuf(&TheRewriter, f);
    writeBuf = &funcBuf;

    std::vector<std::string> par, arg;

    // llvm::errs() << "funcBuffer:\n" << funcBuf.dump() << '\n';

    // cannot rely on getReturnTypeSourceRange() for methods.  Let us not even try,
    // change the whole method here

    bool is_templated_func =
        (f->getTemplatedKind() ==
         FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization);

    int ntemplates = 0;
    std::string template_args = "";
    std::vector<const TemplateArgument *> typeargs = {};

    if (is_templated_func) {
        // Get here the template param->arg mapping for func template
        auto tal = f->getTemplateSpecializationArgs();
        auto tpl = f->getPrimaryTemplate()->getTemplateParameters();
        assert(tal && tpl && tal->size() == tpl->size() &&
               "Method template par/arg error");

        make_mapping_lists(tpl, *tal, par, arg, typeargs, &template_args);
        ntemplates = 1;

        // SourceLocation sl = f->getPointOfInstantiation();
        // llvm::errs() << "Function " << f->getNameAsString() << " instantiated line "
        //              << srcMgr.getSpellingLineNumber(sl) << " file "
        //              << srcMgr.getFilename(f->getBeginLoc()) << '\n';
    }

    // Get template mapping for classes
    // parent is from: CXXRecordDecl * parent = method->getParent();
    if (parent)
        ntemplates += get_param_substitution_list(parent, par, arg, typeargs);

    // llvm::errs() << "Num nesting templates " << ntemplates << '\n';
    // llvm::errs() << "Specializing function " << f->getQualifiedNameAsString()
    //              << " template args: " << template_args << '\n';

    funcBuf.replace_tokens(f->getSourceRange(), par, arg);

    // Check real parameters: default values must be removed, i.e. remove "= value"
    for (unsigned i = 0; i < f->getNumParams(); i++) {
        ParmVarDecl *pvd = f->getParamDecl(i);
        if (pvd->hasDefaultArg() && !pvd->hasInheritedDefaultArg()) {
            // llvm::errs() << "Default arg! " << get_stmt_str(pvd->getDefaultArg())
            //              << " func " << f->getNameAsString() << '\n';

            SourceRange sr = pvd->getDefaultArgRange();

            // If there is a prototype, and the def. parameter is there, it is before in
            // translationunit and nothing needs to be done
            if (srcMgr.isBeforeInTranslationUnit(f->getSourceRange().getBegin(),
                                                 sr.getBegin())) {

                // funcBuf.dump();

                // if default arg is macro, need to read the immediate range
                if (sr.getBegin().isMacroID()) {
                    CharSourceRange CSR =
                        TheRewriter.getSourceMgr().getImmediateExpansionRange(
                            sr.getBegin());
                    sr = CSR.getAsRange();
                }

                SourceLocation b = sr.getBegin();
                SourceLocation m = pvd->getSourceRange().getBegin();

                while (funcBuf.get(b, 1) != "=" && b > m) {
                    b = b.getLocWithOffset(-1);
                }

                sr.setBegin(b);
                funcBuf.remove(sr);
            }
        }
    }


    // if we have special function do not try to explicitly specialize the name
    bool is_special = false;
    if (isa<CXXConstructorDecl>(f) || isa<CXXConversionDecl>(f) ||
        isa<CXXDestructorDecl>(f)) {
        template_args.clear();
        is_special = true;
    }

    // Careful here: for functions which have been declared first, defined later,
    // getNameInfo().getsourceRange() points to the declaration. Thus, it is not in the
    // range of this function definition.
    // TODO: Use this fact to generate specialization declarations?

    SourceRange sr = f->getNameInfo().getSourceRange();

    if (funcBuf.is_in_range(sr)) {
        // remove all to the end of the name
        funcBuf.remove(0, funcBuf.get_index(sr.getBegin()));
        funcBuf.remove(sr);

    } else {

        // now we have to hunt for the function name
        int l = funcBuf.find_original(0, '(');

        // location of first paren - function name should be just before this, right?
        // TODO: what happens with possible keywords with parens, or macro definitions?
        // There could be template arguments after this, but the parameters
        // (template_args) above should take care of this, so kill all
        if (l > 0) {
            // llvm::errs() << "Searching name " << f->getNameAsString() << '\n';
            int j = funcBuf.find_original_word(0, f->getNameAsString());
            if (j < 0 || j > l)
                l = -1; // name not found
        }
        if (l < 0) {
            reportDiag(DiagnosticsEngine::Level::Fatal, f->getSourceRange().getBegin(),
                       "Internal error: Could not locate function name");
            exit(1);
        }
        funcBuf.remove(0, l - 1);
    }

    // FInally produce the function return type and full name + possible templ. args.

    // put right return type and function name
    funcBuf.insert(0, f->getQualifiedNameAsString() + template_args, true, true);
    if (!is_special) {
        // Declarations with a trailing return type behave weirdly, they have empty
        // ReturnTypeSourceRange, but the getDeclaredReturnType is the explicit return
        // type.
        if (TheRewriter.getRewrittenText(f->getReturnTypeSourceRange()) == "") {
            // So this one has a trailing return type. Just add auto.
            funcBuf.insert(0, " auto ", true, true);
        } else {
            // Normal case, just add the declared return type.
            funcBuf.insert(0, f->getDeclaredReturnType().getAsString(PP) + " ", true,
                           true);
        }
    }

    // remove "static" if it is so specified in methods - not needed now
    // if (is_static) {
    //   funcBuf.replace_token(0,
    //                         funcBuf.get_index(f->getNameInfo().getSourceRange().getBegin()),
    //                         "static","");
    // }

    if (!f->isInlineSpecified() && !no_inline)
        funcBuf.insert(0, "inline ", true, true);

    for (int i = 0; i < ntemplates; i++) {
        funcBuf.insert(0, "template <>\n", true, true);
    }

    // Insertion point for the specialization:  it seems that the template
    // specialization "source" is the 1st declaration.  We want to go after

    SourceLocation insertion_point =
        spec_insertion_point(typeargs, global.location.bot, f);

    SourceRange decl_sr = get_func_decl_range(f);
    std::string wheredefined = "";
    if (f->isInlineSpecified() || !no_inline ||
        !in_specialization_db(funcBuf.get(decl_sr), wheredefined)) {
        // Now we should write the spec here

        // llvm::errs() << "new func:\n" << funcBuf.dump() <<'\n';
        // visit the body
        TraverseStmt(f->getBody());

        // llvm::errs() << "new func again:\n" << funcBuf.dump() <<'\n';

        // insert after the current toplevedecl
        std::stringstream sb;
        sb << "\n\n// ++++++++ hilapp generated function/method specialization\n"
           << funcBuf.dump() << "\n// ++++++++\n\n";

        // buffer is not necessarily in toplevelBuf, so:

        srcBuf *filebuf = get_file_srcBuf(insertion_point);
        filebuf->insert(findChar(insertion_point, '\n'), sb.str(), false, true);
    } else {
        // Now the function has been written before (and not inline)
        // just insert declaration, defined on another compilation unit
        toplevelBuf->insert(
            findChar(global.location.bot, '\n'),
            "\n// ++++++++ Generated specialization declaration, defined "
            "in compilation unit " +
                wheredefined + "\n" + funcBuf.get(decl_sr) + ";\n// ++++++++\n",
            false, false);
    }

    writeBuf = writeBuf_saved;
    funcBuf.clear();
    // don't descend again
    parsing_state.skip_children = 1;
}

////////////////////////////////////////////////////////////////////////////////////////
// locate range of specialization "template< ..> .. func<...>( ... )"
// tf is ptr to template, and f to instantiated function

SourceRange TopLevelVisitor::get_func_decl_range(FunctionDecl *f) {

    if (f->hasBody()) {
        SourceLocation a = f->getSourceRange().getBegin();
        SourceLocation b = f->getBody()->getSourceRange().getBegin();

        while (srcMgr.getFileOffset(b) >= srcMgr.getFileOffset(a)) {
            b = b.getLocWithOffset(-1);
            const char *p = srcMgr.getCharacterData(b);
            if (!std::isspace(*p))
                break;
        }
        SourceRange r(a, b);
        return r;
    }

    return f->getSourceRange();
}

////////////////////////////////////////////////////////////////////////////
/// Class template visitor: we check this because we track field and field_storage
/// specializations (not really needed?)
////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitClassTemplateDecl(ClassTemplateDecl *D) {

    // go through with real definitions or as a part of chain
    if (D->isThisDeclarationADefinition()) { // } || state::class_level > 0) {

        // insertion pt for specializations
        //     if (state::class_level == 1) {
        //       global.location.spec_insert =
        //       findChar(D->getSourceRange().getEnd(),'\n');
        //     }

        const TemplateParameterList *tplp = D->getTemplateParameters();
        // save template params in a list, for templates within templates .... ugh!
        // global.class_templ_params.push_back( tplp );

        // this block for debugging
        if (cmdline::funcinfo) {
            std::stringstream SSBefore;
            SSBefore << "// hilapp info:\n"
                     << "//   Begin template class " << D->getNameAsString()
                     << " with template params\n//    ";
            for (unsigned i = 0; i < tplp->size(); i++)
                SSBefore << tplp->getParam(i)->getNameAsString() << " ";
            SourceLocation ST = D->getSourceRange().getBegin();
            SSBefore << '\n';

            writeBuf->insert(ST, SSBefore.str(), true, true);
        }
        // end block

        // global.in_class_template = true;
        // Should go through the template in order to find function templates...
        // Comment out now, let roll through "naturally".
        // TraverseDecl(D->getTemplatedDecl());

        if (D->getNameAsString() == "Field") {
            handle_field_specializations(D);
        } else if (D->getNameAsString() == "field_storage") {
            field_storage_decl = D;
        } else {
        }

        // global.in_class_template = false;

        // Now do traverse the template naturally
        // state::skip_children = 1;
    }

    return true;
}

/////////////////////////////////////////////////////////////////////////////////////
/// Find the element typealias here -- could not work
/// directly with VisitTypeAliasTemplateDecl below, a bug??
/// Also find the 'ast dump' pragma
/////////////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitDecl(Decl *D) {

    if (parsing_state.ast_depth == 1 && has_pragma(D, pragma_hila::AST_DUMP)) {
        ast_dump(D);
    }

    // THis does nothing here

    // auto t = dyn_cast<TypeAliasTemplateDecl>(D);
    // if (t && t->getNameAsString() == "element") {
    //   llvm::errs() << "Got field storage\n";
    // }

    return true;
}

/////////////////////////////////////////////////////////////////////////////////////
/// THis is just to enable ast dump
/////////////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::VisitType(Type *T) {

    auto *recdecl = T->getAsCXXRecordDecl();
    if (recdecl != nullptr) {
        if (has_pragma(recdecl->getInnerLocStart(), pragma_hila::AST_DUMP)) {
            ast_dump_header("type", recdecl->getInnerLocStart(), false);
            recdecl->dumpColor();
        }
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////
/// Check that all template specialization type arguments are defined at the point
/// where the specialization is inserted
///
/// Determine that the candidate insertion point is OK, if not, make new
/// Also do the "kernel insertion point" at the same time, from the default one
///
/////////////////////////////////////////////////////////////////////////////////

SourceLocation
TopLevelVisitor::spec_insertion_point(std::vector<const TemplateArgument *> &typeargs,
                                      SourceLocation ip, FunctionDecl *f) {

    if (f->hasBody() &&
        srcMgr.isBeforeInTranslationUnit(ip, f->getBody()->getSourceRange().getEnd())) {
        // Now the body definition comes after candidate - make new ip
        // This situation arises if func is declared before definition

        SourceLocation sl;

        if (f->isCXXClassMember()) {

            CXXMethodDecl *md = dyn_cast<CXXMethodDecl>(f);
            // method, look at the encompassing classes

            // find the outermost class decl
            CXXRecordDecl *rd, *parent = dyn_cast<CXXRecordDecl>(md->getParent());

            while ((rd = dyn_cast<CXXRecordDecl>(parent->getParent()))) {
                parent = rd;
            }

            // now parent is the outermost class or function, insertion after the end of
            // it
            sl = parent->getEndLoc();

            // skip end chars of definition - sl points to }
            bool error = true;

            char c = getChar(sl);
            if (c == '}') {

                // class defn ends at ;  - there may be variables too
                do {
                    sl = getNextLoc(sl);
                    c = getChar(sl);
                } while (sl.isValid() && c != ';');

                if (c == ';') {
                    sl = getNextLoc(sl);
                    error = false;
                }
            }

            if (error) {

                llvm::errs()
                    << "hilapp internal error: confusion in finding end loc of class\n";
                llvm::errs() << " on line " << srcMgr.getSpellingLineNumber(sl)
                             << " file " << srcMgr.getFilename(f->getBeginLoc())
                             << '\n';
                exit(1);
            }

            // set also the kernel insertion point (if needed at all)
            global.location.kernels = getSourceLocationAtStartOfDecl(parent);

            // It is still possible that the function is defined further down.
            // If that is the case, we insert the
            // specializaton after it.  Must be at the global level (not within
            // a class).
            if (srcMgr.isBeforeInTranslationUnit(
                    sl, f->getBody()->getSourceRange().getBegin())) {

                sl = f->getBody()->getSourceRange().getEnd();
                sl = getNextLoc(sl); // skip }
            }

        } else {

            // Now "std" function, get the location
            f = f->getDefinition();
            sl = f->getBody()->getSourceRange().getEnd();
            sl = getNextLoc(sl); // skip }

            // and the kernel loc too
            global.location.kernels = getSourceLocationAtStartOfDecl(f);
        }

        if (sl.isInvalid() || srcMgr.isBeforeInTranslationUnit(
                                  sl, f->getBody()->getSourceRange().getBegin())) {

            reportDiag(DiagnosticsEngine::Level::Warning,
                       f->getSourceRange().getBegin(),
                       "hilapp internal error: could not resolve the specialization"
                       " insertion point for function  \'%0\'",
                       f->getQualifiedNameAsString().c_str());
        }

        ip = sl;
    }

    // Check if the

    for (const TemplateArgument *tap : typeargs) {
        // llvm::errs() << " - Checking tp type " << tap->getAsType().getAsString() <<
        // '\n';
        const Type *tp = tap->getAsType().getTypePtrOrNull();
        // Builtins are fine too
        if (tp && !tp->isBuiltinType()) {
            RecordDecl *rd = tp->getAsRecordDecl();
            if (rd &&
                srcMgr.isBeforeInTranslationUnit(ip, rd->getSourceRange().getBegin())) {
                reportDiag(DiagnosticsEngine::Level::Warning,
                           f->getSourceRange().getBegin(),
                           "hilapp internal error: specialization insertion point for "
                           "function \'%0\'"
                           " appears to be before the declaration of type \'%1\', code "
                           "might not compile",
                           f->getQualifiedNameAsString().c_str(),
                           tap->getAsType().getAsString().c_str());

                // try to move the insertion point - fails, TODO: more carefully!
                // ip =
                // getRangeWithSemicolon(rd->getSourceRange()).getEnd().getLocWithOffset(1);
            }
        }
    }
    global.location.kernels = ip;
    return ip;
}

////////////////////////////////////////////////////////////////////////////////////
/// Returns the mapping params -> args for class templates, inner first.  Return value
/// is the number of template nestings
////////////////////////////////////////////////////////////////////////////////////

int TopLevelVisitor::get_param_substitution_list(
    CXXRecordDecl *r, std::vector<std::string> &par, std::vector<std::string> &arg,
    std::vector<const TemplateArgument *> &typeargs) {

    if (r == nullptr)
        return 0;

    int level = 0;
    if (r->getTemplateSpecializationKind() ==
        TemplateSpecializationKind::TSK_ImplicitInstantiation) {

        ClassTemplateSpecializationDecl *sp =
            dyn_cast<ClassTemplateSpecializationDecl>(r);
        if (sp) {

            const TemplateArgumentList &tal = sp->getTemplateArgs();
            assert(tal.size() > 0);

            ClassTemplateDecl *ctd = sp->getSpecializedTemplate();
            TemplateParameterList *tpl = ctd->getTemplateParameters();
            assert(tpl && tpl->size() > 0);

            assert(tal.size() == tpl->size());

            make_mapping_lists(tpl, tal, par, arg, typeargs, nullptr);

            level = 1;
        }
    } else {
        // llvm::errs() << "No specialization of class " << r->getNameAsString() <<
        // '\n';
    }

    auto *parent = r->getParent();
    if (parent) {
        if (CXXRecordDecl *pr = dyn_cast<CXXRecordDecl>(parent))
            return level + get_param_substitution_list(pr, par, arg, typeargs);
    }
    return level;
}

///////////////////////////////////////////////////////////////////////////////////
/// mapping of template params <-> args
///////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::make_mapping_lists(
    const TemplateParameterList *tpl, const TemplateArgumentList &tal,
    std::vector<std::string> &par, std::vector<std::string> &arg,
    std::vector<const TemplateArgument *> &typeargs, std::string *argset) {

    if (argset)
        *argset = "< ";

    // Get argument strings without class, struct... qualifiers

    for (int i = 0; i < tal.size(); i++) {
        if (argset && i > 0)
            *argset += ", ";
        switch (tal.get(i).getKind()) {
        case TemplateArgument::ArgKind::Type:
            arg.push_back(tal.get(i).getAsType().getAsString(PP));
            par.push_back(tpl->getParam(i)->getNameAsString());
            if (argset)
                *argset += arg.back();       // write just added arg
            typeargs.push_back(&tal.get(i)); // save type-type arguments
            break;

        case TemplateArgument::ArgKind::Integral:

#if LLVM_VERSION_MAJOR < 13
            arg.push_back(tal.get(i).getAsIntegral().toString(10));
#else
            arg.push_back(llvm::toString(tal.get(i).getAsIntegral(), 10));
#endif

            par.push_back(tpl->getParam(i)->getNameAsString());
            if (argset)
                *argset += arg.back();
            break;

        default:
            llvm::errs() << " debug: ignoring template argument of argument kind "
                         << tal.get(i).getKind() << " with parameter "
                         << tpl->getParam(i)->getNameAsString() << '\n';
            exit(1); // Don't know what to do
        }
    }
    if (argset)
        *argset += " >";

    return;
}

/////////////////////////////////////////////////////////////////////////
/// Hook to set the output buffer

void TopLevelVisitor::set_writeBuf(const FileID fid) {
    writeBuf = get_file_buffer(TheRewriter, fid);
    toplevelBuf = writeBuf;
}
