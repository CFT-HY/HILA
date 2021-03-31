#include "stringops.h"
#include "generalvisitor.h"

//////////////////////////////////////////////////////////////////////////////////
///  Implementation of some methods for generalvisitor
///  For full list of generalvisitor commands, see generalvisitor.h
//////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_duplicate_expr(const Expr *a, const Expr *b) {
    // Use the Profile function in clang, which "fingerprints"
    // statements
    llvm::FoldingSetNodeID IDa, IDb;
    a->Profile(IDa, *Context, true);
    b->Profile(IDb, *Context, true);
    return (IDa == IDb);
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_field_parity_expr(Expr *E) {

    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0 &&
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

/////////////////////////////////////////////////////////////////
/// Checks if E is parity of a field (for example f[X]).
/// Catches both parity and X_plus_direction
bool GeneralVisitor::is_field_with_X_expr(Expr *E) {
    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0 &&
        is_field_expr(OC->getArg(0))) {

        return is_X_index_type(OC->getArg(1));
    }
    return false;
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_field_with_X_and_dir(Expr *E) {
    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0 &&
        is_field_expr(OC->getArg(0))) {

        return is_X_and_dir_type(OC->getArg(1));
    }
    return false;
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_field_with_coordinate(Expr *E) {
    E = E->IgnoreParens();
    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

    if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0 &&
        is_field_expr(OC->getArg(0))) {

        return is_coordinate_type(OC->getArg(1));
    }
    return false;
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_assignment_expr(Stmt *s, std::string *opcodestr,
                                        bool &iscompound) {
    if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
        if (OP->isAssignmentOp()) {

            // TODO: there should be some more elegant way to do this
            const char *sp = getOperatorSpelling(OP->getOperator());
            if (sp[0] != '=' && sp[1] == '=')
                iscompound = true;
            else
                iscompound = false;

            if (opcodestr)
                *opcodestr = getOperatorSpelling(OP->getOperator());

            return true;
        }
    }

    // This is for arithmetic type assignments
    if (BinaryOperator *B = dyn_cast<BinaryOperator>(s)) {
        if (B->isAssignmentOp()) {
            iscompound = B->isCompoundAssignmentOp();
            if (opcodestr)
                *opcodestr = B->getOpcodeStr().str();
            return true;
        }
    }

    return false;
}

///////////////////////////////////////////////////////////////////
/// Check if the RHS of assignment is site dependent
///////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_rhs_site_dependent(Stmt *s, std::vector<var_info *> *vi) {

    if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
        if (OP->isAssignmentOp()) {
            return is_site_dependent(OP->getArg(1), vi);
        }
    }

    if (BinaryOperator *B = dyn_cast<BinaryOperator>(s)) {
        if (B->isAssignmentOp()) {
            return is_site_dependent(B->getRHS(), vi);
        }
    }
    // one of these should have triggered!
    assert(0 && "Internal error in RHS analysis");
}

///////////////////////////////////////////////////////////////////
/// Check the site dependence of the args of "access operators",
/// [] or .e().  These prevent vectorization
///////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_site_dependent_access_op(Expr *e) {
    e = e->IgnoreParens();

    if (ArraySubscriptExpr *ASE = dyn_cast<ArraySubscriptExpr>(e)) {

        return is_site_dependent(ASE->getIdx(), &loop_info.conditional_vars);
    }

    CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(e);
    if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0 &&
        !is_field_expr(OC->getArg(0))) {

        return is_site_dependent(OC->getArg(1), &loop_info.conditional_vars);
    }

    if (CallExpr *CE = dyn_cast<CallExpr>(e)) {
        FunctionDecl *FD = CE->getDirectCallee();
        if (CXXMethodDecl *MD = dyn_cast<CXXMethodDecl>(FD)) {

            std::string method = MD->getNameAsString();
            std::string parent = MD->getParent()->getNameAsString();

            if (method == "e" && (parent == "Matrix" || parent == "Array" ||
                                  parent == "CoordinateVector_t")) {
                bool dep = is_site_dependent(CE->getArg(0), &loop_info.conditional_vars);
                // For matrix or array e may have 1 or 2 args
                if (CE->getNumArgs() > 1) 
                    dep = dep || is_site_dependent(CE->getArg(1), &loop_info.conditional_vars);
                return dep;
            }
        }
    }
    return false;
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_function_call_stmt(Stmt *s) {
    if (auto *Call = dyn_cast<CallExpr>(s)) {
        // llvm::errs() << "Function call found: " << get_stmt_str(s) << '\n';
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////
/// is the stmt pointing now to a member call
bool GeneralVisitor::is_member_call_stmt(Stmt *s) {
    if (auto *Call = dyn_cast<CXXMemberCallExpr>(s)) {
        // llvm::errs() << "Member call found: " << get_stmt_str(s) << '\n';
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////
/// is the stmt pointing now to a constructor
bool GeneralVisitor::is_constructor_stmt(Stmt *s) {
    if (auto *Call = dyn_cast<CXXConstructExpr>(s)) {
        // llvm::errs() << "Constructor found: " << get_stmt_str(s) << '\n';
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////
/// is the stmt pointing now to a cast

bool GeneralVisitor::is_user_cast_stmt(Stmt *s) {

    if (CastExpr *ce = dyn_cast<CastExpr>(s)) {
        if (NamedDecl *p = ce->getConversionFunction()) {
            // now it is a user conversion function
            // llvm::errs() << "GOT USER CONVERSION " << p->getNameAsString() << '\n';
            return true;
        }
    }

    return false;
}

/////////////////////////////////////////////////////////////////
/// Does the statement end with a semicolon
bool GeneralVisitor::isStmtWithSemicolon(Stmt *S) {
    SourceLocation l = Lexer::findLocationAfterToken(S->getEndLoc(), tok::semi,
                                                     TheRewriter.getSourceMgr(),
                                                     Context->getLangOpts(), false);
    if (l.isValid()) {
        //    llvm::errs() << "; found " << get_stmt_str(S) << '\n';
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////

Parity GeneralVisitor::get_parity_val(const Expr *pExpr) {
    SourceLocation SL;
    APValue APV;

    if (pExpr->isCXX11ConstantExpr(*Context, &APV, &SL)) {
        // Parity is now constant
        int64_t val = (APV.getInt().getExtValue());
        Parity p;
        if (0 <= val && val <= (int)Parity::all) {
            p = static_cast<Parity>(val);
        } else {
            reportDiag(DiagnosticsEngine::Level::Fatal,
                       pExpr->getSourceRange().getBegin(),
                       "hilapp internal error, unknown parity");
            exit(1);
        }
        if (p == Parity::none) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       pExpr->getSourceRange().getBegin(),
                       "Parity::none is reserved for internal use");
        }

        return p;
    } else {
        return Parity::none;
    }
}

///////////////////////////////////////////////////////////////////////////////////
/// Source Location utilities
///////////////////////////////////////////////////////////////////////////////////

SourceLocation GeneralVisitor::getSourceLocationAtEndOfRange(SourceRange r) {
    int i = TheRewriter.getRangeSize(r);
    return r.getBegin().getLocWithOffset(i - 1);
}

/// Hunt for SourceLocation at the beginning of a decl, including template<>
/// For templated decls one has to locate the original templatedecl.
/// I am not 100% sure these work for partially specialized templates

SourceLocation GeneralVisitor::getSourceLocationAtStartOfDecl(Decl *d) {

    SourceLocation sl;

    if (ClassTemplateDecl *ctd = dyn_cast<ClassTemplateDecl>(d)) {
        // class template decl

        sl = ctd->getTemplateParameters()->getTemplateLoc();

    } else if (FunctionTemplateDecl *ftd = dyn_cast<FunctionTemplateDecl>(d)) {
        // function template

        sl = ftd->getTemplateParameters()->getTemplateLoc();

    } else if (FunctionDecl *f = dyn_cast<FunctionDecl>(d)) {
        // now it is a function, but is it templated?

        if (FunctionTemplateDecl *ftd = f->getDescribedFunctionTemplate()) {
            sl = ftd->getTemplateParameters()->getTemplateLoc();
        } else {
            sl = f->getSourceRange().getBegin();
        }

    } else if (CXXRecordDecl *rd = dyn_cast<CXXRecordDecl>(d)) {
        // it's a class, is it templated or not?

        if (ClassTemplateDecl *ctd = rd->getDescribedClassTemplate()) {
            sl = ctd->getTemplateParameters()->getTemplateLoc();
        } else {
            sl = rd->getSourceRange().getBegin();
        }

    } else {
        // some other decl -- just give the beginning

        sl = d->getSourceRange().getBegin();
    }

    return sl;
}

/// get next character and sourcelocation, while skipping comments

SourceLocation GeneralVisitor::getNextLoc(SourceLocation sl, bool forward) {
    return ::getNextLoc(TheRewriter.getSourceMgr(), sl, forward);
}

char GeneralVisitor::getChar(SourceLocation sl) {
    return ::getChar(TheRewriter.getSourceMgr(), sl);
}

// Find the location of the next searched for char.
SourceLocation GeneralVisitor::findChar(SourceLocation sloc, char ct) {
    return ::findChar(TheRewriter.getSourceMgr(), sloc, ct);
}

/// Skip paren expression following sl, points after the paren

SourceLocation GeneralVisitor::skipParens(SourceLocation sl) {
    return ::skipParens(TheRewriter.getSourceMgr(), sl);
}

/// Get next word starting from sl

std::string GeneralVisitor::getNextWord(SourceLocation sl, SourceLocation *end) {
    return ::getNextWord(TheRewriter.getSourceMgr(), sl, end);
}

/// Get prev word starting from sl -

std::string GeneralVisitor::getPreviousWord(SourceLocation sl, SourceLocation *start) {
    while (std::isspace(getChar(sl)))
        sl = getNextLoc(sl, false); // skip spaces

    SourceLocation e = sl;
    char c = getChar(sl);
    if (std::isalnum(c) || c == '_') {
        while (sl.isValid() && (std::isalnum(c) || c == '_')) {
            sl = getNextLoc(sl, false);
            c = getChar(sl);
        }
        sl = getNextLoc(sl); // 1 step too much
    }
    if (start != nullptr)
        *start = sl;
    return TheRewriter.getRewrittenText(SourceRange(sl, e));
}

////////////////////////////////////////////////////////////////////////////
/// This processes references to non-field variables within site loops
/// if is_assign==true, this is assigned to with assignop and assign_stmt contains
/// the full assignment op
////////////////////////////////////////////////////////////////////////////

var_info *GeneralVisitor::handle_var_ref(DeclRefExpr *DRE, bool is_assign,
                                         const std::string &assignop, Stmt *assign_stmt,
                                         bool is_raw) {

    if (isa<VarDecl>(DRE->getDecl())) {
        auto decl = dyn_cast<VarDecl>(DRE->getDecl());

        /// we don't want "X" -variable or lattice-> as a kernel parameter
        clang::QualType typ =
            decl->getType().getUnqualifiedType().getNonReferenceType();
        typ.removeLocalConst();

        // if (typ.getAsString(PP) == "lattice_struct *") llvm::errs() << "GOT
        // LATTICE_STRUCT PTR!!!\n";

        if (typ.getAsString(PP) == "X_index_type" ||
            typ.getAsString(PP) == "lattice_struct *")
            return nullptr;

        var_ref vr;
        vr.ref = DRE;
        // vr.ind = writeBuf->markExpr(DRE);
        vr.is_assigned = is_assign;
        if (is_assign)
            vr.assignop = assignop;

        bool foundvar = false;
        var_info *vip = nullptr;
        for (var_info &vi : var_info_list) {
            if (vi.decl == decl && !is_raw) {
                // found already referred to decl
                // check if this particular ref has been handled before
                bool foundref = false;
                for (auto &r : vi.refs)
                    if (r.ref == DRE) {
                        foundref = true;
                        // if old check was not assignment and this is, change status
                        // can happen if var ref is a function "out" argument
                        if (r.is_assigned == false && is_assign == true) {
                            r.is_assigned = true;
                            r.assignop = assignop;
                        }
                        break;
                    }
                if (!foundref) {
                    // a new reference
                    vi.refs.push_back(vr);
                }
                vi.is_assigned |= is_assign;
                if (is_top_level && vi.reduction_type == reduction::NONE) {
                    vi.reduction_type = get_reduction_type(is_assign, assignop, vi);
                }
                vip = &vi;
                foundvar = true;
                break;

            } else if (vi.decl == decl) {
                // now is_raw == true, old variable, do not have to check refs
                foundvar = true;
                break;
            }
        }
        if (!foundvar) {
            // new variable referred to
            vip = new_var_info(decl);

            vip->refs.push_back(vr);
            vip->is_assigned = is_assign;

            if (is_raw) {
                vip->is_raw = true;
                vip->is_site_dependent = true; // probably does not matter
            } else {

                // we know refs contains only 1 element
                if (is_top_level) {
                    vip->reduction_type = get_reduction_type(is_assign, assignop, *vip);
                }
            }
        }

        if (!is_raw && is_assign && assign_stmt != nullptr && !vip->is_site_dependent) {
            vip->is_site_dependent =
                is_rhs_site_dependent(assign_stmt, &vip->dependent_vars);

            // llvm::errs() << "Var " << vip->name << " depends on site: " <<
            // vip->is_site_dependent <<  "\n";
        }
        return vip;

    } else {
        // end of VarDecl - how about other decls, e.g. functions?
        reportDiag(DiagnosticsEngine::Level::Error, DRE->getSourceRange().getBegin(),
                   "Reference to unimplemented (non-variable) type");
    }

    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
///  Insert the new variable info

var_info *GeneralVisitor::new_var_info(VarDecl *decl) {

    var_info vi;
    vi.refs = {};
    vi.decl = decl;
    vi.name = decl->getName().str();
    // Printing policy is somehow needed for printing type without "class" id
    // Unqualified takes away "consts" and namespaces etc. and Canonical typdefs/using.
    // Also need special handling for element type
    clang::QualType type =
        decl->getType().getUnqualifiedType().getCanonicalType().getNonReferenceType();
    type.removeLocalConst();

    vi.type = type.getAsString(PP);
    vi.type = remove_extra_whitespace(vi.type);
    // bool is_elem = (vi.type.find("element<") == 0);
    // vi.type = type.getAsString(PP);
    // if (is_elem)
    //     vi.type = "element<" + vi.type + ">";
    // llvm::errs() << " + Got " << vi.type << '\n';

    // is it loop-local?
    vi.is_loop_local = false;
    for (var_decl &d : var_decl_list) {
        if (d.scope >= 0 && vi.decl == d.decl) {
            // llvm::errs() << "loop local var ref! " << vi.name << '\n';
            vi.is_loop_local = true;
            break;
        }
    }
    vi.is_site_dependent = false; // default case
    vi.reduction_type = reduction::NONE;

    vi.is_special_reduction_type = (vi.type.find("Reduction<") == 0);
    if (vi.is_special_reduction_type) {
        // The signature of the type is Reduction< ... , void>.
        // Remove until last ,

        vi.type = vi.type.substr(10, vi.type.rfind(',') - 10);
    }

    vi.dependent_vars.clear();

    var_info_list.push_back(vi);
    return &(var_info_list.back());
}

//////////////////////////////////////////////////////////////////////////////
/// Now variable is a (loop) local var, add to list
//////////////////////////////////////////////////////////////////////////////

var_info *GeneralVisitor::add_var_to_decl_list(VarDecl *var, int scope_level) {

    // Now it should be automatic local variable decl
    var_decl vd;
    vd.decl = var;
    vd.name = var->getName().str();
    vd.type = var->getType().getAsString();
    vd.scope = scope_level;

    var_decl_list.push_back(vd);

    // insert this to var_info_list too

    var_info *ip = new_var_info(var);
    ip->reduction_type = reduction::NONE;

    // finally, check initialization
    if (var->hasInit()) {
        ip->is_site_dependent = is_site_dependent(var->getInit(), &ip->dependent_vars);
        ip->is_assigned = true;
    } else {
        ip->is_assigned = false;
    }

    //  llvm::errs() << "Local var decl " << vd.name << " of type " << vd.type << '\n';

    return ip;
}

///////////////////////////////////////////////////////////////////////////////////
/// Pragma handling: has_pragma()
///
///
/// Check if the SourceLocation l is preceded by "#pragma hila" on previous line.
/// There cannot be anything except whitespace between l and the beginning of line
/// cannot allow templates because conditionals may contain <> -chars
/// Pragma_args will point to the beginning of arguments of pragma
///////////////////////////////////////////////////////////////////////////////////

bool GeneralVisitor::has_pragma(Stmt *S, const pragma_hila pragma, const char **arg) {

    if (S == nullptr)
        return false;

    return has_pragma(S->getSourceRange().getBegin(), pragma, arg);
}

/// For functiondecl, go through the previous decls (prototype!) too if needed

bool GeneralVisitor::has_pragma(FunctionDecl *F, const pragma_hila pragma,
                                const char **arg) {

    if (F == nullptr)
        return false;
    if (has_pragma(F->getSourceRange().getBegin(), pragma, arg))
        return true;

    if (FunctionDecl *proto = F->getPreviousDecl()) {
        return has_pragma(proto, pragma, arg);
    }
    return false;
}

bool GeneralVisitor::has_pragma(Decl *F, const pragma_hila pragma, const char **arg) {

    if (F == nullptr)
        return false;
    return has_pragma(F->getSourceRange().getBegin(), pragma, arg);
}

/// And here is the main interface to pragma

bool GeneralVisitor::has_pragma(const SourceLocation l, const pragma_hila pragma,
                                const char **pragma_arg) {
    std::string arg;
    SourceLocation pragmaloc, sl = l;

    if (l.isInvalid())
        return false;

    // if macro, get the unexpanded loc
    if (sl.isMacroID()) {
        CharSourceRange CSR = TheRewriter.getSourceMgr().getImmediateExpansionRange(sl);
        sl = CSR.getBegin();
    }

    if (has_pragma_hila(TheRewriter.getSourceMgr(), sl, pragma, pragmaloc,
                        pragma_arg)) {

        // got it, comment out -- check that it has not been commented out before
        // the buffer may not be writeBuf, so be careful

        if (cmdline::comment_pragmas) {
            srcBuf *sb = get_file_srcBuf(pragmaloc);

            int loc;
            if (sb != nullptr) {
                loc = sb->find_original(pragmaloc, '#');
            }

            if (sb == nullptr || loc < 0) {
                llvm::errs() << "internal error in pragma handling\n";
                exit(1);
            }
            std::string s = sb->get(loc, loc + 1);
            if (s.at(0) == '#')
                sb->insert(loc, "//-- ", true, false);
        }

        return true;
    }

    return false;
}
