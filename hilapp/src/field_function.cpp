//------------------------------------------------------------------------------
// Field functions: functions which have a parameter of type X_index_type and
// access fields given as parameters with f[X] or f[X + dir], e.g.
//
//   template <typename T>
//   T gradient(const Field<T> &f, const X_index_type &X, Direction d) {
//       return (f[X + d] - f[X - d]) * 0.5;
//   }
//   ...
//   onsites(ALL) g[X] = gradient(f, X, e_x);
//
// The function definition is rewritten in place, in one of two forms:
//
// Value form, if the directions of all accesses are constants or parameters of the
// function (as above): each access f[X], f[X + dir] becomes a parameter, which gets
// the field value.  The call is rewritten as if it were
//
//   gradient(f[X + e_x], f[X - e_x], X, e_x)
//
// and the accesses are handled as field references of the loop: only the needed
// directions are communicated, and the values are fetched as usual.  The Field
// parameters are removed, X becomes the site index (for X.coordinates() etc.).
//
// Accessor form, if a direction is computed in the function (e.g. foralldir(d) inside):
// the accesses become calls of Field::get_value_at() and Field::get_value_at_nb(), and
// at the call site the fields and X are replaced with the loop variables.  The fields
// are merged into the field info of the loop, all directions are communicated.
//
// On GPUs the function becomes a __device__ function.  In the accessor form the Field
// parameters become field_storage<T> -references (the kernels get field_storage), and
// the accesses use the same expressions as the kernels.  The functions called from it
// are marked as loop functions.
//
// Limitations: const Field parameters, f[X] and f[X + dir] only (no offsets), no
// nested field function calls.  Loops calling field functions are not vectorized.
//------------------------------------------------------------------------------

#include <algorithm>
#include <sstream>
#include <string>

#include "hilapp.h"
#include "toplevelvisitor.h"
#include "stringops.h"

// the name of the site index parameter in the rewritten function
static const std::string ff_index(name_prefix + "_ff_index");

// the names of the value parameters in the value form
static std::string ff_value_name(int k) {
    return name_prefix + "_ff_v" + std::to_string(k);
}

// all field functions in this compilation unit.  Use list, pointers to elements are kept
static std::list<field_function_info> field_functions = {};

void clear_field_functions_in_compilation_unit() {
    field_functions.clear();
}

//////////////////////////////////////////////////////////////////////////////
/// Type of the parameter: Field, X or other
//////////////////////////////////////////////////////////////////////////////

static ff_param get_param_kind(const ParmVarDecl *P, const PrintingPolicy &PP) {
    std::string t =
        P->getType().getNonReferenceType().getCanonicalType().getUnqualifiedType().getAsString(PP);
    if (t == "X_index_type")
        return ff_param::X;
    if (t.find(field_type) == 0)
        return ff_param::FIELD;
    return ff_param::OTHER;
}

/// field function if it has an X_index_type parameter.  Methods and operators (e.g.
/// Field::operator[](X_index_type)) are not field functions
static bool is_field_function_decl(const FunctionDecl *D, const PrintingPolicy &PP) {
    if (D == nullptr || isa<CXXMethodDecl>(D) || D->isOverloadedOperator())
        return false;
    for (const ParmVarDecl *P : D->parameters())
        if (get_param_kind(P, PP) == ff_param::X)
            return true;
    return false;
}

bool TopLevelVisitor::is_field_function(const FunctionDecl *D) {
    return is_field_function_decl(D, PP);
}

/// The X + dir expression of f[X + dir]: return the + or - operator call, nullptr if
/// not of this form
static CXXOperatorCallExpr *get_X_plus_dir_op(Expr *idxE) {
    Expr *e = idxE->IgnoreParens()->IgnoreImplicit();
    auto *Op = dyn_cast<CXXOperatorCallExpr>(e);
    if (!Op)
        if (auto *Ce = dyn_cast<CXXConstructExpr>(e))
            if (Ce->getNumArgs() == 1)
                Op = dyn_cast<CXXOperatorCallExpr>(Ce->getArg(0)->IgnoreImplicit());
    if (Op && Op->getNumArgs() != 2)
        return nullptr;
    return Op;
}

/// f[..] -expression, return the operator[] call or nullptr
static CXXOperatorCallExpr *get_field_subscript(Expr *E, GeneralVisitor &v) {
    auto *OC = dyn_cast<CXXOperatorCallExpr>(E->IgnoreParens()->IgnoreImplicit());
    if (OC && strcmp(getOperatorSpelling(OC->getOperator()), "[]") == 0 &&
        v.is_field_expr(OC->getArg(0)))
        return OC;
    return nullptr;
}

//////////////////////////////////////////////////////////////////////////////
/// Visitor which analyzes and rewrites the body of a field function
//////////////////////////////////////////////////////////////////////////////

class FieldFunctionVisitor : public GeneralVisitor,
                             public RecursiveASTVisitor<FieldFunctionVisitor> {
  public:
    using GeneralVisitor::GeneralVisitor;

    srcBuf *sb;
    FunctionDecl *fd;
    field_function_info *ffi;

    /// f[X] and f[X + dir] in the function, rewritten after all are found
    struct found_access {
        CXXOperatorCallExpr *OC; // the f[..] expression
        int param;               // Field parameter
        Expr *dirE;              // dir, nullptr for f[X]
        bool negate;             // X - dir
    };
    std::vector<found_access> found;

    template <typename visitortype>
    FieldFunctionVisitor(visitortype &v) : GeneralVisitor(v) {}

    /// parameter index of DeclRefExpr to a parameter of fd, -1 if not
    int param_index(Expr *E) {
        // X is often copied with an implicit constructor
        if (auto *DRE = dyn_cast<DeclRefExpr>(E->IgnoreUnlessSpelledInSource())) {
            for (unsigned i = 0; i < fd->getNumParams(); i++)
                if (DRE->getDecl() == fd->getParamDecl(i))
                    return i;
        }
        return -1;
    }

    bool is_X_param(Expr *E) {
        int i = param_index(E);
        return i >= 0 && ffi->param_kind[i] == ff_param::X;
    }

    void error(SourceLocation sl, const char *msg, const char *s = nullptr) {
        ffi->ok = false;
        auto &DE = Context->getDiagnostics();
        auto ID = DE.getCustomDiagID(DiagnosticsEngine::Level::Error, "%0");
        DE.Report(sl, ID) << (std::string(msg) + (s ? s : ""));
    }

    /// f[X] and f[X + dir] in the CPU or GPU code
    std::string value_at_X(const std::string &f) {
        if (target.kernelize)
            return f + ".get(" + ff_index + ", _dev_field_alloc_size)";
        else
            return f + ".get_value_at(" + ff_index + ")";
    }

    std::string value_at_nb(const std::string &f, const std::string &dir) {
        if (!target.kernelize)
            return f + ".get_value_at_nb(" + dir + ", " + ff_index + ")";
        // same as in the kernels, see codegen_gpu.cpp
        if (is_macro_defined("LAYOUT_ON_THE_FLY"))
            return f + ".get(hila_dev_neighbour(" + f + ".neighbours[" + dir +
                   "], hila_dev_site(" + ff_index + "), " + ff_index + ", " + dir +
                   "), _dev_field_alloc_size)";
        else
            return f + ".get(" + f + ".neighbours[" + dir + "][" + ff_index +
                   "], _dev_field_alloc_size)";
    }

    /// Can the access be passed as a value: f[X], or the direction is a constant or a
    /// Direction parameter.  Fill in acc
    bool get_value_access(const found_access &fa, ff_access &acc) {
        acc.param = fa.param;
        acc.atX = (fa.dirE == nullptr);
        acc.is_const = false;
        acc.const_value = 0;
        acc.dir_param = -1;
        acc.negate = fa.negate;
        if (acc.atX)
            return true;
        if (fa.dirE->isIntegerConstantExpr(*Context)) {
            auto v = fa.dirE->getIntegerConstantExpr(*Context);
            acc.is_const = true;
            acc.const_value = v->getExtValue() + (fa.negate ? 50 : 0); // as in field_ref
            acc.const_str = std::string(fa.negate ? "-" : "") + "Direction(" +
                            std::to_string(v->getExtValue()) + ")";
            return true;
        }
        int i = param_index(fa.dirE);
        if (i >= 0 && ffi->param_kind[i] == ff_param::OTHER &&
            get_expr_type(fa.dirE) == "Direction") {
            acc.dir_param = i;
            return true;
        }
        return false;
    }

    static bool same_access(const ff_access &a, const ff_access &b) {
        return a.param == b.param && a.atX == b.atX && a.is_const == b.is_const &&
               (a.atX || (a.is_const ? a.const_value == b.const_value
                                     : a.dir_param == b.dir_param && a.negate == b.negate));
    }

    /// Rewrite the accesses: the value form if possible, otherwise access the fields
    void rewrite_accesses() {
        std::vector<ff_access> acc(found.size());
        ffi->value_form = true;
        for (unsigned i = 0; i < found.size(); i++)
            if (!get_value_access(found[i], acc[i]))
                ffi->value_form = false;

        for (unsigned i = 0; i < found.size(); i++) {
            CXXOperatorCallExpr *OC = found[i].OC;
            if (ffi->value_form) {
                // the same access is the same parameter
                unsigned k;
                for (k = 0; k < ffi->accesses.size(); k++)
                    if (same_access(ffi->accesses[k], acc[i]))
                        break;
                if (k == ffi->accesses.size())
                    ffi->accesses.push_back(acc[i]);
                sb->replace(OC->getSourceRange(), ff_value_name(k));
            } else {
                std::string name = get_stmt_str(OC->getArg(0)->IgnoreImplicit());
                if (found[i].dirE == nullptr) {
                    sb->replace(OC->getSourceRange(), value_at_X(name));
                } else {
                    std::string dir = get_stmt_str(found[i].dirE);
                    if (found[i].negate)
                        dir = "-(" + dir + ")";
                    sb->replace(OC->getSourceRange(), value_at_nb(name, dir));
                }
            }
        }
    }

    /// handle f[X] and f[X + dir].  Returns true if E is a field[] expression
    bool handle_field_ref(Expr *E) {
        CXXOperatorCallExpr *OC = get_field_subscript(E, *this);
        if (!OC)
            return false;

        Expr *nameE = OC->getArg(0)->IgnoreImplicit();
        Expr *idxE = OC->getArg(1)->IgnoreImplicit();

        int p = param_index(nameE);
        if (p < 0 || ffi->param_kind[p] != ff_param::FIELD) {
            error(nameE->getBeginLoc(),
                  "field functions can access only Field parameters of the function");
            return true;
        }

        std::string idxtype = get_expr_type(idxE);

        if (idxtype == "X_index_type") {
            if (!is_X_param(idxE)) {
                error(idxE->getBeginLoc(), "index must be the X parameter of the function");
                return true;
            }
            ffi->read_atX[p] = true;
            found.push_back({OC, p, nullptr, false});

        } else if (idxtype == "X_plus_direction") {
            // X + dir or X - dir, possibly wrapped in a constructor
            CXXOperatorCallExpr *Op = get_X_plus_dir_op(idxE);
            if (!Op || !is_X_param(Op->getArg(0))) {
                error(idxE->getBeginLoc(),
                      "field functions accept only f[X + dir] or f[X - dir] with the X parameter");
                return true;
            }

            bool negate = (strcmp(getOperatorSpelling(Op->getOperator()), "-") == 0);
            ffi->read_nb[p] = true;
            found.push_back({OC, p, Op->getArg(1)->IgnoreImplicit(), negate});

        } else if (idxtype == "X_plus_offset") {
            error(idxE->getBeginLoc(),
                  "offsets f[X + offset] are not (yet) supported in field functions, "
                  "only f[X] and f[X + dir]");
        } else {
            error(idxE->getBeginLoc(), "field functions accept only f[X] and f[X + dir]");
        }
        return true;
    }

    /// handle X.coordinates() etc.  Returns true if Call is a method of X
    bool handle_X_method(CXXMemberCallExpr *MCall) {
        Expr *obj = MCall->getImplicitObjectArgument();
        if (!obj || get_expr_type(obj) != "X_index_type")
            return false;

        if (!is_X_param(obj)) {
            error(obj->getBeginLoc(), "X must be the X parameter of the function");
            return true;
        }

        std::string name = MCall->getMethodDecl()->getNameAsString();
        // the same functions as in the loops, see handle_special_loop_function()
        std::string start = target.kernelize ? "get_dev_" : "lattice.ref().";
        std::string repl;
        if (name == "coordinates")
            repl = start + "coordinates(" + ff_index;
        else if (name == "parity")
            repl = start + "site_parity(" + ff_index;
        else if (name == "coordinate")
            repl = start + "coordinate(" + ff_index + ", ";
        else if (name == "x" || name == "y" || name == "z" || name == "t")
            repl = start + "coordinate(" + ff_index + ", e_" + name;
        else {
            error(MCall->getBeginLoc(), "method not supported in field functions: X.",
                  name.c_str());
            return true;
        }

        // replace up to and including the open paren
        SourceLocation sl = findChar(MCall->getBeginLoc(), '(');
        sb->replace(SourceRange(MCall->getBeginLoc(), sl), repl);

        for (Expr *a : MCall->arguments())
            TraverseStmt(a);
        return true;
    }

    bool TraverseStmt(Stmt *s) {
        if (s == nullptr)
            return true;
        if (Expr *E = dyn_cast<Expr>(s)) {
            if (handle_field_ref(E))
                return true;
            if (auto *MC = dyn_cast<CXXMemberCallExpr>(E))
                if (handle_X_method(MC))
                    return true;
            // hila::global variable g() is the device variable in __device__ functions
            if (target.kernelize)
                if (auto *CE = dyn_cast<CallExpr>(E))
                    if (handle_global_var_method_call(CE, true))
                        return true;
        }
        RecursiveASTVisitor<FieldFunctionVisitor>::TraverseStmt(s);
        return true;
    }

    /// other references to X or Field parameters are errors
    bool VisitDeclRefExpr(DeclRefExpr *E) {
        int p = param_index(E);
        if (p >= 0 && ffi->param_kind[p] != ff_param::OTHER) {
            error(E->getBeginLoc(),
                  "in field functions X and Field parameters can be used only as f[X], "
                  "f[X + dir] or X.method(): ",
                  get_stmt_str(E).c_str());
        } else if (is_field_expr(E) || get_expr_type(E) == "X_index_type") {
            error(E->getBeginLoc(),
                  "field functions can use only the Field and X parameters of the function: ",
                  get_stmt_str(E).c_str());
        }
        return true;
    }

    bool VisitCallExpr(CallExpr *C) {
        FunctionDecl *D = C->getDirectCallee();
        if (D && D != fd) {
            for (const ParmVarDecl *P : D->parameters())
                if (get_param_kind(P, PP) == ff_param::X) {
                    error(C->getBeginLoc(),
                          "calling field functions from field functions is not (yet) supported");
                    break;
                }
        }
        return true;
    }
};

//////////////////////////////////////////////////////////////////////////////
/// Collect the functions and constructors called from a field function.  On GPUs
/// these must be loop functions (__device__).  The field accesses and X methods are
/// rewritten, these are not collected
//////////////////////////////////////////////////////////////////////////////

class FieldFunctionCallees : public GeneralVisitor,
                             public RecursiveASTVisitor<FieldFunctionCallees> {
  public:
    using GeneralVisitor::GeneralVisitor;

    std::vector<Decl *> callees;

    template <typename visitortype>
    FieldFunctionCallees(visitortype &v) : GeneralVisitor(v) {}

    void add(Decl *D) {
        if (D == nullptr || srcMgr.isInSystemHeader(D->getBeginLoc()))
            return;
        if (std::find(callees.begin(), callees.end(), D) == callees.end())
            callees.push_back(D);
    }

    bool TraverseStmt(Stmt *s) {
        if (s == nullptr)
            return true;
        if (Expr *E = dyn_cast<Expr>(s)) {
            if (CXXOperatorCallExpr *OC = get_field_subscript(E, *this)) {
                // only the direction expression remains in the code
                if (CXXOperatorCallExpr *Op = get_X_plus_dir_op(OC->getArg(1)))
                    TraverseStmt(Op->getArg(1));
                return true;
            }
            if (auto *MC = dyn_cast<CXXMemberCallExpr>(E)) {
                Expr *obj = MC->getImplicitObjectArgument();
                if (obj && get_expr_type(obj) == "X_index_type") {
                    for (Expr *a : MC->arguments())
                        TraverseStmt(a);
                    return true;
                }
            }
            // hila::global g() is replaced with the device variable
            if (auto *OC = dyn_cast<CXXOperatorCallExpr>(E))
                if (OC->getOperator() == OverloadedOperatorKind::OO_Call &&
                    get_expr_type(OC->getArg(0)).find("hila::global<") != std::string::npos)
                    return true;
        }
        RecursiveASTVisitor<FieldFunctionCallees>::TraverseStmt(s);
        return true;
    }

    bool VisitCallExpr(CallExpr *C) {
        FunctionDecl *D = C->getDirectCallee();
        if (D && !is_field_function_decl(D, PP))
            add(D);
        return true;
    }

    bool VisitCXXConstructExpr(CXXConstructExpr *C) {
        add(C->getConstructor());
        return true;
    }
};

/// Element type of the Field parameter P as written: Field<double> -> double, and for
/// templates Field<T> -> T.  Empty if not found
static std::string field_elem_type(const ParmVarDecl *P, const PrintingPolicy &PP) {
    QualType t = P->getType().getNonReferenceType();
    std::string elem;
    if (!t->isDependentType()) {
        // Field<elem>
        std::string s = t.getCanonicalType().getUnqualifiedType().getAsString(PP);
        elem = s.substr(field_type.size(), s.rfind('>') - field_type.size());
    } else {
        // template: take the template argument as written, through alias templates
        const TemplateSpecializationType *tst = t->getAs<TemplateSpecializationType>();
        while (tst && tst->isTypeAlias())
            tst = tst->getAliasedType()->getAs<TemplateSpecializationType>();
        if (tst == nullptr || tst->template_arguments().size() != 1)
            return "";
        llvm::raw_string_ostream os(elem);
        tst->template_arguments()[0].print(PP, os, true);
        os.flush();
    }
    return elem;
}

//////////////////////////////////////////////////////////////////////////////
/// Analyze and rewrite the field function D (or the template where D is an instance).
/// This is done once for each function
//////////////////////////////////////////////////////////////////////////////

field_function_info *TopLevelVisitor::get_field_function(FunctionDecl *D) {

    // with templates the source is the template pattern.  The instance has the same
    // source locations, use it in the analysis
    FunctionDecl *pattern = D->getTemplateInstantiationPattern();
    if (pattern == nullptr)
        pattern = D;
    if (pattern->getDefinition() != nullptr)
        pattern = pattern->getDefinition();

    // On GPUs the functions called from field function must be __device__.  Templates
    // may call different functions in each instance, do this for each
    auto mark_callees = [&](field_function_info *ffi) {
        FunctionDecl *def = D->getDefinition();
        if (!target.kernelize || def == nullptr || !def->doesThisDeclarationHaveABody())
            return;
        if (std::find(ffi->instances.begin(), ffi->instances.end(), def) != ffi->instances.end())
            return;
        ffi->instances.push_back(def);

        FieldFunctionCallees fc(*this);
        fc.TraverseStmt(def->getBody());
        for (Decl *c : fc.callees)
            visit_loop_function_decl(c);
    };

    for (field_function_info &f : field_functions)
        if (f.decl == pattern) {
            mark_callees(&f);
            return &f;
        }

    field_functions.emplace_back();
    field_function_info *ffi = &field_functions.back();
    ffi->decl = pattern;
    ffi->ok = true;
    ffi->contains_random = false;
    ffi->value_form = false;

    for (ParmVarDecl *P : D->parameters()) {
        ff_param k = get_param_kind(P, PP);
        ffi->param_kind.push_back(k);
        ffi->read_atX.push_back(false);
        ffi->read_nb.push_back(false);

        if (k == ff_param::FIELD) {
            QualType t = P->getType();
            if (!t->isReferenceType() || !t.getNonReferenceType().isConstQualified()) {
                reportDiag(DiagnosticsEngine::Level::Error, P->getBeginLoc(),
                           "Field parameters of field functions must be const references "
                           "(prototype limitation)");
                ffi->ok = false;
            }
        }
    }

    if (std::count(ffi->param_kind.begin(), ffi->param_kind.end(), ff_param::X) != 1) {
        reportDiag(DiagnosticsEngine::Level::Error, D->getBeginLoc(),
                   "field function must have exactly one X_index_type parameter");
        ffi->ok = false;
    }

    // the declaration with the body (hasBody() is true also for prototypes)
    FunctionDecl *body_decl = D->getDefinition();

    if (body_decl == nullptr || !body_decl->doesThisDeclarationHaveABody()) {
        reportDiag(DiagnosticsEngine::Level::Error, D->getBeginLoc(),
                   "the definition of field function '%0' must be visible",
                   D->getNameAsString().c_str());
        ffi->ok = false;
        return ffi;
    }

    srcBuf *sb = get_file_srcBuf(body_decl->getBeginLoc());
    if (sb == nullptr) {
        reportDiag(DiagnosticsEngine::Level::Error, D->getBeginLoc(),
                   "field function '%0' is in a system file", D->getNameAsString().c_str());
        ffi->ok = false;
        return ffi;
    }

    ffi->contains_random = contains_random(body_decl->getBody());

    // rewrite the body
    FieldFunctionVisitor ffv(*this);
    ffv.sb = sb;
    ffv.fd = body_decl;
    ffv.ffi = ffi;
    ffv.TraverseStmt(body_decl->getBody());
    ffv.rewrite_accesses();

    // and replace the X parameter with the site index in all declarations.  In the value
    // form the Field parameters are replaced with the values of the accesses.  Otherwise
    // on GPUs the function gets field_storage instead of Field
    for (FunctionDecl *R : pattern->redecls()) {
        srcBuf *rb = get_file_srcBuf(R->getBeginLoc());
        if (rb == nullptr)
            continue;
        for (unsigned i = 0; i < R->getNumParams() && i < ffi->param_kind.size(); i++) {
            ParmVarDecl *P = R->getParamDecl(i);
            if (ffi->param_kind[i] == ff_param::X) {
                rb->replace(P->getSourceRange(), "const unsigned " + ff_index);
            } else if (ffi->param_kind[i] == ff_param::FIELD &&
                       (ffi->value_form || target.kernelize)) {
                std::string elem = field_elem_type(P, PP);
                if (elem.empty()) {
                    reportDiag(DiagnosticsEngine::Level::Error, P->getBeginLoc(),
                               "could not resolve the element type of Field parameter");
                    ffi->ok = false;
                } else if (!ffi->value_form) {
                    rb->replace(P->getSourceRange(), "const field_storage<" + elem + "> &" +
                                                         P->getNameAsString());
                } else {
                    std::string params;
                    for (unsigned k = 0; k < ffi->accesses.size(); k++)
                        if (ffi->accesses[k].param == (int)i)
                            params += (params.empty() ? "" : ", ") + std::string("const ") +
                                      elem + " &" + ff_value_name(k);
                    if (params.empty())
                        rb->remove_with_comma(P->getSourceRange()); // field not used
                    else
                        rb->replace(P->getSourceRange(), params);
                }
            }
        }
        // comment before "template <..>" if it is a template, __device__ after it
        std::string comment = ffi->value_form
                                  ? "// hilapp: field function, field values are parameters\n"
                                  : "// hilapp: field function, X parameter is the site index\n";
        std::string device = target.kernelize ? "__device__ " : "";
        if (FunctionTemplateDecl *T = R->getDescribedFunctionTemplate()) {
            rb->insert(T->getBeginLoc(), comment, true, true);
            if (!device.empty())
                rb->insert(R->getInnerLocStart(), device, true, true);
        } else {
            rb->insert(R->getBeginLoc(), comment + device, true, true);
        }
    }

    mark_callees(ffi);

    return ffi;
}

//////////////////////////////////////////////////////////////////////////////
/// Call to a field function inside a site loop.  Returns false if Call is not a
/// field function call
//////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::handle_field_function_call(CallExpr *Call) {

    FunctionDecl *D = Call->getDirectCallee();
    if (!is_field_function(D))
        return false;

    if (isa<CXXMemberCallExpr>(Call) || isa<CXXOperatorCallExpr>(Call)) {
        reportDiag(DiagnosticsEngine::Level::Error, Call->getBeginLoc(),
                   "field functions must be non-member functions (prototype limitation)");
        return true;
    }

    if (target.openacc) {
        reportDiag(DiagnosticsEngine::Level::Error, Call->getBeginLoc(),
                   "field functions are not supported with OpenACC");
        return true;
    }

    field_function_info *ffi = get_field_function(D);

    if (ffi->contains_random)
        loop_info.contains_random = true;

    for (unsigned i = 0; i < Call->getNumArgs(); i++) {
        Expr *arg = Call->getArg(i);
        if (isa<CXXDefaultArgExpr>(arg))
            continue;

        ff_param kind = i < ffi->param_kind.size() ? ffi->param_kind[i] : ff_param::OTHER;

        if (kind == ff_param::FIELD) {
            if (!is_field_expr(arg->IgnoreImplicit())) {
                reportDiag(DiagnosticsEngine::Level::Error, arg->getBeginLoc(),
                           "argument must be a Field variable");
                continue;
            }
            FieldRefChecker frc(*this);
            frc.TraverseStmt(arg);
            if (frc.isLoopLocal()) {
                reportDiag(DiagnosticsEngine::Level::Error, arg->getBeginLoc(),
                           "Field argument cannot depend on variable '%0' defined inside "
                           "onsites-loop",
                           frc.getLocalVarInfo()->name.c_str());
            }
        } else if (kind == ff_param::X) {
            Expr *e = arg->IgnoreParens()->IgnoreImplicit();
            if (!isa<DeclRefExpr>(e) || !is_X_type(e)) {
                reportDiag(DiagnosticsEngine::Level::Error, arg->getBeginLoc(),
                           "argument must be X (prototype limitation)");
            }
        } else {
            TraverseStmt(arg);
        }
    }

    field_function_call fc;
    fc.call = Call;
    fc.ffi = ffi;

    // value form: the accesses are field references of the loop, as if the call were
    // func(f[X + dir], ...).  These have no expressions in the loop, see field_ref
    if (ffi->value_form) {
        for (ff_access &a : ffi->accesses) {
            field_ref r;
            r.nameExpr = Call->getArg(a.param)->IgnoreImplicit();
            r.ff_call = Call;
            r.is_read = true;
            r.is_written = false;
            r.sequence = parsing_state.stmt_sequence;

            if (!a.atX) {
                r.is_direction = true;
                r.dir_negate = a.negate;
                Expr *dirE = nullptr;
                if (!a.is_const) {
                    dirE = Call->getArg(a.dir_param);
                    if (auto *DA = dyn_cast<CXXDefaultArgExpr>(dirE))
                        dirE = DA->getExpr();
                    dirE = dirE->IgnoreImplicit();
                }

                if (a.is_const) {
                    r.is_constant_direction = true;
                    r.constant_value = a.const_value;
                    r.direxpr_s = a.const_str;
                } else if (dirE->isIntegerConstantExpr(*Context)) {
                    // constant argument
                    auto v = dirE->getIntegerConstantExpr(*Context);
                    r.is_constant_direction = true;
                    r.constant_value = v->getExtValue() + (a.negate ? 50 : 0);
                    r.direxpr_s = std::string(a.negate ? "-" : "") + "Direction(" +
                                  std::to_string(v->getExtValue()) + ")";
                } else if (isa<CXXDefaultArgExpr>(Call->getArg(a.dir_param))) {
                    reportDiag(DiagnosticsEngine::Level::Error, Call->getBeginLoc(),
                               "non-constant default Direction argument of field function is "
                               "not supported");
                    continue;
                } else {
                    r.dirE = dirE;
                    r.direxpr_s =
                        std::string(a.negate ? "-" : "") + "(" + get_stmt_str(dirE) + ")";
                    FieldRefChecker frc(*this);
                    frc.TraverseStmt(dirE);
                    r.is_loop_local_dir = frc.isLoopLocal();
                }
            }
            field_ref_list.push_back(r);
            fc.refs.push_back(&field_ref_list.back());
        }
    }

    field_function_call_list.push_back(fc);

    return true;
}

//////////////////////////////////////////////////////////////////////////////
/// Replace the Field and X arguments of field function calls in the loop body.
/// Called after the field names are set
//////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::replace_field_function_args(srcBuf &loopBuf) {
    for (field_function_call &fc : field_function_call_list) {

        if (fc.ffi->value_form) {
            // rewrite the whole call: the field values, the site index and the other
            // arguments as they are in the loop now.  Inner calls in the arguments are
            // before this in the list, these are already rewritten
            if (fc.refs.size() != fc.ffi->accesses.size())
                continue; // there was an error
            std::vector<std::string> args;
            for (unsigned i = 0; i < fc.ffi->param_kind.size(); i++) {
                if (fc.ffi->param_kind[i] == ff_param::FIELD) {
                    for (unsigned k = 0; k < fc.ffi->accesses.size(); k++)
                        if (fc.ffi->accesses[k].param == (int)i)
                            args.push_back(fc.refs[k]->replacement);
                } else if (fc.ffi->param_kind[i] == ff_param::X) {
                    args.push_back(looping_var);
                } else if (i < fc.call->getNumArgs() &&
                           !isa<CXXDefaultArgExpr>(fc.call->getArg(i))) {
                    args.push_back(loopBuf.get(fc.call->getArg(i)->getSourceRange()));
                }
            }
            std::string text = loopBuf.get(fc.call->getCallee()->getSourceRange()) + "(";
            for (unsigned i = 0; i < args.size(); i++)
                text += (i > 0 ? ", " : "") + args[i];
            loopBuf.replace(fc.call->getSourceRange(), text + ")");
            continue;
        }

        for (unsigned i = 0; i < fc.call->getNumArgs() && i < fc.ffi->param_kind.size(); i++) {
            Expr *arg = fc.call->getArg(i);
            if (isa<CXXDefaultArgExpr>(arg))
                continue;
            if (fc.ffi->param_kind[i] == ff_param::FIELD) {
                for (field_info &l : field_info_list)
                    if (is_duplicate_expr(l.nameExpr, arg->IgnoreImplicit())) {
                        loopBuf.replace(arg, l.new_name);
                        break;
                    }
            } else if (fc.ffi->param_kind[i] == ff_param::X) {
                loopBuf.replace(arg, looping_var);
            }
        }
    }
}
