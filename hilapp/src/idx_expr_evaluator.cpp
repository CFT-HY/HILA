#include "idx_expr_evaluator.h"

#include <clang/AST/ASTContext.h>
#include <clang/AST/ParentMapContext.h>

#include <llvm/ADT/Hashing.h>

#include <algorithm>
#include <limits>
#include <unordered_set>

using namespace clang;
using namespace llvm;

static constexpr raw_string_ostream::Colors RED = raw_string_ostream::RED;
static constexpr raw_string_ostream::Colors BLUE = raw_string_ostream::BLUE;
static constexpr raw_string_ostream::Colors GREEN = raw_string_ostream::GREEN;
static constexpr raw_string_ostream::Colors CYAN = raw_string_ostream::CYAN;

static inline raw_ostream &color_msg(raw_ostream &os, raw_ostream::Colors color, const char *msg,
                                     bool bold = true, const char *file = nullptr, int line = -1) {
    os.changeColor(color, bold);
    os << msg;
    os.resetColor();
    if (file)
        os << "(" << file;
    if (line != -1)
        os << ":" << line << ")";
    os << ": ";
    return os;
}

static std::string diag_msg;
static raw_string_ostream diags(diag_msg);

#define DIAGS_TO_STR
#ifdef DIAGS_TO_STR
#define DIAGS_OUT diags
#else
#define DIAGS_OUT llvm::outs()
#endif
#define INTERNAL_NOTE color_msg(DIAGS_OUT, BLUE, "INTERNAL NOTE", true, __FILE__, __LINE__)
#define INTERNAL_ERROR color_msg(DIAGS_OUT, RED, "INTERNAL ERROR", true, __FILE__, __LINE__)
#define DEBUG_OUT color_msg(DIAGS_OUT, CYAN, "INTERNAL DEBUG", true, __FILE__, __LINE__)

// Keep track of for which exprs diags have been reported. So we dont spam them multiple times..
// Should this be moved into Eval_state?
static std::set<const Expr *> reported_diag;
#define INTERNAL_NOTE_IF(E)                                                                        \
    if (reported_diag.find(E) == reported_diag.end() &&                                            \
        reported_diag.insert(E)                                                                    \
            .second /*this should always return true, its here just to insert the Expr*/)          \
    INTERNAL_NOTE

#define loc_str(S) (S)->getSourceRange().getBegin().printToString(ASTctx.getSourceManager())
#define loc_strp(S) (S)->getSourceRange().getBegin().printToString(ASTctx->getSourceManager())

template <unsigned N>
void report_diag(const ASTContext &ctx, DiagnosticsEngine::Level lev, const SourceLocation &SL,
                 const char (&msg)[N], const char *s1 = nullptr, const char *s2 = nullptr,
                 const char *s3 = nullptr) {
    auto &DE = ctx.getDiagnostics();
    auto ID = DE.getCustomDiagID(lev, msg);
    auto DB = DE.Report(SL, ID);
    if (s1 != nullptr)
        DB.AddString(s1);
    if (s2 != nullptr)
        DB.AddString(s2);
    if (s3 != nullptr)
        DB.AddString(s3);
}

bool eval_CFG_block(const CFGBlock &block, Eval_state &state, ASTContext &ASTctx);
bool eval_stmt(const Stmt *S, Eval_state &state, ASTContext &ASTctx);
bool eval_expr(Eval_result &result, const Expr *E, Eval_state &state, ASTContext &ASTctx);
static bool eval_constexpr_function(const FunctionDecl *func, Eval_state &state,
                                    ASTContext &ASTctx);

// The values are represented with int64_t, but follow the semantics of the C++ types of the
// expressions. to_int_type() converts a value to a type as c++ would, int_binary_op() does
// the binop evaluation.
//
// The value `v` converted to type T: truncated to the width of T and sign or zero extended.
// Unknown if T is not an integer type, or if the value does not fit in int64_t.
static std::optional<int64_t> to_int_type(std::optional<int64_t> v, QualType T,
                                          const ASTContext &ASTctx) {
    if (!v.has_value())
        return std::nullopt;
    T = T.getNonReferenceType().getCanonicalType();
    if (T->isBooleanType())
        return v.value() != 0;
    if (!T->isIntegralOrEnumerationType())
        return std::nullopt;
    const uint64_t width = ASTctx.getIntWidth(T);
    const bool is_signed = T->isSignedIntegerOrEnumerationType();
    uint64_t u = static_cast<uint64_t>(v.value());
    if (width < 64) {
        const uint64_t mask = (uint64_t(1) << width) - 1;
        u &= mask;
        if (is_signed && ((u >> (width - 1)) & 1))
            u |= ~mask; // sign extend
        return static_cast<int64_t>(u);
    }
    if (width > 64 || (!is_signed && static_cast<int64_t>(u) < 0))
        return std::nullopt;
    return static_cast<int64_t>(u);
}

// Evaluate binary operator `op` for integers `a` and `b` (also used for the compound assignment
// versions). Gives unknown if the result would overflow int64_t or is undefined (division by zero,
// shift out of range), and then sets `*undefined` to the reason. Gives unknown also
// for operators that do not give an integer value. The result has to be converted to the type of
// the expression with to_int_type().
static std::optional<int64_t> int_binary_op(BinaryOperatorKind op, int64_t a, int64_t b,
                                            const char **undefined = nullptr) {
    auto undefined_result = [&](const char *reason) -> std::optional<int64_t> {
        if (undefined)
            *undefined = reason;
        return std::nullopt;
    };
    constexpr int64_t int64_min = std::numeric_limits<int64_t>::min();
    int64_t r;
    switch (op) {
        // the __builtin_* return true if the operation results into an overflow
    case BO_Mul:
    case BO_MulAssign:
        return __builtin_mul_overflow(a, b, &r) ? undefined_result("overflows int64_t")
                                                : std::optional<int64_t>(r);
    case BO_Add:
    case BO_AddAssign:
        return __builtin_add_overflow(a, b, &r) ? undefined_result("overflows int64_t")
                                                : std::optional<int64_t>(r);
    case BO_Sub:
    case BO_SubAssign:
        return __builtin_sub_overflow(a, b, &r) ? undefined_result("overflows int64_t")
                                                : std::optional<int64_t>(r);
    case BO_Div:
    case BO_DivAssign:
        if (b == 0)
            return undefined_result("division by zero");
        if (a == int64_min && b == -1)
            return undefined_result("overflows int64_t");
        return a / b;
    case BO_Rem:
    case BO_RemAssign:
        if (b == 0)
            return undefined_result("division by zero");
        if (a == int64_min && b == -1)
            return undefined_result("overflows int64_t");
        return a % b;
    case BO_Shl:
    case BO_ShlAssign:
        if (a < 0)
            return undefined_result("left shift of a negative value");
        if (b < 0 || b >= 64)
            return undefined_result("shift amount out of range");
        return static_cast<int64_t>(static_cast<uint64_t>(a) << b);
    case BO_Shr:
    case BO_ShrAssign:
        if (b < 0 || b >= 64)
            return undefined_result("shift amount out of range");
        return a >> b;
    case BO_LT:
        return a < b;
    case BO_GT:
        return a > b;
    case BO_LE:
        return a <= b;
    case BO_GE:
        return a >= b;
    case BO_EQ:
        return a == b;
    case BO_NE:
        return a != b;
    case BO_And:
    case BO_AndAssign:
        return a & b;
    case BO_Xor:
    case BO_XorAssign:
        return a ^ b;
    case BO_Or:
    case BO_OrAssign:
        return a | b;
    default:
        // e.g. pointer to member. Also the C++20 three-way comparison `a <=> b` (BO_Cmp) on
        // purpose: its type is the class std::strong_ordering, so its value would be unknown
        // anyway after to_int_type() ...
        return std::nullopt;
    }
}

// Helper for int_binary_op() that reports a note for expr E, if the result is unknown
static std::optional<int64_t> eval_int_binary_op(BinaryOperatorKind op, int64_t a, int64_t b,
                                                 const Expr *E, ASTContext &ASTctx) {
    const char *undefined = nullptr;
    std::optional<int64_t> r = int_binary_op(op, a, b, &undefined);
    if (undefined) {
        INTERNAL_NOTE_IF(E) << "`" << BinaryOperator::getOpcodeStr(op) << "` with values " << a
                            << " and " << b << ": " << undefined << ", the result is unknown "
                            << loc_str(E) << "\n";
    }
    return r;
}

// Helper that adds the result of Expr E to the state.possible_vals_for_expr
// if we are tracking the Expr E
inline void track_expr(Eval_result &result, const Expr *E, Eval_state &state) {
    auto vals = state.possible_vals_for_expr.find(E);
    if (vals == state.possible_vals_for_expr.end())
        return; // we are not tracking this Expr

    vals->second.insert(result.val);
}

// Reference variables, e.g. `int &r = x;`: a write through `r` changes `x`, and a write to `x`
// changes the value of `r`. A reference can not be bound to another variable later, so when its
// target can be determined from its initializer, every access to `r` is an access to `x` (see
// alias_reference_target()) and nothing about references needs to be in the evaluation state.
// If the target can not be determined, e.g. `int &r = c ? x : z;`, reading the reference gives an
// unknown value and a write through it makes all variables unknown.

// True for a variable of reference type to an integer type (what we track) that refers to some
// other variable. Not for function parameters: their argument values are copied in, see the
// CallExpr case of eval_expr(). Not for references bound to a temporary, e.g.
// `const int &r = x + 1;`, which just hold a value.
static bool is_alias_reference(const VarDecl *VD) {
    if (!VD->getType()->isReferenceType() || isa<ParmVarDecl>(VD) ||
        !VD->getType().getNonReferenceType()->isIntegralOrEnumerationType())
        return false;
    const Expr *init = VD->getInit();
    if (init == nullptr)
        return true; // e.g. `extern int &r;`, refers to something we do not know
    if (const auto *cleanups = dyn_cast<ExprWithCleanups>(init))
        init = cleanups->getSubExpr();
    return !isa<MaterializeTemporaryExpr>(init->IgnoreParens());
}

// The variable the lvalue expression E refers to, nullptr if it is not a variable we can
// determine. E.g. `x`, `(x)`, `++x`, `x = 3`, `x += 3` and `(f(), x)` all refer to `x`.
static const VarDecl *lvalue_var(const Expr *E, ASTContext &ASTctx) {
    while (true) {
        E = E->IgnoreParenNoopCasts(ASTctx);
        if (const auto *cleanups = dyn_cast<ExprWithCleanups>(E)) {
            E = cleanups->getSubExpr();
        } else if (const auto *DRE = dyn_cast<DeclRefExpr>(E)) {
            return dyn_cast<VarDecl>(DRE->getDecl());
        } else if (const auto *UO = dyn_cast<UnaryOperator>(E);
                   UO && (UO->getOpcode() == UO_PreInc || UO->getOpcode() == UO_PreDec)) {
            E = UO->getSubExpr();
        } else if (const auto *BO = dyn_cast<BinaryOperator>(E); BO && BO->isAssignmentOp()) {
            E = BO->getLHS();
        } else if (BO && BO->getOpcode() == BO_Comma) {
            E = BO->getRHS();
        } else {
            return nullptr;
        }
    }
}

#define REFERENCE_CHAIN_MAX_DEPTH 16

// The variable the alias reference VD refers to, nullptr if it can not be determined from the
// initializer, e.g. `int &r = c ? x : z;`
static const VarDecl *alias_reference_target(const VarDecl *VD, ASTContext &ASTctx) {
    // Follow references to references, `int &r2 = r;`. Limit the depth for `int &r = r;`.
    for (int depth = 0;
         depth < REFERENCE_CHAIN_MAX_DEPTH && VD != nullptr && is_alias_reference(VD); ++depth) {
        if (VD->getInit() == nullptr)
            return nullptr;
        VD = lvalue_var(VD->getInit(), ASTctx);
    }
    return (VD != nullptr && !is_alias_reference(VD)) ? VD : nullptr;
}

// A write through an alias reference whose target we do not know: it may have changed any
// variable.
static void write_through_unknown_reference(Eval_state &state) {
    for (auto &[VD, value] : state.vars) {
        value = std::nullopt;
        state.all_possible_vals[VD].insert(std::nullopt);
    }
}

// Reference arguments of calls: an argument bound to a non-const reference to an integer, e.g.
// `x` in f(x) with void f(int &a), can be changed by the call. Such an argument is a non-const
// glvalue, also for int && and forwarding references (std::move(x)), function pointers, member
// functions, operators and constructors. An argument passed by value is a prvalue after the
// lvalue-to-rvalue conversion, one bound to a const reference has a const type.
// The exception is a temporary bound to a reference, eg 5 in f(5) with void f(int &&a): is a
// non-const glvalue too, but no variable refers to it.
static bool is_writable_int_reference_arg(const Expr *arg) {
    return arg->isGLValue() && !arg->getType().isConstQualified() &&
           arg->getType()->isIntegralOrEnumerationType() &&
           !isa<MaterializeTemporaryExpr>(arg->IgnoreParens());
}

// The variable a reference argument refers to, also through reference variables. nullptr if it
// can not be determined, e.g. `c ? x : y`.
static const VarDecl *referenced_var(const Expr *arg, ASTContext &ASTctx) {
    const VarDecl *VD = lvalue_var(arg, ASTctx);
    if (VD != nullptr && is_alias_reference(VD))
        VD = alias_reference_target(VD, ASTctx);
    return VD;
}

// Update the state after a call for a reference argument `arg`, eg `x` in f(x) with
// void f(int &a): the variable that `arg` refers to now has the value `value`.
//  - For an evaluated constexpr call, `value` is the value of the reference parameter `a` at
//    the end of the callee
//  - For a call that is not evaluated, `value` is unknown (std::nullopt).
// If we can not determine which variable `arg` refers to, the call may have changed any of
// them, so all variables become unknown
static void write_reference_arg(const Expr *arg, std::optional<int64_t> value, Eval_state &state,
                                ASTContext &ASTctx) {
    if (const VarDecl *VD = referenced_var(arg, ASTctx)) {
        state.vars[VD] = to_int_type(value, VD->getType(), ASTctx);
        state.all_possible_vals[VD].insert(state.vars[VD]);
    } else {
        write_through_unknown_reference(state);
    }
}

// recursively eval
bool eval_expr(Eval_result &result, const Expr *E, Eval_state &state, ASTContext &ASTctx) {
    result.val = std::nullopt;

    if (state.block_expr_state.find(E) != state.block_expr_state.end()) {
        // we have already evaluated this expression
        result.val = state.block_expr_state[E];
        state.last_res = result;
        return true;
    }

    // Every time we 'successfully' evaluate something we possibly need to save the result of this
    // expression to state.possible_vals_for_expr, and return true: `goto successful_return` does
    // this.

    if (const auto *IL = dyn_cast<IntegerLiteral>(E)) {

        result.val = IL->getValue().getSExtValue();
        goto successful_return;
    } else if (const auto *BL = dyn_cast<CXXBoolLiteralExpr>(E)) {
        result.val = BL->getValue();
        goto successful_return;
    } else if (const auto *CL = dyn_cast<CharacterLiteral>(E)) {
        result.val = CL->getValue();
        goto successful_return;
    } else if (const auto *default_arg = dyn_cast<CXXDefaultArgExpr>(E)) {
        // This takes care of functions with default arguments
        Eval_result default_val{};
        if (!eval_expr(default_val, default_arg->getExpr(), state, ASTctx))
            return false;
        result = default_val;
        goto successful_return;

    } else if (isa<MaterializeTemporaryExpr, ExprWithCleanups>(E)) {
        // A temporary, eg bound to a reference `const int &r = x + 1;`, has the value of its
        // subexpression. ExprWithCleanups marks where temporaries are destroyed
        const Expr *sub = isa<MaterializeTemporaryExpr>(E)
                              ? cast<MaterializeTemporaryExpr>(E)->getSubExpr()
                              : cast<ExprWithCleanups>(E)->getSubExpr();
        if (!eval_expr(result, sub, state, ASTctx))
            return false;
        goto successful_return;

    } else if (const auto *PE = dyn_cast<ParenExpr>(E)) {
        Eval_result paren_expr_val{};
        if (!eval_expr(paren_expr_val, PE->getSubExpr(), state, ASTctx))
            return false;
        result = paren_expr_val;
        goto successful_return;

    } else if (const auto *cast = dyn_cast<CastExpr>(E)) {
        // Implicit and explicit casts. Must come before the DeclRefExpr case below, which looks
        // through implicit casts and would otherwise skip e.g. the int -> bool conversion in
        // `bool t = x;`. Evaluate one cast level at a time so nested casts are all applied.
        if (!eval_expr(result, cast->getSubExpr(), state, ASTctx))
            return false;

        switch (cast->getCastKind()) {
        case CK_NoOp:
        case CK_LValueToRValue:
        case CK_IntegralCast:
        case CK_IntegralToBoolean:
        case CK_UserDefinedConversion: // the value of the conversion function call
        case CK_ConstructorConversion:
            // the value is converted to the type of the cast at successful_return
            break;
        default:
            INTERNAL_NOTE_IF(E) << "Not evaluating cast " << cast->getCastKindName()
                                << ", the value is unknown " << loc_str(E) << "\n";
            result.val = std::nullopt;
        }
        goto successful_return;

    } else if (const auto *DRE = dyn_cast<DeclRefExpr>(E->IgnoreParenImpCasts())) {

        // check if its an enum constant
        if (const auto *EC = dyn_cast<EnumConstantDecl>(DRE->getDecl())) {
            result.val = EC->getInitVal().getSExtValue();
            goto successful_return;
        }

        const auto *VD = dyn_cast<VarDecl>(DRE->getDecl());
        if (!VD) {
            // e.g. a function name, should be nothing we care about
            INTERNAL_NOTE_IF(E) << " Not a variable, the value is unknown " << loc_str(DRE) << "\n";
            goto successful_return;
        }

        if (is_alias_reference(VD)) {
            // Reading the reference reads the variable it refers to
            VD = alias_reference_target(VD, ASTctx);
            if (VD == nullptr) {
                INTERNAL_NOTE_IF(E) << " Reference variable that refers to an unknown variable, "
                                       "its value is unknown "
                                    << loc_str(DRE) << "\n";
                goto successful_return;
            }
        }

        if (state.vars.find(VD) == state.vars.end()) {
            // A constant declared outside the evaluated code, e.g. `constexpr int N = 4;`
            if (VD->isUsableInConstantExpressions(ASTctx)) {
                if (const APValue *value = VD->evaluateValue(); value && value->isInt()) {
                    result.val = value->getInt().tryExtValue();
                    goto successful_return;
                }
            }

            INTERNAL_NOTE_IF(E) << " the variable has an unknown value " << loc_str(DRE) << "\n"
                                << "             : VarDecl " << loc_str(VD) << " " << VD << "\n";
            goto successful_return;
        }
        result.val = state.vars[VD];
        goto successful_return;

    } else if (const auto *BO = dyn_cast<BinaryOperator>(E)) {
        const auto *LHS = BO->getLHS();
        const auto *RHS = BO->getRHS();
        if (BO->isAssignmentOp()) {
            // TODO: for now LHS has to be a simple VarRefDecl
            const auto *dreLHS = dyn_cast<DeclRefExpr>(LHS->IgnoreParenImpCasts());
            if (!dreLHS) {
                // The assigned value is not tracked, but both sides have to be evaluated for
                // their side effects, e.g. `arr[i] = i++`. Since C++17 the RHS is evaluated
                // before the LHS. (previously undefined behaviour)
                Eval_result rhs_res{};
                if (!eval_expr(rhs_res, RHS, state, ASTctx))
                    return false;
                Eval_result lhs_res{};
                if (!eval_expr(lhs_res, LHS, state, ASTctx))
                    return false;

                INTERNAL_NOTE_IF(E)
                    << "Only tracking assignments if LHS is DeclRefExpr: Here it is instead `"
                    << LHS->IgnoreParenBaseCasts()->getStmtClassName()
                    << "` Not tracking the assigned value. " << loc_str(LHS) << "\n";
                // TODO: maybe some cases could be evaluated?
                result.val = std::nullopt;
                goto successful_return;
            }
            const auto *VD = dyn_cast<VarDecl>(dreLHS->getDecl());
            if (!VD) {
                INTERNAL_ERROR << "Expecting binary assignment operators LHS DeclRefExpr to have "
                                  "an associated VarDecl."
                               << loc_str(dreLHS) << "\n";
                return false;
            }

            if (is_alias_reference(VD)) {
                // Assigning through the reference assigns to the variable it refers to
                const VarDecl *target = alias_reference_target(VD, ASTctx);
                if (target == nullptr) {
                    Eval_result rhs_res{};
                    if (!eval_expr(rhs_res, RHS, state, ASTctx))
                        return false;
                    write_through_unknown_reference(state);
                    INTERNAL_NOTE_IF(E) << " Assignment through a reference variable that refers "
                                           "to an unknown variable, all variables are now unknown "
                                        << loc_str(BO) << "\n";
                    goto successful_return;
                }
                VD = target;
            }

            Eval_result rhs_res{};
            if (!eval_expr(rhs_res, RHS, state, ASTctx))
                return false;
            // if rhs has no value should set and return nullopt in all cases
            if (!rhs_res.val.has_value()) {
                state.vars[VD] = std::nullopt;
                result.val = std::nullopt;
                // var changes, record that it may be unknown
                state.all_possible_vals[VD].insert(state.vars[VD]);
                goto successful_return;
            }

            auto op = BO->getOpcode();
            switch (op) {
                //// [C99 6.5.16] Assignment operators.
                // BO_Assign // "="
            case BO_Assign: {
                state.vars[VD] = to_int_type(rhs_res.val, VD->getType(), ASTctx);
                result = rhs_res;
                // var changes
                state.all_possible_vals[VD].insert(state.vars[VD]);
                goto successful_return;
            }
            case BO_MulAssign:
            case BO_DivAssign:
            case BO_RemAssign:
            case BO_AddAssign:
            case BO_SubAssign:
            case BO_ShlAssign:
            case BO_ShrAssign:
            case BO_AndAssign:
            case BO_XorAssign:
            case BO_OrAssign: {
                // the result is converted to the type of the variable, eg unsigned wraps
                if (state.vars[VD].has_value())
                    state.vars[VD] = to_int_type(eval_int_binary_op(op, state.vars[VD].value(),
                                                                    rhs_res.val.value(), E, ASTctx),
                                                 VD->getType(), ASTctx);
                result.val = state.vars[VD];
                state.all_possible_vals[VD].insert(state.vars[VD]);
                goto successful_return;
            }
            default:
                INTERNAL_ERROR << " eval_expr BO opcode not implemented `" << BO->getOpcodeStr()
                               << "` " << loc_str(BO) << "\n";
                return false;
            }
        } else {
            auto op = BO->getOpcode();
            // special cases || and && are handeled on the level of CFG blocks,
            // not evaluated recursively.
            // ASSUMPTION: assumes the CFGBlock's are layed out as in LLVM(21.1.8), i.e. LAnd/LOr
            // are preceded by the last relevant boolean expression.
            //// [C99 6.5.13] Logical AND operator.
            //// [C99 6.5.14] Logical OR operator.
            if (op == BO_LAnd || op == BO_LOr) {
                // We can reach this with an unknown last_res, so dont try to get a concrete value
                // out of the std::optional here.
                result.val = state.last_res.val;
                goto successful_return;
            }

            // LHS evaluated first. Evaluate both before checking for unknown values: the side
            // effects of the RHS (e.g. `y + (i++)`) happen even if the LHS value is unknown.
            Eval_result lhs_res{};
            if (!eval_expr(lhs_res, LHS, state, ASTctx))
                return false;

            Eval_result rhs_res{};
            if (!eval_expr(rhs_res, RHS, state, ASTctx))
                return false;

            //// [C99 6.5.17] Comma operator.
            // BO_Comma // ","
            if (op == BO_Comma) {
                // comma operator evals both and returns the rhs, whether the lhs is known or not
                // This is before the !lhs_res.val.has_value() check below so that exprs like
                // (unknown, 123) is still evaluated to the known value 123.
                result = rhs_res;
                goto successful_return;
            }

            if (!lhs_res.val.has_value() || !rhs_res.val.has_value())
                goto successful_return; // We dont know the values just return the nullopt

            // the result is converted to the type of the expression at successful_return
            result.val =
                eval_int_binary_op(op, lhs_res.val.value(), rhs_res.val.value(), E, ASTctx);
            goto successful_return;
        }
    } else if (const auto *UO = dyn_cast<UnaryOperator>(E)) {
        const auto *operand = UO->getSubExpr();
        auto op = UO->getOpcode();

        if (op == UO_PostInc || op == UO_PostDec || op == UO_PreInc || op == UO_PreDec) {
            // TODO: for now operand has to be a simple VarRefDecl
            const auto *dre = dyn_cast<DeclRefExpr>(operand->IgnoreParenImpCasts());
            if (!dre) {
                INTERNAL_ERROR
                    << "Expecting unary operator operands to be DeclRefExpr: insted it is `"
                    << operand->IgnoreParenBaseCasts()->getStmtClassName() << "` "
                    << loc_str(operand) << "\n";
                return false;
            }
            const auto *VD = dyn_cast<VarDecl>(dre->getDecl());
            if (!VD) {
                INTERNAL_ERROR << "Expecting unary operator operands DeclRefExpr to have an "
                                  "associated VarDecl."
                               << loc_str(operand) << "\n";
                return false;
            }

            if (is_alias_reference(VD)) {
                // Incrementing through the reference increments the variable it refers to
                const VarDecl *target = alias_reference_target(VD, ASTctx);
                if (target == nullptr) {
                    write_through_unknown_reference(state);
                    INTERNAL_NOTE_IF(E) << " Increment / decrement through a reference variable "
                                           "that refers to an unknown variable, all variables "
                                           "are now unknown "
                                        << loc_str(UO) << "\n";
                    goto successful_return;
                }
                VD = target;
            }

            {
                // the new value is converted to the type of the variable
                const std::optional<int64_t> before = state.vars[VD];
                std::optional<int64_t> after;
                if (before.has_value()) {
                    after = to_int_type(
                        eval_int_binary_op((op == UO_PostInc || op == UO_PreInc) ? BO_Add : BO_Sub,
                                           before.value(), 1, E, ASTctx),
                        VD->getType(), ASTctx);
                }
                state.vars[VD] = after;
                // postfix gives the value before, prefix the value after
                result.val = (op == UO_PostInc || op == UO_PostDec) ? before : after;
                state.all_possible_vals[VD].insert(after);
                goto successful_return;
            }

        } else {

            if (op == UO_AddrOf) {
                // var we track could be changed through pointer, we do not follow these for now..
                // needs to be checked also when the value of the variable is unknown
                if (const VarDecl *VD = lvalue_var(operand, ASTctx);
                    VD && VD->getType().getNonReferenceType()->isIntegralOrEnumerationType()) {
                    INTERNAL_ERROR << " Taking the address of an integer variable is not "
                                      "supported "
                                   << loc_str(UO) << "\n";
                    return false;
                }
            }

            // eval operand
            Eval_result operand_res{};
            if (!eval_expr(operand_res, operand, state, ASTctx))
                return false;

            // now we dont know, so return the nullopt
            if (!operand_res.val.has_value())
                goto successful_return;

            // clang-format off
            switch (op) {
                // [C99 6.5.3.3] Unary arithmetic
            case UO_Plus:  { result.val = +operand_res.val.value(); goto successful_return; } break;
            case UO_Minus: { result.val = eval_int_binary_op(BO_Sub, 0, operand_res.val.value(), E, ASTctx); goto successful_return; } break; // We can reuse the int_binary_op here.. (0 - a)
            case UO_Not:   { result.val = ~operand_res.val.value(); goto successful_return; } break;
            case UO_LNot:  { result.val = !operand_res.val.value(); goto successful_return; } break;
                // __extension__ marker.
                // UO_Extension "__extension__"
            case UO_Extension: { result = operand_res; goto successful_return; } break; // the __extension__ marker is just for suppressing compiler messages, so it is just a no-op
                // [C++ Coroutines] co_await operator
                // UO_Coawait "co_await"
                // [C99 6.5.3.2] Address and indirection
                // UO_AddrOf "&"
            case UO_AddrOf: // of something we do not track, checked above
                // UO_Deref "*"
            case UO_Deref:
                // "__real expr"/"__imag expr" Extension.
                // UO_Real "__real"
            case UO_Real:
                // UO_Imag "__imag"
            case UO_Imag:
                result.val = std::nullopt;
                goto successful_return;
                // clang-format on
            default:
                INTERNAL_ERROR << " TODO: eval_expr UO opcode not implemented `"
                               << UO->getOpcodeStr(op) << "` " << loc_str(UO) << "\n";
                return false;
            }
        }
    } else if (dyn_cast<ConditionalOperator>(E)) {
        // this is the ternary op (x?a:b)
        // Similarly as the short circuiting for LAnd/LOr, the result is the last evaluated
        // expression. We can reach this with an unknown last_res.
        result.val = state.last_res.val;
        goto successful_return;

    } else if (const auto *CE = dyn_cast<CallExpr>(E)) {
        const auto *func = CE->getDirectCallee();

        bool do_eval = true;
        // func is nullptr for calls through a function pointer
        if (!func || !func->isConstexpr()) {

            INTERNAL_NOTE_IF(E) << "Call to a function that is not a constexpr " << loc_str(CE)
                                << "\n";
            do_eval = false;
        }

        // the call needs a separate Eval_state
        Eval_state call_state;

        // eval the argument values to the function call
        // CE->arguments() map to func->parameters() one to one
        // NOTE edge cases:
        //  - call to a member operator, e.g. `obj(a)` or `obj + a` (CXXOperatorCallExpr), has
        //    the object as its first argument. A member function call `obj.f(a)`
        //    (CXXMemberCallExpr) does not, and the parameters never include the implicit 'this'
        //    pointer.
        //    TODO: C++23 explicit object parameters `this Self &self` are not handled
        //  - C-style variadic function has more arguments than parameters.
        // Arguments without a parameter are only evaluated for their side effects.
        // TODO: the object is not passed to member functions, so its members are unknown there
        const auto *method = dyn_cast_or_null<CXXMethodDecl>(func);
        const unsigned first_param_arg =
            (isa<CXXOperatorCallExpr>(CE) && method != nullptr && !method->isStatic()) ? 1 : 0;

        // arguments bound to non-const references to integers, which the call can change, with
        // the parameter (nullptr if none or the call is not evaluated)
        std::vector<std::pair<const Expr *, const ParmVarDecl *>> reference_args;

        unsigned i = 0;
        for (const Expr *arg_expr : CE->arguments()) {
            // the ->arguments() contains CXXDefaultArgExpr which are now taken care of the
            // eval_expr. So default arguments are taken care of.
            Eval_result arg_val{};
            // Even if we dont want to eval the call we want to try eval the argument expressions
            // which can affect the result..
            if (!eval_expr(arg_val, arg_expr, state, ASTctx))
                return false;
            const ParmVarDecl *param = nullptr;
            if (do_eval && i >= first_param_arg && i - first_param_arg < func->getNumParams()) {
                param = func->getParamDecl(i - first_param_arg);
                // A worklist path may reach a constexpr call with an unknown
                // argument. Preserve that fact in the callee instead of
                // throwing std::bad_optional_access here.
                call_state.vars[param] = arg_val.val;
            }
            if (is_writable_int_reference_arg(arg_expr))
                reference_args.push_back({arg_expr, param});
            i++;
        }

        if (do_eval) {
            auto *func_body = func->getBody();
            if (!func_body) {
                INTERNAL_ERROR << " no body found for function " << loc_str(func) << "\n";
                return false;
            }

            if (state.call_depth >= MAX_CALL_DEPTH) {
                INTERNAL_ERROR << " constexpr function calls nested deeper than " << MAX_CALL_DEPTH
                               << ", giving up. " << loc_str(CE) << "\n";
                return false;
            }
            call_state.call_depth = state.call_depth + 1;

            // The callee is evaluated with the same worklist traversal as exec()
            call_state.max_block_states = state.max_block_states;
            call_state.func_cache = state.func_cache;
            if (!eval_constexpr_function(func, call_state, ASTctx))
                return false;

            // Write back the values of reference parameters, e.g. `dir` in
            // constexpr Direction &operator++(Direction &dir). The callee state holds the values
            for (const auto &[arg, param] : reference_args) {
                std::optional<int64_t> value;
                if (param != nullptr) {
                    if (const auto it = call_state.vars.find(param); it != call_state.vars.end())
                        value = it->second;
                }
                write_reference_arg(arg, value, state, ASTctx);
            }

            // The callee can execute successfully while still having an
            // unknown return value (when argument is unknown)
            result.val = call_state.return_val;
            goto successful_return;
        }

        // A function we do not evaluate may change the variables passed by reference
        for (const auto &[arg, param] : reference_args) {
            INTERNAL_NOTE_IF(arg) << "Variable passed by reference to a function that is not "
                                     "evaluated, its value is unknown after the call "
                                  << loc_str(arg) << "\n";
            write_reference_arg(arg, std::nullopt, state, ASTctx);
        }
        goto successful_return;

    } else if (const auto *LAM = dyn_cast<LambdaExpr>(E)) {
        // TODO: Some lambda functions could be analysed
        // But things get a bit tricky
        // clang-format off
        //     int dir = 0; // = the variable we wish to track
        //     auto sneaky = [&]() { dir = dir + 1; };
        //     for (int i = 0; i < 10; ++i) {
        //         // If we encounter an function call that does not take any of the tracked variables by ref
        //         // or pointer we can safely assume they do not affect the variables we care about.
        //         // HOWEVER: if the call is lambda function call: CXX_something_something then it may take
        //         // the variables we care about by ref and change them. We need to detect this. Need to
        //         // analyze the function call and see if it modifies any of the variables we need to track.
        //         // For simplicity we can just not allow these kind of non constexpr lambdas for now. TODO.
        //         sneaky();
        //     }
        // clang-format on
        INTERNAL_ERROR << " TODO: lambda functions not implemented yet" << loc_str(LAM) << "\n";
        return false;
    } else if (const auto *constructor_expr = dyn_cast<CXXConstructExpr>(E)) {
        // The possible constructor argument have to be evaluated.
        // For example if we encounter something crazy like `some_constructor( x=1 );` where `x` is
        // a variable we care about.
        for (auto arg_expr : constructor_expr->arguments()) {
            Eval_result arg_val{};
            if (!eval_expr(arg_val, arg_expr, state, ASTctx))
                return false;
        }
        // the constructor is not evaluated, it may change the variables passed by reference
        for (auto arg_expr : constructor_expr->arguments()) {
            if (is_writable_int_reference_arg(arg_expr)) {
                INTERNAL_NOTE_IF(arg_expr) << "Variable passed by reference to a constructor, its "
                                              "value is unknown after the call "
                                           << loc_str(arg_expr) << "\n";
                write_reference_arg(arg_expr, std::nullopt, state, ASTctx);
            }
        }

        // TODO(CXXConstructExpr):  ... not trying to eval these for now

        INTERNAL_NOTE_IF(E) << " Not trying to eval CXXConstructExpr. Ignoring for now "
                            << loc_str(constructor_expr) << "\n";
        goto successful_return;
    } else if (const auto *array_subscript = dyn_cast<ArraySubscriptExpr>(E)) {
        // the LHS[RHS] expressions have to be evaluated since they might affect variables we care
        // about E.g. `array[x++]`
        auto *LHS = array_subscript->getLHS();
        Eval_result lhs_res{};
        if (!eval_expr(lhs_res, LHS, state, ASTctx))
            return false;
        Eval_result rhs_res{};
        auto *RHS = array_subscript->getRHS();
        if (!eval_expr(rhs_res, RHS, state, ASTctx))
            return false;

        // TODO: ... not trying to eval these for now

        INTERNAL_NOTE_IF(E) << " Not trying to eval ArraySubscriptExpr. Ignoring for now "
                            << loc_str(array_subscript) << "\n";
        goto successful_return;

    } else if (const auto *member = dyn_cast<MemberExpr>(E)) {
        // The base has to be evaluated since it might affect variables we care about, e.g.
        // `arr[x++].m`. `this` has no side effects.
        const Expr *base = member->getBase();
        if (!isa<CXXThisExpr>(base->IgnoreParenImpCasts())) {
            Eval_result base_res{};
            if (!eval_expr(base_res, base, state, ASTctx))
                return false;
        }

        // TODO: ... not trying to eval these for now

        INTERNAL_NOTE_IF(E) << " Not trying to eval MemberExpr. Ignoring for now "
                            << loc_str(member) << "\n";
        goto successful_return;
    }

    // Other kinds of expressions.
    { // this catches e.g. `sizeof(int)`
        Expr::EvalResult constant;
        if (E->EvaluateAsInt(constant, ASTctx)) {
            result.val = constant.Val.getInt().tryExtValue();
            goto successful_return;
        }
    }
    // cases where the expression is a type we dont care about, but its subexpressions may still
    // change tracked vars, e.g. `S{x++}`
    if (!isa<CXXNewExpr, CXXDeleteExpr, CXXThrowExpr, StmtExpr, BlockExpr, AtomicExpr,
             CoroutineSuspendExpr, CXXTypeidExpr, GenericSelectionExpr, UnaryExprOrTypeTraitExpr,
             CXXNoexceptExpr>(E) &&
        std::all_of(E->child_begin(), E->child_end(),
                    [](const Stmt *child) { return child == nullptr || isa<Expr>(child); })) {
        for (const Stmt *child : E->children()) {
            if (const auto *child_expr = dyn_cast_or_null<Expr>(child)) {
                Eval_result child_res{};
                if (!eval_expr(child_res, child_expr, state, ASTctx))
                    return false;
            }
        }
        INTERNAL_NOTE_IF(E) << " Not trying to eval " << E->getStmtClassName()
                            << ", the value is unknown " << loc_str(E) << "\n";
        goto successful_return;
    }

    INTERNAL_ERROR << "TODO: UNHANDELED EXPR: eval_expr() " << loc_str(E) << "\n";
    E->dump(diags, ASTctx);
    return false;

successful_return:
    // the value as the type of the expression, unknown if not an integer type
    result.val = to_int_type(result.val, E->getType(), ASTctx);
    track_expr(result, E, state);
    state.block_expr_state[E] = result.val;
    state.last_res = result;

    // record the values of the tuple based on E, see add_tuple_to_tracking_list(). The child_exprs
    // have been evaluated before in this block, so their values are in block_expr_state.
    if (!state.tracked_tuples.empty()) {
        if (const auto tuple = state.tracked_tuples.find(E); tuple != state.tracked_tuples.end()) {
            Value_tuple values;
            for (const Expr *member : tuple->second) {
                const auto it = state.block_expr_state.find(member);
                values.push_back(it != state.block_expr_state.end() ? it->second : std::nullopt);
            }
            state.possible_tuples_for_expr[E].insert(std::move(values));
        }
    }
    return true;
}

bool eval_stmt(const Stmt *S, Eval_state &state, ASTContext &ASTctx) {

    if (const auto *CS = dyn_cast<CompoundStmt>(S)) {
        for (const auto s : CS->children()) {
            if (!eval_stmt(s, state, ASTctx))
                return false;
        }
        return true;

    } else if (const auto *DS = dyn_cast<DeclStmt>(S)) {
        // TODO: we should ignore all Decls that have a type thats underlying primitive
        // representation is not an integer

        // DeclStmt may have more than one decl:  int x,y=0,z;
        for (const Decl *decl : DS->decls()) {
            if (const auto *VD = dyn_cast<VarDecl>(decl)) {
                if (is_alias_reference(VD)) {
                    // Evaluate the initializer for its side effects, e.g. `int &r = ++x;`. The
                    // reference has no value of its own, accesses to it go to the variable it
                    // refers to. If that is not known, the reference has an unknown value.
                    if (VD->hasInit()) {
                        Eval_result res{};
                        if (!eval_expr(res, VD->getInit(), state, ASTctx))
                            return false;
                    }
                    if (alias_reference_target(VD, ASTctx) == nullptr)
                        state.all_possible_vals[VD].insert(std::nullopt);
                    return true;
                }
                if (VD->hasInit()) {
                    const Expr *init = VD->getInit();
                    Eval_result res{};
                    if (!eval_expr(res, init, state, ASTctx))
                        return false;
                    state.vars[VD] = to_int_type(res.val, VD->getType(), ASTctx);
                    // var changes
                    state.all_possible_vals[VD].insert(state.vars[VD]);
                    return true;
                } else {
                    // if no init it is std::nullopt
                    state.vars[VD] = std::nullopt;
                    // var changes
                    state.all_possible_vals[VD].insert(state.vars[VD]);
                    return true;
                }
            }
        }
    } else if (const auto *RE = dyn_cast<ReturnStmt>(S)) {
        const auto *ret_expr = RE->getRetValue();
        if (!ret_expr) {
            return true;
        }
        Eval_result ret_val{};
        if (!eval_expr(ret_val, ret_expr, state, ASTctx))
            return false;
        state.return_val = ret_val.val;

        return true;

    } else if (const auto *E = dyn_cast<Expr>(S)) {
        Eval_result res{};
        if (!eval_expr(res, E, state, ASTctx))
            return false;
        return true;
    }
    INTERNAL_ERROR << " UNHANDELED STMT eval_stmt():\n" << loc_str(S) << "\n";
    S->dump(diags, ASTctx);
    return false;
}

// old stuff
// We depend on the successor ordering defined by CFGBlock, see comment above class CFGBlock in
// CFG.h in llvm sources
#define MAX_BLOCK_VISITS 100 // for the old eval_CFG_block() stuff
bool eval_CFG_block(const CFGBlock &block, Eval_state &state, ASTContext &ASTctx) {
    // reset block specific state
    state.block_expr_state.clear();

    if (state.block_visited.find(&block) == state.block_visited.end())
        state.block_visited[&block] = 0;
    if (state.block_visited[&block] > MAX_BLOCK_VISITS)
        return true;
    state.block_visited[&block]++;

    // the exit_block is just an empty block signaling that we are done
    if (&block == state.exit_block)
        return true;

    const auto *trm_stmt = block.getTerminatorStmt();
    Eval_result cond_res{};
    const CFGElement &last_element = *block.rbegin();
    // Eval all the elements in the block
    for (const CFGElement &element : block) {
        if (element.getKind() == CFGElement::Statement) {
            const auto *stmt = element.castAs<CFGStmt>().getStmt();
            // If there is a terminator statement the last statement is in most cases the condition
            // expression that is used to choose which successor to go to
            if (trm_stmt && &element == &last_element) {
                if (const auto *E = dyn_cast<Expr>(stmt)) {
                    if (!eval_expr(cond_res, E, state, ASTctx))
                        return false;
                } else {
                    INTERNAL_ERROR << " Expecting the last element of an block with terminator "
                                      "stmt to be an expression "
                                   << loc_str(stmt) << "\n";
                    return false;
                }
            } else {
                if (!eval_stmt(stmt, state, ASTctx))
                    return false;
            }
        } else {
            INTERNAL_ERROR << " TODO: Unhandled CFGElement Kind: " << element.getKind() << "\n";
            return false;
        }
    }


    // Now based on the terminatorStmt result pick which block we go to
    auto succ = block.succ_begin();
    if (trm_stmt == nullptr // if no terminator statement or terminator is one of these:
        || dyn_cast<BreakStmt>(trm_stmt) // break, continue, goto we should have only one successor
        || dyn_cast<ContinueStmt>(trm_stmt) //
        || dyn_cast<GotoStmt>(trm_stmt)     //
    ) {
        if (block.succ_empty()) {
            // we have reached an early return
            return true;
        }
        if (block.succ_size() != 1) {
            INTERNAL_ERROR << " Expected block[" << block.getBlockID()
                           << "] to have 1 successors but it has succ_size " << block.succ_size()
                           << " terminator StmtClass: `" << trm_stmt->getStmtClassName() << "` \n";
            if (block.Elements.size() > 0)
                block.Elements.begin()->dump();
        }
        CFGBlock *s0 = succ[0]; // this is a cast, succ[0] is actually an AdjacentBlock
        return eval_CFG_block(*s0, state, ASTctx);
    } else if (dyn_cast<IfStmt>(trm_stmt)                 //
               || dyn_cast<ForStmt>(trm_stmt)             //
               || dyn_cast<WhileStmt>(trm_stmt)           //
               || dyn_cast<DoStmt>(trm_stmt)              //
               || dyn_cast<BinaryOperator>(trm_stmt)      //
               || dyn_cast<ConditionalOperator>(trm_stmt) //
    ) {
        if (block.succ_size() != 2) {
            INTERNAL_ERROR << "ERROR: Expected block[" << block.getBlockID()
                           << "] to have 2 successors but it has succ_size " << block.succ_size()
                           << " terminator StmtClass: `" << trm_stmt->getStmtClassName() << "` \n";
            return false;
        }

        // ASSUMPTION: This assumes this order is quaranteed by LLVM
        CFGBlock *true_branch = succ[0];
        CFGBlock *false_branch = succ[1];

        if (!cond_res.val.has_value()) {
            INTERNAL_ERROR << "Could not determine the value of the condition. "
                              "Cant follow the control flow correctly. Giving up. "
                           << loc_str(trm_stmt) << "\n";
            // Report relevant error for the user
            report_diag(ASTctx, DiagnosticsEngine::Level::Error,
                        trm_stmt->getSourceRange().getBegin(),
                        "Could not determine the value of the condition. Cant follow the control "
                        "flow correctly. Giving up.");
            return false;
        }

        if (const auto *BO = dyn_cast<BinaryOperator>(trm_stmt)) {
            if (!(BO->getOpcode() == BO_LAnd || BO->getOpcode() == BO_LOr)) {
                // Unreachable
                INTERNAL_ERROR << "Expected terminator StmtClass: `BinaryOperator` to have Opcode"
                                  "BO_LAnd or BO_LOr but found "
                               << BO->getOpcodeStr() << loc_str(trm_stmt) << "\n";
                return false;
            }
        }

        if (cond_res.val.value()) {
            return eval_CFG_block(*true_branch, state, ASTctx);
        } else {
            return eval_CFG_block(*false_branch, state, ASTctx);
        }

    } else if (dyn_cast<SwitchStmt>(trm_stmt)) {
        // switch statement needs special handling
        // has as many successor as there are cases (+1 if there is no default case present)

        CFGBlock *default_block = block.succ_end()[-1];
        for (auto it = block.succ_begin(); it != block.succ_end(); it++) {
            const CFGBlock *case_block = *it;
            auto case_label_stmt = case_block->getLabel();
            if (case_label_stmt == nullptr) {
                // this is hit if there is no `default:` in the source code
            } else if (auto cs = dyn_cast<CaseStmt>(case_label_stmt)) {

                auto case_val_expr = cs->getLHS();
                Expr::EvalResult result, result_rhs;
                if (!case_val_expr->EvaluateAsConstantExpr(result, ASTctx) ||
                    (cs->getRHS() && !cs->getRHS()->EvaluateAsConstantExpr(result_rhs, ASTctx))) {
                    // probably not valid c++ and should not happen
                    INTERNAL_ERROR << " the case label is not an constant expression"
                                   << loc_str(case_val_expr) << "\n";
                    return false;
                }
                // GNU case ranges `case 1 ... 3:` have a RHS
                const int64_t case_low = result.Val.getInt().getExtValue();
                const int64_t case_high =
                    cs->getRHS() ? result_rhs.Val.getInt().getExtValue() : case_low;
                if (!cond_res.val.has_value()) {
                    INTERNAL_ERROR << " Could not determine the value of the switch variable. "
                                      "Cant follow the control flow correctly. Giving up. "
                                   << loc_str(trm_stmt) << "\n";
                    // Report relevant error for the user
                    report_diag(ASTctx, DiagnosticsEngine::Level::Error,
                                trm_stmt->getSourceRange().getBegin(),
                                " Could not determine the value of the switch variable. "
                                "Cant follow the control flow correctly. Giving up. ");
                    return false;
                }
                if (case_low <= cond_res.val.value() && cond_res.val.value() <= case_high) {
                    return eval_CFG_block(*case_block, state, ASTctx);
                }

            } else if (auto _ = dyn_cast<DefaultStmt>(case_label_stmt)) {
                // do nothing
            } else {
                INTERNAL_ERROR
                    << " case is not CaseStmt or DefaultStmt or null, should be unreachable\n ";
                return false;
            }
        }
        //  Nothing has matched: eval default
        return eval_CFG_block(*default_block, state, ASTctx);

    } else {
        INTERNAL_ERROR << " TODO: Unhandeled terminatorStmt class " << loc_str(trm_stmt) << " : "
                       << trm_stmt->getStmtClassName() << "\n";
        return false;
    }

    // Should be unreachable
    INTERNAL_ERROR << " Something unhandeled eval_CFG_block fell through\n";
    return false;
}

bool Idx_expr_evaluator::init(const FunctionDecl *func, Stmt *statement_in_the_func,
                              ASTContext *ctx) {
    assert(func != nullptr);
    assert(statement_in_the_func != nullptr);

    CFG::BuildOptions CFG_opts;
    cfg = CFG::buildCFG(func, statement_in_the_func, ctx, CFG_opts);

    if (!cfg) {
        ASTContext &ASTctx = *ctx;
        INTERNAL_ERROR << " Could not build a CFG for statement starting at "
                       << loc_str(statement_in_the_func) << "\n"
                       << "               : The supplied FunctionDecl is at " << loc_str(func)
                       << "\n";
        return false;
    }
    ASTctx = ctx;
    main_entry_block = &cfg->getEntry();
    main_exit_block = &cfg->getExit();
    // main_state.exit_block = main_exit_block; // only used by the old eval_CFG_block()

    diags.enable_colors(true);
    reported_diag.clear();

    return true;
}

void Idx_expr_evaluator::add_expr_to_tracking_list(const Expr *E) {
    main_state.possible_vals_for_expr[E] = {}; // init as empty set
}

void Idx_expr_evaluator::add_tuple_to_tracking_list(const Expr *main_expr,
                                                    std::vector<const Expr *> child_expr) {
    main_state.tracked_tuples[main_expr] = std::move(child_expr);
    main_state.possible_tuples_for_expr[main_expr] = {}; // init as empty set
}

std::pair<const Tuple_set *, bool>
Idx_expr_evaluator::get_tuples_for_expr(const Expr *main_expr) const {
    const auto it = main_state.possible_tuples_for_expr.find(main_expr);
    if (it == main_state.possible_tuples_for_expr.end())
        return {nullptr, false};

    for (const Value_tuple &tuple : it->second)
        for (const std::optional<int64_t> &value : tuple)
            if (!value.has_value())
                return {&it->second, false};

    return {&it->second, true};
}

std::pair<const Value_set *, bool> Idx_expr_evaluator::get_vals_for_expr(Expr *E) const {

    auto it = main_state.possible_vals_for_expr.find(E);
    if (it == main_state.possible_vals_for_expr.end())
        return {nullptr, false};

    const Value_set *ret;
    ret = &(it->second);

    if (it->second.find(std::nullopt) != it->second.end())
        return {ret, false};

    return {ret, true};
}

void Idx_expr_evaluator::print_vals_for_expr(llvm::raw_ostream &os, Expr *E) const {
    const auto pair = get_vals_for_expr(E);
    const auto *vals = pair.first;
    bool ok = pair.second;
    if (!vals) {
        os << "NO EXPR FOUND: " << "print_vals_for_expr: no expr in possible_vals_for_expr: "
           << loc_strp(E) << "\n";
        return;
    }
    if (!ok) {
        os << "NOTE: Could not track the Expr: ";
    }
    os << loc_strp(E) << " Expr = [ ";
    for (auto val : *vals) {
        os << val << ", ";
    }
    os << " ]\n";
}

void Idx_expr_evaluator::print_all_vals_vars(llvm::raw_ostream &os) const {
    for (auto var_vals : main_state.all_possible_vals) {
        os << loc_strp(var_vals.first) << " `" << var_vals.first->getNameAsString() << "` = [ ";
        for (auto val : var_vals.second) {
            os << val << ", ";
        }
        os << " ]\n";
    }
}

void Idx_expr_evaluator::print_all_vals_exprs(llvm::raw_ostream &os) const {
    for (auto expr_vals : main_state.possible_vals_for_expr) {
        os << loc_strp(expr_vals.first) << " = [ ";
        for (auto val : expr_vals.second) {
            os << val << ", ";
        }
        os << " ]\n";
    }
}

void Idx_expr_evaluator::print_eval_diags(llvm::raw_ostream &os) const {
    os << diag_msg;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Misc debug helpers
#define BRED_STR(str) "\033[41m" str "\033[0m"
#define RED_STR(str) "\033[0;31m" str "\033[0m"
bool debug_print_CFG_block(const CFGBlock &block, Eval_state &state, ASTContext &ASTctx) {
    bool ret = true; // set this to false if we encounter something we don't handle

    if (state.block_visited.find(&block) == state.block_visited.end())
        state.block_visited[&block] = 0;
    if (state.block_visited[&block] > 1 /*visit each block only once*/)
        return true;
    state.block_visited[&block]++;

    const auto *trm_stmt = block.getTerminatorStmt();
    const Stmt *cond = nullptr;
    bool is_binaryOperator_LAnd_Lor = false;
    if (trm_stmt) {
        if (const auto *IF = dyn_cast<IfStmt>(trm_stmt)) {
            cond = IF->getCond();
        } else if (const auto *FOR = dyn_cast<ForStmt>(trm_stmt)) {
            cond = FOR->getCond();
        } else if (const auto *WHILE = dyn_cast<WhileStmt>(trm_stmt)) {
            cond = WHILE->getCond();
        } else if (const auto *DO = dyn_cast<DoStmt>(trm_stmt)) {
            cond = DO->getCond();
        } else if (const auto *SW = dyn_cast<SwitchStmt>(trm_stmt)) {
            cond = SW->getCond();

            // TODO:
            //} else if (const auto *CS = dyn_cast<ConditionalOperator>(trm_stmt)) {
            //    cond = CS->getCond();
        } else if (const auto *BO = dyn_cast<BinaryOperator>(trm_stmt)) {
            is_binaryOperator_LAnd_Lor = true;
            (void)BO; // Unused
        } else if (const auto *BRK = dyn_cast<BreakStmt>(trm_stmt)) {
            // BreakStmt does not have any associated condition, so we dont need anything
            (void)BRK; // Unused
        } else if (const auto *CON = dyn_cast<ContinueStmt>(trm_stmt)) {
            // ContinueStmt does not have any associated condition, so we dont need anything
            (void)CON; // Unused
        } else if (const auto *GOTO = dyn_cast<GotoStmt>(trm_stmt)) {
            // ContinueStmt does not have any associated condition, so we dont need anything
            (void)GOTO; // Unused
        } else {
            llvm::outs() << BRED_STR("TODO:") " Unhandeled terminatorStmt class "
                         << loc_str(trm_stmt) << " : " << trm_stmt->getStmtClassName() << "\n";
            ret = false;
        }
    }

    llvm::outs() << "------------------------------\n";
#define INDENT0 "            : "
    //          "- block[n]  :
    llvm::outs() << "- block[" << block.getBlockID() << "] " << " : \n";

    const auto &last_element = *block.rbegin();
    int eidx = 0;
    for (const CFGElement &element : block) {

        if (&last_element == &element)
            llvm::outs() << "LAST ELEMENT\n";
        llvm::outs() << INDENT0 "element[" << eidx << "]:  ";
        if (element.getKind() == CFGElement::Statement) {
            const auto *stmt = element.castAs<CFGStmt>().getStmt();
            // stmt->dumpColor();
            llvm::outs() << stmt->getStmtClassName() << " `";
            element.dump();
#define INDENT1 "            "
            llvm::outs() << INDENT1 INDENT0 << loc_str(stmt);

            if (cond && cond == stmt) {
                llvm::outs() << " cond " << cond << " stmt " << stmt << "\n";
                llvm::outs() << RED_STR(" [Terminator Expr] : ") << stmt->getStmtClassName();
            } else if (cond && is_binaryOperator_LAnd_Lor && &last_element == &element) {
                llvm::outs() << " cond " << cond << " stmt " << stmt << "\n";
                llvm::outs() << RED_STR(" [Terminator Expr] : is_binaryOperator_LAnd_Lor=true ")
                             << stmt->getStmtClassName();
            }
            llvm::outs() << "\n";

            if (const auto *E = dyn_cast<Expr>(stmt)) {
                if (const auto *_ = dyn_cast<CallExpr>(E)) {

                    llvm::outs() << INDENT1 INDENT0 << "Expr ";
                    auto paren = ASTctx.getParents(*E);
                    if (!paren.empty()) {
                        llvm::outs() << "paren len " << paren.size() << "  paren[0] "
                                     << (paren[0].getSourceRange().getBegin().printToString(
                                            ASTctx.getSourceManager()));
                        const auto *pe = paren[0].get<Expr>();
                        llvm::outs() << " Expr: " << E << " paren: " << pe;
                        if (pe != nullptr) {
                            llvm::outs() << RED_STR(" Extra control flow call ");
                        }
                    }
                    llvm::outs() << "\n";
                }
            }

        } else {
            element.dump();
            llvm::outs() << INDENT1 INDENT0 " Unhandled CFGElement Kind: " << element.getKind()
                         << "\n";
            ret = false;
        }
        eidx++;
    }

    // auto term = block.getTerminator();
    if (trm_stmt) {
        llvm::outs() << INDENT0 RED_STR(" [Terminator Stmt] :  ") << loc_str(trm_stmt) << " : "
                     << trm_stmt->getStmtClassName() << " " << trm_stmt << "\n";
        // term.getStmt()->dumpColor();
    }

    uint64_t i = 0;
    for (auto it = block.succ_begin(); it < block.succ_end(); it++, i++) {
        llvm::outs() << INDENT0 " -> succ[" << i << "] = block[" << (*it)->getBlockID() << "]\n";
    }

    for (const CFGBlock *succ : block.succs()) {
        bool r = debug_print_CFG_block(*succ, state, ASTctx);
        if (!r)
            ret = false;
    }
    return ret;
}


//////////////////////////////////////////////////////////////////////////////////////////
// Worklist fixpoint CFG traversal
//
// The worklist holds (block, state) pairs. Popping a pair that
// has not been seen at that block before evaluates the block from that state and pushes the
// resulting state to each successor the terminator can take: a known condition or switch value
// selects one successor, an unknown one all of them. A pair already seen at its block is
// skipped, so each (block, state) pair is evaluated exactly once.
//
// This reaches the least fixpoint of the dataflow equations
//     states(entry) = { initial state }
//     states(B) = union over predecessors P of transfer(P -> B, states(P))
// over the lattice of sets of states (join = set union, states are never merged). It is
// Kildall's algorithm, but as each state is evaluated independently (the transfer distributes
// over the union), only newly reached states are evaluated instead of re-evaluating a block
// for all its states on every visit. The result does not depend on the worklist order. Loops
// end when their iterations stop producing new states. At the end states(B) holds every state
// reaching B, and possible_vals_for_expr / all_possible_vals every value observed on the way.
//
// We do no widening, since merging states would make the results not exact.
// To terminate eg unknown loop bounds there are fixed limits on how many times CFG blocks are
// visited. Exceeding a limit stops the evaluation with a failure
//  - Eval_state::max_block_states (default DEFAULT_MAX_BLOCK_STATES): maximum number of
//    distinct states entering one CFG block.
//  - MAX_CALL_DEPTH: maximum depth of nested constexpr function calls

// The state of one execution path at a CFG block boundary
struct Path_eval_state {
    std::unordered_map<const VarDecl *, std::optional<int64_t>> vars;
    // The value of the last evaluated expression. Needed when a block starts with an `&&`, `||` or
    // `?:` whose value was decided in the predecessor block, see eval_expr().
    std::optional<int64_t> last_res;
    std::optional<int64_t> return_val;

    bool operator==(const Path_eval_state &other) const {
        return vars == other.vars && last_res == other.last_res && return_val == other.return_val;
    }
};

// hashing Path_eval_state so we can use it in unordered_set
// and quickly check if a state has been encountered
struct Path_eval_state_hash {
    size_t operator()(const Path_eval_state &p) const {
        // vars is unordered, so combine its entries with an order independent sum.
        size_t vars_hash = 0;
        for (const auto &[decl, value] : p.vars)
            vars_hash += hash_combine(decl, value.has_value(), value.value_or(0));
        return hash_combine(vars_hash, p.last_res.has_value(), p.last_res.value_or(0),
                            p.return_val.has_value(), p.return_val.value_or(0));
    }
};

using Path_set = std::unordered_set<Path_eval_state, Path_eval_state_hash>;

static Path_eval_state merge_paths(const Path_eval_state &a, const Path_eval_state &b) {
    Path_eval_state merged;
    for (const auto &[decl, value] : a.vars) {
        const auto it = b.vars.find(decl);
        if (it != b.vars.end() && it->second == value)
            merged.vars.emplace(decl, value);
        else
            merged.vars.emplace(
                decl, std::nullopt); // variables (and the return value) that differ become unknown
    }
    for (const auto &[decl, value] : b.vars) {
        if (a.vars.find(decl) == a.vars.end())
            merged.vars.emplace(decl, std::nullopt);
    }
    if (a.return_val == b.return_val)
        merged.return_val = a.return_val;
    return merged;
}

static std::string block_loc_str(const CFGBlock &block, ASTContext &ASTctx) {
    if (const Stmt *terminator = block.getTerminatorStmt())
        return loc_str(terminator);
    for (const CFGElement &element : block) {
        if (const auto cfg_stmt = element.getAs<CFGStmt>())
            return loc_str(cfg_stmt->getStmt());
    }
    return "<unknown location>";
}

// Evaluate all execution paths through `cfg`, starting with the variable values in state.vars.
// Observed values are recorded in state.possible_vals_for_expr and state.all_possible_vals.
// On success state.vars and state.return_val are the merge of all paths reaching the exit.
// Returns false (and stops) if something can not be evaluated, or if more than
// state.max_block_states distinct states reach a block.
static bool eval_CFG_worklist(const CFG &cfg, Eval_state &state, ASTContext &ASTctx) {
    std::unordered_map<const CFGBlock *, Path_set> seen;
    Path_set exit_paths;
    // Used as a stack: finishes one path before starting the others, which keeps it short.
    std::vector<std::pair<const CFGBlock *, Path_eval_state>> worklist;

    const auto push = [&](const CFGBlock *succ, Path_eval_state path) {
        if (succ != nullptr) // nullptr is an unreachable successor
            worklist.emplace_back(succ, std::move(path));
    };

    Path_eval_state initial;
    initial.vars = std::move(state.vars);
    push(&cfg.getEntry(), std::move(initial));

    while (!worklist.empty()) {
        const CFGBlock *block = worklist.back().first;
        Path_eval_state path = std::move(worklist.back().second);
        worklist.pop_back();

        Path_set &block_seen = seen[block];
        if (block_seen.find(path) != block_seen.end())
            continue; // this path has already been explored from here on
        if (block_seen.size() >= state.max_block_states) {
            INTERNAL_ERROR << "More than " << state.max_block_states
                           << " different states reach the same CFG block (e.g. a loop with an "
                              "unknown or too large bound), giving up. "
                           << block_loc_str(*block, ASTctx) << "\n";
            return false;
        }
        block_seen.insert(path);

        if (block == &cfg.getExit()) {
            exit_paths.insert(std::move(path));
            continue;
        }

        // Evaluate the block. The block_expr_state cache is only valid within one block
        state.vars = std::move(path.vars);
        state.last_res.val = path.last_res;
        state.return_val = path.return_val;
        state.block_expr_state.clear();
        for (const CFGElement &element : *block) {
            if (const auto cfg_stmt = element.getAs<CFGStmt>()) {
                if (!eval_stmt(cfg_stmt->getStmt(), state, ASTctx))
                    return false;
            } else if (element.getAs<CFGScopeBegin>() || element.getAs<CFGScopeEnd>() ||
                       element.getAs<CFGLifetimeEnds>() || element.getAs<CFGLoopExit>()) {
                // skip blocks (markers) that do not execute any code
            } else {
                INTERNAL_ERROR << " Unhandled CFGElement Kind " << element.getKind()
                               << " in CFG block at " << block_loc_str(*block, ASTctx) << "\n";
                return false;
            }
        }

        Path_eval_state next;
        next.vars = std::move(state.vars);
        next.last_res = state.last_res.val;
        next.return_val = state.return_val;

        // Get the value of the terminator condition. Prioritise the conditions own cached value
        // over last_res, which may come from a later element in the block. The condition is not
        // cached when it is not evaluated in this block, eg the whole `a && b` of `if (a && b)`:
        // the block branching on it evaluates `b`, so last_res is the right value there
        const auto condition_value = [&]() -> std::optional<int64_t> {
            if (const auto *cond = dyn_cast_or_null<Expr>(block->getTerminatorCondition())) {
                const auto it = state.block_expr_state.find(cond);
                if (it != state.block_expr_state.end())
                    return it->second;
            }
            return state.last_res.val;
        };

        // Pick the successors. We depend on the successor ordering defined by CFGBlock, see the
        // comment above class CFGBlock in CFG.h in llvm sources.
        const Stmt *trm_stmt = block->getTerminatorStmt();
        if (trm_stmt == nullptr || isa<BreakStmt, ContinueStmt, GotoStmt>(trm_stmt)) {
            // Only one successor, or none at the end of the function
            for (const CFGBlock *succ : block->succs())
                push(succ, next);

        } else if (isa<IfStmt, ForStmt, WhileStmt, DoStmt, BinaryOperator, ConditionalOperator>(
                       trm_stmt)) {
            if (block->succ_size() != 2) {
                INTERNAL_ERROR << "Expected block[" << block->getBlockID()
                               << "] to have 2 successors but it has succ_size "
                               << block->succ_size() << " terminator StmtClass: `"
                               << trm_stmt->getStmtClassName() << "` \n";
                return false;
            }
            if (const auto *BO = dyn_cast<BinaryOperator>(trm_stmt)) {
                if (!(BO->getOpcode() == BO_LAnd || BO->getOpcode() == BO_LOr)) {
                    // Unreachable
                    INTERNAL_ERROR << "Expected terminator StmtClass: `BinaryOperator` to have "
                                      "Opcode BO_LAnd or BO_LOr but found "
                                   << BO->getOpcodeStr() << loc_str(trm_stmt) << "\n";
                    return false;
                }
            }
            // succ[0] is the true branch, succ[1] the false branch
            // Unknown condition value follows both branches
            const std::optional<int64_t> cond = condition_value();
            for (const bool branch : {true, false}) {
                if (cond.has_value() && static_cast<bool>(cond.value()) != branch)
                    continue;
                Path_eval_state taken = next;
                // Taking this branch means the condition evaluated to `branch`. This is also the
                // value of a short circuited `&&` / `||` used as a value in the successor block
                taken.last_res = branch;
                push(block->succ_begin()[branch ? 0 : 1], std::move(taken));
            }

        } else if (isa<SwitchStmt>(trm_stmt)) {
            // One successor for each case label. The last successor is the `default:` block, or
            // the block after the switch if there is no default label
            const std::optional<int64_t> value = condition_value();
            if (!value.has_value()) {
                // Any case is possible
                for (const CFGBlock *succ : block->succs())
                    push(succ, next);
            } else {
                const CFGBlock *selected = block->succ_rbegin()->getReachableBlock();
                for (const CFGBlock::AdjacentBlock &succ : block->succs()) {
                    // Look at the label also if the block is unreachable, then only
                    // getPossiblyUnreachableBlock() is set
                    const CFGBlock *succ_block = succ.getReachableBlock()
                                                     ? succ.getReachableBlock()
                                                     : succ.getPossiblyUnreachableBlock();
                    const auto *case_stmt =
                        succ_block ? dyn_cast_or_null<CaseStmt>(succ_block->getLabel()) : nullptr;
                    if (!case_stmt)
                        continue;
                    Expr::EvalResult lhs, rhs;
                    if (!case_stmt->getLHS()->EvaluateAsInt(lhs, ASTctx) ||
                        (case_stmt->getRHS() && !case_stmt->getRHS()->EvaluateAsInt(rhs, ASTctx))) {
                        // probably not valid c++ and should not happen
                        INTERNAL_ERROR << " the case label is not an constant expression "
                                       << loc_str(case_stmt) << "\n";
                        return false;
                    }
                    // support GNU extension for fun, GNU case ranges `case 1 ... 3:` have a RHS
                    const int64_t low = lhs.Val.getInt().getExtValue();
                    const int64_t high = case_stmt->getRHS() ? rhs.Val.getInt().getExtValue() : low;
                    if (low <= value.value() && value.value() <= high) {
                        selected = succ.getReachableBlock();
                        break;
                    }
                }
                push(selected, std::move(next));
            }

        } else {
            INTERNAL_ERROR << " TODO: Unhandeled terminatorStmt class " << loc_str(trm_stmt)
                           << " : " << trm_stmt->getStmtClassName() << "\n";
            return false;
        }
    }

    state.block_expr_state.clear();
    state.vars.clear();
    state.return_val = std::nullopt;
    if (!exit_paths.empty()) {
        // Publish the merged state at the exit: the API has one value per variable, so differing
        // exit paths become unknown there. All observed values remain in all_possible_vals.
        auto it = exit_paths.begin();
        Path_eval_state merged = *it;
        for (++it; it != exit_paths.end(); ++it)
            merged = merge_paths(merged, *it);
        state.vars = std::move(merged.vars);
        state.return_val = merged.return_val;
    }
    return true;
}

// Evaluate a constexpr function called during exec()
// state.vars holds the argument values, the result is in state.return_val.
static bool eval_constexpr_function(const FunctionDecl *func, Eval_state &state,
                                    ASTContext &ASTctx) {
    Function_eval_cache *cache = state.func_cache;
    // The cache holds only the return value, not the values of reference parameters
    for (const ParmVarDecl *param : func->parameters()) {
        const QualType T = param->getType();
        if (T->isReferenceType() && !T.getNonReferenceType().isConstQualified())
            cache = nullptr;
    }

    // The callee does not track any expressions and its variables are not reported, so the
    // return value is all we need from an earlier call with the same arguments.
    std::vector<std::optional<int64_t>> args;
    if (cache) {
        for (const ParmVarDecl *param : func->parameters()) {
            const auto it = state.vars.find(param);
            args.push_back(it != state.vars.end() ? it->second : std::nullopt);
        }
        const auto memo = cache->return_vals.find({func, args});
        if (memo != cache->return_vals.end()) {
            state.return_val = memo->second;
            return true;
        }
    }

    std::unique_ptr<CFG> uncached_cfg;
    std::unique_ptr<CFG> &func_cfg = cache ? cache->cfgs[func] : uncached_cfg;
    if (!func_cfg) {
        CFG::BuildOptions CFG_opts;
        func_cfg = CFG::buildCFG(func, func->getBody(), &ASTctx, CFG_opts);
        if (!func_cfg) {
            INTERNAL_ERROR << " Could not build a CFG for constexpr function " << loc_str(func)
                           << "\n";
            return false;
        }
    }

    if (!eval_CFG_worklist(*func_cfg, state, ASTctx))
        return false;

    if (cache)
        cache->return_vals.emplace(std::make_pair(func, std::move(args)), state.return_val);
    return true;
}

bool Idx_expr_evaluator::exec() {
    assert(cfg != nullptr);
    assert(ASTctx != nullptr);

    // Clear prior results while retaining the registered expression keys.
    main_state.vars.clear();
    main_state.all_possible_vals.clear();
    for (auto &[expr, values] : main_state.possible_vals_for_expr)
        values.clear();
    for (auto &[main_expr, tuples] : main_state.possible_tuples_for_expr)
        tuples.clear();

    main_state.func_cache = &func_cache;
    main_state.eval_succeeded = eval_CFG_worklist(*cfg, main_state, *ASTctx);
    return main_state.eval_succeeded;
}
