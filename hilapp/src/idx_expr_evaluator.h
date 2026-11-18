#ifndef IDX_EXPR_EVALUATOR_H_
#define IDX_EXPR_EVALUATOR_H_

/**
 * @file idx_expr_evaluator.h
 * @brief Utility for analysing integer 'index like' expressions during AST analysis. Evaluates
 * integer expressions (as much as it can) and accumulates possible different values variables and
 * expressions obtain during evaluation.
 * @details
 * @note Should be tested if there are major changes to llvm. Unit tests for Idx_expr_evaluator
 * are in ../unit_tests/test_idx_expr_evaluator. Tested with llvm 18, 23.1.4
 *
 */

// Limitations:
//
//  - Only integer, enum and bool values are tracked (as int64_t) 
//  - Values that do not fit in int64_t are unknown
//
// Function calls:
//  - Only constexpr functions with found body are evaluated (nested calls at most MAX_CALL_DEPTH)
//  - NOT IMPLEMENTED: lambda functions: evaluation fails
//
// Pointers and references:
//  - Writes are followed through the variable itself, through local reference variables
//    (`int &r = x;`) and through reference arguments at the call.
//  - A variable whose address or non-const reference is stored, so that it can be changed later
//    in a way the evaluation does not follow, "escapes", e.g. `int *p = &x;`, `Holder h(x);`
//    with a reference member, `std::tie(x, y) = ...` or `pick(x) = 1` with `int &pick(int &)`.
//    These are found before the evaluation (find_escaped_vars()), and the values of escaped
//    variables are always unknown.
//  - ASSUMPTION: a function does not store the address of its reference argument in a global
//    variable, or in an object passed to it by value
//
// Control flow:
//  - Statements that are not handled make the evaluation fail (most things are implemented)
//  - Range based for loops run an unknown number of times, the loop variable is unknown
//
//  Over-approximation:
//  - A condition with an unknown value follows both branches, and the values are not narrowed
//    by the condition (inside `if (i < 3)` with i unknown, i is still unknown). The results are
//    an over-approximation: the set of possible values contains all the values that
//    occur, but can still contain values that never occur at run time... TODO: distinquish between certain and uncertain values?
//
//  - The paths are enumerated with concrete values (no intervals or widening): each distinct
//    state entering a CFG block is evaluated separately. More than max_block_states in
//    a block (a loop with an unknown bound, or very many combinations of values) stops
//    the evaluation and it fails
//
// Relevant variables (only_relevant_vars):
//  - Found from the AST of the statement: the variables in the tracked expressions, in the
//    assignments to relevant variables, and in the conditions of if, switch, ?:, && and || that
//    decide whether a relevant expression is evaluated. Not the conditions of loops, unless they
//    contain a relevant variable: an unrelated loop just runs an unknown number of times.
//    break and continue are relevant if their loop contains a relevant expression.
//
// Tuples (add_tuple_to_tracking_list()):
//  - The child_exprs have to be evaluated before the main_expr in the same CFG block,
//    otherwise an unknown value is recorded for them...


#include <map>
#include <memory>
#include <set>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "clang/Analysis/CFG.h"

using namespace clang;


#define DEFAULT_MAX_BLOCK_STATES 10000 // Default for Eval_state::max_block_states used by exec()
#define MAX_CALL_DEPTH 64 // Maximum depth of nested constexpr function calls we evaluate

using Value_set = std::set<std::optional<int64_t>>;

// values of multiple exprs tracked together
using Value_tuple = std::vector<std::optional<int64_t>>;
using Tuple_set = std::set<Value_tuple>;

struct Eval_result {
    std::optional<int64_t> val = std::nullopt;
};

// Variables that may be changed through a stored reference or pointer, see find_escaped_vars().
// Their values are always unknown.
struct Escaped_vars {
    std::unordered_set<const VarDecl *> vars;
    // a reference to a variable that could not be determined escapes: any variable may change
    bool all = false;

    bool contains(const VarDecl *VD) const {
        return all || vars.find(VD) != vars.end();
    }
};

// exec(): results reused between calls to the same constexpr function. Owned by
// Idx_expr_evaluator, as the cached pointers are only valid for one ASTContext.
struct Function_eval_cache {
    std::unordered_map<const FunctionDecl *, std::unique_ptr<CFG>> cfgs;
    // the escaped variables (parameters and locals) of each function
    std::unordered_map<const FunctionDecl *, Escaped_vars> escaped;
    // (function, argument values) -> return value
    std::map<std::pair<const FunctionDecl *, std::vector<std::optional<int64_t>>>,
             std::optional<int64_t>>
        return_vals;
};

struct Eval_state {

    // The maximum number of distinct states (variable values) allowed to enter a single
    // CFG block. More means e.g. a loop with an unknown or a too large trip count, in which case
    // the evaluation is stopped and does not succeed.
    uint64_t max_block_states = DEFAULT_MAX_BLOCK_STATES;

    // Shared with the Eval_states of called constexpr functions, may be nullptr.
    Function_eval_cache *func_cache = nullptr;

    // Depth of nested constexpr function calls
    uint64_t call_depth = 0;

    // The escaped variables of the evaluated code, their values are unknown. May be nullptr.
    const Escaped_vars *escaped = nullptr;

    // If not nullptr, only these variables are tracked, the values of others are unknown. See
    // Idx_expr_evaluator::only_relevant_vars
    const std::unordered_set<const VarDecl *> *relevant = nullptr;

    ////////////////////////////////////////////
    // The state during evaluation:

    // Map of all variables we are tracking ´VarDecl*` to their current value
    // ´std::optional<int64_t>´ during the evaluation.
    // A value of std::nullopt means we dont know the value.
    std::unordered_map<const VarDecl *, std::optional<int64_t>> vars;
    // When the evaluation ends in a return statement the return value is put here.
    // If the value could not be evaluated, or the function does not return anything
    // return_val is nullopt.
    std::optional<int64_t> return_val = std::nullopt;

    // Store the last evaluated expr used for short circuiting in (x?a:b) , || , && expressions.
    Eval_result last_res{};

    std::unordered_map<const Expr *, std::optional<int64_t>> block_expr_state;

    // Only used by the old eval_CFG_block(), which is not used anymore:
    // Track how many times we have visited each block in the control flow graph `CFG`.
    // This needs to be done so we can limit the amount of times we allow the evaluation to evaluate
    // e.g. loops.
    std::unordered_map<const CFGBlock *, uint64_t> block_visited;
    //
    CFGBlock *exit_block;

    ////////////////////////////////////////////
    // The results that we care about in the end

    // maps the Expr* to all its possible values it obtains during evaluation
    std::unordered_map<const Expr *, Value_set> possible_vals_for_expr;

    // Exprs whose values are tracked together: main_expr -> child_exprs, and main_expr -> the
    // tuples of the values of the child_exprs recorded when the anchor was evaluated
    std::unordered_map<const Expr *, std::vector<const Expr *>> tracked_tuples;
    std::unordered_map<const Expr *, Tuple_set> possible_tuples_for_expr;
    // set to true if `possible_vals_for_expr` has reliable results
    bool eval_succeeded = false;

    // list all possible values of all variables we encountered during evaluation
    std::unordered_map<const VarDecl *, Value_set> all_possible_vals;
};


struct Idx_expr_evaluator {

    bool init(const FunctionDecl *func, Stmt *statement_in_the_func, ASTContext *ASTctx);

    void add_expr_to_tracking_list(const Expr *E);

    /**
     * @brief Track the values of the expressions `child_exprs` together: each time `main_expr` is
     * evaluated, the values the child_exprs have in that state are saved as one tuple.
     * @details  this keeps only the combinations of values that occur. E.g. for (field index,
     * direction) in U[mu][X - mu] inside foralldir(mu) the tuples are (e_x, e_x), (e_y, e_y) ...
     * and not all pairs. The child_exprs have to be evaluated before the main_expr in the same CFG
     * block, typically they are subexpressions of it. Otherwise an unknown value is recorded for
     * them.
     */
    void add_tuple_to_tracking_list(const Expr *main_expr, std::vector<const Expr *> child_exprs);

    // Track only the variables that the values of the tracked expressions, or whether they
    // are evaluated, depend on (find_relevant_vars()). Others are unknown, so that e.g.
    //     for (int i = 0; i < n; i++) t += a[i];
    // with an unknown n does not stop the evaluation. As unknown values are followed both ways,
    // the results are still an over-approximation. Set before exec().
    bool only_relevant_vars = false;

    // Is VD tracked in exec(), see only_relevant_vars
    bool is_relevant(const VarDecl *VD) const;

    // Values of variables defined outside the statement given to init() when it starts, e.g.
    // a Direction defined outside a site loop. Variables not here have unknown values.
    std::unordered_map<const VarDecl *, std::optional<int64_t>> initial_vars;

    // Evaluate all execution paths through the statement given to init() with a worklist: a
    // condition with an unknown value follows both branches. Sets main_state.eval_succeeded.
    // Can be called again, e.g. with different initial_vars.
    // Returns false if we encountered something we dont handle or more than
    // main_state.max_block_states distinct states reach a CFG block, the results are not reliable
    bool exec();

    /**
     * @param E the Expr* we want to get the possible values for.
     * @return pair.first = const Value_set *, pointer to the set of possible values, nullptr if the
     * Expr *E is not in the tracked expressions. pair.second = bool, false if E not found or if the
     * set of possible values contains unknown values (std::nullopt), meaning that we could not
     * track the value.
     * note: Empty set is a valid return.
     */
    std::pair<const Value_set *, bool> get_vals_for_expr(Expr *E) const;

    /**
     * @param main_expr given to add_tuple_to_tracking_list()
     * @return pair.first = the set of recorded tuples of the values of the child_exprs, nullptr if
     * main_expr is not tracked. pair.second = false if main_expr is not tracked or if a tuple
     * contains an unknown value (std::nullopt). note: Empty set is a valid return, the main_expr
     * was never evaluated.
     */
    std::pair<const Tuple_set *, bool> get_tuples_for_expr(const Expr *main_expr) const;


    void print_vals_for_expr(llvm::raw_ostream &os, Expr *E) const;
    void print_all_vals_vars(llvm::raw_ostream &os) const;
    void print_all_vals_exprs(llvm::raw_ostream &os) const;
    void print_eval_diags(llvm::raw_ostream &os) const;
    std::string &get_eval_diags_string();

    Eval_state main_state{};

    //  private:
    std::unique_ptr<CFG> cfg;
    ASTContext *ASTctx;
    Function_eval_cache func_cache;
    // the escaped variables of the function given to init()
    Escaped_vars escaped;
    // the statement given to init()
    const Stmt *statement = nullptr;
    // the variables tracked with only_relevant_vars, all if relevant_all
    std::unordered_set<const VarDecl *> relevant;
    bool relevant_all = false;
    // relevant is up to date with the tracked expressions
    bool relevant_valid = false;
    CFGBlock *main_entry_block = nullptr;
    CFGBlock *main_exit_block = nullptr;
};

#endif // IDX_EXPR_EVALUATOR_H_
