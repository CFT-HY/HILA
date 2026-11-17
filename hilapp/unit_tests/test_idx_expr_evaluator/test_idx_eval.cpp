
// Some unit tests for Idx_expr_evaluator

#include "clang/ASTMatchers/ASTMatchFinder.h"
#include "clang/ASTMatchers/ASTMatchers.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Tooling/CompilationDatabase.h"

using namespace clang;
using namespace clang::tooling;
using namespace clang::ast_matchers;

#include "tracked_var.h"

#include <cstdlib>
#include <set>
#include <functional>

#include "../../src/idx_expr_evaluator.h"
#define DIAGS_TO_STR
#include "../../src/idx_expr_evaluator.cpp" // yes this is intentional, diag_msg

/////////////////////////////////////////////////////////////////
// Setup stuff for tests

int g_verbosity = 0;
bool g_make_dot_graph = false;
bool g_print_all_possible_vals = false;
// If non zero, overrides Eval_state::max_block_states
uint64_t g_max_block_states = 0;

std::string g_test_output;

inline llvm::raw_ostream &test_out() {
    static llvm::raw_string_ostream out(g_test_output);
    // the output ends up in llvm::outs()
    static const bool colors_set = (out.enable_colors(llvm::outs().has_colors()), true);
    (void)colors_set;
    return out;
}

// Print the output of a test case collected in test_out()
inline void print_test_output(bool test_pass, unsigned indent = 8) {
    test_out().flush();
    if (g_verbosity >= 2 || (g_verbosity == 1 && !test_pass)) {
        StringRef rest(g_test_output);
        while (!rest.empty()) {
            auto [line, tail] = rest.split('\n');
            llvm::outs().indent(indent) << line << "\n";
            rest = tail;
        }
    }
    g_test_output.clear();
}

// Prints the compiler diagnostics for the test input to test_out(), one line each
class Test_diag_consumer : public DiagnosticConsumer {
  public:
    void HandleDiagnostic(DiagnosticsEngine::Level level, const Diagnostic &info) override {
        DiagnosticConsumer::HandleDiagnostic(level, info);
        SmallString<128> msg;
        info.FormatDiagnostic(msg);
        llvm::raw_ostream &os = test_out();
        if (info.hasSourceManager() && info.getLocation().isValid())
            os << info.getLocation().printToString(info.getSourceManager()) << ": ";
        switch (level) {
        case DiagnosticsEngine::Note:
            os << "note: ";
            break;
        case DiagnosticsEngine::Remark:
            os << "remark: ";
            break;
        case DiagnosticsEngine::Warning:
            os << "warning: ";
            break;
        case DiagnosticsEngine::Error:
        case DiagnosticsEngine::Fatal:
            os << "error: ";
            break;
        default:
            break;
        }
        os << msg << "\n";
    }
};

class TestFunctionAnalyzer : public MatchFinder::MatchCallback {
  public:
    Idx_expr_evaluator evaluator{};

    std::unordered_map<uintptr_t, std::string> names;

    bool exec_ret;
    // set if the test input itself does not compile
    bool input_has_errors = false;

    // The values of the arguments of track_pair(a, b) calls in the test function, tracked
    // together with Idx_expr_evaluator::add_tuple_to_tracking_list(): the union of the tuples of
    // all calls. pairs_reliable is false if an unknown value was recorded.
    Tuple_set pairs;
    bool pairs_reliable = true;

  public:
    void run(const MatchFinder::MatchResult &Result) override {
        diags.enable_colors(true);
        diag_msg.clear();
        ASTContext &ASTctx = *Result.Context;

        // Parsing is done at this point, so errors so far come from the test input. Errors
        // reported later may come from the evaluator (report_diag).
        input_has_errors = ASTctx.getDiagnostics().hasErrorOccurred();

        // the function we analyse for testing:
        const FunctionDecl *func = Result.Nodes.getNodeAs<FunctionDecl>("func");
        if (!func || !func->hasBody())
            return;
        auto *FB = func->getBody();

        evaluator = {};
        evaluator.init(func, FB, &ASTctx);

        if (g_make_dot_graph)
            evaluator.cfg->viewCFG(LangOptions());

        if (g_max_block_states != 0)
            evaluator.main_state.max_block_states = g_max_block_states;

        // track the arguments of track_pair(a, b) calls together
        std::vector<const CallExpr *> track_pair_calls;
        std::function<void(const Stmt *)> find_track_pair_calls = [&](const Stmt *S) {
            if (S == nullptr)
                return;
            if (const auto *CE = dyn_cast<CallExpr>(S)) {
                if (const FunctionDecl *callee = CE->getDirectCallee();
                    callee && callee->getNameAsString() == "track_pair" && CE->getNumArgs() == 2)
                    track_pair_calls.push_back(CE);
            }
            for (const Stmt *child : S->children())
                find_track_pair_calls(child);
        };
        find_track_pair_calls(FB);
        for (const CallExpr *CE : track_pair_calls)
            evaluator.add_tuple_to_tracking_list(CE, {CE->getArg(0), CE->getArg(1)});

        exec_ret = evaluator.exec();

        for (const CallExpr *CE : track_pair_calls) {
            const auto [tuples, reliable] = evaluator.get_tuples_for_expr(CE);
            if (tuples)
                pairs.insert(tuples->begin(), tuples->end());
            pairs_reliable = pairs_reliable && reliable;
        }

        // After ToolInvocation the Compiler context may get freed so store the names now
        // Use the pointers as keys
        for (auto &all_vals : evaluator.main_state.all_possible_vals) {
            names[(uintptr_t)all_vals.first] = all_vals.first->getNameAsString();
        }

        if (g_print_all_possible_vals) {
            test_out() << "------- state.all_possible_vals -------\n";
            for (auto var_vals : evaluator.main_state.all_possible_vals) {
                test_out() << loc_str(var_vals.first) << " `" << var_vals.first->getNameAsString()
                           << "` = [ ";
                for (auto val : var_vals.second) {
                    test_out() << val << ", ";
                }
                test_out() << " ]\n";
            }
        }

        if (!exec_ret) {
            test_out() << diag_msg;
        }
    }
};

#define TEST_EVAL_COMMON_                                                                          \
    std::string fname = "input.cpp";                                                               \
    /* -fno-caret-diagnostics: no "N warnings generated." to stderr, Test_diag_consumer prints */  \
    std::vector<std::string> args = {"clang-tool", "-std=c++17", "-fno-caret-diagnostics", fname}; \
    FixedCompilationDatabase Compilations(".", args);                                              \
    IntrusiveRefCntPtr<llvm::vfs::InMemoryFileSystem> InMemoryFS(                                  \
        new llvm::vfs::InMemoryFileSystem());                                                      \
    InMemoryFS->addFile(fname, 0, llvm::MemoryBuffer::getMemBuffer(code));                         \
    FileSystemOptions FileSysOpts;                                                                 \
    FileManager *Files = new FileManager(FileSysOpts, InMemoryFS);                                 \
    TestFunctionAnalyzer test_analyzer;                                                            \
    MatchFinder Finder;                                                                            \
    Finder.addMatcher(functionDecl(hasName("test")).bind("func"), &test_analyzer);                 \
    auto factory = newFrontendActionFactory(&Finder);                                              \
    std::unique_ptr<FrontendAction> action(factory->create());                                     \
    ToolInvocation Invocation(args, std::move(action),                                             \
                              Files); /*this seems to free Files, .. eventhoug a comment says it   \
                                         does not take ownership..*/                               \
    Test_diag_consumer test_diag_consumer;                                                         \
    Invocation.setDiagnosticConsumer(&test_diag_consumer);                                         \
    const auto &test_state = test_analyzer.evaluator.main_state;                                   \
    const auto &test_names = test_analyzer.names;                                                  \
    (void)test_state; /* not all tests use these */                                               \
    (void)test_names;                                                                              \
    Invocation.run();                                                                              \
    if (test_analyzer.input_has_errors) {                                                          \
        test_pass = false;                                                                         \
        color_msg(test_out(), RED, "Test input does not compile", true, __FILE__, __LINE__)        \
            << "\n";                                                                               \
    }

// Eval the compound statement given in string
#define TEST_EVAL_STRING_CSTMT(str)                                                                \
    std::string code = "void test() {" str "}";                                                    \
    TEST_EVAL_COMMON_

// Eval the compound statement given in string
#define TEST_EVAL_STRING_TOPLVL(str)                                                               \
    std::string code = str;                                                                        \
    TEST_EVAL_COMMON_

const std::set<std::optional<int64_t>> *
get_possible_values_for_var(std::string var_name, const Eval_state &state,
                            const std::unordered_map<uintptr_t, std::string> &names) {
    for (auto &all_vals : state.all_possible_vals) {
        auto name = names.find((uintptr_t)all_vals.first);
        if (name != names.end()) {
            if (var_name == name->second) {
                return &all_vals.second;
            }
        }
    }
    return nullptr;
}

template <typename T>
inline void print_set(raw_ostream &os, const std::set<T> &set) {
    os << "[ ";
    for (auto val : set) {
        // print char types as numbers, not as (possibly unprintable) characters
        if constexpr (std::is_integral_v<T> && sizeof(T) == 1)
            os << static_cast<int>(val) << ", ";
        else
            os << val << ", ";
    }
    os << "]\n";
}

template <typename T>
inline bool compare_var_expect(std::string var_name, const std::set<T> &expected_values,
                               const Eval_state &state,
                               const std::unordered_map<uintptr_t, std::string> &names) {
    const auto *var_vals = get_possible_values_for_var(var_name, state, names);
    if (!var_vals) {
        test_out() << "Did not find var in test_state with name `" << var_name << "`\n";
        return false;
    }
    // quick hack, convert std::optionals to int
    std::set<T> vars;
    for (auto v : *var_vals) {
        if (!v.has_value()) {
            test_out() << "Found a nullopt in possible values for var `" << var_name
                       << "` in test_state \n";
            return false;
        }
        vars.insert(v.value());
    }

    if (!(vars == expected_values)) {
        return false;
    }
    return true;
}

#define print_var_result(VAR)                                                                      \
    do {                                                                                           \
        auto *a_vals = get_possible_values_for_var(VAR, test_state, test_names);                   \
        if (!a_vals) {                                                                             \
            test_out() << "Var result: Did not find var `" VAR "`\n";                              \
        } else {                                                                                   \
            test_out() << "Var result: `" VAR "` = ";                                              \
            print_set(test_out(), *a_vals);                                                        \
        }                                                                                          \
    } while (0)
#define print_var_expect(VAR, VALUES)                                                              \
    do {                                                                                           \
        test_out() << "Var expect: `" #VAR "` = ";                                                 \
        print_set(test_out(), VALUES);                                                             \
    } while (0)

// Compare the possible values of the variable `VAR` in the evaluated code against the set VALUES.
// Use when the expected values are collected from several runs of the test code, e.g. for
// different values of an unknown variable.
#define COMPARE_VALUES_AND_SET_TEST_PASS(VAR, VALUES)                                              \
    do {                                                                                           \
        if (!compare_var_expect(#VAR, VALUES, test_state, test_names)) {                           \
            test_pass = false;                                                                     \
            color_msg(test_out(), RED, "Results differ for var ", true, __FILE__, __LINE__)        \
                << "`" #VAR "`\n";                                                                 \
            print_var_result(#VAR);                                                                \
            print_var_expect(VAR, VALUES);                                                         \
        } else if (g_verbosity >= 2) {                                                             \
            print_var_result(#VAR);                                                                \
            print_var_expect(VAR, VALUES);                                                         \
        }                                                                                          \
    } while (0)

#define COMPARE_AND_SET_TEST_PASS(VAR) COMPARE_VALUES_AND_SET_TEST_PASS(VAR, (VAR).values)

// For variables that can not be tracked: check that the possible values of `VAR` in the
// evaluated code contain an unknown value (std::nullopt), so that they are not reported as
// reliable.
#define EXPECT_UNKNOWN_AND_SET_TEST_PASS(VAR)                                                      \
    do {                                                                                           \
        const auto *u_vals = get_possible_values_for_var(#VAR, test_state, test_names);            \
        if (!u_vals || u_vals->count(std::nullopt) == 0) {                                         \
            test_pass = false;                                                                     \
            color_msg(test_out(), RED, "Expected an unknown value for var ", true, __FILE__,       \
                      __LINE__)                                                                    \
                << "`" #VAR "`\n";                                                                 \
            print_var_result(#VAR);                                                                \
        } else if (g_verbosity >= 2) {                                                             \
            print_var_result(#VAR);                                                                \
        }                                                                                          \
    } while (0)

inline void print_tuples(raw_ostream &os, const Tuple_set &tuples) {
    os << "[ ";
    for (const Value_tuple &t : tuples) {
        os << "(";
        for (size_t i = 0; i < t.size(); ++i)
            os << (i ? ", " : "") << t[i];
        os << ") ";
    }
    os << "]\n";
}

// Compare the pairs recorded for the track_pair(a, b) calls in the evaluated code (see
// TestFunctionAnalyzer::pairs) against PAIRS, a std::set<std::vector<int64_t>> recorded by the
// normally run code.
#define COMPARE_PAIRS_AND_SET_TEST_PASS(PAIRS)                                                     \
    do {                                                                                           \
        Tuple_set expected_pairs;                                                                  \
        for (const auto &p : PAIRS)                                                                \
            expected_pairs.insert(Value_tuple(p.begin(), p.end()));                                \
        if (!test_analyzer.pairs_reliable || test_analyzer.pairs != expected_pairs) {              \
            test_pass = false;                                                                     \
            color_msg(test_out(), RED, "Results differ for track_pair", true, __FILE__, __LINE__)  \
                << "\n";                                                                           \
            test_out() << "Pairs result: ";                                                        \
            print_tuples(test_out(), test_analyzer.pairs);                                         \
            test_out() << "Pairs expect: ";                                                        \
            print_tuples(test_out(), expected_pairs);                                              \
        } else if (g_verbosity >= 2) {                                                             \
            test_out() << "Pairs result: ";                                                        \
            print_tuples(test_out(), test_analyzer.pairs);                                         \
        }                                                                                          \
    } while (0)

#define TEST_CASE_BEGIN                                                                            \
    bool test_pass = true;                                                                         \
    g_test_output.clear()

#define TEST_CASE_END                                                                              \
    print_test_output(test_pass);                                                                  \
    if (!test_pass) {                                                                              \
        color_msg(llvm::outs(), RED, " TEST FAILED: ", true, __FILE__, __LINE__)                   \
            << __func__ << "\n";                                                                   \
    } else {                                                                                       \
        color_msg(llvm::outs(), GREEN, " TEST PASSED: ", true, __FILE__, __LINE__)                 \
            << __func__ << "\n";                                                                   \
    }                                                                                              \
    return test_pass

/////////////////////////////////////////////////////////////////
// test cases

bool test_basic_if() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = 0;
    if( i < 10) a = 10;
    if( i < -10 ) a = -10;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = 0;
    if( i < 10) a = 10;
    if( i < -10 ) a = -10;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_if_else() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = 0;
    if( i < 10) {a = 10;}
    else        {a = -10;}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = 0;
    if( i < 10) {a = 10;}
    else        {a = -10;}
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_if_elseifs() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = 0;
    if( i < 10) {a = 10;}
    else if(i < -10) {a = -10;}
    else             {a = -999;}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = 0;
    if( i < 10) {a = 10;}
    else if(i < -10) {a = -10;}
    else             {a = -999;}
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_if_decl() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = 0;
    if( int x = 0 ) {a = 10;}
    else {a = 1234;}
    )")
    /////////////////////////////////////////////
    // clang-format off
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
// your code for which the warning gets suppressed 
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = 0;
    if( int x = 0 ) {a = 10;}
    else {a = 1234;}
    // clang-format on
    /////////////////////////////////////////////
#pragma clang diagnostic pop
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_for_loop() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = -1;
    for(i = 0; i < 10; ++i){
        a += i;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = -1;
    for(i = 0; i < 10; ++i){
        a += i;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_for_loop_break() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = -1;
    for(i = 0; i < 10; ++i){
        a += i;
        if(i > 5) break;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = -1;
    for(i = 0; i < 10; ++i){
        a += i;
        if(i > 5) break;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_for_loop_continue() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 0;
    int i = -1;
    for(i = 0; i < 10; ++i){
        a += i;
        if(i > 5) continue;
        a++;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 0;
    TrackedVar<int> i = -1;
    for(i = 0; i < 10; ++i){
        a += i;
        if(i > 5) continue;
        a++;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_while() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    while(i < 10){
        a += i + 2;
        ++i;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    while(i < 10){
        a += i + 2;
        ++i;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_while_break() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    while(i < 10){
        a += i + 2;
        ++i;
        if(i < 6) break;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    while(i < 10){
        a += i + 2;
        ++i;
        if(i < 6) break;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_while_continue() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    while(i < 10){
        a += i + 2;
        ++i;
        if(i < 3) continue;
        a += 100;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    while(i < 10){
        a += i + 2;
        ++i;
        if(i < 3) continue;
        a += 100;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_do_while() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    do {
        a *= ++i;
    } while(i < 10);
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    do {
        a *= ++i;
    } while(i < 10);
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_do_while_break() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    do {
        a *= ++i;
        if(i > 4) break;
    } while(i < 10);
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    do {
        a *= ++i;
        if(i > 4) break;
    } while(i < 10);
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_do_while_continue() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    do {
        a *= ++i;
        if(i > 4) continue;
        a += 1;
    } while(i < 10);
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    do {
        a *= ++i;
        if(i > 4) continue;
        a += 1;
    } while(i < 10);
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_switch() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = -999;
    int i = 0;
    switch(i){
        case 1: a = 1; break;
        case 2: a = 2; break;
        case 0: a = 0; break;
        default: a = 999; break;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = -999;
    TrackedVar<int> i = 0;
    switch(i){
        case 1: a = 1; break;
        case 2: a = 2; break;
        case 0: a = 0; break;
        default: a = 999; break;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_switch_no_default() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = -999;
    int i = 0;
    switch(i){
        case 1: a = 1; break;
        case 2: a = 2; break;
        case 0: a = 0; break;
    }
    a = 123;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = -999;
    TrackedVar<int> i = 0;
    switch(i){
        case 1: a = 1; break;
        case 2: a = 2; break;
        case 0: a = 0; break;
    }
    a = 123;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_switch_no_default_fallthrough() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = -999;
    int i = 0;
    switch(i){
        case 1: a = 1; break;
        case 2: a = 2; break;
    }
    a = 123;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = -999;
    TrackedVar<int> i = 0;
    switch(i){
        case 1: a = 1; break;
        case 2: a = 2; break;
    }
    a = 123;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_switch_case_range() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // GNU case ranges `case low ... high:`, also falling through into the next case
    TEST_EVAL_STRING_CSTMT(R"(
    int a = -999;
    int i = 0;
    for(i = 0; i < 8; ++i){
        switch(i){
            case 0: a = 0; break;
            case 1 ... 3: a = 100 + i; break;
            case 4 ... 5: a = 200 + i;
            case 6: a += 1000; break;
            default: a = 999; break;
        }
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-case-range"
    TrackedVar<int> a = -999;
    TrackedVar<int> i = 0;
    for(i = 0; i < 8; ++i){
        switch(i){
            case 0: a = 0; break;
            case 1 ... 3: a = 100 + i; break;
            case 4 ... 5: a = 200 + i;
            case 6: a += 1000; break;
            default: a = 999; break;
        }
    }
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_goto() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 1;
    int i = 0;
    if( i == 0 ) goto test_basic_goto_label_1;
    a++;
test_basic_goto_label_1:
    a = 123;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    if( i == 0 ) goto test_basic_goto_label_1;
    a++;
test_basic_goto_label_1:
    a = 123;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_basic_LOr() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 11;
    int i = 0;
    int c = 0;
    a = i || ++i;
    if(a) c++;
    a = i || ++i;
    if(a) c++;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 11;
    TrackedVar<int> i = 0;
    TrackedVar<int> c = 0;
    a = i || ++i;
    if(a) c++;
    a = i || ++i;
    if(a) c++;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(c);

    TEST_CASE_END;
}

bool test_basic_LAnd() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int a = 123;
    int i = 0;
    int c = 0;
    a = (++i && ++i);
    if(a) c++;
    a = ( (i-2) && ++i);
    if(a) c++;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> a = 123;
    TrackedVar<int> i = 0;
    TrackedVar<int> c = 0;
    a = (++i && ++i);
    if(a) c++;
    a = ( (i-2) && ++i);
    if(a) c++;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(c);

    TEST_CASE_END;
}

bool test_basic_LAnd_LOr() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 123;
    int a = 10;
    int i = 0;
    int c = 0;
    a = ( !(x && ++i) || (++i && a) ) && ++i;
    if(a) c++;
    a = (x++ && (x-123)) && ( (x=0) && ++i );
    if(a) c++;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> x = 123;
    TrackedVar<int> a = 10;
    TrackedVar<int> i = 0;
    TrackedVar<int> c = 0;
    a = ( !(x && ++i) || (++i && a) ) && ++i;
    if(a) c++;
    a = (x++ && (x-123)) && ( (x=0) && ++i );
    if(a) c++;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(c);

    TEST_CASE_END;
}


bool test_basic_ternary() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    int y = 0;
    int a = x ? y++ : y--;
    a = (x+1) ? ++x + 1 : --x - 1;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> x = 0;
    TrackedVar<int> y = 0;
    TrackedVar<int> a = x ? y++ : y--;
    a = (x+1) ? ++x + 1 : --x - 1;
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(y);
    COMPARE_AND_SET_TEST_PASS(a);

    TEST_CASE_END;
}

bool test_int_to_bool_conversion() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 5;
    bool t = x;
    int z = 0;
    z = t + 1;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    int x = 5;
    bool t = x;
    TrackedVar<int> z = 0;
    z = t + 1;
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(z);

    TEST_CASE_END;
}

bool test_side_effects_with_unknown_operand() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The side effects of the RHS must happen even if the LHS value is unknown.
    // The comma operator's value is its RHS even if the LHS value is unknown.
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int i = 0;
    int z = 0;
    z = y + (i++);
    z = (i++) + y;
    int w = 0;
    w = (z = y, 5);
    )")
    /////////////////////////////////////////////
    // clang-format off
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
    int y = 0;
    TrackedVar<int> i = 0;
    int z = 0;
    z = y + (i++);
    z = (i++) + y;
    TrackedVar<int> w = 0;
    w = (z = y, 5);
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(w);

    TEST_CASE_END;
}

bool test_side_effects_in_store_to_non_variable() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The side effects on both sides of an assignment to something that is not a plain
    // variable (array element, struct member) must happen.
    TEST_EVAL_STRING_CSTMT(R"(
    int arr[4];
    struct S { int m; } s;
    int i = 0;
    arr[0] = i++;
    arr[i++] = 5;
    arr[1] = 0;
    arr[1] += i++;
    s.m = i++;
    int y = 0;
    y = i;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    int arr[4];
    struct S { int m; } s;
    TrackedVar<int> i = 0;
    arr[0] = i++;
    arr[i++] = 5;
    arr[1] = 0;
    arr[1] += i++;
    s.m = i++;
    TrackedVar<int> y = 0;
    y = i;
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(y);

    TEST_CASE_END;
}

bool test_ref_arg_unknown_function() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A function we can not evaluate may change `x` through the reference parameter, so `x`
    // and everything computed from it after the call is unknown.
    TEST_EVAL_STRING_TOPLVL(R"(
void add10(int &a);
void test(){
    int x = 0;
    add10(x);
    int y = 0;
    y = x;
}
    )")
    /////////////////////////////////////////////
    // Can not be evaluated normally, the body of add10() is not known.
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(x);
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(y);

    TEST_CASE_END;
}

bool test_ref_arg_constexpr_function() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A constexpr function changes `x` through the reference parameter
    TEST_EVAL_STRING_TOPLVL(R"(
constexpr void add10(int &a){
    a += 10;
}
void test(){
    int x = 0;
    add10(x);
    int y = 0;
    y = x;
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// use a lambda so we can keep it contained in the test function, `auto &` so that it also
// takes a TrackedVar ...
// constexpr void add10(int &a){
constexpr auto add10 = [](auto &a){
    a += 10;
};
//void test(){
    TrackedVar<int> x = 0;
    add10(x);
    TrackedVar<int> y = 0;
    y = x;
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(y);

    TEST_CASE_END;
}

// For test_ref_arg_constexpr_increment(): an operator can not be defined inside the test function
namespace ref_arg_increment_test {
enum Dir : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRS };
constexpr Dir &operator++(Dir &d) {
    return d = static_cast<Dir>(static_cast<unsigned>(d) + 1);
}
} // namespace ref_arg_increment_test

bool test_ref_arg_constexpr_increment() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A loop stepped with a constexpr operator++ taking the loop variable by reference, like
    // HILA's operator++(Direction &). Without writing back the parameter, d would stay e_x.
    TEST_EVAL_STRING_TOPLVL(R"(
enum Dir : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRS };
constexpr Dir &operator++(Dir &d){
    return d = static_cast<Dir>(static_cast<unsigned>(d) + 1);
}
void test(){
    int sum = 0;
    for (Dir d = e_x; d < NDIRS; ++d)
        sum += d;
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, with the enum and operator above
    using namespace ref_arg_increment_test;
//void test(){
    TrackedVar<int> sum = 0;
    for (Dir d = e_x; d < NDIRS; ++d)
        sum += d;
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(sum);

    TEST_CASE_END;
}

bool test_ref_arg_temporary() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A temporary bound to a reference parameter, e.g. 5 in consume(5) with consume(int &&),
    // can not change any variable: the variables stay known.
    TEST_EVAL_STRING_TOPLVL(R"(
void consume(int &&a);
void test(){
    int x = 0;
    int y = 2;
    consume(5);
    consume(x + 1);
    y = x;
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// void consume(int &&a);
auto consume = [](int &&a) { (void)a; };
//void test(){
    TrackedVar<int> x = 0;
    TrackedVar<int> y = 2;
    consume(5);
    consume(x + 1);
    y = x;
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(y);

    TEST_CASE_END;
}

bool test_ref_arg_constructor() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A constructor we do not evaluate may change `x` through the reference parameter
    TEST_EVAL_STRING_TOPLVL(R"(
struct Set7 {
    Set7(int &a) { a = 7; }
};
void test(){
    int x = 0;
    Set7 s(x);
    int y = 0;
    y = x;
}
    )")
    /////////////////////////////////////////////
    // The evaluator does not evaluate constructors, so the values can not be compared.
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(x);
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(y);

    TEST_CASE_END;
}

bool test_ref_variable() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Assigning through a reference variable changes `x`
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    int &r = x;
    r = 5;
    int y = 0;
    y = x;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    TrackedVar<int> x = 0;
    TrackedVar<int> &r = x;
    r = 5;
    TrackedVar<int> y = 0;
    y = x;
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(y);

    TEST_CASE_END;
}

bool test_ref_variable_chains() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A reference to a reference, a write to `x` read through `r`, a reference bound to the result
    // of `++x`, and a const reference bound to a temporary, which keeps its value.
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    int &r = x;
    int &r2 = r;
    r2 += 2;
    x = 7;
    int y = 0;
    y = r;
    int &r3 = ++x;
    r3 *= 2;
    const int &t = x + 1;
    x = 0;
    int z = 0;
    z = t;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    TrackedVar<int> x = 0;
    TrackedVar<int> &r = x;
    TrackedVar<int> &r2 = r;
    r2 += 2;
    x = 7;
    TrackedVar<int> y = 0;
    y = r;
    TrackedVar<int> &r3 = ++x;
    r3 *= 2;
    const int &t = x + 1;
    x = 0;
    TrackedVar<int> z = 0;
    z = t;
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(y);
    COMPARE_AND_SET_TEST_PASS(z);

    TEST_CASE_END;
}

bool test_ref_variable_unknown_target() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The variable `r` refers to is not determined from its initializer, so a write through it may
    // change any variable: they all become unknown.
    TEST_EVAL_STRING_CSTMT(R"(
    int c = 1;
    int x = 0;
    int z = 0;
    int &r = c ? x : z;
    r = 5;
    int a = 0;
    a = x;
    )")
    /////////////////////////////////////////////
    // The evaluator does not track this precisely, so the values can not be compared.
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(x);
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(z);
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(a);

    TEST_CASE_END;
}

bool test_other_expression_kinds() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Expressions the evaluator does not model (floating point literals, initializer lists) give
    // unknown values but do not stop the evaluation. Constants like sizeof are evaluated.
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    double d = 0.5;
    float f = 1.5f;
    int arr[3] = {1, 2, 3};
    x = sizeof(int);
    x = 'a';
    bool b = true;
    x += b;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    TrackedVar<int> x = 0;
    double d = 0.5;
    float f = 1.5f;
    int arr[3] = {1, 2, 3};
    x = sizeof(int);
    x = 'a';
    bool b = true;
    x += b;
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);

    TEST_CASE_END;
}

bool test_bool_and_char_literals() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // bool literals in arithmetic, conversions and conditions, and character literals of the
    // different kinds, also with escapes and arithmetic that wraps in char.
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    x = true + true;
    x = false;
    x = !false;
    bool b = 5;
    x = b + 10;
    int n = 0;
    bool flag = true;
    while (flag) {
        n++;
        if (n == 3) flag = false;
    }
    int c = 0;
    c = 'z' - 'a';
    c = '\n';
    c = '\0';
    char ch = 'a';
    ch += 2;
    c = ch - 'a';
    char wrap = 127;
    wrap += 1;
    char sc = '\xff';
    unsigned char uc = '\xff';
    c = L'A';
    c = u'ä';
    c = U'\U0001F600';
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconstant-conversion"
    TrackedVar<int> x = 0;
    x = true + true;
    x = false;
    x = !false;
    TrackedVar<bool> b = 5;
    x = b + 10;
    TrackedVar<int> n = 0;
    TrackedVar<bool> flag = true;
    while (flag) {
        n++;
        if (n == 3) flag = false;
    }
    TrackedVar<int> c = 0;
    c = 'z' - 'a';
    c = '\n';
    c = '\0';
    TrackedVar<char> ch = 'a';
    ch += 2;
    c = ch - 'a';
    TrackedVar<char> wrap = 127;
    wrap += 1;
    TrackedVar<char> sc = '\xff';
    TrackedVar<unsigned char> uc = '\xff';
    c = L'A';
    c = u'ä';
    c = U'\U0001F600';
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(b);
    COMPARE_AND_SET_TEST_PASS(n);
    COMPARE_AND_SET_TEST_PASS(flag);
    COMPARE_AND_SET_TEST_PASS(c);
    COMPARE_AND_SET_TEST_PASS(ch);
    COMPARE_AND_SET_TEST_PASS(wrap);
    COMPARE_AND_SET_TEST_PASS(sc);
    COMPARE_AND_SET_TEST_PASS(uc);

    TEST_CASE_END;
}

bool test_unevaluated_operands() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The operands of sizeof and decltype are not evaluated: i++ does not increment i
    TEST_EVAL_STRING_CSTMT(R"(
    int i = 0;
    int s = 0;
    s = sizeof(i++);
    s = sizeof i++;
    decltype(i++) t = 0;
    i += 1;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunevaluated-expression"
#pragma clang diagnostic ignored "-Wunused-variable"
    TrackedVar<int> i = 0;
    TrackedVar<int> s = 0;
    s = sizeof(i++);
    s = sizeof i++;
    decltype(i++) t = 0;
    i += 1;
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(s);

    TEST_CASE_END;
}

bool test_floating_point_is_unknown() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A floating point value can not be tracked as an integer: with integer arithmetic k would be
    // (3 / 2) * 2 = 2, but it is (1.5) * 2 = 3.
    TEST_EVAL_STRING_CSTMT(R"(
    int i = 3;
    double d = i;
    int k = 0;
    k = (d / 2) * 2;
    )")
    /////////////////////////////////////////////
    // The evaluator does not track floating point values, so k can not be compared.
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    EXPECT_UNKNOWN_AND_SET_TEST_PASS(k);

    TEST_CASE_END;
}

bool test_unsigned_and_narrowing() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Values follow the types: unsigned arithmetic wraps, narrowing conversions truncate, also
    // with explicit casts.
    TEST_EVAL_STRING_CSTMT(R"(
    unsigned u = 0;
    u--;
    int y = 0;
    if (u > 0) y = 1;
    else y = 2;
    int x = 300;
    unsigned char c = 0;
    c = x;
    int n = 0;
    n = static_cast<unsigned char>(x);
    n = (short)(x * 1000);
    n = int(u / 2);
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    TrackedVar<unsigned> u = 0;
    u--;
    TrackedVar<int> y = 0;
    if (u > 0) y = 1;
    else y = 2;
    int x = 300;
    TrackedVar<unsigned char> c = 0;
    c = x;
    TrackedVar<int> n = 0;
    n = static_cast<unsigned char>(x);
    n = (short)(x * 1000);
    n = int(u / 2);
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(u);
    COMPARE_AND_SET_TEST_PASS(y);
    COMPARE_AND_SET_TEST_PASS(c);
    COMPARE_AND_SET_TEST_PASS(n);

    TEST_CASE_END;
}

bool test_constants_declared_outside() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Constants declared outside the evaluated code have their values
    TEST_EVAL_STRING_TOPLVL(R"(
constexpr int N = 3;
const int M = 2;
void test(){
    int x = 0;
    for(int i = 0; i < N; ++i){
        x += M;
    }
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
constexpr int N = 3;
const int M = 2;
//void test(){
    TrackedVar<int> x = 0;
    for(int i = 0; i < N; ++i){
        x += M;
    }
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);

    TEST_CASE_END;
}

bool test_enum_direction_loop() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A loop over an unsigned enum like HILA's foralldir(d)
    TEST_EVAL_STRING_TOPLVL(R"(
enum Dir : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRS };
constexpr Dir next_dir(Dir d){
    return static_cast<Dir>(static_cast<unsigned>(d) + 1);
}
constexpr Dir opp_dir(Dir d){
    return static_cast<Dir>(NDIRS - 1 - static_cast<int>(d));
}
void test(){
    int sum = 0;
    int opp = 0;
    for(Dir d = e_x; d < 3; d = next_dir(d)){
        sum += d;
        opp = opp_dir(d);
    }
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// local enum and lambdas so we can keep them contained in the test function ...
enum Dir : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRS };
// constexpr Dir next_dir(Dir d){
constexpr auto next_dir = [](Dir d){
    return static_cast<Dir>(static_cast<unsigned>(d) + 1);
};
// constexpr Dir opp_dir(Dir d){
constexpr auto opp_dir = [](Dir d){
    return static_cast<Dir>(NDIRS - 1 - static_cast<int>(d));
};
//void test(){
    TrackedVar<int> sum = 0;
    TrackedVar<int> opp = 0;
    for(Dir d = e_x; d < 3; d = next_dir(d)){
        sum += d;
        opp = opp_dir(d);
    }
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(sum);
    COMPARE_AND_SET_TEST_PASS(opp);

    TEST_CASE_END;
}

bool test_pair_tracking_correlated() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The values of i and j tracked together: only the 6 combinations that occur, not all 9
    TEST_EVAL_STRING_TOPLVL(R"(
void track_pair(long a, long b);
void test(){
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (i != j) track_pair(i, j);
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
std::set<std::vector<int64_t>> pairs;
// void track_pair(long a, long b);
auto track_pair = [&](long a, long b) { pairs.insert({a, b}); };
//void test(){
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (i != j) track_pair(i, j);
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_PAIRS_AND_SET_TEST_PASS(pairs);

    TEST_CASE_END;
}

bool test_pair_tracking_direction_like() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Like the field index and direction in U[mu][X - mu] inside foralldir(mu)
    TEST_EVAL_STRING_TOPLVL(R"(
enum Dir : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRS };
constexpr Dir next_dir(Dir d){
    return static_cast<Dir>(static_cast<unsigned>(d) + 1);
}
constexpr Dir opp_dir(Dir d){
    return static_cast<Dir>(NDIRS - 1 - static_cast<int>(d));
}
void track_pair(long a, long b);
void test(){
    for (Dir d = e_x; d < 3; d = next_dir(d))
        track_pair(d, opp_dir(d));
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// local enum and lambdas so we can keep them contained in the test function ...
enum Dir : unsigned { e_x = 0, e_y, e_z, e_z_down, e_y_down, e_x_down, NDIRS };
// constexpr Dir next_dir(Dir d){
constexpr auto next_dir = [](Dir d){
    return static_cast<Dir>(static_cast<unsigned>(d) + 1);
};
// constexpr Dir opp_dir(Dir d){
constexpr auto opp_dir = [](Dir d){
    return static_cast<Dir>(NDIRS - 1 - static_cast<int>(d));
};
std::set<std::vector<int64_t>> pairs;
// void track_pair(long a, long b);
auto track_pair = [&](long a, long b) { pairs.insert({a, b}); };
//void test(){
    for (Dir d = e_x; d < 3; d = next_dir(d))
        track_pair(d, opp_dir(d));
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_PAIRS_AND_SET_TEST_PASS(pairs);

    TEST_CASE_END;
}

bool test_pair_tracking_unknown_branch() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Both branches of the unknown condition, but a and b stay correlated: (0, 0) and (1, 10),
    // not (0, 10) or (1, 0)
    TEST_EVAL_STRING_TOPLVL(R"(
void track_pair(long a, long b);
extern int y;
void test(){
    int a = 0;
    if (y > 0) a = 1;
    int b = a * 10;
    track_pair(a, b);
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that take every branch
std::set<std::vector<int64_t>> pairs;
// void track_pair(long a, long b);
auto track_pair = [&](long a, long b) { pairs.insert({a, b}); };
// extern int y;
for (int y : {0, 1}) {
//void test(){
    int a = 0;
    if (y > 0) a = 1;
    int b = a * 10;
    track_pair(a, b);
//}
}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_PAIRS_AND_SET_TEST_PASS(pairs);

    TEST_CASE_END;
}

bool test_pair_tracking_unknown_value() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A member with an unknown value makes the tuples unreliable
    TEST_EVAL_STRING_TOPLVL(R"(
void track_pair(long a, long b);
extern int y;
void test(){
    int i = 0;
    track_pair(i, y);
}
    )")
    /////////////////////////////////////////////
    // The value of y is not known, so the pairs can not be compared.
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret || test_analyzer.pairs_reliable || test_analyzer.pairs.size() != 1)
        test_pass = false;

    TEST_CASE_END;
}

bool test_lots_of_control_flow_things_1() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int i = 0;
    int sum = 0;
    int i2 = 0;
    for(i = 1; i < 100; ++i){
        if(i % 2 == 0) sum += i;
        if(sum == 3){
            while(i2 < 10 + i){
                sum -= i2;
                if(i2 > 8) break;
                i2++;
            }
        }
        if(i%3 == 0) continue;
        else sum += 1;
        sum += 1;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> i = 0;
    TrackedVar<int> sum = 0;
    TrackedVar<int> i2 = 0;
    for(i = 1; i < 100; ++i){
        if(i % 2 == 0) sum += i;
        if(sum == 3){
            while(i2 < 10 + i){
                sum -= i2;
                if(i2 > 8) break;
                i2++;
            }
        }
        if(i%3 == 0) continue;
        else sum += 1;
        sum += 1;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(sum);
    COMPARE_AND_SET_TEST_PASS(i2);

    TEST_CASE_END;
}

bool test_lots_of_control_flow_things_2() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    int i = 0;
    int sum = 0;
    int i2 = 0;
    for(i = 1; i < 100; ++i){
        if(i % 2 == 0 || i % 3 == 0) sum += i;
        if(sum == 3){
            while(i2 < 10 + i){
                sum -= i2;
                if(i2 > 8 && i > 10) break;
                i2++;
            }
        }
        if(i%3 == 0) continue;
        else sum += 1;
        sum += 1;
        if(sum > 2000 && i%5==0) sum += 10000;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
    TrackedVar<int> i = 0;
    TrackedVar<int> sum = 0;
    TrackedVar<int> i2 = 0;
    for(i = 1; i < 100; ++i){
        if(i % 2 == 0 || i % 3 == 0) sum += i;
        if(sum == 3){
            while(i2 < 10 + i){
                sum -= i2;
                if(i2 > 8 && i > 10) break;
                i2++;
            }
        }
        if(i%3 == 0) continue;
        else sum += 1;
        sum += 1;
        if(sum > 2000 && i%5==0) sum += 10000;
    }
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(sum);
    COMPARE_AND_SET_TEST_PASS(i2);

    TEST_CASE_END;
}

bool test_unknown_branch() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_CSTMT(R"(
    // y = some unknown runtime value
    extern int y;
    int x=0,i=0;
    for(i=0;i<10;++i){
        x++;
        if(y > 0) {1+1+y;}
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Execute the same statements concretely and compare every observed value.
    // y = some unknown runtime value
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
    int y=1;
    TrackedVar<int> x=0,i=0;
    for(i=0;i<10;++i){
        x++;
        if(y > 0) {1+1+y;}
    }
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_unknown_ternary() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Unknown condition: both arms are possible. Known condition with an unknown
    // arm: the result is unknown (used to throw std::bad_optional_access).
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int x = 0;
    x = (y > 0) ? 5 : 7;
    int k = 1;
    int z = 0;
    z = (k > 0) ? y : 1;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that take every branch
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
    std::set<int> x_values;
    for (int y : {0, 1}) {
        TrackedVar<int> x = 0;
        x = (y > 0) ? 5 : 7;
        int k = 1;
        int z = 0;
        z = (k > 0) ? y : 1;
        x_values.insert(x.values.begin(), x.values.end());
    }
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    COMPARE_VALUES_AND_SET_TEST_PASS(x, x_values);

    TEST_CASE_END;
}

bool test_unknown_LOr_LAnd_value() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // `||` / `&&` used as values with an unknown LHS (used to throw
    // std::bad_optional_access).
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int i = 0;
    int b = 0;
    b = (y > 0) || (i > 3);
    int c = 0;
    c = (y > 0) && (i > 3);
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that take every branch
    std::set<int> b_values, c_values;
    for (int y : {0, 1}) {
        int i = 0;
        TrackedVar<int> b = 0;
        b = (y > 0) || (i > 3);
        TrackedVar<int> c = 0;
        c = (y > 0) && (i > 3);
        b_values.insert(b.values.begin(), b.values.end());
        c_values.insert(c.values.begin(), c.values.end());
    }
    // clang-format on
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    COMPARE_VALUES_AND_SET_TEST_PASS(b, b_values);
    COMPARE_VALUES_AND_SET_TEST_PASS(c, c_values);

    TEST_CASE_END;
}

bool test_unknown_int_to_bool_cond() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Unknown int used as a condition (CK_IntegralToBoolean on an unknown
    // value used to throw std::bad_optional_access).
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int x = 0;
    if (y + 1) x = 1;
    else x = 2;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that take every branch
    std::set<int> x_values;
    for (int y : {-1, 0}) {
        TrackedVar<int> x = 0;
        if (y + 1) x = 1;
        else x = 2;
        x_values.insert(x.values.begin(), x.values.end());
    }
    // clang-format on
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    COMPARE_VALUES_AND_SET_TEST_PASS(x, x_values);

    TEST_CASE_END;
}

bool test_constexpr_func() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_TOPLVL(R"(
constexpr int func1(int a){
    return a + 2;
}

void test(){
    int a = 1;
    int i = 0;
    for(i = 0; i < 10; ++i){
        a = func1(a);
    }
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it 
// use a lambda so we can keep it contained in the test_constexpr_func() ...
// constexpr int func1(int a){
constexpr auto func1 = [](int a){
    return a + 2;
};
//void test(){
    TrackedVar<int> a = 1;
    TrackedVar<int> i = 0;
    for(i = 0; i < 10; ++i){
        a = func1(a);
    }
//}
    // clang-format on
    /////////////////////////////////////////////
    // Compare results
    COMPARE_AND_SET_TEST_PASS(a);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_division_by_zero_on_unknown_path() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // `y > 0` and `y <= 0` are both unknown, so the path with d = 0 also takes the division,
    // which never happens in the program. It must not crash the evaluator.
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int d = 1;
    if (y > 0) d = 0;
    int x = 0;
    if (y <= 0) x = 10 / d;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that take every branch.
    // `x` is not compared: the evaluator also follows the path that never happens.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
    std::set<int> d_values;
    for (int y : {0, 1}) {
        TrackedVar<int> d = 1;
        if (y > 0) d = 0;
        int x = 0;
        if (y <= 0) x = 10 / d;
        d_values.insert(d.values.begin(), d.values.end());
    }
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    COMPARE_VALUES_AND_SET_TEST_PASS(d, d_values);

    TEST_CASE_END;
}

bool test_switch_fallthrough() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The case 1 block is reached from the switch (s == 1) and by falling through from case 0.
    TEST_EVAL_STRING_CSTMT(R"(
    int s = 0;
    int x = 0;
    switch(s){
        case 0: x = 10;
        // fall through
        case 1: x += 1; break;
        default: x = 30;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    int s = 0;
    TrackedVar<int> x = 0;
    switch(s){
        case 0: x = 10;
        // fall through
        case 1: x += 1; break;
        default: x = 30;
    }
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);

    TEST_CASE_END;
}

bool test_switch_unknown_and_case_range() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // Unknown switch value: every case is possible, including no match.
    // Known value inside a GNU case range.
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int x = 0;
    switch(y){
        case 0: x = 10; break;
        case 1 ... 3: x = 20; break;
    }
    int s = 2;
    int z = 0;
    switch(s){
        case 0: z = 10; break;
        case 1 ... 3: z = 20; break;
        default: z = 30;
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that take every case
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-case-range"
    std::set<int> x_values, z_values;
    for (int y : {0, 2, 5}) {
        TrackedVar<int> x = 0;
        switch(y){
            case 0: x = 10; break;
            case 1 ... 3: x = 20; break;
        }
        int s = 2;
        TrackedVar<int> z = 0;
        switch(s){
            case 0: z = 10; break;
            case 1 ... 3: z = 20; break;
            default: z = 30;
        }
        x_values.insert(x.values.begin(), x.values.end());
        z_values.insert(z.values.begin(), z.values.end());
    }
#pragma clang diagnostic pop
    // clang-format on
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret)
        test_pass = false;
    COMPARE_VALUES_AND_SET_TEST_PASS(x, x_values);
    COMPARE_VALUES_AND_SET_TEST_PASS(z, z_values);

    TEST_CASE_END;
}

bool test_switch_in_loop() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    TEST_EVAL_STRING_CSTMT(R"(
    int k = 0;
    int d = 0;
    for(d = 0; d < 4; d++){
        switch(d){
            case 0: k = 5; break;
            case 2: k += 100; break;
            case 7: k = -1; break;
            default: k += 1;
        }
    }
    )")
    /////////////////////////////////////////////
    // clang-format off
    TrackedVar<int> k = 0;
    TrackedVar<int> d = 0;
    for(d = 0; d < 4; d++){
        switch(d){
            case 0: k = 5; break;
            case 2: k += 100; break;
            case 7: k = -1; break;
            default: k += 1;
        }
    }
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(k);
    COMPARE_AND_SET_TEST_PASS(d);

    TEST_CASE_END;
}

bool test_long_loop() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0, i = 0;
    for(i = 0; i < 3000; ++i){ x += 2; }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    TrackedVar<int> x = 0, i = 0;
    for(i = 0; i < 3000; ++i){ x += 2; }
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_nested_loops() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0, i = 0, j = 0;
    for(i = 0; i < 60; ++i)
        for(j = 0; j < 60; ++j)
            x++;
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
    TrackedVar<int> x = 0, i = 0, j = 0;
    for(i = 0; i < 60; ++i)
        for(j = 0; j < 60; ++j)
            x++;
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(i);
    COMPARE_AND_SET_TEST_PASS(j);

    TEST_CASE_END;
}

bool test_unknown_loop_converges() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The trip count is unknown, but the loop stops producing new states.
    TEST_EVAL_STRING_CSTMT(R"(
    extern int y;
    int x = 0;
    while (y > 0) { x = 1; y--; }
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it, for values of the unknown `y` that run the loop zero and many times
    std::set<int> x_values;
    for (int y : {0, 3}) {
        TrackedVar<int> x = 0;
        while (y > 0) { x = 1; y--; }
        x_values.insert(x.values.begin(), x.values.end());
    }
    // clang-format on
    /////////////////////////////////////////////
    if (!test_analyzer.exec_ret || !test_state.eval_succeeded)
        test_pass = false;
    COMPARE_VALUES_AND_SET_TEST_PASS(x, x_values);

    TEST_CASE_END;
}

bool test_constexpr_member_operator_call() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A member function call `obj.f(a)` does not have the object in its arguments, a member
    // operator call `obj(a, b)` (CXXOperatorCallExpr) has it as the first argument. The
    // parameters of both are without the implicit object.
    TEST_EVAL_STRING_TOPLVL(R"(
struct Add2 {
    constexpr int f(int a) const { return a + 2; }
};
struct Sub {
    constexpr int operator()(int a, int b) const { return a - b; }
};
void test(){
    Add2 add2;
    Sub sub;
    int x = 0;
    x = add2.f(5);
    x = sub(5, 1);
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// local classes so we can keep them contained in the test function
struct Add2 {
    constexpr int f(int a) const { return a + 2; }
};
struct Sub {
    constexpr int operator()(int a, int b) const { return a - b; }
};
//void test(){
    Add2 add2;
    Sub sub;
    TrackedVar<int> x = 0;
    x = add2.f(5);
    x = sub(5, 1);
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);

    TEST_CASE_END;
}

bool test_constexpr_variadic_func() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // A C-style variadic function has more arguments than parameters. The extra arguments are
    // only evaluated for their side effects.
    TEST_EVAL_STRING_TOPLVL(R"(
constexpr int f(int a, ...){
    return a + 2;
}
void test(){
    int x = 0;
    int i = 0;
    x = f(5, i++, 3);
    x = f(i);
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// use a lambda so we can keep it contained in the test function ...
// constexpr int f(int a, ...){
constexpr auto f = [](int a, ...){
    return a + 2;
};
//void test(){
    TrackedVar<int> x = 0;
    TrackedVar<int> i = 0;
    x = f(5, i++, 3);
    x = f(i);
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);
    COMPARE_AND_SET_TEST_PASS(i);

    TEST_CASE_END;
}

bool test_constexpr_func_unknown_branch() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // The callee is evaluated with the worklist too: the branch on the unknown argument `b`
    // returns the same value on both paths.
    TEST_EVAL_STRING_TOPLVL(R"(
constexpr int func1(int a, int b){
    int r = a;
    if (b > 0) r = a;
    return r + 2;
}
extern int y;
void test(){
    int x = 0;
    for(int i = 0; i < 3; ++i){
        x = func1(i, y);
    }
}
    )")
    /////////////////////////////////////////////
    // clang-format off
    // Really evaluate it
// use a lambda so we can keep it contained in the test function ...
// constexpr int func1(int a, int b){
constexpr auto func1 = [](int a, int b){
    int r = a;
    if (b > 0) r = a;
    return r + 2;
};
// extern int y;
int y = 0;
//void test(){
    TrackedVar<int> x = 0;
    for(int i = 0; i < 3; ++i){
        x = func1(i, y);
    }
//}
    // clang-format on
    /////////////////////////////////////////////
    COMPARE_AND_SET_TEST_PASS(x);

    TEST_CASE_END;
}

////////////////////////////////////////////////////////
// test things that should fail exec: limits

bool test_should_fail_address_of_variable() {
    TEST_CASE_BEGIN;
    /////////////////////////////////////////////
    // `x` could be changed through the pointer, which is not followed
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    int *p = &x;
    *p = 5;
    int y = 0;
    y = x;
    )")
#pragma clang diagnostic pop
    /////////////////////////////////////////////
    if (test_analyzer.exec_ret || test_state.eval_succeeded)
        test_pass = false;

    TEST_CASE_END;
}

bool test_should_fail_max_block_states() {
    TEST_CASE_BEGIN;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    /////////////////////////////////////////////
    // A known loop with more iterations than the configured limit
    g_max_block_states = 20;
    TEST_EVAL_STRING_CSTMT(R"(
    int x = 0;
    for(int i = 0; i < 50; ++i){ x++; }
    )")
    g_max_block_states = 0;
#pragma clang diagnostic pop

    if (test_analyzer.exec_ret || test_state.eval_succeeded)
        test_pass = false;

    TEST_CASE_END;
}

////////////////////////////////////////////////////////
// test things that should fail exec

bool test_should_fail_1() {
    TEST_CASE_BEGIN;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    /////////////////////////////////////////////
    // Our ast_eval it
    TEST_EVAL_STRING_TOPLVL(R"(

void test(int unknown){
    int a = 0;
    for(int i = 0; i < unknown; ++i){
        a++;
    }
}
    )")
#pragma clang diagnostic pop

    if (test_analyzer.exec_ret || test_state.eval_succeeded)
        test_pass = false;

    TEST_CASE_END;
}

////////////////////////////////////////////////////////
// test things that should succeed to exec but fails to track some value

// TODO

typedef bool (*test_func)();
static std::vector<test_func> test_cases = {
    // clang-format off
    // for
    test_basic_for_loop,
    test_basic_for_loop_break,
    test_basic_for_loop_continue,

    // if
    test_basic_if,
    test_basic_if_else,
    test_basic_if_elseifs,
    test_basic_if_decl,

    // while
    test_basic_while,
    test_basic_while_break,
    test_basic_while_continue,

    // do while
    test_basic_do_while,
    test_basic_do_while_break,
    test_basic_do_while_continue,

    // switch
    test_basic_switch,
    test_basic_switch_no_default,
    test_basic_switch_no_default_fallthrough,
    test_basic_switch_case_range,
    test_switch_fallthrough,
    test_switch_in_loop,

    // goto
    test_basic_goto,

    // LOr, LAnd
    test_basic_LOr,
    test_basic_LAnd,
    test_basic_LAnd_LOr,

    // ternary
    test_basic_ternary,

    // side effects
    test_side_effects_with_unknown_operand,
    test_side_effects_in_store_to_non_variable,

    // changes through references
    test_ref_arg_unknown_function,
    test_ref_arg_constexpr_function,
    test_ref_arg_constexpr_increment,
    test_ref_arg_constructor,
    test_ref_arg_temporary,
    test_ref_variable,
    test_ref_variable_chains,
    test_ref_variable_unknown_target,

    // implicit conversions
    test_int_to_bool_conversion,

    // types, constants and other kinds of expressions
    test_other_expression_kinds,
    test_bool_and_char_literals,
    test_unevaluated_operands,
    test_floating_point_is_unknown,
    test_unsigned_and_narrowing,
    test_constants_declared_outside,
    test_enum_direction_loop,

    //
    test_lots_of_control_flow_things_1,
    test_lots_of_control_flow_things_2,

    // tracking the values of several expressions together
    test_pair_tracking_correlated,
    test_pair_tracking_direction_like,
    test_pair_tracking_unknown_branch,
    test_pair_tracking_unknown_value,

    // unknown conditions
    test_unknown_branch,
    test_unknown_ternary,
    test_unknown_LOr_LAnd_value,
    test_unknown_int_to_bool_cond,
    test_division_by_zero_on_unknown_path,
    test_switch_unknown_and_case_range,

    // long loops
    test_long_loop,
    test_nested_loops,
    test_unknown_loop_converges,

    // calling constexpr functions
    test_constexpr_func,
    test_constexpr_member_operator_call,
    test_constexpr_variadic_func,
    test_constexpr_func_unknown_branch,

    // should fail
    test_should_fail_1,
    test_should_fail_max_block_states,
    test_should_fail_address_of_variable,

    // weird stuff
    // test_Lor_in_switch,
    // clang-format on
};

#if TEST_IDX_EVAL
// Run all the test cases above
int main(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        char *end = nullptr;
        if (std::string(argv[i]) == "-v" && i + 1 < argc &&
            (g_verbosity = std::strtol(argv[i + 1], &end, 10),
             *end == '\0' && end != argv[i + 1])) {
            ++i;
        } else {
            llvm::errs() << "Usage: " << argv[0] << " [-v <number>]\n"
                         << "  -v 0   print only the test case results (TEST PASSED / FAILED)\n"
                         << "  -v 1   also print the output of the failed test cases (default)\n"
                         << "  -v 2   print the output of all test cases\n";
            return 1;
        }
    }

    int n_failed = 0;
    const size_t n_runs = test_cases.size();
    for (auto test_case : test_cases) {
        if (!test_case())
            n_failed++;
    }

    if (n_failed != 0) {
        llvm::outs() << "Test statistics: " << n_failed << " / " << n_runs;
        llvm::outs().changeColor(RED, true) << " test cases failed.";
        llvm::outs().resetColor() << "\n";
    } else {
        llvm::outs() << "Test statistics: " << n_failed << " / " << n_runs
                     << " test cases failed.\n";
        llvm::outs().changeColor(GREEN, true) << " All test cases passed.";
        llvm::outs().resetColor() << "\n";
    }

    return n_failed;
}

#elif IDX_EVAL_FILE
// For manual testing.
// Runs the idx_expr_evaluator on an input file and analyse a function with name `test`.

#include "clang/Tooling/CommonOptionsParser.h"
static cl::OptionCategory LoopVarTrackerOpts("LoopVarTrackerOpts options");

static cl::opt<bool> make_dot_graph("dot",                                 // The flag name
                                    cl::desc("Make .dot graph from CFG."), // Description for --help
                                    cl::init(false),                       // Default value
                                    cl::cat(LoopVarTrackerOpts) // Associate with your category
);

int main(int argc, const char **argv) {
    auto OptionsParser = CommonOptionsParser::create(argc, argv, LoopVarTrackerOpts);
    if (!OptionsParser) {
        exit(1);
    }
    auto &optParser = OptionsParser.get();
    ClangTool Tool(optParser.getCompilations(), optParser.getSourcePathList());

    TestFunctionAnalyzer Analyzer;
    MatchFinder Finder;

    Finder.addMatcher(functionDecl(hasName("test")).bind("func"), &Analyzer);

    g_print_all_possible_vals = true;
    if (make_dot_graph)
        g_make_dot_graph = true;

    const int ret = Tool.run(newFrontendActionFactory(&Finder).get());
    // the output of TestFunctionAnalyzer is collected in test_out()
    test_out().flush();
    llvm::outs() << g_test_output;
    return ret;
}
#endif
