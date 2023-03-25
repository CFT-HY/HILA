#ifndef GENERALVISITOR_H
#define GENERALVISITOR_H

#include <string>
#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"

#include "srcbuf.h" //srcbuf class interface
#include "hilapp.h" //global vars needed


/////////////////////////////////////////////////////////////////////////////////////
/// This implements a "general visitor", a placeholder to many utility functions
/// which can be used to visit AST
///
/// Example:
///
/// class exampleVisitor : public GeneralVisitor, public
/// RecursiveASTVisitor<exampleVisitor> {
///
///    // if you need a special constructor, do it like this:
///   public:
///    exampleVisitor(Rewriter & R, ASTContext *C, ... other vars... ) :
///    GeneralVisitor(R,C) {
///       ... custom constructor stmts ...
///    }
///    ... whatever definitions ...
///
///    // if constructor is not needed, use the GeneralVisitor constructor:
///    using GeneralVisitor::GeneralVisitor;
///
/// };
///
/// exampleVisitor ev(TheRewriter,Context, .... );
///    // here TheRewriter is of type Rewriter &, Context of type ASTContext *
///    // these are defined in TopLevelVisitor, so this works directly there
/// ev.TraverseStmt( .. some Stmt * variable ..);   // or TraverseDecl or some other
/// Traverse
///
/////////////////////////////////////////////////////////////////////////////////////
class GeneralVisitor {

  private:
    // put stubs of these variables here.  These are either unused or
    // _lists below will point to these

    std::list<var_info> _var_info_list;
    std::list<var_decl> _var_decl_list;

  protected:
    /// Reference to the Clang rewriter object
    Rewriter &TheRewriter;

    /// Store the context pointer
    ASTContext *Context;

    /// SourceManager is used all over the place
    SourceManager &srcMgr;

    /// Store a printing policy. It is required quite often
    PrintingPolicy PP;

    /// are we on TopLevelVisitor?  False otherwise
    bool is_top_level = false;

    /// store the variables at this level in var_info_list
    /// use list because there will be pointers to elements and list elems
    /// do not move

    std::list<var_info> &var_info_list;

    std::list<var_decl> &var_decl_list;

    // calls to loop functions done in this loop or function.  Always cleared for
    // every new visitor
    std::vector<call_info_struct> loop_function_calls = {};

  public:
    /// Construct from previous visitor - we inherit the lists of the caller
    template <typename Visitor_type>
    GeneralVisitor(Visitor_type &caller)
        : TheRewriter(caller.TheRewriter), Context(caller.Context),
          PP(caller.Context->getLangOpts()), srcMgr(caller.srcMgr),
          var_info_list(caller.var_info_list), var_decl_list(caller.var_decl_list) {}

    /// Construct with rewriter and context
    GeneralVisitor(Rewriter &R, ASTContext *C)
        : TheRewriter(R), Context(C), srcMgr(R.getSourceMgr()), PP(C->getLangOpts()),
          var_info_list(_var_info_list), var_decl_list(_var_decl_list) {}

    Rewriter &getRewriter() {
        return TheRewriter;
    }
    ASTContext *getASTContext() {
        return Context;
    }

    /// Report diagnostic info
    template <unsigned N>
    void reportDiag(DiagnosticsEngine::Level lev, const SourceLocation &SL, const char (&msg)[N],
                    const char *s1 = nullptr, const char *s2 = nullptr, const char *s3 = nullptr) {
        // we'll do reporting only when output is on, avoid double reports
        auto &DE = Context->getDiagnostics();
        auto ID = DE.getCustomDiagID(lev, msg);
        auto DB = DE.Report(SL, ID);
        if (s1 != nullptr)
            DB.AddString(s1);
        if (s2 != nullptr)
            DB.AddString(s2);
        if (s3 != nullptr)
            DB.AddString(s3);
    }

    /// Return the unexpanded range (without macro expansion) from the source
    SourceRange get_real_range(SourceRange r) {
        if (r.getBegin().isMacroID()) {
            CharSourceRange CSR =
                TheRewriter.getSourceMgr().getImmediateExpansionRange(r.getBegin());
            r = SourceRange(CSR.getAsRange().getBegin(), r.getEnd());
        }
        if (r.getEnd().isMacroID()) {
            CharSourceRange CSR = TheRewriter.getSourceMgr().getImmediateExpansionRange(r.getEnd());
            r = SourceRange(r.getBegin(), CSR.getAsRange().getEnd());
        }
        return r;
    }


    /// Get the written expression of a C++ statement
    std::string get_stmt_str(const Stmt *s) {
        return TheRewriter.getRewrittenText(get_real_range(s->getSourceRange()));
    }

    /// Get the type of an expression (i.e. int, double...) as a string
    std::string get_expr_type(const Expr *e) {
        // This is somehow needed for printing type without "class" id
        // TODO: perhaps reanalyse and make more robust?
        return e->getType().getCanonicalType().getUnqualifiedType().getAsString(PP);
    }

    /// check if stmt contains random number generator
    bool contains_random(Stmt *s);

    // shofthand for obtaining file buffer within this class
    srcBuf *get_file_srcBuf(SourceLocation sl) {
        return get_file_buffer(TheRewriter, TheRewriter.getSourceMgr().getFileID(sl));
    }

    /// Get FileID for SourceLocation
    FileID get_FileId(SourceLocation sl) {
        return TheRewriter.getSourceMgr().getFileID(sl);
    }

    /// a list of utility inspection functions
    /// getCanonicalType takes away typedefs, getUnqualifiedType() qualifiers, pp just
    /// in case the type string needs to begin with the string
    bool is_field_storage_expr(Expr *E) {
        return (E && E->getType().getCanonicalType().getUnqualifiedType().getAsString(PP).find(
                         field_storage_type) == 0);
    }

    /// Check the expression of a field expression
    bool is_field_expr(Expr *E) {
        return (E && E->getType().getCanonicalType().getUnqualifiedType().getAsString(PP).find(
                         field_type) == 0);
    }

    /// Check if declaration declate a field
    bool is_field_decl(ValueDecl *D) {
        return (D && D->getType().getCanonicalType().getUnqualifiedType().getAsString(PP).find(
                         field_type) == 0);
    }

    /// try to figure out whether expressions are duplicates
    bool is_duplicate_expr(const Expr *a, const Expr *b);

    /// Just check that the expression is of type Parity
    bool is_parity_index_type(Expr *E) {
        return (get_expr_type(E) == "Parity");
    }

    /// Checks if E is of type Field[Parity]  Parity=EVEN,ODD,ALL
    bool is_field_parity_expr(Expr *E);

    /// true if X, X+dir, X+offset -type expression
    bool is_X_type(Expr *E) {
        std::string s = get_expr_type(E);
        return (s == "X_index_type");
    }

    /// true if X, X+dir, X+offset -type expression
    bool is_X_index_type(Expr *E) {
        std::string s = get_expr_type(E);
        return (s == "X_index_type" || s == "X_plus_direction" || s == "X_plus_offset");
    }

    /// true if X+dir - type
    bool is_X_and_dir_type(Expr *E) {
        std::string s = get_expr_type(E);
        return (s == "X_plus_direction" || s == "X_plus_offset");
    }

    bool is_coordinate_type(Expr *E) {
        std::string s = get_expr_type(E);
        return (s.find("CoordinateVector") == 0);
    }

    bool is_field_with_X_expr(Expr *E);

    /// Checks if E is Parity plus Direction of a field (for example f[X+dir]).
    bool is_field_with_X_and_dir(Expr *E);

    bool is_field_with_coordinate(Expr *E);

    bool is_assignment_expr(Stmt *s, std::string *opcodestr, bool &iscompound,
                            Expr **assignee = nullptr, Expr **assigned_expr = nullptr);

    bool is_simple_reduction(const std::string &opcode, Expr *assignee);

    bool is_increment_expr(Stmt *s, Expr **assignee = nullptr);

    bool is_site_dependent(Expr *e, std::vector<var_info *> *dependent_var);

    bool is_loop_constant(Expr *e);

    bool class_has_base_type(CXXRecordDecl *d);

    bool is_rhs_site_dependent(Stmt *s, std::vector<var_info *> *vi);

    bool is_site_dependent_access_op(Expr *e);

    bool contains_loop_local_var(Expr *e, std::vector<var_info *> *loop_var = nullptr);

    bool contains_special_reduction_var(Expr *e);

    /// is the stmt pointing now to a function call
    bool is_function_call_stmt(Stmt *s);

    /// is the stmt pointing now to a member call
    bool is_member_call_stmt(Stmt *s);

    /// is the stmt pointing now to a constructor
    bool is_constructor_stmt(Stmt *s);

    /// is the stmt pointing now to a constructor
    bool is_user_cast_stmt(Stmt *s);

    /// inspect if the type name is vectorizable
    /// returns the vectorized type in vectorized_type, if it is
    bool is_vectorizable_type(const std::string &type_name, vectorization_info &vi);
    bool is_vectorizable_type(const QualType &QT, vectorization_info &vi);

    /// Does the statement end with a semicolon
    bool isStmtWithSemicolon(Stmt *S);

    Parity get_parity_val(const Expr *pExpr);

    /// utility used in inserting stuff after new line in buffer
    // SourceLocation getSourceLocationAtEndOfLine( SourceLocation l );
    /// another utility (cannot trust r.getEnd())
    SourceLocation getSourceLocationAtEndOfRange(SourceRange r);

    /// Find the start loc of declaration, taking into account all template<>'s
    SourceLocation getSourceLocationAtStartOfDecl(Decl *d);

    // get next char and loc, while skipping comments
    SourceLocation getNextLoc(SourceLocation sl, bool forward = true);

    char getChar(SourceLocation sl);
    SourceLocation findChar(SourceLocation sl, char ch);

    // get next word or symbol, if it is not a legal name symbol
    std::string getNextWord(SourceLocation sl, SourceLocation *end = nullptr);
    std::string getPreviousWord(SourceLocation sl, SourceLocation *start = nullptr);

    /// jump over following () expr
    SourceLocation skipParens(SourceLocation sl);

    var_info *handle_var_ref(DeclRefExpr *E, bool is_assign, const std::string &op,
                             Stmt *assign_stmt = nullptr, bool is_raw = false);

    var_info *new_var_info(VarDecl *decl);

    var_info *add_var_to_decl_list(VarDecl *var, int scope);

    void handle_constructor_in_loop(Stmt *s);

    bool handle_loop_function_if_needed(call_info_struct &ci);

    call_info_struct handle_loop_function_args(FunctionDecl *D, CallExpr *Call, bool sitedep,
                                               bool is_assignment = false);

    bool handle_call_argument(Expr *E, ParmVarDecl *pv, bool sitedep,
                              std::vector<var_info *> *out_variables,
                              std::vector<var_info *> *dep_variables, argument_info &ai);

    bool attach_dependent_vars(std::vector<var_info *> &variables, bool sitedep,
                               std::vector<var_info *> &dep_variables);

    void backend_handle_loop_function(call_info_struct &ci);
    void backend_handle_loop_constructor(call_info_struct &ci);

    /// Handle functions called in a loop
    void handle_loop_function_gpu(call_info_struct &ci);
    void handle_loop_function_openacc(FunctionDecl *fd);
    void handle_loop_function_avx(call_info_struct &ci);

    /// Handle functions called in a loop
    void handle_loop_constructor_gpu(call_info_struct &ci);
    void handle_loop_constructor_openacc(CXXConstructorDecl *fd);
    void handle_loop_constructor_avx(call_info_struct &ci);

    /// True if the decl is preceded by "#pragma hila <string>" where s is the string
    bool has_pragma(Decl *d, const pragma_hila p, const char **arg = nullptr);
    bool has_pragma(FunctionDecl *fd, const pragma_hila p, const char **arg = nullptr);
    bool has_pragma(Stmt *S, const pragma_hila p, const char **arg = nullptr);
    bool has_pragma(const SourceLocation sl, const pragma_hila p, const char **arg = nullptr);
};

#endif // ifdef GENERALVISITOR_H