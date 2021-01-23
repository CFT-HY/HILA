#include <sstream>
#include <iostream>
#include <string>

#include "stringops.h"
#include "hilapp.h"
#include "toplevelvisitor.h"


// helper function to finally resolve site dependence of calls

void check_site_dependence(call_info_struct & ci) {

  for (auto & arg : ci.arguments) {
    if (!arg.is_site_dependent) {
      for (auto & dv : arg.dependent_vars) {
        arg.is_site_dependent |= dv->is_site_dependent;
      }
    }
    ci.is_site_dependent |= arg.is_site_dependent;
  }

  // Was it a method?
  if (ci.is_method) {
    if (ci.object.is_lvalue) { 
      ci.object.is_site_dependent |= ci.is_site_dependent;
    } else if (!ci.is_site_dependent) {
      for (auto & dv : ci.object.dependent_vars) {
        ci.is_site_dependent |= dv->is_site_dependent;
      }
    }
  }
}


// This variable keeps track of visited decls - needed in order to avoid infinite recursion
static std::vector<Stmt *> visited_decls;


//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking if functions can be called from loops
/// 
/// Logic: find if it contains X (X_index_type), or field variables, or global variables
///  - not allowed in loop functions
/// Check also if they can be vectorized
//////////////////////////////////////////////////////////////////////////////


class loopFunctionVisitor : public GeneralVisitor, public RecursiveASTVisitor<loopFunctionVisitor> {

public:
  using GeneralVisitor::GeneralVisitor;

  bool contains_field;

  call_info_struct * this_ci;

  // store calls in this function
  std::vector<call_info_struct> calls;
  std::vector<var_info *> conditional_vars;

  std::string assignment_op;
  bool is_assginment, is_compound_assign;
  Stmt * assign_stmt;


  loopFunctionVisitor(Rewriter &R, ASTContext *C) : GeneralVisitor(R,C) {

    contains_field = false;
    is_assignment = false;
    var_info_list = {};
    var_decl_list = {};
    calls = {};
    conditional_vars = {};
  }


  
  bool VisitStmt(Stmt *s) { 

    // flag assignments from Stmt
    if (is_assignment_expr(s, &assignment_op, is_compound_assign)) {
      // This checks the "element<> -style assigns which we do not want now!
      assign_stmt = s;
      is_assignment = true;
      // next visit to declrefexpr will be to the assigned to variable
      return true;
    }
    
    // further function calls
    if( is_function_call_stmt(s) ){


      handle_function_call_in_loop(s);


    // And conditional stmts
    if (!this_ci->has_site_dependent_conditional) {
      Expr * condexpr = nullptr;
      if      (IfStmt    * IS = dyn_cast<IfStmt>(s))     condexpr = IS->getCond();
      else if (ForStmt   * FS = dyn_cast<ForStmt>(s))    condexpr = FS->getCond();
      else if (WhileStmt * WS = dyn_cast<WhileStmt>(s))  condexpr = WS->getCond();
      else if (DoStmt    * DS = dyn_cast<DoStmt>(s))     condexpr = DS->getCond();
      else if (SwitchStmt* SS = dyn_cast<SwitchStmt>(s)) condexpr = SS->getCond();
      else if (ConditionalOperator * CO = dyn_cast<ConditionalOperator>(s)) 
                                                         condexpr = CO->getCond();

      if (condexpr != nullptr) {
        this_ci->has_site_dependent_conditional = 
          is_site_dependent(condexpr, &conditional_vars);
        if (this_ci->has_site_dependent_conditional) {
          this_ci->condExpr = condexpr;
        }
      }

      return true;
    }



    return true;
  }

  //////////////////////////////////////////////////////////////////////////
  /// variable references are found here
  //////////////////////////////////////////////////////////////////////////


  bool VisitDeclRefExpr(DeclRefExpr * e) {
    /// if we see X or field, not good for loop function
    /// there are of course many others, I/O, Vectors, memory allocation...
    ///
    if (is_X_index_type(e) || is_field_expr(e)) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 e->getSourceRange().getBegin(),
                 "Field references are not allowed in functions called from site loops.");
      return false;  // stop here for this function

    }

    if (VarDecl * vdecl = dyn_cast<VarDecl>(e->getDecl())) {
      if (vdecl->hasExternalStorage() || vdecl->hasGlobalStorage()) {
        if (!cmdline::allow_func_globals) {
          reportDiag(DiagnosticsEngine::Level::Error,
                     e->getSourceRange().getBegin(),
                     "global or extern variable references in functions called from site loops are not allowed."
                     "\nThis can be enable in non-kernelized code with option '-allow-func-globals'" );
          return false;
        } else {
          if (e->isLValue()) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       e->getSourceRange().getBegin(),
                       "modification of global or extern variables in functions called from site loops is not allowed.");
            return false;
          }
          reportDiag(DiagnosticsEngine::Level::Warning,
                     e->getSourceRange().getBegin(),
                     "global or extern variable references in site loop functions make "
                     "code non-portable to kernelized code (e.g. GPU code).");
          // just continue after this warning
        }
      }
      
      // handle the local var ref.
      handle_var_ref(e,vdecl);

      
    }
    return true;
  }


  //////////////////////////////////////////////////////////////////////////////////
  /// and variable refs - pretty much copied from myastvisitor
  //////////////////////////////////////////////////////////////////////////////////

  void handle_var_ref(DeclRefExpr * DRE, VarDecl *decl) {

    var_ref vr;
    vr.ref = DRE;
    //vr.ind = writeBuf->markExpr(DRE);
    vr.is_assigned = is_assignment;
    if (is_assign) vr.assignop = assignment_op;


    bool foundvar = false;
    var_info *vip = nullptr;
    for (var_info & vi : var_info_list) {
      if (vi.decl == decl) {
        // found already referred to decl
        // check if this particular ref has been handled before
        bool foundref = false;
        for (auto & r : vi.refs) if (r.ref == DRE) {
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
        vi.is_assigned |= is_assignment;
        vi.reduction_type = reduction::NONE
        vip = &vi;
        foundvar = true;
        break;
      }
    }
    if (!foundvar) {
      // new variable referred to
      vip = new_var_info(decl);

      vip->refs.push_back(vr);
      vip->is_assigned = is_assignment;
    }

    if (is_assignment && assign_stmt != nullptr && !vip->is_site_dependent) {
      vip->is_site_dependent = is_rhs_site_dependent(assign_stmt, &vip->dependent_vars );
      
      // llvm::errs() << "Var " << vip->name << " depends on site: " << vip->is_site_dependent <<  "\n";
    }

    is_assignment = false;  // not an assignment any more

    return;
  }

  //////////////////////////////////////////////////////////////////////////////////////

  var_info * new_var_info(VarDecl *decl) {

    var_info vi;
    vi.refs = {};
    vi.decl = decl;
    vi.name = decl->getName().str();
    // Printing policy is somehow needed for printing type without "class" id
    // Unqualified takes away "consts" etc and Canonical typdefs/using.
    // Also need special handling for element type
    clang::QualType type = decl->getType().getUnqualifiedType().getNonReferenceType();
    vi.type = type.removeLocalConst().getAsString(PP);
    vi.type = remove_all_whitespace(vi.type);
    bool is_elem = (vi.type.find("element<") == 0);
    vi.type = type.getAsString(PP);
    if (is_elem) vi.type = "element<" + vi.type + ">";
    // llvm::errs() << " + Got " << vi.type << '\n';

    vi.is_loop_local = true;

    vi.is_site_dependent = false;  // default case

    // if this is a parameter var decl, mark the site dep.
    if (decl->isLocalVarDeclOrParm() && !decl->isLocalVarDecl()) {
      ParmVarDecl * pv = dyn_cast<ParmVarDecl>(decl);
      // decl is now parameter var
      llvm::errs() << "FOUND PARAMETER VAR " << vi.name << " of type " << vi.type;

      bool found = false;
      for (auto & arg : this_ci->arguments) {
        if (arg.PV == pv) {
          vi.is_site_dependent = arg.is_site_dependent;
          llvm::errs() << " site dep " << vi.is_site_dependent;
          found = true;
          break;
        }
      }
      llvm::errs() << '\n';

      if (!found) {
        llvm::errs() << " ERROR: PARAMETER NOT FOUND IN CALLINFO!  EXIT HERE\n";
        exit(1);
      }      
    }

    vi.dependent_vars.clear();

    var_info_list.push_back(vi);
    return &(var_info_list.back());
  }


  ///////////////////////////////////////////////////////////////////////////////////
  /// Check non-trivial declarations
  ///////////////////////////////////////////////////////////////////////////////////

  bool VisitDecl(Decl *D) {

    if (VarDecl * V = dyn_cast<VarDecl>(D)) {
      // it's a variable decl inside function
      if (V->getStorageClass() == StorageClass::SC_Extern ||
          V->getStorageClass() == StorageClass::SC_Static ||
          V->getStorageClass() == StorageClass::SC_PrivateExtern) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   D->getSourceRange().getBegin(),
                   "cannot declare static or extern variables in functions called from site loops.");
        return false;
      }

    }
    return true;
  }

  void handle_function_call(Stmt * s) {

    // Get the call expression
    CallExpr *Call = dyn_cast<CallExpr>(s);
    
    // Handle special loop functions
    if( handle_special_functions(Call) ){
      return;
    }

    // Get the declaration of the function
    Decl* decl = Call->getCalleeDecl();
    FunctionDecl* D = (FunctionDecl*) llvm::dyn_cast<FunctionDecl>(decl);

    bool contains_rng = false;
    if (D->hasBody()) {
      // trivially site dep if it has random
      contains_rng = contains_random(D->getBody());
    } else {
      // TODO - these functions are at least not vectorizable ...

      SourceManager &SM = TheRewriter.getSourceMgr();
      llvm::errs() << "FUNC DECL WITHOUT BODY IN LOOP FUNC - " << D->getNameAsString() << '\n';
      llvm::errs() << "  Call appears on line " << SM.getSpellingLineNumber(Call->getBeginLoc())
          << " in file " << SM.getFilename(Call->getBeginLoc()) << '\n';

    }

    // check the arg list
    call_info_struct ci = handle_loop_function_args(D,Call,contains_rng);
    ci.call = Call;
    ci.decl = D;
    ci.contains_random = contains_rng;
    
    /// add to function calls to be checked ...
    calls.push_back(ci);
  }
  


  ///////////////////////////////////////////////////////////////////////////////////
  /// Loop through functions seen
  ///////////////////////////////////////////////////////////////////////////////////


  void visit_calls( std::vector<call_info_struct> & calls ) {

    for (auto & ci : calls ) {
      check_site_dependence(ci);

      // spin new visitor for all calls here

      loopFunctionVisitor visitor(TheRewriter,Context);
      visitor.start_visit(ci);
    }
  }


  void start_visit( call_info_struct & ci ) {

    // TODO: Should one here check first the vectorization of parameters?

    this_ci = &ci;

    Stmt * decl_body = nullptr;
    if (ci.call != nullptr) {
      // a function now, does it have a body?
      if (ci.decl->hasBody()) decl_body = ci.decl->getBody();
    } else if (ci.constructor != nullptr) {
      // same stuff for constructor
      if (ci.ctordecl->hasBody()) decl_body = ci.ctordecl->getBody();
    }

    if (decl_body) {
      for (auto d : visited_decls) if (d == decl_body) {
        // it was checked, return and continue
        return true;
      }
    } else {
      // now decl has no body: TODO: handling!!!

      if (ci.call != nullptr)
        llvm::errs() << "Loop func decl has no body: " << ci.decl->getNameAsString() << '\n';
      else if (ci.constructor != nullptr)
        llvm::errs() << "Loop constructor decl has no body: " << ci.ctordecl->getNameAsString() << '\n';

      return true;
    }

    // mark this as visited
    visited_decls.push_back(decl_body);

    // push the param vars to var list
    for (auto & arg : ci.arguments) {
      var_info vi;
      vi.
      arg.
    }

    // start the visit here
    TraverseStmt(decl_body);

    return true;

  }


};






void TopLevelVisitor::visit_loop_functions( std::vector<call_info_struct> & calls ) {


  visited_decls.clear();

  for (auto & ci : calls ) {
    // verify first the site dep of args and the function
    // probably this stage is needed only in vanishingly obscure loops

    check_site_dependence(ci);

    // spin new visitor for all calls here

    loopFunctionVisitor visitor(TheRewriter,Context);
    visitor.start_visit(ci);

  }

  visited_decls.clear();

  // Now the calls should contain full info about the calls
  // Visit all calls, and functions inside them hiearchially

}
