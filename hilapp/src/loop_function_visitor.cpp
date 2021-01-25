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
  std::vector<call_info_struct> function_calls;
  std::vector<var_info *> conditional_vars;

  std::vector<special_function_call> special_call_list;

  std::string assignment_op;
  bool is_assginment, is_compound_assign;
  Stmt * assign_stmt;

  // use the constructor which does not inherit var_info_list and var_decl_list

  loopFunctionVisitor(Rewriter &R, ASTContext *C) : GeneralVisitor(R,C) {

    contains_field = false;
    is_assignment = false;
    var_info_list = {};
    var_decl_list = {};
    special_call_list = {};
    function_calls = {};
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
      handle_function_call(s);
      return true;
    }

    // further function calls
    if( is_constructor_stmt(s) ){
      handle_constructor_in_loop(s);
      return true;
    }


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
        this_ci->condExpr = condexpr;
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



  ///////////////////////////////////////////////////////////////////////////////////
  /// Check non-trivial declarations
  ///////////////////////////////////////////////////////////////////////////////////

  bool VisitVarDecl(VarDecl *V) {

    // it's a variable decl inside function
    if (V->getStorageClass() == StorageClass::SC_Extern ||
        V->getStorageClass() == StorageClass::SC_Static ||
        V->getStorageClass() == StorageClass::SC_PrivateExtern) {
      reportDiag(DiagnosticsEngine::Level::Error,
                  D->getSourceRange().getBegin(),
                  "cannot declare static or extern variables in functions called from site loops.");
      return false;
    }
    // Now it should be automatic local variable decl

    add_var_to_decl_list(V,0);

    return true;
  } 

  //////////////////////////////////////////////////////////////////////////
  /// log function calls in loop functions
  //////////////////////////////////////////////////////////////////////////


  void handle_function_call(Stmt * s) {

    // Get the call expression
    CallExpr *Call = dyn_cast<CallExpr>(s);
    
    // Handle special loop functions
    if( handle_special_function(Call) ){
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

      llvm::errs() << "FUNC DECL WITHOUT BODY IN LOOP FUNC - " << D->getNameAsString() << '\n';
      llvm::errs() << "  Call appears on line " << srcMgr.getSpellingLineNumber(Call->getBeginLoc())
          << " in file " << srcMgr.getFilename(Call->getBeginLoc()) << '\n';

    }

    // check the arg list
    call_info_struct ci = handle_loop_function_args(D,Call,contains_rng);
    ci.call = Call;
    ci.decl = D;
    ci.contains_random = contains_rng;
    
    /// add to function calls to be checked ...
    function_calls.push_back(ci);
  }
  



  /////////////////////////////////////////////////////////////////////////////////
  /// Special functions:  methods X.method(), and hila_random()
  /////////////////////////////////////////////////////////////////////////////////

  bool handle_special_function(CallExpr *Call) {
    // If the function is in a list of defined loop functions, add it to a list
    // Return true if the expression is a special function

    std::string name = Call->getDirectCallee()->getNameInfo().getAsString();

    if (CXXMemberCallExpr *MCall = dyn_cast<CXXMemberCallExpr>(Call)) {
      // llvm::errs() << "It's a member call, name " << name << " objarg "
      //       << MCall->getImplicitObjectArgument()->getType().getAsString() << "\n";
      //    std::string objtype = MCall->getImplicitObjectArgument()->getType().getAsString();
      std::string objtype = get_expr_type(MCall->getImplicitObjectArgument());  
      if (objtype == "lattice_struct *") {

        llvm::errs() << "CALL in LOOP FUNC: " << get_stmt_str(Call) << '\n';
        special_function_call sfc;
        sfc.fullExpr = Call;
        sfc.name = name;
        sfc.argsExpr = nullptr;

        SourceLocation sl = findChar(Call->getSourceRange().getBegin(),'(');
        if (sl.isInvalid()) {
          reportDiag(DiagnosticsEngine::Level::Fatal,
                    Call->getSourceRange().getBegin(),
                    "Open parens '(' not found, internal error");
          exit(1);
        }
        sfc.replace_range = SourceRange(sfc.fullExpr->getSourceRange().getBegin(), sl);

        bool replace_this = true;   // for non-cuda code replace only cases which are needed
        if (name == "size") {
          sfc.replace_expression = "loop_lattice_size(";
          sfc.add_loop_var = false;
          replace_this = target.CUDA;
        } else if (name == "volume") {
          sfc.replace_expression = "loop_lattice_volume(";
          sfc.add_loop_var = false;
          replace_this = target.CUDA;
        } else {
          reportDiag(DiagnosticsEngine::Level::Error,
            Call->getSourceRange().getBegin(),
          "Method 'lattice->.%0()' not allowed inside site loops", name.c_str() );
        }

        if (replace_this) special_call_list.push_back(sfc);

        return true;

      } else {
        
        // other method calls?
        return false;
      }

    } else {
      if( name == "hila_random" ){

        // no need to do anything (yet)

        // special_function_call sfc;
        // sfc.fullExpr = Call;
        // sfc.argsExpr = nullptr;
        // sfc.scope = parsing_state.scope_level;
        // sfc.name = name;
        // sfc.replace_expression = "hila_random()";
        // sfc.replace_range = Call->getSourceRange();  // replace full range
        // sfc.add_loop_var = false;
        // special_function_call_list.push_back(sfc);

        return true;
      }
    }

    return true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /// Check now that the references to variables are as required
  /////////////////////////////////////////////////////////////////////////////

  void check_var_info() {

    // iterate through var_info_list until no more is_site_dependent -relations found
    // this should not leave any corner cases behind

    int found;
    do {
      found = 0;
      for (var_info & vi : var_info_list) {
        if (vi.is_site_dependent == false) {
          for (var_info * d : vi.dependent_vars) if (d->is_site_dependent) {
            vi.is_site_dependent = true;
            found++;
            break;  // go to next var
          }
        }
      } 
    } while (found > 0);

    // and also get the vectorized type for them, to be prepared...

    if (target.vectorize) {
      for (var_info & vi : var_info_list) {
        vi.vecinfo.is_vectorizable = is_vectorizable_type(vi.decl->getType(),vi.vecinfo);
      }
    }

  }



  
  ///////////////////////////////////////////////////////////////////////////////////
  /// Visit each function seen
  /// return false if things went wrong
  ///////////////////////////////////////////////////////////////////////////////////

  bool do_visit( call_info_struct & ci ) {


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

      return false;
    }

    // mark this as visited
    visited_decls.push_back(decl_body);

    // push the param vars to var list
    for (auto & arg : ci.arguments) {
      add_var_to_decl_list( arg.PV, 0 );     // no need to worry about scope levels
    }

    // start the visit here
    TraverseStmt(decl_body);

    // recheck variables
    check_var_info();

    // go and visit all of the calls inside (recursively)
    visit_calls();

    // post-process vectorization info
    if (target.vectorize) {
      ci.is_vectorizable = !ci.contains_random;
      if (ci.is_vectorizable) {
        for (auto & func : function_calls) {
          ci.is_vectorizable &= ( func.is_vectorizable || !func.is_site_dependent);
        }
      }
      if (ci.is_vectorizable) {
        for (auto & vi : var_info_list) {
          ci.is_vectorizable &= ( vi.is_vectorizable || !vi.is_site_dependent);
        }
      }
    }

    

    return true;

  }

  ///////////////////////////////////////////////////////////////////////////////////
  /// Loop through functions seen - almost copy of the visit_loop_functions below, but
  /// callable from the visitor itself
  ///////////////////////////////////////////////////////////////////////////////////


  void visit_calls() {

    for (auto & ci : function_calls ) {
      check_site_dependence(ci);

      // spin new visitor for all calls here

      loopFunctionVisitor visitor(TheRewriter,Context);
      visitor.do_visit(ci);
    }
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
    visitor.do_visit(ci);

  }

  visited_decls.clear();

  // Now the calls should contain full info about the calls
  // Visit all calls, and functions inside them hiearchially

}
