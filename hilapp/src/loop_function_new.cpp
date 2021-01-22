#include <sstream>
#include <iostream>
#include <string>
#include "clang/Analysis/CallGraph.h"

#include "myastvisitor.h"
#include "hilapp.h"
#include "stringops.h"

// collect here all loop functions in a compilation unit

static std::vector<FunctionDecl *> loop_functions = {};
static std::vector<CXXConstructorDecl *> loop_constructors = {};

// and these hold the loop-local function calls, for revisit


static std::vector<call_info_struct> loop_function_calls = {};


// clear function ptrs in whole compilation unit
void clear_loop_functions_in_compilation_unit() {
  loop_functions.clear();
  loop_constructors.clear();
}

// and clear calls for this loop
void MyASTVisitor::clear_loop_function_calls() {
  loop_function_calls.clear();
}


////////////////////////////////////////////////////////////////////////////
/// Entry for function calls inside loops.  The call requires a bit
/// special treatment, the arguments can be field[X] etc. elements
////////////////////////////////////////////////////////////////////////////

void MyASTVisitor::handle_function_call_in_loop(Stmt * s) {

  // Get the call expression
  CallExpr *Call = dyn_cast<CallExpr>(s);

  assert(Call && "Loop function call not valid");

  // Handle special loop functions
  if( handle_special_loop_function(Call) ){
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

    llvm::errs() << "FUNC DECL WITHOUT BODY - " << D->getNameAsString() << '\n';
    llvm::errs() << "  Call appears on line " << SM.getSpellingLineNumber(Call->getBeginLoc())
        << " in file " << SM.getFilename(Call->getBeginLoc()) << '\n';

  }

  // check the arg list
  call_info_struct ci = handle_loop_function_args(D,Call,contains_rng);
  ci.call = Call;
  ci.decl = D;
  ci.contains_random = contains_rng;
  
  /// add to function calls to be checked ...
  loop_function_calls.push_back(ci);

  

  // Store functions used in loops, recursively...
  loop_function_check(decl);

}

////////////////////////////////////////////////////////////////////////////
/// And entry point for constructors inside loops.
////////////////////////////////////////////////////////////////////////////


void MyASTVisitor::handle_constructor_in_loop(Stmt * s) {

  // Get the call expression
  CXXConstructExpr *CtorE = dyn_cast<CXXConstructExpr>(s);

  assert(CtorE && "Constructor call in loop not valid");

  // Get the declaration of the constructor
  CXXConstructorDecl* decl = CtorE->getConstructor();

  //llvm::errs() << " callee:\n";
  //decl->dump();

  llvm::errs() << "  Constructor " << decl->getNameAsString() << '\n';

  // Store functions used in loops, recursively...

  // handle args if any

  bool contains_rng = false;
  if (decl->hasBody()) {
    // trivially site dep if it has random
    contains_rng = contains_random(decl->getBody());
  } else {
    llvm::errs() << "CONSTRUCTOR DECL WITHOUT BODY - TODO HANDLING\n";
  }

  // Handle parameters
  if (decl->getNumParams() != CtorE->getNumArgs()) {
    reportDiag(DiagnosticsEngine::Level::Fatal,
               CtorE->getSourceRange().getBegin(),
               "Internal error: #params %0, #args %1 in constructor", 
               std::to_string(decl->getNumParams()).c_str(),
               std::to_string(CtorE->getNumArgs()).c_str() );
    return;
  }
  
  std::vector<var_info *> out_variables;
  std::vector<var_info *> dep_variables;

  call_info_struct ci;

  ci.constructor = CtorE;
  ci.ctordecl = decl;
  ci.contains_random = contains_rng;

  // go through the args - If contains field[X] or site dep. vars 
  // whole call is site dep.
  // Lvalue vars can change
  bool is_site_dependent = false;
  for( int i=0; i<CtorE->getNumArgs(); i++) {

    Expr * E = CtorE->getArg(i);
    const ParmVarDecl * pv = decl->getParamDecl(i);

    argument_info ai;
    is_site_dependent |= handle_call_argument(E, pv, is_site_dependent, 
                                             &out_variables, &dep_variables, ai);
    ci.is_site_dependent |= ai.is_site_dependent;
    ci.arguments.push_back(ai);

  } // end of arg loop.

  // and check the var dependency
  is_site_dependent = attach_dependent_vars(out_variables, is_site_dependent, dep_variables);
  ci.is_site_dependent = is_site_dependent;

  // and add the call  to check-up list
  loop_function_calls.push_back(ci);

  SourceManager &SM = Context->getSourceManager();
  bool handle_decl = !SM.isInSystemHeader(decl->getBeginLoc());
  
  // check if we already have this declaration - either the pointer is the same
  // or the source location (actually, source location should do all, no need for 
  // CXXConstructorDecl *, but it does not hurt)
  for (int i=0; handle_decl && i<loop_constructors.size(); i++) { 
    if (decl == loop_constructors[i] || 
        decl->getSourceRange().getBegin() == 
          loop_constructors[i]->getSourceRange().getBegin() ) {
      handle_decl = false;
    }
  }
  if (handle_decl) {
    if (decl->isTrivial()) {
      llvm::errs() << "TRIVIAL CONSTRUCTOR " << decl->getNameAsString() << " DO NOTHING\n";
    } else {
      loop_constructors.push_back(decl);
        llvm::errs() << "NEW LOOP CONSTRUCTOR " << decl->getNameAsString() << '\n';
        // " parameters ";
        //  for (int i=0; i<fd->getNumParams(); i++) 
        //  llvm::errs() << fd->getParamDecl(i)->getOriginalType().getAsString();
        // llvm::errs() << '\n';
    
      backend_handle_loop_constructor(decl);
    }
  }
}



/////////////////////////////////////////////////////////////////////////////////
/// This handles the args of a "first level" function, i.e. function
/// called directly from site loop
/// Need to find out:
///  - are args site dependent -> function site dependent
///  - lvalue args inherit site dependence
///  - lvalue field element args are assumed changed
///  - output_only flagged field elements need not be input
/////////////////////////////////////////////////////////////////////////////////


call_info_struct MyASTVisitor::handle_loop_function_args(FunctionDecl *D, CallExpr *Call,
                                                         bool sitedep) {

  call_info_struct cinfo;


  #define LOOP_FUNC_DEBUG
  #ifdef LOOP_FUNC_DEBUG
  llvm::errs() << "LOOP FUNC " << D->getNameAsString() << " with "
  << D->getNumParams() << " parameters and " << Call->getNumArgs() << " arguments\n";

  llvm::errs() << "Is it a method? " << isa<CXXMemberCallExpr>(Call) << '\n';

  llvm::errs() << "   Func args: ";
  for (Expr * E : Call->arguments()) {
    llvm::errs() << get_stmt_str(E);
    if (E->isLValue()) llvm::errs() << "-LVALUE";
    llvm::errs() << ", ";
  }
  llvm::errs() << "\n   Func params: ";
  for (int i=0; i<D->getNumParams(); i++) llvm::errs() << D->getParamDecl(i)->getNameAsString() << ", ";
  llvm::errs() << "\n";

  #endif

  cinfo.is_site_dependent = sitedep;

  // for operators, the dependency in args is handled separately in myastvisitor 
  // (Really, should treat everything here but started with that)
  if (isa<CXXOperatorCallExpr>(Call)) {
    cinfo.is_operator = true;
    return cinfo;
  }

  if (D->getNumParams() != Call->getNumArgs()) {
      SourceManager &SM = Context->getSourceManager();

      llvm::errs() << "Internal error: #params != #args, function " << D->getNameAsString() << '\n';
      llvm::errs() << "  Call appears on line " << SM.getSpellingLineNumber(Call->getBeginLoc())
                   << " in file " << SM.getFilename(Call->getBeginLoc()) << '\n';
      llvm::errs() << "  Function is defined on line " << SM.getSpellingLineNumber(D->getBeginLoc())
                   << " in file " << SM.getFilename(D->getBeginLoc()) << '\n';

      exit(1);
  }
  
  std::vector<var_info *> out_variables;
  std::vector<var_info *> dep_variables;

  // go through the args - If contains field[X] or site dep. vars 
  // whole call is site dep.
  // Lvalue vars can change
  for( int i=0; i<Call->getNumArgs(); i++) {

    Expr * E = Call->getArg(i);
    const ParmVarDecl * pv = D->getParamDecl(i);

    argument_info ai;
    sitedep = handle_call_argument(E, pv, sitedep, &out_variables, &dep_variables, ai);

    cinfo.is_site_dependent |= ai.is_site_dependent;
    cinfo.arguments.push_back(ai);

  } // end of arg loop.

  // now the arguments have been handled
  // If the function is a method, check the member call arg too 

  if ( CXXMemberCallExpr * MCE = dyn_cast<CXXMemberCallExpr>(Call) ) {
    Expr * E = MCE->getImplicitObjectArgument();
    E = E->IgnoreParens();
    E = E->IgnoreImplicit();

    CXXMethodDecl * MD = MCE->getMethodDecl();
    bool is_const = MD->isConst();

    // try this method...
    SourceLocation sl = MD->getNameInfo().getEndLoc();
    // scan parens after name
    bool output_only = false;
    // llvm::errs() << "METHOD WORD AFTER PARENS " << getNextWord(skipParens(sl)) 
    //              << " is const? " << is_const << '\n';
    if (getNextWord(skipParens(sl)) == output_only_keyword) {
      output_only = true;
    }

    if (output_only && is_const) {
       reportDiag(DiagnosticsEngine::Level::Error,
                 sl,
                 "'output_only' cannot be used with 'const'");
       reportDiag(DiagnosticsEngine::Level::Note,
                 Call->getSourceRange().getBegin(),
                 "Called from here");
    }

    cinfo.is_method = true;
    cinfo.object.E = E;

    if (is_field_with_X_expr(E)) {
      handle_field_X_expr(E,!is_const,(!is_const && !output_only),true);
      sitedep = true;

      cinfo.object.is_site_dependent = true;
      cinfo.object.is_lvalue = !is_const;

    } else if (isa<DeclRefExpr>(E)) {
      // some other variable reference
      DeclRefExpr * DRE = dyn_cast<DeclRefExpr>(E);
      var_info * vip = handle_var_ref(DRE,!is_const,"method",nullptr);

      if (vip != nullptr) {
        // site dep is additive
        sitedep = (sitedep || vip->is_site_dependent);

        if (!is_const) out_variables.push_back(vip);

        cinfo.object.is_lvalue = !is_const;

        if (is_const) {
          cinfo.object.is_site_dependent = vip->is_site_dependent;
        } else {
          cinfo.object.is_site_dependent = sitedep;
        }

      }
    } else {
      // now some other expression -- is it site dependent

      cinfo.object.is_lvalue = false;

      sitedep |= is_site_dependent(E,&dep_variables);
    }


  }

  sitedep = attach_dependent_vars( out_variables, sitedep, dep_variables );
  cinfo.is_site_dependent = sitedep;

  return cinfo;
}


///////////////////////////////////////////////////////////////////////////////////
/// single argument in function or constructor call
/// out_variables gives a list of variables which are lvalue refs (output),
/// and dep_variables a list of vars whose site dependence may affect the function
/// and out variables
///////////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::handle_call_argument( Expr *E, const ParmVarDecl * pv, bool sitedep, 
                                         std::vector<var_info *> * out_variables, 
                                         std::vector<var_info *> * dep_variables,
                                         argument_info & ai) {

  // keep filling argument info
  ai.E = E;

  // check parameter type

  // check if we have output_only qualifier
  bool out_only = false;

  if (pv != nullptr) {
    QualType q = pv->getOriginalType();

    if ( getPreviousWord(pv->getSourceRange().getBegin().getLocWithOffset(-1)) 
          == output_only_keyword) {
      out_only = true;
    }
  }

  bool is_lvalue = E->isLValue();
  ai.is_lvalue = is_lvalue;


  if (!is_lvalue && out_only) {
    reportDiag(DiagnosticsEngine::Level::Error,
                E->getSourceRange().getBegin(),
                "'output_only' can be used only with lvalue reference");
  }

  if (is_lvalue) {

    // "output" vars
    if (is_field_with_X_expr(E)) {
      // Mark it as changed
      // llvm::errs() << "FIELD CAN CHANGE HERE!\n";
      sitedep = true;
      handle_field_X_expr(E, is_lvalue, (is_lvalue && !out_only), true, true);

      ai.is_site_dependent = true;

    } else {

      // Now it must be a non-field var reference
      if (DeclRefExpr * DRE = dyn_cast<DeclRefExpr>(E)) {
        var_info * vip = handle_var_ref(DRE,true,"ref",nullptr);

        if (vip == nullptr) return false;  // special var, nothing to do

        // site dep is additive
        sitedep = (sitedep || vip->is_site_dependent);

        if (out_variables != nullptr) {
          out_variables->push_back(vip);
        }

        ai.is_site_dependent |= sitedep;
      }
    }
  } else {

    // now arg is not lvalue - it is input argument
    // sitedep needs to be checked for argument

    ai.is_site_dependent = is_site_dependent(E, &ai.dependent_vars );

    sitedep |= ai.is_site_dependent;

    if (!sitedep && dep_variables != nullptr) {
      // append ai.dependent_vars to dep_variables
      dep_variables->insert(dep_variables->end(), ai.dependent_vars.begin(), 
                                                  ai.dependent_vars.end() );
    }

  }

  return sitedep;

}


////////////////////////////////////////////////////////////////////////////////
/// add the site dependency chain to variables
/// here variables contain pointers to var_info for some "real" variables,
/// and dep_variables a list of var_infos for variables which may bring site dep
/// to all variables in "variables"
/// if "sitedep" is true, mark all unconditionally site dep.
////////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::attach_dependent_vars( std::vector<var_info *> & variables, 
                                          bool sitedep,
                                          std::vector<var_info *> & dep_variables ) {

  // if sitedep == true or one of the dep_variables is sitedep, we can mark
  // all vars in "variables" as site dep.
  if (!sitedep) {
    for (auto d : dep_variables) {
      if (d->is_site_dependent) {
        sitedep = true;
        break;
      }
    }
  }
  if (sitedep) {
    for (auto v : variables) 
      v->is_site_dependent = true;

  } else {
    for (auto v : variables) {

      // if not site dep. seen so far, variables appearing in args may bring it later.
      // add dep var, if not there already
      for (auto d : dep_variables) {
        bool found = false;
        for (auto vdep : v->dependent_vars ) {
          if (d == vdep) {
            found = true;
            break;
          }
        }
        if (!found) v->dependent_vars.push_back(d);
      }
    }
  }
  return sitedep;
}


/////////////////////////////////////////////////////////////////////////////////
/// Special functions:  methods X.method(), and hila_random()
/////////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::handle_special_loop_function(CallExpr *Call) {
  // If the function is in a list of defined loop functions, add it to a list
  // Return true if the expression is a special function

  std::string name = Call->getDirectCallee()->getNameInfo().getAsString();

  // Check if one of the arguments is 'X' (X_index_type, X_plus_direction, X_plus_offset)
  // in that case do nothing, this is handled in code generation
  for (int i=0; i<Call->getNumArgs(); i++) {
    if (is_X_index_type(Call->getArg(0))) {
      return true;
    }
  }

  // check here if this is a special method call (X.method() )
  if (CXXMemberCallExpr *MCall = dyn_cast<CXXMemberCallExpr>(Call)) {
    // llvm::errs() << "It's a member call, name " << name << " objarg "
    //       << MCall->getImplicitObjectArgument()->getType().getAsString() << "\n";
    //    std::string objtype = MCall->getImplicitObjectArgument()->getType().getAsString();
    std::string objtype = get_expr_type(MCall->getImplicitObjectArgument());  
    if (objtype == "X_index_type" || objtype == "lattice_struct *") {
      // now it is a method of X
      // llvm::errs() << " X-method name " << get_stmt_str(Call) << '\n';

      llvm::errs() << "CALL: " << get_stmt_str(Call) << '\n';
      special_function_call sfc;
      sfc.fullExpr = Call;
      sfc.scope = parsing_state.scope_level;
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
      if (name == "coordinates") {
        sfc.replace_expression = "loop_lattice->coordinates(";
        sfc.add_loop_var = true;
      } else if (name == "parity") {
        sfc.replace_expression = "loop_lattice->site_parity(";
        sfc.add_loop_var = true;
      } else if (name == "coordinate") {
        sfc.replace_expression = "loop_lattice->coordinate(";
        sfc.argsExpr = MCall->getArg(0);
        sfc.add_loop_var = true;
      } else if (name == "random") {
        sfc.replace_expression = "hila_random(";
        sfc.add_loop_var = false;
      } else if (name == "size") {
        sfc.replace_expression = "loop_lattice_size(";
        sfc.add_loop_var = false;
        replace_this = target.CUDA;
      } else if (name == "volume") {
        sfc.replace_expression = "loop_lattice_volume(";
        sfc.add_loop_var = false;
        replace_this = target.CUDA;
      } else {
        if (objtype == "const class X_index_type") {
          reportDiag(DiagnosticsEngine::Level::Error,
            Call->getSourceRange().getBegin(),
          "Unknown method X.%0()", name.c_str() );
        } else {
          reportDiag(DiagnosticsEngine::Level::Error,
            Call->getSourceRange().getBegin(),
          "Method 'lattice->.%0()' not allowed inside site loops", name.c_str() );
        }
      }

      if (replace_this) special_function_call_list.push_back(sfc);
      return true;

    } else {
      
      // other method calls?
      return false;
    }

  } else {
    if( name == "hila_random" ){
      llvm::errs() << get_stmt_str(Call) << '\n';
      special_function_call sfc;
      sfc.fullExpr = Call;
      sfc.argsExpr = nullptr;
      sfc.scope = parsing_state.scope_level;
      sfc.name = name;
      sfc.replace_expression = "hila_random()";
      sfc.replace_range = Call->getSourceRange();  // replace full range
      sfc.add_loop_var = false;
      special_function_call_list.push_back(sfc);
      return true;
    }
  }
  return false;
}



///////////////////////////////////////////////////////////////////////////////
/// Start running through loop functions again; call a special visitor
///////////////////////////////////////////////////////////////////////////////

void MyASTVisitor::process_loop_functions() {

  // spin off to a new visitor
  // visit_loop_functions( loop_function_calls );
}



///////////////////////////////////////////////////////////////////////////////
/// Utility for checking if need to handle decls and do it
///////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::handle_loop_function_if_needed(FunctionDecl *fd) {
  // Check if it is in a system header. If so, skip
  SourceManager &SM = Context->getSourceManager();
  bool handle_decl = !SM.isInSystemHeader(fd->getBeginLoc());

  // check if we already have this declaration - either the pointer is the same
  // or the source location (actually, source location should do all, no need for 
  // FunctionDecl *, but it does not hurt)
  for (int i=0; handle_decl && i<loop_functions.size(); i++) { 
    if (fd == loop_functions[i] || 
        fd->getSourceRange().getBegin() == loop_functions[i]->getSourceRange().getBegin() )
      handle_decl = false;
  }
  if (handle_decl) {
    loop_functions.push_back(fd);
    // llvm::errs() << "NEW LOOP FUNCTION " << fd->getNameAsString() << 
    //   " parameters ";
    // for (int i=0; i<fd->getNumParams(); i++) 
    //   llvm::errs() << fd->getParamDecl(i)->getOriginalType().getAsString() << '\n';
    
    backend_handle_loop_function(fd);
  }
  return handle_decl;
}



////////////////////////////////////////////////////////////////////
/// Check if the function is allowed to be within site loops.
/// Returns true if OK to be included; false (and flags error) if not
////////////////////////////////////////////////////////////////////

bool MyASTVisitor::loop_function_check(Decl *d) {
  assert(d != nullptr);
  
  FunctionDecl *fd = dyn_cast<FunctionDecl>(d);
  if (fd) {
    // fd may point to declaration (prototype) without a body.
    // First handle it here in either case

    bool is_new_func = handle_loop_function_if_needed(fd);

    // Now find the declaration of the function body
    // Argument of hasBody becomes the pointer to definition if it is in this compilation unit
    // needs to be const FunctionDecl *
    const FunctionDecl * cfd;
    if (fd->hasBody(cfd)) {
  
      // take away const
      FunctionDecl *fbd = const_cast<FunctionDecl *>(cfd);
      
      if (fbd != fd) {
        is_new_func = handle_loop_function_if_needed(fbd);
      }

      // Now is_new_func is true if the function body has not been scanned before
      if (is_new_func) {
        // get params of the function
        
        // And check also functions called by this func
        CallGraph CG;
        // addToCallGraph takes Decl *: cast 
        // llvm::errs() << " ++ callgraph for " << fbd->getNameAsString() << '\n';

        CG.addToCallGraph( dyn_cast<Decl>(fbd) );
        // CG.dump();
        int i = 0;
        for (auto iter = CG.begin(); iter != CG.end(); ++iter, ++i) {
          // loop through the nodes - iter is of type map<Decl *, CallGraphNode *>
          // root i==0 is "null function", skip
          if (i > 0) {
            Decl * nd = iter->second->getDecl();
            assert(nd != nullptr);
            if (nd != fd) {
              loop_function_check(nd);
            }
          }
          // llvm::errs() << "   ++ loop_function loop " << i << '\n';
        }
      }
      return true;
    } else {
      // Now function has no body - could be in other compilation unit or in system library.
      // TODO: should we handle these?
      // llvm::errs() << "   Function has no body!\n";
    }
  } else {
    // now not a function - should not happen
  }
  return false;
}



