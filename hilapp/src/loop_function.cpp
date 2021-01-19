#include "myastvisitor.h"
#include "hilapp.h"
#include "clang/Analysis/CallGraph.h"
#include <sstream>
#include <iostream>
#include <string>


std::vector<FunctionDecl *> loop_functions = {};


// Go through each parameter of function calls and handle
// any field references.
// Assume non-const references can be assigned to.

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

  // llvm::errs() << " callee:\n";
  // decl->dumpColor();

  FunctionDecl* D = (FunctionDecl*) llvm::dyn_cast<FunctionDecl>(decl);

  // Store functions used in loops, recursively...
  loop_function_check(decl);

  // Now inspect function call parameters/args

  // don't inspect args for operator calls, assume it is done separately
  if (isa<CXXOperatorCallExpr>(Call)) return;

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

  // Handle parameters
  if (D->getNumParams() != Call->getNumArgs()) {
    llvm::errs() << "Internal error: #params != #args, function " << D->getNameAsString() << '\n';
    exit(1);
  }
#endif
  

  for( int i=0; i<Call->getNumArgs(); i++) {
    Expr * E = Call->getArg(i);
    
    const ParmVarDecl * pv = D->getParamDecl(i);
    QualType q = pv->getOriginalType ();


    // check if we have output_only qualifier
    bool output_only = false;
 
    if ( getPreviousWord(pv->getSourceRange().getBegin().getLocWithOffset(-1)) 
         == output_only_keyword) {
      output_only = true;
    }
    
    bool is_lvalue = E->isLValue();

    // if (output_only) llvm::errs() << "It is an out parameter!\n";
    // if (is_lvalue) llvm::errs() << " LVALUE\n";

    if (!is_lvalue && output_only) 
      reportDiag(DiagnosticsEngine::Level::Error,
                 Call->getSourceRange().getBegin(),
                 "'output_only' can be used only with lvalue reference");
     
    if (is_field_with_X_expr(E)) {
      // Mark it as changed
      // llvm::errs() << "FIELD CAN CHANGE HERE!\n";
      handle_field_X_expr(E, is_lvalue, (is_lvalue && !output_only), true, true);
    }
  }

  // If the function is a method, check the member call arg too 

  if ( CXXMemberCallExpr * MCE = dyn_cast<CXXMemberCallExpr>(s) ) {
    Expr * E = MCE->getImplicitObjectArgument();
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

    if (is_field_with_X_expr(E)) {
      handle_field_X_expr(E,!is_const,(!is_const && !output_only),true);
    }
  }

}


// this one not used now...

void MyASTVisitor::handle_member_call_in_loop(Stmt * s) {

  // Get the call expression
  CXXMemberCallExpr *Call = dyn_cast<CXXMemberCallExpr>(s);

  assert(Call && "Loop method call not valid");

  // Get the declaration of the function
  Decl* decl = Call->getCalleeDecl();
  FunctionDecl* D = (FunctionDecl*) llvm::dyn_cast<FunctionDecl>(decl);

  llvm::errs() << "  Member call " << D->getNameAsString() << '\n';

  // Store functions used in loops, recursively...
  loop_function_check(decl);
  
  // Handle parameters
  int i=0;
  for( Expr * E : Call->arguments() ){
    if( is_field_with_X_expr(E) ) {
      const ParmVarDecl * pv = D->getParamDecl(i);
      QualType q = pv->getOriginalType ();
        
      // Check for const qualifier
      if( !q.isConstQualified ()) {
        // Mark it as changed
        handle_field_X_expr(E, true, true, true);
      }
    }
    i++;
  }

  // Handle the object itself
  Expr * E = Call->getImplicitObjectArgument();
  E = E->IgnoreImplicit();
  if( is_field_with_X_expr(E) ) {
    CXXMethodDecl * Decl = Call->getMethodDecl();

    // Check for const qualifier
    if( !Decl->isConst() ) {
      // Mark it as changed
      handle_field_X_expr(E, true, true, true);
    }
  }
}



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

  static std::vector<CXXConstructorDecl *> loop_constructors;
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


bool MyASTVisitor::handle_special_loop_function(CallExpr *Call) {
  // If the function is in a list of defined loop functions, add it to a list
  // Return true if the expression is a special function

  std::string name = Call->getDirectCallee()->getNameInfo().getAsString();

  // check here if this is a special method call (X.method() )
  if (CXXMemberCallExpr *MCall = dyn_cast<CXXMemberCallExpr>(Call)) {
    // llvm::errs() << "It's a member call, name " << name << " objarg "
    //       << MCall->getImplicitObjectArgument()->getType().getAsString() << "\n";
    std::string objtype = MCall->getImplicitObjectArgument()->getType().getAsString();
    if (objtype == "const class X_index_type" || objtype == "const class lattice_struct *") {
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
        fd->getSourceRange().getBegin() == 
          loop_functions[i]->getSourceRange().getBegin() ) {
      handle_decl = false;
    }
  }
  if (handle_decl) {
    if (fd->isTrivial()) {
      llvm::errs() << "TRIVIAL LOOP FUNCTION " << fd->getNameAsString() << " NO HANDLING\n";
    } else {
      loop_functions.push_back(fd);
          llvm::errs() << "NEW LOOP FUNCTION " << fd->getNameAsString() << 
                " parameters ";
          for (int i=0; i<fd->getNumParams(); i++) 
             llvm::errs() << fd->getParamDecl(i)->getOriginalType().getAsString();
             llvm::errs() << '\n';

          llvm::errs() << "decl: " << getRangeText(SM,fd->getSourceRange().getBegin(), 
                                findChar(fd->getSourceRange().getBegin(),'\n')) << '\n';
    
      backend_handle_loop_function(fd);
    }
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
