#include "myastvisitor.h"
#include "transformer.h"
#include "clang/Analysis/CallGraph.h"
#include <sstream>
#include <iostream>
#include <string>


std::vector<FunctionDecl *> loop_functions = {};


// Go through each parameter of function calls and handle
// any field references.
// Assume non-const references can be assigned to.

void MyASTVisitor::handle_function_call_in_loop(Stmt * s, bool is_assignment, bool is_compound) {

  // Get the call expression
  CallExpr *Call = dyn_cast<CallExpr>(s);

  assert(Call && "Loop function call not valid");

  // Handle special loop functions
  if( handle_special_loop_function(Call) ){
    return;
  }

  // Get the declaration of the function
  Decl* decl = Call->getCalleeDecl();

  //llvm::errs() << " callee:\n";
  //decl->dump();

  FunctionDecl* D = (FunctionDecl*) llvm::dyn_cast<FunctionDecl>(decl);

  // Store functions used in loops, recursively...
  loop_function_check(decl);

  // Handle parameters
  int i=0;
  for( Expr * E : Call->arguments() ){
    if( is_field_parity_expr(E) ) {
      if(i < D->getNumParams()){
        const ParmVarDecl * pv = D->getParamDecl(i);
        QualType q = pv->getOriginalType ();

        // Check for const qualifier
        if( !q.isConstQualified ()) {
          // Mark it as changed
          handle_field_parity_expr(E, is_assignment, is_compound);
        }
      }
    }
    i++;
  }
}

void MyASTVisitor::handle_function_call_in_loop(Stmt * s) {
  handle_function_call_in_loop(s,true,true);
}



void MyASTVisitor::handle_member_call_in_loop(Stmt * s) {

  // Get the call expression
  CXXMemberCallExpr *Call = dyn_cast<CXXMemberCallExpr>(s);

  assert(Call && "Loop method call not valid");

  // Get the declaration of the function
  Decl* decl = Call->getCalleeDecl();
  FunctionDecl* D = (FunctionDecl*) llvm::dyn_cast<FunctionDecl>(decl);

  // Store functions used in loops, recursively...
  loop_function_check(decl);
  
  // Handle parameters
  int i=0;
  for( Expr * E : Call->arguments() ){
    if( is_field_parity_expr(E) ) {
      const ParmVarDecl * pv = D->getParamDecl(i);
      QualType q = pv->getOriginalType ();
        
      // Check for const qualifier
      if( !q.isConstQualified ()) {
        // Mark it as changed
        handle_field_parity_expr(E, true, true);
      }
    }
    i++;
  }

  // Handle the object itself
  Expr * E = Call->getImplicitObjectArgument();
  if( is_field_parity_expr(E) ) {
    QualType q = E->getType();

    // Check for const qualifier
    if( !q.isConstQualified ()) {
      // Mark it as changed
      handle_field_parity_expr(E, true, true);
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

  // Store functions used in loops, recursively...
  loop_function_check(decl);
}



bool MyASTVisitor::handle_special_loop_function(CallExpr *Call) {
  // If the function is in a list of defined loop functions, add it to a list
  // Return true if the expression is a special function
  std::string name = Call->getDirectCallee()->getNameInfo().getAsString();
  if( name == "coordinates" ){
    llvm::errs() << get_stmt_str(Call) << '\n';
    special_function_call sfc;
    sfc.fullExpr = Call;
    sfc.scope = parsing_state.scope_level;
    sfc.name = name;
    sfc.replace_expression = "loop_lattice->coordinates";
    sfc.add_loop_var = true;
    special_function_call_list.push_back(sfc);
    return true;
  }
  if( name == "hila_random" ){
    llvm::errs() << get_stmt_str(Call) << '\n';
    special_function_call sfc;
    sfc.fullExpr = Call;
    sfc.scope = parsing_state.scope_level;
    sfc.name = name;
    sfc.replace_expression = "hila_random";
    sfc.add_loop_var = false;
    special_function_call_list.push_back(sfc);
    return true;
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
        fd->getSourceRange().getBegin() == loop_functions[i]->getSourceRange().getBegin() )
      handle_decl = false;
  }
  if (handle_decl) {
    loop_functions.push_back(fd);
    backend_handle_loop_function(fd);
  }
  return handle_decl;
}



////////////////////////////////////////////////////////////////////
/// Check if the function is allowed to be within field loops.
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
