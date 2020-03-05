#include "myastvisitor.h"
#include "transformer.h"
#include "stringops.h"
#include "specialization_db.h"
#include "clang/Analysis/CallGraph.h"
#include <sstream>
#include <iostream>
#include <string>

/////////
/// Implementation of most myastvisitor methods
/////////

//function used for development
std::string print_TemplatedKind(const enum FunctionDecl::TemplatedKind kind) {
  switch (kind) {
    case FunctionDecl::TemplatedKind::TK_NonTemplate:  
      return "TK_NonTemplate";
    case FunctionDecl::TemplatedKind::TK_FunctionTemplate:
      return "TK_FunctionTemplate";
    case FunctionDecl::TemplatedKind::TK_MemberSpecialization:
      return "TK_MemberSpecialization";
    case FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization:  
      return "TK_FunctionTemplateSpecialization";
    case FunctionDecl::TemplatedKind::TK_DependentFunctionTemplateSpecialization:  
      return "TK_DependentFunctionTemplateSpecialization";
  }
}

/// -- Identifier utility functions --

bool MyASTVisitor::is_field_storage_expr(Expr *E) {
  return( E && E->getType().getAsString().find(field_storage_type) != std::string::npos);
}

bool MyASTVisitor::is_field_expr(Expr *E) {
  return( E && E->getType().getAsString().find(field_type) != std::string::npos);
}

bool MyASTVisitor::is_field_decl(ValueDecl *D) {
  return( D && D->getType().getAsString().find(field_type) != std::string::npos);
}

bool MyASTVisitor::is_duplicate_expr(const Expr * a, const Expr * b) {
  // Use the Profile function in clang, which "fingerprints"
  // statements
  llvm::FoldingSetNodeID IDa, IDb;
  a->Profile(IDa, *Context, true);
  b->Profile(IDb, *Context, true);
  return ( IDa == IDb );
}

// Checks if E is a parity Expr. Catches both parity and parity_plus_direction 
bool MyASTVisitor::is_field_parity_expr(Expr *E) {
  E = E->IgnoreParens();
  CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);
  if (OC &&
      strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
      is_field_expr(OC->getArg(0))) {
    std::string s = get_expr_type(OC->getArg(1)); 
    if (s == "parity" || s == "parity_plus_direction") {
      // llvm::errs() << " <<<Parity type " << get_expr_type(OC->getArg(1)) << '\n';
      return true;
    }
  } else {
    // This is for templated expressions
    // for some reason, expr a[X] "getBase() gives X, getIdx() a...
    if (ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(E)) {
      Expr * lhs = ASE->getLHS()->IgnoreParens();
      
      if (is_field_expr(ASE->getLHS()->IgnoreParens())) {
        // llvm::errs() << " FP: and field\n";
        std::string s = get_expr_type(ASE->getRHS());
        if (s == "parity" || s == "parity_plus_direction") {
          // llvm::errs() << " <<<Parity type " << get_expr_type(ASE->getRHS()) << '\n';
          return true;
        }
      }
    }
  }
  return false;   
}

/// is the stmt pointing now to assignment
bool MyASTVisitor::is_assignment_expr(Stmt * s, std::string * opcodestr, bool &iscompound) {
  if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
    if (OP->isAssignmentOp()) {
      // TODO: there should be some more elegant way to do this
      const char *sp = getOperatorSpelling(OP->getOperator());
      if ((sp[0] == '+' || sp[0] == '-' || sp[0] == '*' || sp[0] == '/')
          && sp[1] == '=') iscompound = true;
      else iscompound = false;
      if (opcodestr)
        *opcodestr = getOperatorSpelling(OP->getOperator());

      // Need to mark/handle the assignment method if necessary
      if( is_function_call_stmt(s) ){
        handle_function_call_in_loop(s, true, iscompound);
      } else if ( is_constructor_stmt(s) ){
        handle_constructor_in_loop(s);
      }

      return true;
    }
  } 

  // TODO: this is for templated expr, I think -- should be removed (STILL USED; WHY)
  if (BinaryOperator *B = dyn_cast<BinaryOperator>(s)) {
    if (B->isAssignmentOp()) {
      iscompound = B->isCompoundAssignmentOp();
      if (opcodestr)
        *opcodestr = B->getOpcodeStr();
      return true;
    }
  }

  return false;
}

// is the stmt pointing now to a function call
bool MyASTVisitor::is_function_call_stmt(Stmt * s) {
  if (auto *Call = dyn_cast<CallExpr>(s)){
    llvm::errs() << "Function call found: " << get_stmt_str(s) << '\n';
    return true;
  }
  return false;
}
bool MyASTVisitor::is_constructor_stmt(Stmt * s) {
  if (auto *Call = dyn_cast<CXXConstructExpr>(s)){
    llvm::errs() << "Constructor found: " << get_stmt_str(s) << '\n';
    return true;
  }
  return false;
}

bool MyASTVisitor::isStmtWithSemi(Stmt * S) {
  SourceLocation l = Lexer::findLocationAfterToken(S->getEndLoc(),
                                                   tok::semi,
                                                   TheRewriter.getSourceMgr(),
                                                   Context->getLangOpts(),
                                                   false);
  if (l.isValid()) {
    //    llvm::errs() << "; found " << get_stmt_str(S) << '\n';
    return true;
  }
  return false;
}



/// Check the validity a variable reference in a loop
bool FieldRefChecker::VisitDeclRefExpr(DeclRefExpr *e) {
  // It must be declared already. Get the declaration and check
  // the variable list. (If it's not in the list, it's not local)
  // llvm::errs() << "LPC variable reference: " <<  get_stmt_str(e) << "\n" ;
  for ( var_info & vi : var_info_list ) if( vi.is_loop_local ){
    // Checks that an expression is does not refer to loop local variables
    if( vi.decl == dyn_cast<VarDecl>(e->getDecl()) ){
      // It is local! Generate a warning
      reportDiag(DiagnosticsEngine::Level::Error,
           e->getSourceRange().getBegin(),
          "Field reference depends on loop-local variable");
      break;
    }
  }
  return true;
}


// Walk the tree recursively 
bool FieldRefChecker::TraverseStmt(Stmt *s) {
  RecursiveASTVisitor<FieldRefChecker>::TraverseStmt(s);
  return true;
}

// Walk the tree recursively 
bool LoopAssignChecker::TraverseStmt(Stmt *s) {
  RecursiveASTVisitor<LoopAssignChecker>::TraverseStmt(s);
  return true;
}

/// Check the validity a variable reference in a loop
bool LoopAssignChecker::VisitDeclRefExpr(DeclRefExpr *e) {
  std::string type = e->getType().getAsString();
  type = remove_all_whitespace(type);
  if(type.rfind("element<",0) != std::string::npos){
    reportDiag(DiagnosticsEngine::Level::Error,
      e->getSourceRange().getBegin(),
      "Cannot assing an element to a non-element type");
  }
  return true;
}

/// Check if an assignment is allowed
void MyASTVisitor::check_allowed_assignment(Stmt * s) {
  if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
    if(OP->getNumArgs() == 2){
      // Walk the right hand side to check for element types. None are allowed.
      std::string type = OP->getArg(0)->getType().getAsString();
      type = remove_all_whitespace(type);
      if(type.rfind("element<",0) == std::string::npos){
        LoopAssignChecker lac(TheRewriter, Context);
        lac.TraverseStmt(OP->getArg(1));
      } else {
        llvm::errs() << " ** Element type : " << type << '\n';
        PrintingPolicy pp(Context->getLangOpts());
        llvm::errs() << " ** Canonical type without keywords: " << OP->getArg(0)->getType().getCanonicalType().getAsString(pp) << '\n';
      }
    }
  }
}




/// -- Handler utility functions -- 

//////////////////////////////////////////////////////////////////////////////
/// Go through one field reference within parity loop and store relevant info
//////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::handle_field_parity_expr(Expr *e, bool is_assign, bool is_compound) {
    
  e = e->IgnoreParens();
  field_ref lfe;
  if (CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(e)) {
    lfe.fullExpr   = OC;
    // take name 
    lfe.nameExpr   = OC->getArg(0);
    lfe.parityExpr = OC->getArg(1);
  } else if (ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(e)) {
    // In template definition TODO: should be removed?

    lfe.fullExpr   = ASE;
    lfe.nameExpr   = ASE->getLHS();
    lfe.parityExpr = ASE->getRHS();
    //llvm::errs() << lfe.fullExpr << " " << lfe.nameExpr << " " << lfe.parityExpr << "\n";
  } else {
    llvm::errs() << "Should not happen! Error in field parity\n";
    exit(1);
  }

  // Check if the expression is already handled
  for( field_ref r : field_ref_list)
    if( r.fullExpr == lfe.fullExpr  ){
      return(true);
  }

  // Check that there are no local variable references up the AST
  FieldRefChecker frc(TheRewriter, Context);
  frc.TraverseStmt(lfe.fullExpr);

  //lfe.nameInd    = writeBuf->markExpr(lfe.nameExpr); 
  //lfe.parityInd  = writeBuf->markExpr(lfe.parityExpr);
  
  lfe.dirExpr  = nullptr;  // no neighb expression
    
  if (get_expr_type(lfe.parityExpr) == "parity") {
    if (parsing_state.accept_field_parity) {
      // 1st parity statement on a single line lattice loop
      loop_parity.expr  = lfe.parityExpr;
      loop_parity.value = get_parity_val(loop_parity.expr);
      loop_parity.text  = get_stmt_str(loop_parity.expr);
    } else {
      require_parity_X(lfe.parityExpr);
    }
  }

  lfe.is_written = is_assign;
  lfe.is_read = (is_compound || !is_assign);
  
  // next ref must have wildcard parity
  parsing_state.accept_field_parity = false;
        
  if (get_expr_type(lfe.parityExpr) == "parity_plus_direction") {

    if (is_assign) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 lfe.parityExpr->getSourceRange().getBegin(),
                 "Neighbour offset not allowed on the LHS of an assignment");
    }

    // Now need to split the expr to parity and dir-bits
    // Need to descent quite deeply into the expr chain
    Expr* e = lfe.parityExpr->IgnoreParens();
    // The next line is necessary with the clang version on Puhti
    e = e->IgnoreImplicit();
    CXXOperatorCallExpr* Op = dyn_cast<CXXOperatorCallExpr>(e);
    // descent into expression
    // TODO: must allow for arbitrary offset!

    if (!Op) {
      CXXConstructExpr * Ce = dyn_cast<CXXConstructExpr>(e);
      if (Ce) {
        // llvm::errs() << " ---- got Ce, args " << Ce->getNumArgs() << '\n';
        if (Ce->getNumArgs() == 1) {
          e = Ce->getArg(0)->IgnoreImplicit();
          Op = dyn_cast<CXXOperatorCallExpr>(e);
        }
      }
    }
    if (!Op) {
      reportDiag(DiagnosticsEngine::Level::Fatal,
                 lfe.parityExpr->getSourceRange().getBegin(),
                 "Internal error: could not decipher parity + dir statement" );
      exit(1);
    }
      
    // if (!Op) {
    //   e = e->IgnoreImplicit();
    //   Op = dyn_cast<CXXOperatorCallExpr>(e);        
    // }

    if (Op &&
        (strcmp(getOperatorSpelling(Op->getOperator()),"+") == 0 ||
         strcmp(getOperatorSpelling(Op->getOperator()),"-") == 0) &&
        get_expr_type(Op->getArg(0)) == "parity") {
        llvm::errs() << " ++++++ found parity + dir\n";

        require_parity_X(Op->getArg(0));
        lfe.dirExpr = Op->getArg(1)->IgnoreImplicit();
        lfe.dirname = get_stmt_str(lfe.dirExpr);

        // If the direction is a variable, add it to the list
        DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(lfe.dirExpr);
        static std::string assignop;
        if(DRE && isa<VarDecl>(DRE->getDecl())) {
          handle_var_ref(DRE, false, assignop);
        }
    }
  }
    
  // llvm::errs() << "field expr " << get_stmt_str(lfe.nameExpr)
  //              << " parity " << get_stmt_str(lfe.parityExpr)
  //              << "\n";

   
  field_ref_list.push_back(lfe);
      
  return(true);
}

// This processes references to non-field variables within field loops
void MyASTVisitor::handle_var_ref(DeclRefExpr *DRE,
                                  bool is_assign,
                                  std::string &assignop) {

  
  if (isa<VarDecl>(DRE->getDecl())) {
    auto decl = dyn_cast<VarDecl>(DRE->getDecl());
    var_ref vr;
    vr.ref = DRE;
    //vr.ind = writeBuf->markExpr(DRE);
    vr.is_assigned = is_assign;
    if (is_assign) vr.assignop = assignop;
    
    bool found = false;
    var_info *vip = nullptr;
    for (var_info & vi : var_info_list) {
      if (vi.decl == decl) {
        // found already referred to decl
        vi.refs.push_back(vr);
        vi.is_assigned |= is_assign;
        vi.reduction_type = get_reduction_type(is_assign, assignop, vi);
        
        vip = &vi;
        found = true;
        break;
      }
    }
    if (!found) {
      // new variable referred to
      var_info vi;
      vi.refs = {};
      vi.refs.push_back(vr);
      vi.decl = decl;
      vi.name = decl->getName();
      // Printing policy is somehow needed for printing type without "class" id
      // Unqualified takes away "consts" etc and Canonical typdefs/using.
      // Also need special handling for element type
      PrintingPolicy pp(Context->getLangOpts());
      vi.type = DRE->getType().getUnqualifiedType().getAsString(pp);
      vi.type = remove_all_whitespace(vi.type);
      bool is_elem = (vi.type.find("element<") == 0);
      vi.type = DRE->getType().getUnqualifiedType().getCanonicalType().getAsString(pp);
      if (is_elem) vi.type = "element<" + vi.type + ">";
      llvm::errs() << " + Got " << vi.type << '\n';

      // is it loop-local?
      vi.is_loop_local = false;
      for (var_decl & d : var_decl_list) {
        if (d.scope >= 0 && vi.decl == d.decl) {
          llvm::errs() << "loop local var ref! " << vi.name << '\n';
          vi.is_loop_local = true;
          vi.var_declp = &d;
          break;
        }
      }
      vi.is_assigned = is_assign;
      // we know refs contains only 1 element
      vi.reduction_type = get_reduction_type(is_assign, assignop, vi);
      
      var_info_list.push_back(vi);
      vip = &(var_info_list.back());
    }
  } else { 
    // end of VarDecl - how about other decls, e.g. functions?
    reportDiag(DiagnosticsEngine::Level::Error,
               DRE->getSourceRange().getBegin(),
               "Reference to unimplemented (non-variable) type");
  }
}


// handle an array subscript expression
void MyASTVisitor::handle_array_var_ref(ArraySubscriptExpr *E) {
  PrintingPolicy pp(Context->getLangOpts());

  // Check if it's local
  bool is_loop_local = false;
  Expr * RE = E->getBase()->IgnoreImplicit();
  while(!dyn_cast<DeclRefExpr>(RE)){
    if(dyn_cast<ArraySubscriptExpr>(RE)){
      ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(RE);
      RE = ASE->getBase()->IgnoreImplicit();
    } else if(dyn_cast<MemberExpr>(RE)){
      MemberExpr * ME = dyn_cast<MemberExpr>(RE);
      RE = ME->getBase()->IgnoreImplicit();
    } else if(dyn_cast<CXXOperatorCallExpr>(RE)) {
      CXXOperatorCallExpr * OCE = dyn_cast<CXXOperatorCallExpr>(RE);
      if(strcmp(getOperatorSpelling(OCE->getOperator()),"[]") == 0){
        RE = OCE->getArg(0)->IgnoreImplicit();
      } else {
        // It's not a variable
        return;
      }
    } else {
      // It's not a variable
      return;
    }
  }
  DeclRefExpr * DRE = dyn_cast<DeclRefExpr>(RE);
  VarDecl * decl = dyn_cast<VarDecl>(DRE->getDecl());
  for (var_decl & d : var_decl_list) {
    if (d.scope >= 0 && decl == d.decl) {
      llvm::errs() << "loop local var ref! \n";
      is_loop_local = true;
      break;
    }
  }

  // If it's defined outside the loop, it may need to be handled
  // (communicated to kernels, for example). So record it here.
  if(!is_loop_local){
    array_ref ar;
    ar.ref = E;
    
    // Find the type of the full expression (that is an element of the array)
    const Type *t = E->getType().getTypePtr();
    llvm::errs() << E->getType().getCanonicalType().getUnqualifiedType().getAsString(pp) << "\n";
    ar.type = E->getType().getCanonicalType().getUnqualifiedType().getAsString(pp);
  
    array_ref_list.push_back(ar);
  }
}




///////////////////////////////////////////////////////////////////////////////
/// handle_full_loop_stmt() is the starting point for the analysis of all
/// "parity" -loops
///////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::handle_full_loop_stmt(Stmt *ls, bool field_parity_ok ) {
  // init edit buffer
  // Buf.create( &TheRewriter, ls );
          
  field_ref_list.clear();
  special_function_call_list.clear();
  var_info_list.clear();
  var_decl_list.clear();
  array_ref_list.clear();
  remove_expr_list.clear();
  global.location.loop = ls->getSourceRange().getBegin();
  
  parsing_state.accept_field_parity = field_parity_ok;
    
  // the following is for taking the parity from next elem
  parsing_state.scope_level = 0;
  parsing_state.in_loop_body = true;
  TraverseStmt(ls);
  parsing_state.in_loop_body = false;

  // Remove exprs which we do not want
  for (Expr * e : remove_expr_list) writeBuf->remove(e);
  
  // check and analyze the field expressions
  check_field_ref_list();
  check_var_info_list();
  // check that loop_parity is not X
  if (loop_parity.value == parity::x) {
    reportDiag(DiagnosticsEngine::Level::Error,
               loop_parity.expr->getSourceRange().getBegin(),
               "Parity of the full loop cannot be \'X\'");
  }
  
  generate_code(ls, target);
  
  // Buf.clear();
          
  // Emit the original command as a commented line
  writeBuf->insert(ls->getSourceRange().getBegin(),
                   comment_string(global.full_loop_text) + "\n",true,true);
  
  global.full_loop_text = "";

  // don't go again through the arguments
  parsing_state.skip_children = 1;

  state::loop_found = true;
  
  return true;
}


////////////////////////////////////////////////////////////////////////////////
///  act on statements within the parity loops.  This is called 
///  from VisitStmt() if the status state::in_loop_body is true
////////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::handle_loop_body_stmt(Stmt * s) {

  // This keeps track of the assignment to field
  // must remember the set value across calls
  static bool is_assignment = false;
  static bool is_compound = false;
  static std::string assignop;
 
  // Need to recognize assignments lf[X] =  or lf[X] += etc.
  // And also assignments to other vars: t += norm2(lf[X]) etc.
  if (is_assignment_expr(s,&assignop,is_compound)) {
    check_allowed_assignment(s);
    is_assignment = true;
    // next visit here will be to the assigned to variable
    return true;
  }

  // Check for function calls parameters. We need to determine if the 
  // function can assign to the a field parameter (is not const).
  if( is_function_call_stmt(s) ){
    handle_function_call_in_loop(s);
    // let this ripple trough, for now ...
    // return true;
  }

  if ( is_constructor_stmt(s) ){
    handle_constructor_in_loop(s);
    // return true;
  }
   
  // catch then expressions
      
  if (Expr *E = dyn_cast<Expr>(s)) {
    
    // Avoid treating constexprs as variables
     if (E->isCXX11ConstantExpr(*Context, nullptr, nullptr)) {
       parsing_state.skip_children = 1;   // nothing to be done
       return true;
     }
    
    //if (is_field_element_expr(E)) {
      // run this expr type up until we find field variable refs
    if (is_field_parity_expr(E)) {
      // Now we know it is a field parity reference
      // get the expression for field name
          
      handle_field_parity_expr(E, is_assignment, is_compound);
      is_assignment = false;  // next will not be assignment
      // (unless it is a[] = b[] = c[], which is OK)

      parsing_state.skip_children = 1;
      return true;
    }

    if (is_field_expr(E)) {
      // field without [parity], bad usually (TODO: allow  scalar func(field)-type?)
      reportDiag(DiagnosticsEngine::Level::Error,
                 E->getSourceRange().getBegin(),
                 "Field expressions without [..] not allowed within field loop");
      parsing_state.skip_children = 1;  // once is enough
      return true;
    }

    if (DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(E)) {
      if (isa<VarDecl>(DRE->getDecl())) {
        // now it should be var ref non-field
      
        handle_var_ref(DRE,is_assignment,assignop);
        is_assignment = false;
      
        parsing_state.skip_children = 1;
        llvm::errs() << "Variable ref: "
                     << TheRewriter.getRewrittenText(E->getSourceRange()) << '\n';

        parsing_state.skip_children = 1;
        return true;
      }
      // TODO: function ref?
    }



#if 1

    if (isa<ArraySubscriptExpr>(E)) {
      llvm::errs() << "  It's array expr "
                   << TheRewriter.getRewrittenText(E->getSourceRange()) << "\n";
      //parsing_state.skip_children = 1;
      auto a = dyn_cast<ArraySubscriptExpr>(E);

      // At this point this should be an allowed expression?
      handle_array_var_ref(a);
      
      return true;
    }
#endif          
        
    if (0){

      // not field type non-const expr
      llvm::errs() << "Non-const other Expr: " << get_stmt_str(E) << '\n';
      // loop-local variable refs inside? If so, we cannot evaluate this as "whole"

      // check_local_loop_var_refs = 1;
        
      // TODO: find really uniq variable references
      //var_ref_list.push_back( handle_var_ref(E) );

      parsing_state.skip_children = 1;          
      return true;
    }
  } // Expr checking branch - now others...

  // This reached only if s is not Expr

 
  // start {...} -block or other compound
  if (isa<CompoundStmt>(s) || isa<ForStmt>(s) || isa<IfStmt>(s)
      || isa<WhileStmt>(s)) {

    static bool passthrough = false;
    // traverse each stmt - use passthrough trick if needed
    if (passthrough) {
      passthrough = false;
      return true;
    }
    
    parsing_state.scope_level++;
    passthrough = true;     // next visit will be to the same node, skip
    TraverseStmt(s);
    parsing_state.scope_level--;
    remove_vars_out_of_scope(parsing_state.scope_level);
    parsing_state.skip_children = 1;
    return true;
  }
    
  return true;    
}

int MyASTVisitor::handle_field_specializations(ClassTemplateDecl *D) {
  // save global, perhaps needed (perhaps not)
  field_decl = D;

  // llvm::errs() << "+++++++\n Specializations of field\n";

  int count = 0;
  for (auto spec = D->spec_begin(); spec != D->spec_end(); spec++ ) {
    count++;
    auto & args = spec->getTemplateArgs();

    if (args.size() != 1) {
      llvm::errs() << " *** Fatal: More than one type arg for field<>\n";
      exit(1);
    }
    if (TemplateArgument::ArgKind::Type != args.get(0).getKind()) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 D->getSourceRange().getBegin(),
                 "Expect type argument in \'field\' template" );
      return(0);
    }

    // Get typename without class, struct... qualifiers
    PrintingPolicy pp(Context->getLangOpts());
    std::string typestr = args.get(0).getAsType().getAsString(pp);
    llvm::errs() << "arg type " << typestr << "\n";

    // Type of field<> can never be field?  This always is true
    if( typestr.find("field<") ){ // Skip for field templates
      if (spec->isExplicitSpecialization()) llvm::errs() << " explicit\n";
    }
    
  }
  return(count);
      
} // end of "field"

/// Source Location utilities

SourceLocation MyASTVisitor::getSourceLocationAtEndOfLine( SourceLocation l ) {
  SourceManager &SM = TheRewriter.getSourceMgr();
  for (int i=0; i<10000; i++) {
    bool invalid = false;
    const char * c = SM.getCharacterData(l.getLocWithOffset(i),&invalid);
    if (invalid) {
      // no new line found in buffer.  return previous loc, could be false!
      llvm::errs() << program_name + ": no new line found in buffer, internal error\n";
      return( l.getLocWithOffset(i-1) );
    }
    if (*c == '\n') return( l.getLocWithOffset(i) );
  }
  return l;
}

SourceLocation MyASTVisitor::getSourceLocationAtEndOfRange( SourceRange r ) {
  int i = TheRewriter.getRangeSize(r);
  return r.getBegin().getLocWithOffset(i-1);
}


/// Find the source range of the previous line from source location
SourceRange MyASTVisitor::getSourceRangeAtPreviousLine( SourceLocation l ){
  SourceManager &SM = TheRewriter.getSourceMgr();
  SourceRange r = SourceRange(l,l);
  bool found_line_break = false;
  
  // Move backward from location l and find two linebreaks
  // These are the beginning and the end of the previous line
  for (int i=0; i<1000; i++) {
    bool invalid = false;
    const char * c = SM.getCharacterData(l,&invalid);
    if (invalid) {
      // Invalid character found. This is probably the beginning of a file
      r.setBegin(l.getLocWithOffset(1));
      return(r);
    }
    if (*c == '\n'){
      if(!found_line_break) {
        // First line break
        r.setEnd(l);
        found_line_break = true;
      } else {
        // Second line break, set as start. Don't include the line break
        r.setBegin(l.getLocWithOffset(1));
        return r;
      }
    }
    l=l.getLocWithOffset(-1);
  }
  return r;
}


// By overriding these methods in MyASTVisitor we can control which nodes are visited. 

bool MyASTVisitor::TraverseStmt(Stmt *S) {

  if (parsing_state.check_loop && state::loop_found) return true;
    
  // if state::skip_children > 0 we'll skip all until return to level up
  if (parsing_state.skip_children > 0) parsing_state.skip_children++;
    
  // go via the original routine...
  if (!parsing_state.skip_children) RecursiveASTVisitor<MyASTVisitor>::TraverseStmt(S);

  if (parsing_state.skip_children > 0) parsing_state.skip_children--;
      
  return true;
}

bool MyASTVisitor::TraverseDecl(Decl *D) {

  if (parsing_state.check_loop && state::loop_found) return true;

  // if state::skip_children > 0 we'll skip all until return to level up
  if (parsing_state.skip_children > 0) parsing_state.skip_children++;
    
  // go via the original routine...
  if (!parsing_state.skip_children) RecursiveASTVisitor<MyASTVisitor>::TraverseDecl(D);

  if (parsing_state.skip_children > 0) parsing_state.skip_children--;

  return true;
}

template <unsigned N>
void GeneralVisitor::reportDiag(DiagnosticsEngine::Level lev, const SourceLocation & SL,
                              const char (&msg)[N],
                              const char *s1,
                              const char *s2,
                              const char *s3 ) {
  // we'll do reporting only when output is on, avoid double reports
  auto & DE = Context->getDiagnostics();    
  auto ID = DE.getCustomDiagID(lev, msg );
  auto DB = DE.Report(SL, ID);
  if (s1 != nullptr) DB.AddString(s1);
  if (s2 != nullptr) DB.AddString(s2);
  if (s3 != nullptr) DB.AddString(s3);
}

parity MyASTVisitor::get_parity_val(const Expr *pExpr) {
  SourceLocation SL;
  APValue APV;

  if (pExpr->isCXX11ConstantExpr( *Context, &APV, &SL )) {
    // Parity is now constant
    int64_t val = (APV.getInt().getExtValue());
    parity p;
    if (0 <= val && val <= (int)parity::x) {
      p = static_cast<parity>(val);
    } else {
      reportDiag(DiagnosticsEngine::Level::Fatal,
                 pExpr->getSourceRange().getBegin(),
                 "Transformer internal error, unknown parity" );
      exit(1);
    }
    if (p == parity::none) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 pExpr->getSourceRange().getBegin(),
                 "parity::none is reserved for internal use" );
    }
        
    return p;
  } else {
    return parity::none;
  }
}

void MyASTVisitor::require_parity_X(Expr * pExpr) {
  // Now parity has to be X (or the same as before?)
  if (get_parity_val(pExpr) != parity::x) {
    reportDiag(DiagnosticsEngine::Level::Error,
               pExpr->getSourceRange().getBegin(),
               "Use wildcard parity \"X\" or \"parity::x\"" );
  }
}

// finish the field_ref_list, and
// construct the field_info_list

bool MyASTVisitor::check_field_ref_list() {

  bool no_errors = true;
  
  global.assert_loop_parity = false;

  field_info_list.clear();
    
  for( field_ref & p : field_ref_list ) {

    p.direction = -1;  // reset the direction
    
    std::string name = get_stmt_str(p.nameExpr);
      
    field_info * lfip = nullptr;

    // search for duplicates: if found, lfip is non-null

    for (field_info & li : field_info_list) {
      if (name.compare(li.old_name) == 0) {
        lfip = &li;
        break;
      }
    }

    if (lfip == nullptr) {
      field_info lfv;
      lfv.old_name = name;
      lfv.type_template = get_expr_type(p.nameExpr);
      if (lfv.type_template.find("field",0) != 0) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   p.nameExpr->getSourceRange().getBegin(),
                   "Confused: type of field expression?");
        no_errors = false;
      }
      lfv.type_template.erase(0,5);  // Remove "field"  from field<T>
      lfv.is_written = p.is_written;
      lfv.is_read    = p.is_read;
      
      field_info_list.push_back(lfv);
      lfip = & field_info_list.back();
    }
    // now lfip points to the right info element
    // copy that to lf reference
    p.info = lfip;

    if (p.is_written) lfip->is_written = true;
    if (p.is_read)    lfip->is_read    = true;
      
    // save expr record
    lfip->ref_list.push_back(&p);

    if (p.dirExpr != nullptr) {

      if (p.is_written) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   p.parityExpr->getSourceRange().getBegin(),
                   "Neighbour offset not allowed on the LHS of an assignment");
        no_errors = false;
      }

      // does this dir with this name exist before?
      unsigned i = 0;
      bool found = false;
      for (dir_ptr & d : lfip->dir_list) {
        if (is_duplicate_expr(d.e, p.dirExpr)) {
          d.count++;
          p.direction = i;
          found = true;
          break;
        }
        i++;
      }
        
      if (!found) {
        dir_ptr dp;
        dp.e = p.dirExpr;
        dp.count = 1;
        p.direction = lfip->dir_list.size();

        lfip->dir_list.push_back(dp);
      }
    } // dirExpr
  } // p-loop
  
  // check for f[ALL] = f[X+dir] -type use, which is undefined
  
  for (field_info & l : field_info_list) {
    if (l.is_written && l.dir_list.size() > 0) {
      if (loop_parity.value == parity::all) {
        // There's error, find culprits
        for (field_ref * p : l.ref_list) {
          if (p->dirExpr != nullptr && !p->is_written) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       p->parityExpr->getSourceRange().getBegin(),
                       "Accessing field '%0' undefined when assigning to '%1' with parity ALL, flagging as error",
                       get_stmt_str(p->fullExpr).c_str(),
                       l.old_name.c_str());
            no_errors = false;
          }
        }

        for (field_ref * p : l.ref_list) {
          if (p->is_written && p->dirExpr == nullptr) {
            reportDiag(DiagnosticsEngine::Level::Note,
                       p->fullExpr->getSourceRange().getBegin(),
                       "Location of assignment");
              
          }
        }
      } else if (loop_parity.value == parity::none) {
        // not sure if there's an error, emit an assertion
        global.assert_loop_parity = true;
        reportDiag(DiagnosticsEngine::Level::Note,
                   l.ref_list.front()->fullExpr->getSourceRange().getBegin(),
                   "Assign to '%0' and access with offset may be undefined with parity '%1', inserting assertion",
                   l.old_name.c_str(),
                   loop_parity.text.c_str());
      }
    }
  }
  return no_errors;
}

/// Check now that the references to variables are according to rules
void MyASTVisitor::check_var_info_list() {
  for (var_info & vi : var_info_list) {
    if (!vi.is_loop_local) {
      if (vi.reduction_type != reduction::NONE) {
        if (vi.refs.size() > 1) {
          // reduction only once
          int i=0;
          for (auto & vr : vi.refs) {
            if (vr.assignop == "+=" || vr.assignop == "*=") {
              reportDiag(DiagnosticsEngine::Level::Error,
                         vr.ref->getSourceRange().getBegin(),
                         "Reduction variable \'%0\' used more than once within one field loop",
                         vi.name.c_str());
              break;
            }
            i++;
          }
          int j=0;
          for (auto & vr : vi.refs) {
            if (j!=i) reportDiag(DiagnosticsEngine::Level::Note,
                                 vr.ref->getSourceRange().getBegin(),
                                 "Other reference to \'%0\'", vi.name.c_str());
            j++;
          }
        }
      } else if (vi.is_assigned) {
        // now not reduction
        for (auto & vr : vi.refs) {
          if (vr.is_assigned) 
            reportDiag(DiagnosticsEngine::Level::Error,
                       vr.ref->getSourceRange().getBegin(),
                       "Cannot assign to variable defined outside field loop (unless reduction \'+=\' or \'*=\')");
        }
      }
    }
  }
}


/// flag_error = true by default in myastvisitor.h
SourceRange MyASTVisitor::getRangeWithSemi(Stmt * S, bool flag_error) {
  SourceRange range(S->getBeginLoc(),
                    Lexer::findLocationAfterToken(S->getEndLoc(),
                                                  tok::semi,
                                                  TheRewriter.getSourceMgr(),
                                                  Context->getLangOpts(),
                                                  false));
  if (!range.isValid()) {
    if (flag_error) {
      reportDiag(DiagnosticsEngine::Level::Fatal,
                 S->getEndLoc(),
                 "Expecting ';' after expression");
    }
    // put a valid value in any case
    range = S->getSourceRange();        
  }
    
  // llvm::errs() << "Range w semi: " << TheRewriter.getRewrittenText(range) << '\n';
  return range;
}

bool MyASTVisitor::control_command(VarDecl *var) {
  std::string n = var->getNameAsString();
  if (n.find("_transformer_ctl_",0) == std::string::npos) return false;
  
  if (n == "_transformer_ctl_dump_ast") {
    parsing_state.dump_ast_next = true;
  } else if(n == "_transformer_ctl_loop_function") {
    parsing_state.loop_function_next = true;
  } else {
    reportDiag(DiagnosticsEngine::Level::Warning,
               var->getSourceRange().getBegin(),
               "Unknown command for transformer_ctl(), ignoring");
  }
  // remove the command
  SourceRange sr = var->getSourceRange();
  if (sr.getBegin().isMacroID()){
    CharSourceRange csr = TheRewriter.getSourceMgr().getImmediateExpansionRange( sr.getBegin() );
    sr = csr.getAsRange();
  }
  writeBuf->remove(sr);
  return true;
}


bool MyASTVisitor::VisitVarDecl(VarDecl *var) {

  // catch the transformer_ctl -commands here
  if (control_command(var)) return true;
  
  if (parsing_state.check_loop && state::loop_found) return true;

  if (parsing_state.dump_ast_next) {
    // llvm::errs() << "**** Dumping declaration:\n" + get_stmt_str(s)+'\n';
    var->dump();
    parsing_state.dump_ast_next = false;
  }

  
  if (parsing_state.in_loop_body) {
    // for now care only loop body variable declarations

    if (!var->hasLocalStorage()) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 var->getSourceRange().getBegin(),
                 "Static or external variable declarations not allowed within field loops");
      return true;
    }

    if (var->isStaticLocal()) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 var->getSourceRange().getBegin(),
                 "Cannot declare static variables inside field loops");
      return true;
    }

    if (is_field_decl(var)) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 var->getSourceRange().getBegin(),
                 "Cannot declare field variables within field loops");
      parsing_state.skip_children = 1;
      return true;
    }

    // Now it should be automatic local variable decl
    var_decl vd;
    vd.decl = var;
    vd.name = var->getName();
    vd.type = var->getType().getAsString();
    vd.scope = parsing_state.scope_level;
    var_decl_list.push_back(vd);
    
    llvm::errs() << "Local var decl " << vd.name << " of type " << vd.type << '\n';
    return true;
  } 

  // if (is_field_decl(var)) {
  //   llvm::errs() << "FIELD DECL \'" << var->getName() << "\' of type "
  //                << var->getType().getAsString() << '\n';
  //   if (var->isTemplated()) llvm::errs() << " .. was templated\n";
  // }
  
  return true;
}


void MyASTVisitor::remove_vars_out_of_scope(unsigned level) {
  while (var_decl_list.size() > 0 && var_decl_list.back().scope > level)
    var_decl_list.pop_back();
}

///////////////////////////////////////////////////////////////////////////////
/// VisitStmt is called for each statement in AST.  Thus, when traversing the
/// AST or part of it we start here, and branch off depending on the statements
/// and parsing_state.lags
///////////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::VisitStmt(Stmt *s) {

  if (parsing_state.check_loop && state::loop_found) return true;
  
  if (parsing_state.dump_ast_next) {
    llvm::errs() << "**** Dumping statement:\n" + get_stmt_str(s)+'\n';
    s->dump();
    parsing_state.dump_ast_next = false;
  }

  // Entry point when inside field[par] = .... body
  if (parsing_state.in_loop_body) {
    return handle_loop_body_stmt(s);
  }
    
  // loop of type "onsites(p)"
  // Defined as a macro, needs special macro handling
  if (isa<ForStmt>(s)) {

    ForStmt *f = cast<ForStmt>(s);
    SourceLocation startloc = f->getSourceRange().getBegin();

    if (startloc.isMacroID()) {
      Preprocessor &pp = myCompilerInstance->getPreprocessor();
      static std::string loop_call("onsites");
      if (pp.getImmediateMacroName(startloc) == loop_call) {
        // Now we know it is onsites-macro

        if (parsing_state.check_loop) {
          state::loop_found = true;
          return true;
        }

        
        
        CharSourceRange CSR = TheRewriter.getSourceMgr().getImmediateExpansionRange( startloc );
        std::string macro = TheRewriter.getRewrittenText( CSR.getAsRange() );
        bool internal_error = true;

        // llvm::errs() << "macro str " << macro << '\n';
        
        DeclStmt * init = dyn_cast<DeclStmt>(f->getInit());
        if (init && init->isSingleDecl() ) {
          VarDecl * vd = dyn_cast<VarDecl>(init->getSingleDecl());
          if (vd) {
            const Expr * ie = vd->getInit();
            if (ie) {
              loop_parity.expr  = ie;
              loop_parity.value = get_parity_val(loop_parity.expr);
              loop_parity.text  = remove_initial_whitespace(macro.substr(loop_call.length(),
                                                                         std::string::npos));
                
              global.full_loop_text = macro + " " + get_stmt_str(f->getBody());

              // Delete "onsites()" -text

              // TheRewriter.RemoveText(CSR);
              writeBuf->remove(CSR);
              
              handle_full_loop_stmt(f->getBody(), false);
              internal_error = false;
            }
          }
        }
        if (internal_error) {
          reportDiag(DiagnosticsEngine::Level::Error,
                     f->getSourceRange().getBegin(),
                     "\'onsites\'-macro: not a parity type argument" );
          return true;
        }
      }
    }        
    return true;      
  }

                                   
  //  Starting point for fundamental operation
  //  field[par] = ....  version with field<class>
  
  // isStmtWithSemi(s);

  CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s);
  
  if (OP && OP->isAssignmentOp() && is_field_parity_expr(OP->getArg(0))) {
    // now we have a[par] += ...  -stmt.  Arg(0) is
    // the lhs of the assignment
    
    if (parsing_state.check_loop) {
      state::loop_found = true;
      return true;
    }

    SourceRange full_range = getRangeWithSemi(OP,false);
    global.full_loop_text = TheRewriter.getRewrittenText(full_range);
        
    handle_full_loop_stmt(OP, true);
    return true;
  } 
  // now the above when type is field<double> or some other non-class element
  
  BinaryOperator *BO = dyn_cast<BinaryOperator>(s);
  if (BO && BO->isAssignmentOp() && is_field_parity_expr(BO->getLHS())) {
    if (parsing_state.check_loop) {
      state::loop_found = true;
      return true;
    }
    
    SourceRange full_range = getRangeWithSemi(BO,false);
    global.full_loop_text = TheRewriter.getRewrittenText(full_range);        
  
    handle_full_loop_stmt(BO, true);
    return true;    
  }

  return true;
}

//////// Functiondecl and templates below

bool MyASTVisitor::does_function_contain_loop( FunctionDecl *f ) {
  // Currently simple: buffer the function and traverse through it

  srcBuf buf(&TheRewriter,f);
  srcBuf *bp = writeBuf;
  writeBuf = &buf;
  
  buf.off();
  
  bool lf = state::loop_found;
  state::loop_found = false;  // use this to flag

  bool retval;
  if (f->hasBody()) {
    
    // llvm::errs() << "About to check function " << f->getNameAsString() << '\n';
    // llvm::errs() << buf.dump() << '\n';
    
    parsing_state.check_loop = true;
    TraverseStmt(f->getBody());
    parsing_state.check_loop = false;
    
    // llvm::errs() << "Func check done\n";
    
    retval = state::loop_found;
  } else {
    retval = false;
  }
  state::loop_found = lf;
  writeBuf = bp;
  
  buf.clear();
  
  return retval;
}


bool MyASTVisitor::has_pragma(Decl *f, const char * s) {
  SourceRange sr = getSourceRangeAtPreviousLine(f->getSourceRange().getBegin());
  std::string line = TheRewriter.getRewrittenText( sr );

  const std::vector<std::string> a{"#pragma", "transformer", s};
  if (contains_word_list(line,a)) {

    // got it, comment out -- check that it has not been commented out before
    int loc = writeBuf->find_original(sr.getBegin(),'#');
    std::string s = writeBuf->get(loc,loc+1);
    if (s.at(0) == '#') writeBuf->insert(loc ,"//-- ",true,false);
    return true;
  }

  return false;
}


bool MyASTVisitor::VisitFunctionDecl(FunctionDecl *f) {
  // Only function definitions (with bodies), not declarations.
  // also only non-templated functions
  // this does not really do anything

  if (parsing_state.dump_ast_next) {
    llvm::errs() << "**** Dumping funcdecl:\n";
    f->dump();
    parsing_state.dump_ast_next = false;
  }
  if(!parsing_state.check_loop && (parsing_state.loop_function_next || has_pragma(f,"loop_function"))) {
    // This function can be called from a loop,
    // handle as if it was called from one
    loop_function_check(f);
    parsing_state.loop_function_next = false;
  }

  // Check if the function can be called from a loop
  bool loop_callable = true;
  // llvm::errs() << "Function " << f->getNameInfo().getName() << "\n";
  
  if (f->isThisDeclarationADefinition() && f->hasBody()) {
    global.currentFunctionDecl = f;
    
    Stmt *FuncBody = f->getBody();

    // Type name as string
    QualType QT = f->getReturnType();
    std::string TypeStr = QT.getAsString();

    // Function name
    DeclarationName DeclName = f->getNameInfo().getName();
    std::string FuncName = DeclName.getAsString();

    // llvm::errs() << " - Function "<< FuncName << "\n";

      if (does_function_contain_loop(f)) {
        loop_callable = false;
      }

     
    switch (f->getTemplatedKind()) {
      case FunctionDecl::TemplatedKind::TK_NonTemplate:
        // Normal, non-templated class method -- nothing here
        break;
        
      case FunctionDecl::TemplatedKind::TK_FunctionTemplate:
        // not descent inside templates
        parsing_state.skip_children = 1;
        break;
        
      case FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization:

        if (does_function_contain_loop(f)) {
          specialize_function_or_method(f);
        } else {
          parsing_state.skip_children = 1;  // no reason to look at it further
        }
        break;
        
      default:
        // do nothing
        break;
    }

    SourceLocation ST = f->getSourceRange().getBegin();
    global.location.function = ST;

    if (cmdline::funcinfo) {
      // Add comment before
      std::stringstream SSBefore;
      SSBefore << "// Begin function " << FuncName << " returning " << TypeStr
               << " of template type " << print_TemplatedKind(f->getTemplatedKind())
               << "\n";
      writeBuf->insert(ST, SSBefore.str(), true,true);
    }
    
  }

  return true;
}



void MyASTVisitor::specialize_function_or_method( FunctionDecl *f ) {
  // This handles all functions and methods. Parent is non-null for methods,
  // and then is_static gives the static flag
  
  bool no_inline;
  bool is_static = false;
  CXXRecordDecl * parent = nullptr;

  /* Check if the function is a class method */
  if(f->isCXXClassMember()){
    // method is defined inside template class.  Could be a chain of classes!
    CXXMethodDecl *method = dyn_cast<CXXMethodDecl>(f);
    parent = method->getParent();
    is_static = method->isStatic();
    no_inline = cmdline::method_spec_no_inline;
  } else {
    no_inline = cmdline::function_spec_no_inline;
  }
  
  srcBuf * writeBuf_saved = writeBuf;
  srcBuf funcBuf(&TheRewriter,f);
  writeBuf = &funcBuf;
  
  std::vector<std::string> par, arg;

  // llvm::errs() << "funcBuffer:\n" << funcBuf.dump() << '\n';

  // cannot rely on getReturnTypeSourceRange() for methods.  Let us not even try,
  // change the whole method here
  
  bool is_templated = ( f->getTemplatedKind() ==
                        FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization );
  
  int ntemplates = 0;
  std::string template_args = "";
  std::vector<const TemplateArgument *> typeargs = {};

  if (is_templated) {
    // Get here the template param->arg mapping for func template
    auto tal = f->getTemplateSpecializationArgs();
    auto tpl = f->getPrimaryTemplate()->getTemplateParameters();
    assert( tal && tpl && tal->size() == tpl->size() && "Method template par/arg error");

    make_mapping_lists(tpl, *tal, par, arg, typeargs, &template_args);
    ntemplates = 1;
  }

  // Get template mapping for classes
  // parent is from: CXXRecordDecl * parent = method->getParent();   
  if (parent) ntemplates += get_param_substitution_list( parent, par, arg, typeargs );
  llvm::errs() << "Num nesting templates " << ntemplates << '\n';

  funcBuf.replace_tokens(f->getSourceRange(), par, arg );

  // template_args adds template specialization args after the name, name<args>(..)
  funcBuf.replace(f->getNameInfo().getSourceRange(),
                  f->getQualifiedNameAsString() + template_args);
  
// #define use_ast_type
#ifdef use_ast_type
  // replace type as written with the type given in ast (?qualifiers)
  // we could also leave the "written" type as is.  Problems with array types?
  int i = funcBuf.get_index(f->getNameInfo().getSourceRange().getBegin());
  if (i > 0)
    funcBuf.replace(0,i-1,remove_class_from_type(f->getReturnType().getAsString()) + " ");
  else 
    funcBuf.insert(0,remove_class_from_type(f->getReturnType().getAsString()) + " ",true,false);
  
#else 
  
  // remove "static" if it is so specified in methods
  if (is_static) { 
    funcBuf.replace_token(0,
                          funcBuf.get_index(f->getNameInfo().getSourceRange().getBegin()),
                          "static","");
  }

#endif

  if (!f->isInlineSpecified() && !no_inline)
    funcBuf.insert(0, "inline ", true, true);

  for (int i=0; i<ntemplates; i++) {
    funcBuf.insert(0,"template <>\n",true,true);
  }

  check_spec_insertion_point(typeargs, global.location.bot, f);

  SourceRange decl_sr = get_func_decl_range(f);
  std::string wheredefined = "";
  if (f->isInlineSpecified() || !no_inline ||
      !in_specialization_db(funcBuf.get(decl_sr), wheredefined)) {
    // Now we should write the spec here
      
    // llvm::errs() << "new func:\n" << funcBuf.dump() <<'\n';
    // visit the body
    TraverseStmt(f->getBody());

    // llvm::errs() << "new func again:\n" << funcBuf.dump() <<'\n';

    // insert after the current toplevedecl
    std::stringstream sb;
    sb << "\n// ++++++++ Generated function/method specialization\n"
       << funcBuf.dump() 
       << "\n// ++++++++\n";
    toplevelBuf->insert( getSourceLocationAtEndOfLine(global.location.bot),
                         sb.str(), false, true );
  } else { 
    // Now the function has been written before (and not inline)
    // just insert declaration, defined on another compilation unit
    toplevelBuf->insert( getSourceLocationAtEndOfLine(global.location.bot),
            "\n// ++++++++ Generated specialization declaration, defined in compilation unit "
                         + wheredefined + "\n"
                         + funcBuf.get(decl_sr)
                         + ";\n// ++++++++\n",
                         false, false);
  }
    
  writeBuf = writeBuf_saved;
  funcBuf.clear();
  // don't descend again
  parsing_state.skip_children = 1;
}


// locate range of specialization "template< ..> .. func<...>( ... )"
// tf is ptr to template, and f to instantiated function
SourceRange MyASTVisitor::get_func_decl_range(FunctionDecl *f) {

  if (f->hasBody()) {
    SourceLocation a = f->getSourceRange().getBegin();
    SourceLocation b = f->getBody()->getSourceRange().getBegin();
    SourceManager &SM = TheRewriter.getSourceMgr();
    while (SM.getFileOffset(b) >= SM.getFileOffset(a)) {
      b = b.getLocWithOffset(-1);
      const char * p = SM.getCharacterData(b);
      if (!std::isspace(*p)) break;
    }
    SourceRange r(a,b);
    return r;
  }
  
  return f->getSourceRange();
}

bool MyASTVisitor::VisitClassTemplateDecl(ClassTemplateDecl *D) {
  
  if (parsing_state.dump_ast_next) {
    llvm::errs() << "**** Dumping class template declaration: \'" << D->getNameAsString() << "\'\n";
    D->dump();
    parsing_state.dump_ast_next = false;
  }

  // go through with real definitions or as a part of chain
  if (D->isThisDeclarationADefinition()) { // } || state::class_level > 0) {

    // insertion pt for specializations
//     if (state::class_level == 1) {
//       global.location.spec_insert = getSourceLocationAtEndOfLine(D->getSourceRange().getEnd());
//     }

    const TemplateParameterList * tplp = D->getTemplateParameters();
    // save template params in a list, for templates within templates .... ugh!
    // global.class_templ_params.push_back( tplp );
    
    // this block for debugging
    if (cmdline::funcinfo) {
      std::stringstream SSBefore;
      SSBefore << "// Begin template class "
               << D->getNameAsString()
               << " with template params " ;
      for (unsigned i = 0; i < tplp->size(); i++) 
        SSBefore << tplp->getParam(i)->getNameAsString() << " ";
      SourceLocation ST = D->getSourceRange().getBegin();
      SSBefore << '\n';
    
      writeBuf->insert(ST, SSBefore.str(), true, true);
    }
    // end block
    
    // global.in_class_template = true;
    // Should go through the template in order to find function templates...
    // Comment out now, let roll through "naturally".
    // TraverseDecl(D->getTemplatedDecl());

    if (D->getNameAsString() == "field") {
      handle_field_specializations(D);
    } else if (D->getNameAsString() == "field_storage") {
      field_storage_decl = D;
    } else {
    }

    // global.in_class_template = false;

    // Now do traverse the template naturally
    // state::skip_children = 1;
    
  }    
  
  return true;
}

// Find the element typealias here -- could not work
// directly with VisitTypeAliasTemplateDecl below, a bug??
bool MyASTVisitor::VisitDecl( Decl * D) {
  if (parsing_state.check_loop && state::loop_found) return true;

  if (parsing_state.dump_ast_next) {
    llvm::errs() << "**** Dumping declaration:\n";
    D->dump();
    parsing_state.dump_ast_next = false;
  }

  auto t = dyn_cast<TypeAliasTemplateDecl>(D);
  if (t && t->getNameAsString() == "element") {
    llvm::errs() << "Got field storage\n";
  }
  
  return true;
}

#if 0
bool MyASTVisitor::
VisitClassTemplateSpecalializationDecl(ClassTemplateSpecializationDecl *D) {
  if (D->getNameAsString() == "field") {    
    const TemplateArgumentList & tal = D->getTemplateArgs();
    llvm::errs() << " *** field with args ";
    for (unsigned i = 0; i < tal.size(); i++) 
      llvm::errs() << TheRewriter.getRewrittenText(tal.get(i).getAsExpr()->getSourceRange())
                   << " ";
    llvm::errs() << "\n";
  }
  return true;
}
#endif




/////////////////////////////////////////////////////////////////////////////////
/// Check that all template specialization type arguments are defined at the point
/// where the specialization is inserted
/// TODO: change the insertion point
/////////////////////////////////////////////////////////////////////////////////

void MyASTVisitor::check_spec_insertion_point(std::vector<const TemplateArgument *> & typeargs,
                                              SourceLocation ip, 
                                              FunctionDecl *f) 
{
  SourceManager &SM = TheRewriter.getSourceMgr();

  for (const TemplateArgument * tap : typeargs) {
    llvm::errs() << " - Checking tp type " << tap->getAsType().getAsString() << '\n';
    const Type * tp = tap->getAsType().getTypePtrOrNull();
    // Builtins are fine too
    if (tp && !tp->isBuiltinType()) {
      RecordDecl * rd = tp->getAsRecordDecl();
      if (rd && SM.isBeforeInTranslationUnit( ip, rd->getSourceRange().getBegin() )) {
        reportDiag(DiagnosticsEngine::Level::Warning,
                   f->getSourceRange().getBegin(),
    "Specialization point for function appears to be before the declaration of type \'%0\', code might not compile",
                   tap->getAsType().getAsString().c_str());
      } 
    }
  }
}

/// Returns the mapping params -> args for class templates, inner first.  Return value
/// the number of template nestings
int MyASTVisitor::get_param_substitution_list( CXXRecordDecl * r,
                                               std::vector<std::string> & par,
                                               std::vector<std::string> & arg,
                                               std::vector<const TemplateArgument *> & typeargs ) {
  
  if (r == nullptr) return 0;

  int level = 0;
  if (r->getTemplateSpecializationKind() == TemplateSpecializationKind::TSK_ImplicitInstantiation) {

    ClassTemplateSpecializationDecl * sp = dyn_cast<ClassTemplateSpecializationDecl>(r);
    if (sp) {
      llvm::errs() << "Got specialization of " << sp->getNameAsString() << '\n';
      const TemplateArgumentList & tal = sp->getTemplateArgs();
      assert(tal.size() > 0);
    
      ClassTemplateDecl * ctd = sp->getSpecializedTemplate();
      TemplateParameterList * tpl = ctd->getTemplateParameters();
      assert(tpl && tpl->size() > 0);

      assert(tal.size() == tpl->size());
    
      make_mapping_lists(tpl, tal, par, arg, typeargs, nullptr);
    
      level = 1;
    }
  } else {
    llvm::errs() << "No specialization of class " << r->getNameAsString() << '\n';
  }
  
  auto * parent = r->getParent();
  if (parent) {
    if (CXXRecordDecl * pr = dyn_cast<CXXRecordDecl>(parent))
      return level + get_param_substitution_list(pr, par, arg, typeargs);
  }
  return level;
}

void MyASTVisitor::make_mapping_lists( const TemplateParameterList * tpl, 
                                       const TemplateArgumentList & tal,
                                       std::vector<std::string> & par,
                                       std::vector<std::string> & arg,
                                       std::vector<const TemplateArgument *> & typeargs,
                                       std::string * argset ) {

  if (argset) *argset = "< ";

  // Get argument strings without class, struct... qualifiers
  PrintingPolicy pp(Context->getLangOpts()); 

  
  for (int i=0; i<tal.size(); i++) {
    if (argset && i>0) *argset += ", ";
    switch (tal.get(i).getKind()) {
      case TemplateArgument::ArgKind::Type:
        arg.push_back( tal.get(i).getAsType().getAsString(pp) );
        par.push_back( tpl->getParam(i)->getNameAsString() );
        if (argset) *argset += arg.back();  // write just added arg
        typeargs.push_back( &tal.get(i) );  // save type-type arguments
        break;
        
      case TemplateArgument::ArgKind::Integral:
        arg.push_back( tal.get(i).getAsIntegral().toString(10) );
        par.push_back( tpl->getParam(i)->getNameAsString() );
        if (argset) *argset += arg.back();
        break;
        
      default:
        llvm::errs() << " debug: ignoring template argument of argument kind " 
                     << tal.get(i).getKind() 
                     << " with parameter "
                     << tpl->getParam(i)->getNameAsString() << '\n';
        exit(1);  // Don't know what to do
    }
  }
  if (argset) *argset += " >";

  return;

}

void MyASTVisitor::set_writeBuf(const FileID fid) {
  writeBuf = get_file_buffer(TheRewriter, fid);
  toplevelBuf = writeBuf;
}




