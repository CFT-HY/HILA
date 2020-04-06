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


bool MyASTVisitor::is_parity_index_type(Expr *E) {
  return (get_expr_type(E) == "parity");
}

// Checks if E is a parity Expr. Catches both parity and X_plus_direction 
bool MyASTVisitor::is_field_parity_expr(Expr *E) {
  E = E->IgnoreParens();
  CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

  if (OC &&
      strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
      is_field_expr(OC->getArg(0))) {

    return is_parity_index_type(OC->getArg(1));

  } else {
    // DON'T DO TEMPLATES NOW!  ONLY SPECIALIZATIONS
    #if 0
    // This is for templated expressions
    // for some reason, expr a[X] "getBase() gives X, getIdx() a...
    if (ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(E)) {
      Expr * lhs = ASE->getLHS()->IgnoreParens();
      
      if (is_field_expr(ASE->getLHS()->IgnoreParens())) {
        // llvm::errs() << " FP: and field\n";
        return is_parity_index_type(ASE->getRHS());
      }
    }
    #endif
  }
  return false;   
}

bool MyASTVisitor::is_X_index_type(Expr *E) {
  std::string s = get_expr_type(E);
  if (s == "X_index_type" || s == "X_plus_direction" || s == "X_plus_offset") 
    return true;
  else 
    return false;
}

// Checks if E is a parity Expr. Catches both parity and X_plus_direction 
bool MyASTVisitor::is_field_with_X_expr(Expr *E) {
  E = E->IgnoreParens();
  CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

  if (OC &&
      strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
      is_field_expr(OC->getArg(0))) {

    return is_X_index_type(OC->getArg(1));

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
bool MyASTVisitor::is_member_call_stmt(Stmt * s) {
  if (auto *Call = dyn_cast<CXXMemberCallExpr>(s)){
    llvm::errs() << "Member call found: " << get_stmt_str(s) << '\n';
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

bool MyASTVisitor::handle_field_parity_X_expr(Expr *e, bool is_assign, bool is_compound, bool is_X) {
    
  e = e->IgnoreParens();
  field_ref lfe;

  // we know here that Expr is of field-parity type
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


  //lfe.nameInd    = writeBuf->markExpr(lfe.nameExpr); 
  //lfe.parityInd  = writeBuf->markExpr(lfe.parityExpr);
  
  lfe.is_written = is_assign;
  lfe.is_read = (is_compound || !is_assign);
  lfe.sequence = parsing_state.stmt_sequence;


  std::string parity_expr_type = get_expr_type(lfe.parityExpr);

  if (parity_expr_type == "parity") {
    if (is_X) {
      llvm::errs() << "Internal error in handle_loop_parity\n";
      exit(-1);
    }
    if (parsing_state.accept_field_parity) {
      // 1st parity statement on a single line lattice loop
      loop_parity.expr  = lfe.parityExpr;
      loop_parity.value = get_parity_val(loop_parity.expr);
      loop_parity.text  = get_stmt_str(loop_parity.expr);
    } else {
      reportDiag(DiagnosticsEngine::Level::Error,
                 lfe.parityExpr->getSourceRange().getBegin(),
                 "field[parity] not allowed here, use field[X] -type instead" );
    }
  }
  
  // next ref must have wildcard parity
  parsing_state.accept_field_parity = false;
        
  if (parity_expr_type == "X_plus_direction" || 
      parity_expr_type == "X_plus_offset") {

    if (is_assign) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 lfe.parityExpr->getSourceRange().getBegin(),
                 "X + dir -type reference not allowed on the LHS of an assignment");
    }

    // Now need to split the expr to parity and dir-bits
    // Because of offsets this is pretty complicated to do in AST.
    // We now know that the expr is of type
    // [X+direction]  or  [X+coordinate_vector] -- just
    // use the textual form of the expression!

    bool has_X;
    lfe.direxpr_s = remove_X( get_stmt_str(lfe.parityExpr), &has_X );

    if (!has_X) {
      reportDiag(DiagnosticsEngine::Level::Fatal,
                 lfe.parityExpr->getSourceRange().getBegin(),
                 "Internal error: index should have been X" );
      exit(-1);
    }

    // llvm::errs() << "Direxpr " << lfe.direxpr_s << '\n';

    lfe.is_direction = true;

    if (parity_expr_type == "X_plus_offset") {

      // It's an offset, no checking here to be done
      lfe.is_offset = true;

    } else {

      // Now make a check if the reference is just constant (XUP etc.)
      // Need to descent quite deeply into the expr chain
      Expr* e = lfe.parityExpr->IgnoreParens()->IgnoreImplicit();
      CXXOperatorCallExpr* Op = dyn_cast<CXXOperatorCallExpr>(e);
      if (!Op) {
        if (CXXConstructExpr * Ce = dyn_cast<CXXConstructExpr>(e)) {
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
                   "Internal error: could not parse X + direction/offset -statement" );
        exit(1);
      }

      Expr *dirE = Op->getArg(1)->IgnoreImplicit();
      llvm::APSInt result;
      if (dirE->isIntegerConstantExpr(result, *Context)) {
        // Got constant
        lfe.is_constant_direction = true;
        lfe.constant_value = result.getExtValue();
        // llvm::errs() << " GOT DIR CONST, value " << lfe.constant_value << "  expr " << lfe.direxpr_s << '\n';
      } else {
        lfe.is_constant_direction = false;
        // llvm::errs() << "GOT DIR NOT-CONST " << lfe.direxpr_s << '\n';
      
        // If the direction is a variable, add it to the list
        // DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(lfe.dirExpr);
        // static std::string assignop;
        // if(DRE && isa<VarDecl>(DRE->getDecl())) {
        //   handle_var_ref(DRE, false, assignop);
        // }

        // traverse the dir-expression to find var-references etc.
        TraverseStmt(lfe.parityExpr);
      }
    }
  } // end of "direction"-branch
    
  // llvm::errs() << "field expr " << get_stmt_str(lfe.nameExpr)
  //              << " parity " << get_stmt_str(lfe.parityExpr)
  //              << "\n";


  // Check that there are no local variable references up the AST
  FieldRefChecker frc(TheRewriter, Context);
  frc.TraverseStmt(lfe.fullExpr);
   
  field_ref_list.push_back(lfe);
      
  return(true);
}


///  Utility to find the reduction typ

reduction get_reduction_type(bool is_assign,
                             std::string & assignop,
                             var_info & vi) {
  if (is_assign && (!vi.is_loop_local)) {
    if (assignop == "+=") return reduction::SUM;
    if (assignop == "*=") return reduction::PRODUCT;
  }
  return reduction::NONE;
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
      // llvm::errs() << " + Got " << vi.type << '\n';

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



/// Find the the base of a compound variable expression
DeclRefExpr * find_base_variable(Expr * E){
  Expr * RE = E;

  while(!dyn_cast<DeclRefExpr>(RE)){
    // RE may be a compound expression. We want the base variable.
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
        return nullptr;
      }
    } else {
      // It's not a variable
      return nullptr;
    }
  }
  return dyn_cast<DeclRefExpr>(RE);
}


bool is_variable_loop_local(VarDecl * decl){
  for (var_decl & d : var_decl_list) {
    if (d.scope >= 0 && decl == d.decl) {
      llvm::errs() << "loop local var ref! \n";
      return true;
    }
  }
  return false;
}



// handle an array subscript expression
void MyASTVisitor::handle_array_var_ref(ArraySubscriptExpr *E,
                                        bool is_assign,
                                        std::string &assignop) {
  PrintingPolicy pp(Context->getLangOpts());

  // array refs are OK if they're inside the field element type,
  // for example  f[X].c[i][j]
  // Try to find these

  // Check if it's local
  DeclRefExpr * DRE = find_base_variable(E);
  VarDecl * decl = dyn_cast<VarDecl>(DRE->getDecl());
  bool array_local = is_variable_loop_local(decl);

  // Also check the index
  bool index_local;
  DRE = dyn_cast<DeclRefExpr>(E->getIdx()->IgnoreImplicit());
  if(DRE){
    decl = dyn_cast<VarDecl>(DRE->getDecl());
    index_local = is_variable_loop_local(decl);
  } else {
    // This happens when the index is not a variable.
    // It's probably a compile time constant
    index_local = true;
  }

  if(!array_local){
    llvm::errs() << "Non-local array\n";
    if(!index_local){
      llvm::errs() << "Non-local index\n";
      // It's defined completely outside the loop. Can be replaced with a
      // temporary variable
      array_ref ar;
      ar.ref = E;
      
      // Find the type of the full expression (that is an element of the array)
      auto type = E->getType().getCanonicalType().getUnqualifiedType();
      ar.type = type.getAsString(pp);
    
      array_ref_list.push_back(ar);
    } else {
      llvm::errs() << "Local index\n";

      // The array is defined outside, but the index is local. This is 
      // the most problematic case.
      // For now, only allow this if it's a
      // histogram reduction

      return;

      reportDiag(DiagnosticsEngine::Level::Error,
                 E->getSourceRange().getBegin(),
                 "Cannot mix loop local index and predefined array. You must use std::vector in a vector reduction." );
    }
  } else {
    llvm::errs() << "Local array\n";
    if( !index_local ){
      llvm::errs() << "Local index\n";
      // The index needs to be communicated to the loop. It's a normal variable,
      // so we can handle it as such
      handle_var_ref(DRE, is_assign, assignop);
    } else {
      llvm::errs() << "Local index\n";
      // Array and index are local. This does not require any action.
    }
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
  parsing_state.ast_depth = 0;   // renormalize to the beginning of loop
  parsing_state.stmt_sequence = 0;
  TraverseStmt(ls);
  parsing_state.in_loop_body = false;
  parsing_state.ast_depth = 0;

  // Remove exprs which we do not want
  for (Expr * e : remove_expr_list) writeBuf->remove(e);
  
  // check and analyze the field expressions
  check_field_ref_list();
  check_var_info_list();

  // check that loop_parity is not X  -- impossible now
  // if (loop_parity.value == parity::x) {
  //   reportDiag(DiagnosticsEngine::Level::Error,
  //              loop_parity.expr->getSourceRange().getBegin(),
  //              "Parity of the full loop cannot be \'X\'");
  // }
  
  generate_code(ls);
  
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

  // depth = 1 is the "top level" statement, should give fully formed
  // c++ statements separated by ';'.  These act as sequencing points
  // This is used to obtain assignment and read ordering
  if (parsing_state.ast_depth == 1) parsing_state.stmt_sequence++;

 
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

  // Check c++ methods  -- HMM: it seems above function call stmt catches these first
  if( is_member_call_stmt(s) ){
    handle_member_call_in_loop(s);
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
    if (is_field_with_X_expr(E)) {
      // It is field[X] reference
      // get the expression for field name
          
      handle_field_parity_X_expr(E, is_assignment, is_compound, true);
      is_assignment = false;  // next will not be assignment
      // (unless it is a[] = b[] = c[], which is OK)

      parsing_state.skip_children = 1;
      return true;
    }


    if (is_field_parity_expr(E)) {
      // Now we know it is a field parity reference
      // get the expression for field name
          
      handle_field_parity_X_expr(E, is_assignment, is_compound, false);
      is_assignment = false;  // next will not be assignment
      // (unless it is a[] = b[] = c[], which is OK)

      parsing_state.skip_children = 1;
      return true;
    }

    if (is_field_expr(E)) {
      // field without [X], bad usually (TODO: allow  scalar func(field)-type?)
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
      handle_array_var_ref(a, is_assignment, assignop);
      
      // We don't want to handle the array variable or the index separately
      //parsing_state.skip_children = 1;
      return true;
    }

    // Check for a vector reduction
    if( is_assignment && isa<CXXOperatorCallExpr>(s) ){
      CXXOperatorCallExpr * OC = dyn_cast<CXXOperatorCallExpr>(s);
      std::string type = OC->getArg(0)->getType().getAsString();
      if( type.rfind("std::vector<",0) != std::string::npos ){
        // It's an assignment to a vector element
        // Still need to check if it's a reduction
        DeclRefExpr * DRE = dyn_cast<DeclRefExpr>(OC->getArg(0)->IgnoreImplicit());
        VarDecl * vector_decl = dyn_cast<VarDecl>(DRE->getDecl());
        bool array_local = is_variable_loop_local(vector_decl);
        
        DRE = dyn_cast<DeclRefExpr>(OC->getArg(1)->IgnoreImplicit());
        VarDecl * index_decl = dyn_cast<VarDecl>(DRE->getDecl());
        bool index_local = is_variable_loop_local(index_decl);

        // Handle the index as a variable (it's local, so the name won't change)
        handle_var_ref(DRE,false,assignop);

        if( !array_local && index_local ) {
          llvm::errs() << "Found a vector reduction\n";
          vector_reduction_ref vrf;
          vrf.ref = OC;
          vrf.vector_name = vector_decl->getName();
          vrf.index_name = index_decl->getName();
          if( type.rfind("float",0) != std::string::npos ){
            vrf.type = "float";
          } else {
            vrf.type = "double";
          }
          if (assignop == "+=") {
            vrf.reduction_type = reduction::SUM;
          } else if (assignop == "*="){
            vrf.reduction_type = reduction::PRODUCT;
          } else {
            vrf.reduction_type = reduction::NONE;
          }
          vector_reduction_ref_list.push_back(vrf);
          parsing_state.skip_children = 1;
        }

      }
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

    // Reset ast_depth, so that depth == 0 again for the block.
    if (isa<CompoundStmt>(s)) parsing_state.ast_depth = -1; 

    TraverseStmt(s);

    parsing_state.ast_depth = 0;

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



bool MyASTVisitor::has_pragma(Stmt *S, const char * n) {
  return has_pragma( S->getSourceRange().getBegin(), n);
}

bool MyASTVisitor::has_pragma(Decl *F, const char * n) {
  return has_pragma( F->getSourceRange().getBegin(), n);
}

bool MyASTVisitor::has_pragma(const SourceLocation l, const char * n) {
  std::string arg;
  SourceLocation pragmaloc,sl = l;

  // if macro, get the unexpanded loc
  if (sl.isMacroID()) {
    CharSourceRange CSR = TheRewriter.getSourceMgr().getImmediateExpansionRange( sl );
    sl = CSR.getBegin();
  }

  if (is_preceded_by_pragma(sl, arg, pragmaloc) && (arg.find(n) != std::string::npos) ) {

    // got it, comment out -- check that it has not been commented out before
    // the buffer may not be writeBuf, so be careful

    srcBuf * sb = get_file_srcBuf(pragmaloc);

    int loc = sb->find_original(pragmaloc,'#');
    if (loc < 0) {
      llvm::errs() << "internal error in pragma handling\n";
      exit(1);
    }
    std::string s = sb->get(loc,loc+1);
    if (s.at(0) == '#') sb->insert(loc ,"//-- ",true,false);

    return true;
  }

  return false;
}


/// Check if the SourceLocation l is preceded by "#pragma transformer" on previous line.
/// There cannot be anything except whitespace between l and the beginning of line
/// Pragma_args will point to the beginning of arguments of pragma
bool MyASTVisitor::is_preceded_by_pragma( SourceLocation l0 , std::string & arguments, 
                                          SourceLocation & pragmaloc ) {
  SourceLocation l = l0;
  SourceLocation lend;
  SourceManager &SM = TheRewriter.getSourceMgr();
  bool found_line_break = false;
  bool got_non_space_chars = false;
  
  // Move backward from location l and find two linebreaks
  // These are the beginning and the end of the previous line
  constexpr int maxiter = 5000;
  for (int i=0; i<maxiter; i++) {
    l=l.getLocWithOffset(-1);

    bool invalid = false;
    const char * c = SM.getCharacterData(l,&invalid);
    if (invalid) {
      // Invalid character found. This is probably the beginning of a file
      return false;
    }
    if (*c == '\n'){
      if(!found_line_break) {
        // this is the 1st line break
        found_line_break = true;
      } else {
        // Second line break, exit here if this line was not empty
        if (got_non_space_chars) break;
      }
    } else {
      if (!got_non_space_chars && !isspace(*c)) {
        // non-space chars before the beginning of line where l was
        if (!found_line_break) return false;

        // now we got non-empty line above
        got_non_space_chars = true;
        lend = l;    // end loc of non-trivial txt
      } 
    }
    if (i == maxiter-1) {
      llvm::errs() << "Search error in is_preceded_by_pragma\n";
      return false;
    }
  }

  // l points to \n, skip
  //l = l.getLocWithOffset(1);
  // Now l points to the beginning of prospective #pragma transformer -line
  // Get first the source text 
  std::string txt = TheRewriter.getRewrittenText(SourceRange(l,lend));

  txt = remove_extra_whitespace(txt);

  std::string comp = "#pragma transformer";
  if (txt.compare(0,comp.length(),comp) == 0) {
    // found it, set return value
    arguments = txt.substr(comp.length()+1, std::string::npos);
    pragmaloc = l;
    return true;
  }
  // there could be space between # and pragma
  comp = "# pragma transformer";
  if (txt.compare(0,comp.length(),comp) == 0) {
    // found it, set return value
    arguments = txt.substr(comp.length()+1, std::string::npos);
    pragmaloc = l;
    return true;
  }
  
  return false; 
}


/// These are the main traverse methods
/// By overriding these methods in MyASTVisitor we can control which nodes are visited.
/// These are control points for the depth of the traversal;
///  check_loop, skip_children,  ast_depth

bool MyASTVisitor::TraverseStmt(Stmt *S) {

  if (parsing_state.check_loop && state::loop_found) return true;
    
  // if state::skip_children > 0 we'll skip all until return to level up
  if (parsing_state.skip_children > 0) parsing_state.skip_children++;
    
  // go via the original routine...
  if (!parsing_state.skip_children) {
    parsing_state.ast_depth++;
    RecursiveASTVisitor<MyASTVisitor>::TraverseStmt(S);
    if (parsing_state.ast_depth > 0) parsing_state.ast_depth--;
  }

  if (parsing_state.skip_children > 0) parsing_state.skip_children--;
      
  return true;
}

bool MyASTVisitor::TraverseDecl(Decl *D) {

  if (parsing_state.check_loop && state::loop_found) return true;

  // if state::skip_children > 0 we'll skip all until return to level up
  if (parsing_state.skip_children > 0) parsing_state.skip_children++;
    
  // go via the original routine...
  if (!parsing_state.skip_children) {
    parsing_state.ast_depth++;
    RecursiveASTVisitor<MyASTVisitor>::TraverseDecl(D);
    if (parsing_state.ast_depth > 0) parsing_state.ast_depth--;
  }

  if (parsing_state.skip_children > 0) parsing_state.skip_children--;

  return true;
}



parity MyASTVisitor::get_parity_val(const Expr *pExpr) {
  SourceLocation SL;
  APValue APV;

  if (pExpr->isCXX11ConstantExpr( *Context, &APV, &SL )) {
    // Parity is now constant
    int64_t val = (APV.getInt().getExtValue());
    parity p;
    if (0 <= val && val <= (int)parity::all) {
      p = static_cast<parity>(val);
    } else {
      reportDiag(DiagnosticsEngine::Level::Fatal,
                 pExpr->getSourceRange().getBegin(),
                 "Transformer internal error, unknown parity" );
      exit(-1);
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

//  Obsolete when X is new type
// void MyASTVisitor::require_parity_X(Expr * pExpr) {
//   // Now parity has to be X (or the same as before?)
//   if (get_parity_val(pExpr) != parity::x) {
//     reportDiag(DiagnosticsEngine::Level::Error,
//                pExpr->getSourceRange().getBegin(),
//                "Use wildcard parity \"X\" or \"parity::x\"" );
//   }
// }

// finish the field_ref_list, and
// construct the field_info_list

bool MyASTVisitor::check_field_ref_list() {

  bool no_errors = true;
  
  global.assert_loop_parity = false;

  field_info_list.clear();
    
  for( field_ref & p : field_ref_list ) {

    std::string name = get_stmt_str(p.nameExpr);
      
    field_info * fip = nullptr;

    // search for duplicates: if found, lfip is non-null

    for (field_info & li : field_info_list) {
      if (name.compare(li.old_name) == 0) {
        fip = &li;
        break;
      }
    }

    if (fip == nullptr) {
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
      
      field_info_list.push_back(lfv);
      fip = & field_info_list.back();
    }
    // now lfip points to the right info element
    // copy that to lf reference
    p.info = fip;

    if (p.is_written && !fip->is_written) {
      // first write to this field var
      fip->first_assign_seq = p.sequence;
      fip->is_written = true;
    }

    // note special treatment of file[X] -read: it is real read only if it 
    // comes before or at the same time than assign
    if (p.is_read) {
      if (p.is_direction) {
        fip->is_read_nb = true;
      } else if ( !fip->is_written || fip->first_assign_seq >= p.sequence) {
        fip->is_read_atX = true;
      }
    }

    if (p.is_offset) fip->is_read_offset = true;
      
    // save expr record
    fip->ref_list.push_back(&p);

    if (p.is_direction) {

      if (p.is_written) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   p.parityExpr->getSourceRange().getBegin(),
                   "Neighbour offset not allowed on the LHS of an assignment");
        no_errors = false;
      }

      // does this dir with this field name exist before?
      // Use is_duplicate_expr() to resolve the ptr, it has (some) intelligence to
      // find equivalent constant expressions
      // TODO: use better method?
      bool found = false;
      for (dir_ptr & dp : fip->dir_list) {
        if (p.is_constant_direction) {
          found = (dp.is_constant_direction && dp.constant_value == p.constant_value);
        } else {
          found = is_duplicate_expr(dp.parityExpr, p.parityExpr);
        }

        if (found) {
          dp.count += (p.is_offset == false); // only nn in count
          dp.ref_list.push_back(&p);
          break;
        }
      }
      
      if (!found) {
        dir_ptr dp;
        dp.parityExpr = p.parityExpr;
        dp.count = (p.is_offset == false);
        dp.is_offset = p.is_offset;
        dp.is_constant_direction = p.is_constant_direction;
        dp.constant_value = p.constant_value;
        dp.direxpr_s = p.direxpr_s;     // copy the string expr of direction

        dp.ref_list.push_back(&p);

        fip->dir_list.push_back(dp);
      }
    } // direction
  } // p-loop
  
  // check for f[ALL] = f[X+dir] -type use, which is undefined
  
  for (field_info & l : field_info_list) {
    if (l.is_written && l.dir_list.size() > 0) {
 
      // There may be error, find culprits
      bool found_error = false;
      for (field_ref * p : l.ref_list) {
        if (p->is_direction && !p->is_written && !p->is_offset) {
          if (loop_parity.value == parity::all) {

            reportDiag(DiagnosticsEngine::Level::Error,
                       p->parityExpr->getSourceRange().getBegin(),
                       "Simultaneous access '%0' and assignment '%1' not allowed with parity ALL",
                       get_stmt_str(p->fullExpr).c_str(),
                       l.old_name.c_str());
            no_errors = false;
            found_error = true;

          } else if (loop_parity.value == parity::none) {
            reportDiag(DiagnosticsEngine::Level::Remark,
                       p->parityExpr->getSourceRange().getBegin(),
                       "Simultaneous access '%0' and assignment '%1' is allowed only with parity %2 is EVEN or ODD.  Inserting assertion",
                       get_stmt_str(p->fullExpr).c_str(),
                       l.old_name.c_str(),
                       loop_parity.text.c_str());
            found_error = true;
          }
        }
      }

      if (found_error) {
        for (field_ref * p : l.ref_list) {
          if (p->is_written && p->is_direction) {
            reportDiag(DiagnosticsEngine::Level::Remark,
                       p->fullExpr->getSourceRange().getBegin(),
                       "Location of assignment");
          }
        }
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
            if (j!=i) reportDiag(DiagnosticsEngine::Level::Remark,
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
SourceRange MyASTVisitor::getRangeWithSemicolon(Stmt * S, bool flag_error) {
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


bool MyASTVisitor::VisitVarDecl(VarDecl *var) {
  
  if (parsing_state.check_loop && state::loop_found) return true;
  
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

void MyASTVisitor::ast_dump_header(const char *s, const SourceRange sr_in) {
  SourceManager &SM = TheRewriter.getSourceMgr();
  SourceRange sr = sr_in;
  unsigned linenumber = SM.getSpellingLineNumber(sr.getBegin());

  // check if it is macro
  if (sr.getBegin().isMacroID()) {
    CharSourceRange CSR = TheRewriter.getSourceMgr().getImmediateExpansionRange( sr.getBegin() );
    sr = CSR.getAsRange();
  }

  std::string source = TheRewriter.getRewrittenText(sr);
  auto n = source.find('\n');

  if (n == std::string::npos) {
    llvm::errs() << "**** AST dump of " << s << " \'" << source << "\' on line "
                 << linenumber << '\n';
  } else {
    llvm::errs() << "**** AST dump of " << s << " starting with \'" 
                 << source.substr(0,n) << "\' on line " << linenumber << '\n';
  }
}


void MyASTVisitor::ast_dump(const Stmt *S) {
  ast_dump_header("statement", S->getSourceRange());
  S->dumpColor();
  llvm::errs() << "*****************************\n";
}


void MyASTVisitor::ast_dump(const Decl *D) {
  ast_dump_header("declaration", D->getSourceRange());
  D->dumpColor();
  llvm::errs() << "*****************************\n";
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
  
 
  if ( !parsing_state.check_loop && parsing_state.ast_depth == 1 &&
       has_pragma(s,"ast dump") ) {
    ast_dump(s);
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
  //  Arg(0)  is the LHS of assignment
  
  CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s);
  bool found = false;
  if (OP && OP->isAssignmentOp() && is_field_parity_expr(OP->getArg(0))) found = true;
  else {
    // check also field<double> or some other non-class var
    BinaryOperator *BO = dyn_cast<BinaryOperator>(s);
    if (BO && BO->isAssignmentOp() && is_field_parity_expr(BO->getLHS())) found = true;
  }

  if (found) {
    
    if (parsing_state.check_loop) {
      state::loop_found = true;
      return true;
    }

    SourceRange full_range = getRangeWithSemicolon(s,false);
    global.full_loop_text = TheRewriter.getRewrittenText(full_range);
        
    handle_full_loop_stmt(s, true);
    return true;
  }

  // And, for correct level for pragma handling - turns to 0 for stmts inside
  if (isa<CompoundStmt>(s)) parsing_state.ast_depth = -1;

  //  Finally, if we get to a field[parity] -expression without a loop or assignment flag error
  if (!parsing_state.check_loop) {
    Expr * E = dyn_cast<Expr>(s);
    if (E && is_field_parity_expr(E)) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 E->getSourceRange().getBegin(),
                 "field[parity] -expression is allowed only in LHS of field assignment statements (field[par] = ...)");
    } else if (E && is_field_with_X_expr(E)) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 E->getSourceRange().getBegin(),
                 "field[X] -expressions allowed only in RHS of field assignment statements or in \"onsites()\" blocks");

    }

  }

  // and add special handling for special function calls here

  if (CallExpr * CE = dyn_cast<CallExpr>(s)) {
    if (FunctionDecl * FD = CE->getDirectCallee()) {

      // is it memalloc(size) -call -> substitute with
      // memalloc( size, filename, linenumber)
      // don't know if this is really useful

      if (FD->getNameAsString() == "memalloc" && CE->getNumArgs() == 1) {
        SourceLocation sl = CE->getRParenLoc();
        SourceManager &SM = TheRewriter.getSourceMgr();
        // generate new args
        std::string name(SM.getFilename(sl));
        std::size_t i = name.rfind('/');
        if (i != std::string::npos) name = name.substr(i);

        std::string args(", \"");
        args.append(name).append( "\", " ).append(
            std::to_string( SM.getSpellingLineNumber(sl)));

        // if (!writeBuf->is_edited(sl)) writeBuf->insert(sl,args);
      }
    } 
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


bool MyASTVisitor::VisitFunctionDecl(FunctionDecl *f) {
  // Only function definitions (with bodies), not declarations.
  // also only non-templated functions
  // this does not really do anything

  if (!parsing_state.check_loop && has_pragma(f,"loop_function")) {
    // This function can be called from a loop,
    // handle as if it was called from one
    loop_function_check(f);
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

  if ( !parsing_state.check_loop && parsing_state.ast_depth == 1 &&
       has_pragma(D,"ast dump") ) {
    ast_dump(D);
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




