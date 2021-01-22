#include "generalvisitor.h"

//////////////////////////////////////////////////////////////////////////////////
///  Implementation of some methods for generalvisitor
///  For full list of generalvisitor commands, see generalvisitor.h
//////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_duplicate_expr(const Expr * a, const Expr * b) {
  // Use the Profile function in clang, which "fingerprints"
  // statements
  llvm::FoldingSetNodeID IDa, IDb;
  a->Profile(IDa, *Context, true);
  b->Profile(IDb, *Context, true);
  return ( IDa == IDb );
}


/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_field_parity_expr(Expr *E) {

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

/////////////////////////////////////////////////////////////////
/// Checks if E is parity of a field (for example f[X]).
/// Catches both parity and X_plus_direction 
bool GeneralVisitor::is_field_with_X_expr(Expr *E) {
  E = E->IgnoreParens();
  CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

  if (OC &&
      strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
      is_field_expr(OC->getArg(0))) {

    return is_X_index_type(OC->getArg(1));

  }
  return false;   
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_field_with_X_and_dir(Expr *E) {
  E = E->IgnoreParens();
  CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);

  if (OC &&
      strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
      is_field_expr(OC->getArg(0))) {
    
    return is_X_and_dir_type(OC->getArg(1));

  }
  return false;   
}

/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_assignment_expr(Stmt * s, std::string * opcodestr, bool &iscompound) {
  if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
    if (OP->isAssignmentOp()) {

      // TODO: there should be some more elegant way to do this
      const char *sp = getOperatorSpelling(OP->getOperator());
      if ((sp[0] == '+' || sp[0] == '-' || sp[0] == '*' || sp[0] == '/')
          && sp[1] == '=') iscompound = true;
      else iscompound = false;
      if (opcodestr)
        *opcodestr = getOperatorSpelling(OP->getOperator());

      return true;
    }
  }

  // This is for arithmetic type assignments
  if (BinaryOperator *B = dyn_cast<BinaryOperator>(s)) {
    if (B->isAssignmentOp()) {
      iscompound = B->isCompoundAssignmentOp();
      if (opcodestr)
        *opcodestr = B->getOpcodeStr().str();
      return true;
    }
  }
  
  return false;
}

///////////////////////////////////////////////////////////////////
/// Check if the RHS of assignment is site dependent
///////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_rhs_site_dependent(Stmt *s, std::vector<var_info *> * vi) {

  if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {
    if (OP->isAssignmentOp()) {
      return is_site_dependent(OP->getArg(1),vi);
    }
  }

  if (BinaryOperator *B = dyn_cast<BinaryOperator>(s)) {
    if (B->isAssignmentOp()) {
      return is_site_dependent(B->getRHS(),vi);
    }
  }
  // one of these should have triggered!  
  assert(0 && "Internal error in RHS analysis");
}


/////////////////////////////////////////////////////////////////

bool GeneralVisitor::is_function_call_stmt(Stmt * s) {
  if (auto *Call = dyn_cast<CallExpr>(s)){
    // llvm::errs() << "Function call found: " << get_stmt_str(s) << '\n';
    return true;
  }
  return false;
}

/////////////////////////////////////////////////////////////////
/// is the stmt pointing now to a member call
bool GeneralVisitor::is_member_call_stmt(Stmt * s) {
  if (auto *Call = dyn_cast<CXXMemberCallExpr>(s)){
    // llvm::errs() << "Member call found: " << get_stmt_str(s) << '\n';
    return true;
  }
  return false;
}

/////////////////////////////////////////////////////////////////
/// is the stmt pointing now to a constructor
bool GeneralVisitor::is_constructor_stmt(Stmt * s) {
  if (auto *Call = dyn_cast<CXXConstructExpr>(s)){
    // llvm::errs() << "Constructor found: " << get_stmt_str(s) << '\n';
    return true;
  }
  return false;
}

/////////////////////////////////////////////////////////////////
/// is the stmt pointing now to a constructor
bool GeneralVisitor::is_user_cast_stmt(Stmt * s) {
  if (auto *ce = dyn_cast<ImplicitCastExpr>(s)) {
  // if (CastExpr *ce = dyn_cast<CastExpr>(s)) {
    if (ce->getCastKind() == CK_UserDefinedConversion) {
    // llvm::errs() << "Constructor found: " << get_stmt_str(s) << '\n';
      return true;
    }
  }
  return false;
}

/////////////////////////////////////////////////////////////////
/// Does the statement end with a semicolon
bool GeneralVisitor::isStmtWithSemicolon(Stmt * S) {
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

/////////////////////////////////////////////////////////////////

parity GeneralVisitor::get_parity_val(const Expr *pExpr) {
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
                 "hilapp internal error, unknown parity" );
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
