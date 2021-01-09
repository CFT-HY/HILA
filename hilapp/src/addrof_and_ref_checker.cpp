#include <sstream>
#include <iostream>
#include <string>

#include "myastvisitor.h"
#include "hilapp.h"


//////////////////////////////////////////////////////////////////////////////
/// An AST Visitor for checking address of operators and references inside field loops
///   - Flag "& f[X]""  -type operations as errors
///   - Taking a reference to field:  
///        "const auto & ref = f[X]"  : field is read (has already been flagged, nothing here)
///        "auto & ref = f[X]" : field is read and written, disallow [X+..]
///   - If "var" is an loop extern variable:
///       - Do not allow  "& var" (could allow const pointers?)
///       - Do not allow "auto & r = var", but allow "const auto & r = var"
///       
//////////////////////////////////////////////////////////////////////////////
class addrOfAndRefChecker : public GeneralVisitor, public RecursiveASTVisitor<addrOfAndRefChecker> {

public:
  /// Use general visitor constructor
  using GeneralVisitor::GeneralVisitor;

  /// Used to skip child nodes while walking the syntax tree
  int skip_children;
  /// Buffer for the generated loop text, which replaces parity-loop expressions
  srcBuf *loopBuf;

  /// Constructor
  addrOfAndRefChecker(Rewriter &R,ASTContext *C, srcBuf *wb) : GeneralVisitor(R,C) {
    skip_children = 0;
    loopBuf = wb;
  }

  /// implement the "skip_children" method from MyASTVisitor also here
  bool TraverseStmt(Stmt *S) {  
    if (skip_children > 0) skip_children++;    
    if (!skip_children) RecursiveASTVisitor<addrOfAndRefChecker>::TraverseStmt(S);
    if (skip_children > 0) skip_children--;  
    return true;
  }

  /// implement the "skip_children" method from MyASTVisitor also here
  bool TraverseDecl(Decl *D) {
    if (skip_children > 0) skip_children++;    
    if (!skip_children) RecursiveASTVisitor<addrOfAndRefChecker>::TraverseDecl(D);
    if (skip_children > 0) skip_children--;  
    return true;
  }

  ///////////////////////////////////////////////////////////////////////
  /// there does not appear to be unaryop visitor, so go through Expr visitor
  ///////////////////////////////////////////////////////////////////////
  bool VisitExpr(Expr *E) {
    if (UnaryOperator * UO = dyn_cast<UnaryOperator>(E)) {
      if (UO->getOpcode() == UnaryOperatorKind::UO_AddrOf) {
        // got addrof op, check the expression

        Expr * E = UO->getSubExpr();

        E = E->IgnoreParens();
        E = E->IgnoreImplicit();
        E = E->IgnoreCasts();

        /// Check here the declrefexpr

        if (DeclRefExpr * DRE = dyn_cast<DeclRefExpr>(E)) {

          // is it X-expression?  Highly unlikely
          if (is_X_type(DRE)) {
            // now trying to take address of f[X]
            reportDiag(DiagnosticsEngine::Level::Error,
                       UO->getSourceRange().getBegin(),
                       "address of expression '%0' is not allowed",
                       get_stmt_str(UO->getSubExpr()).c_str() );
            skip_children = 1;
            return true;
          }  // end of X_type

          // then variable refs

          if (VarDecl * decl = dyn_cast<VarDecl>(DRE->getDecl())) {
            // if we are inspecting address of operator
            for (var_info & vi : var_info_list) if (vi.decl == decl) {
              // found decl
              if (!vi.is_loop_local) {
                reportDiag(DiagnosticsEngine::Level::Error,
                           UO->getSourceRange().getBegin(),
                           "taking address of field loop external variable '%0' is not allowed, suggest using references. "
                           "If a pointer is necessary, copy first: 'auto v = %1; auto *p = &v;'",
                           get_stmt_str(DRE).c_str(),
                           get_stmt_str(UO->getSubExpr()).c_str() );
              }
              return true;  // found the var, can return
            }
            // we get here only if this was unknown var
            llvm::errs() << "Mystery var found!\n";
            exit(1);
          }
          return true;
        }  // end of declRefExpr

        // operator, esp. field refs

        if (CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E)) {
          if (is_field_parity_expr(OC) || is_field_with_X_expr(OC)) {
            // now trying to take address of f[X]
            reportDiag(DiagnosticsEngine::Level::Error,
                       UO->getSourceRange().getBegin(),
                       "taking address of Field variable '%0' is not allowed, suggest using references. "
                       "If a pointer is necessary, copy first: 'auto v = %1; auto *p = &v;'",
                       get_stmt_str(OC).c_str(),
                       get_stmt_str(UO->getSubExpr()).c_str() );

            skip_children = 1;
            return true;
          }
          return true;
        }
      } // end of &-op
    } // end of unary


    return true;
  } // end of visitexpr

  //////////////////////////////////////////////////////////////////////////
  /// Visit separately function call arguments
  //////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////
  /// Variable declarator visitor
  /// this should catch normal declarations
  //////////////////////////////////////////////////////////////////////////
  bool VisitDecl(Decl *D) {
    if (VarDecl * V = dyn_cast<VarDecl>(D)) {
      // it's a variable decl inside field loop
      if (V->getStorageClass() == StorageClass::SC_Extern ||
          V->getStorageClass() == StorageClass::SC_Static ||
          V->getStorageClass() == StorageClass::SC_PrivateExtern) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   D->getSourceRange().getBegin(),
                   "cannot declare static or extern variables in field loops.");
        return true;
      }

      // Check now the reference var declarations
      if (V->getType().getTypePtr()->isReferenceType()) {
        // a reference, check to what
        if (!V->hasInit()) return true; // this error must have been flagged already
        
        // For reference vars:  const double & d = dval;  it's the referred to val which is const qualified!
        // bool is_const_qualified = V->getType().isConstQualified() || V->getType().isLocalConstQualified();

        bool pure_out_ref = false;
        if ( loopBuf->get_previous_original_word( V->getSourceRange().getBegin() ) ==
             output_only_keyword) {
          pure_out_ref = true;
        }  

        Expr * E = V->getInit();
        bool is_const_qualified = E->getType().isConstQualified();

        E = E->IgnoreCasts();
        E = E->IgnoreParens();
        E = E->IgnoreImplicit();

        // reference to field var
        if (is_field_with_X_expr(E)) {
          // this must have been scanned before
          for( field_ref & r : field_ref_list) if( r.fullExpr == E ){
            r.is_read = !pure_out_ref;
            if (!is_const_qualified) {
              r.is_written = true;
              if (is_field_with_X_and_dir(E)) {
                reportDiag(DiagnosticsEngine::Level::Error,
                            E->getSourceRange().getBegin(),
                            "cannot take non-const reference of a Field variable with X+direction");
                return true;
              }
            }
            return true;
          }
          // we should never get here
          llvm::errs() << " -- Unexplored field! \n"; 
          exit(1);
        }

        // check non-const qualified references to variables
        if (!is_const_qualified) {
          if (DeclRefExpr * DRE = dyn_cast<DeclRefExpr>(E)) {
            if (VarDecl * decl = dyn_cast<VarDecl>(DRE->getDecl())) {
              // it's variable ref.  
              for (var_info & vi : var_info_list) if (vi.decl == decl) {
                if (!vi.is_loop_local) {
                  reportDiag(DiagnosticsEngine::Level::Error,
                             DRE->getSourceRange().getBegin(),
                             "cannot take non-const reference of a variable declared outside field loop");
                  return true;
                }
                // now var is loop local - nothing special to do
                
                return true;
              }
              // again should never get here
              llvm::errs() << " -- Mystery var!\n";
            }
          }
        }
      } 
    } // end of varDecl

    return true;
  } // end of visitDecl

};

///////////////////////////////////////////////////////////////////////////////////
/// And loop checker interface here
///////////////////////////////////////////////////////////////////////////////////
void MyASTVisitor::check_addrofops_and_refs(Stmt * S) {

  addrOfAndRefChecker arf(TheRewriter,Context,writeBuf);
  arf.TraverseStmt(S);
}

