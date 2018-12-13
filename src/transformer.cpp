//------------------------------------------------------------------------------
// Transformer tools to convert "lattice loops" into
// hardware-dependent "kernels".
//
// Uses Clang RecursiveASTVisitor and Rewriter 
// interfaces
//
// Kari Rummukainen 2017-18
// 
//------------------------------------------------------------------------------
#include <sstream>
#include <string>

#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
//#include "llvm/Support/raw_ostream.h"

using namespace clang;
//using namespace clang::driver;
using namespace clang::tooling;

#include "transformer.h"
#include "stringops.h"
#include "srcbuf.h"
#include "myastvisitor.h"


// collection of variables holding the state of parsing
namespace state {
  unsigned skip_children = 0;
  unsigned scope_level = 0;
  int skip_next = 0;
  bool in_loop_body = false;
  bool accept_field_parity = false;
  bool is_file_modified = false;
}


static llvm::cl::OptionCategory ToolingSampleCategory("Transformer");

// command line options
static llvm::cl::opt<bool> dump_ast("d", llvm::cl::desc("dump AST tree"));
// static llvm::cl::extrahelp CommonHelp("  -d dump AST");

// local global vars
global_state global;
// and global loop_parity
loop_parity_struct loop_parity;

static const std::string field_element_type = "field_element<";
static const std::string field_type = "field<";

// global lists for functions
std::list<field_ref> field_ref_list = {};
std::list<field_info> field_info_list = {};
std::list<var_expr> var_expr_list = {};
std::list<var_decl> var_decl_list = {};

// take global CI just in case
CompilerInstance *myCompilerInstance;


// utility functions


bool MyASTVisitor::is_field_element_expr(Expr *E) {
  return( E && E->getType().getAsString().find(field_element_type) != std::string::npos);
}

bool MyASTVisitor::is_field_expr(Expr *E) {
  return( E && E->getType().getAsString().find(field_type) != std::string::npos);
}

bool MyASTVisitor::is_field_decl(ValueDecl *D) {
  return( D && D->getType().getAsString().find(field_type) != std::string::npos);
}
  
  
// // Define a pragma handler for #pragma example_pragma
// class heLppPragmaHandler : public PragmaHandler {
// public:
//   heLppPragmaHandler() : PragmaHandler("heLpp") { }
//   void HandlePragma(Preprocessor &PP, PragmaIntroducerKind Introducer,
//                     Token &PragmaTok) {
//     // Handle the pragma
//     llvm::errs() << "Got the pragma! name " << getName() << '\n';
//     llvm::errs() << TheRewriter.getRewrittenText(SourceRange(PragmaTok.getLocation(),
//                                                              PragmaTok.getLastLoc())) << '\n';
    
//   }
// };

// static PragmaHandlerRegistry::Add<heLppPragmaHandler> Y("heLpp","heL pragma description");

  
// By implementing RecursiveASTVisitor, we can specify which AST nodes
// we're interested in by overriding relevant methods.

bool MyASTVisitor::TraverseStmt(Stmt *S) {

  if (state::skip_next > 0) {
    state::skip_next--;
    return true;
  }
    
  if (global.in_func_template > 0) {
    global.in_func_template ++;
  }

  // if state::skip_children > 0 we'll skip all until return to level up
  if (state::skip_children > 0) state::skip_children++;
    
  // go via the original routine...
  if (!state::skip_children) RecursiveASTVisitor<MyASTVisitor>::TraverseStmt(S);

  if (state::skip_children > 0) state::skip_children--;
  
  if (global.in_func_template > 0) {
    global.in_func_template--;
    if (global.in_func_template == 1) {
      // Now we finished loop
      //llvm::errs() << "out of func template, last stmt " <<
      //  TheRewriter.getRewrittenText(S->getSourceRange()) << "\n";
      global.in_func_template = 0;
    }
  }
    
  return true;
}

template <unsigned N>
void MyASTVisitor::reportDiag(DiagnosticsEngine::Level lev, const SourceLocation & SL,
                              const char (&msg)[N],
                              const char *s1,
                              const char *s2,
                              const char *s3 ) {
  auto & DE = Context->getDiagnostics();    
  auto ID = DE.getCustomDiagID(lev, msg );
  auto DB = DE.Report(SL, ID);
  if (s1 != nullptr) DB.AddString(s1);
  if (s2 != nullptr) DB.AddString(s2);
  if (s3 != nullptr) DB.AddString(s3);
    
}
  
bool MyASTVisitor::is_duplicate_expr(const Expr * a, const Expr * b) {
  // Use the Profile function in clang, which "fingerprints"
  // statements
  llvm::FoldingSetNodeID IDa, IDb;
  a->Profile(IDa, *Context, true);
  b->Profile(IDb, *Context, true);
  return ( IDa == IDb );
}

  
// catches both parity and parity_plus_direction 
bool MyASTVisitor::is_lf_parity_expr(Expr *e) {
  e = e->IgnoreParens();
  CXXOperatorCallExpr *O = dyn_cast<CXXOperatorCallExpr>(e);
  if (O &&
      strcmp(getOperatorSpelling(O->getOperator()),"[]") == 0 && 
      is_field_expr(O->getArg(0))) {
    std::string s = get_expr_type(O->getArg(1)); 
    if (s == "parity" || s == "parity_plus_direction") {
      // llvm::errs() << " <<<Parity type " << get_expr_type(O->getArg(1)) << '\n';
      return true;
    } 
    // llvm::errs() << " <<<Parity type " << get_expr_type(O->getArg(1)) << '\n';
  }
  return false;
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
      exit(0);
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
  
bool MyASTVisitor::check_lf_ref_list() {

  bool no_errors = true;
  
  global.assert_loop_parity = false;

  field_info_list.clear();
    
  for( field_ref & p : field_ref_list ) {

    p.dirInd = -1;  // reset the direction
    
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
      lfv.type = get_expr_type(p.nameExpr);
      lfv.is_changed = p.is_changed;
        
      field_info_list.push_back(lfv);
      lfip = & field_info_list.back();
    }
    // now lfip points to the right info element
    // copy that to lf reference
    p.info = lfip;

    if (p.is_changed) lfip->is_changed = true;
      
    // save expr record
    lfip->ref_list.push_back(&p);

    if (p.dirExpr != nullptr) {

      if (p.is_changed) {
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
          p.dirInd = i;
          found = true;
          break;
        }
        i++;
      }
        
      if (!found) {
        dir_ptr dp;
        dp.e = p.dirExpr;
        dp.count = 1;
        p.dirInd = lfip->dir_list.size();

        lfip->dir_list.push_back(dp);
      }
    } // dirExpr
  } // p-loop
  
  // check for f[ALL] = f[X+dir] -type use, which is undefined
  
  for (field_info & l : field_info_list) {
    if (l.is_changed && l.dir_list.size() > 0) {
      if (loop_parity.value == parity::all) {
        // There's error, find culprits
        for (field_ref * p : l.ref_list) {
          if (p->dirExpr != NULL && !p->is_changed) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       p->parityExpr->getSourceRange().getBegin(),
                       "Accessing field '%0' undefined when assigning to '%1' with parity ALL, flagging as error",
                       get_stmt_str(p->fullExpr).c_str(),
                       l.old_name.c_str());
            no_errors = false;
          }
        }

        for (field_ref * p : l.ref_list) {
          if (p->is_changed && p->dirExpr == NULL) {
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
          
        
// This routine goes through one field reference and
// pushes the info to lists
  
bool MyASTVisitor::handle_lf_parity_expr(Expr *e, bool is_assign) {
    
  e = e->IgnoreParens();
  CXXOperatorCallExpr *O = dyn_cast<CXXOperatorCallExpr>(e);
  field_ref lfe;
  bool no_errors = true;
    
  lfe.fullExpr   = O;
  // take name 
  lfe.nameExpr   = O->getArg(0);
  lfe.nameInd    = Buf.markExpr(lfe.nameExpr); 
  
  lfe.parityExpr = O->getArg(1);
  lfe.parityInd  = Buf.markExpr(lfe.parityExpr);

  lfe.dirExpr  = nullptr;  // no neighb expression
    
  if (get_expr_type(lfe.parityExpr) == "parity") {
    if (state::accept_field_parity) {
      // 1st parity statement on a single line lattice loop
      loop_parity.expr  = lfe.parityExpr;
      loop_parity.value = get_parity_val(loop_parity.expr);
      loop_parity.text  = get_stmt_str(loop_parity.expr);
    } else {
      require_parity_X(lfe.parityExpr);
    }
  }

  lfe.is_changed = is_assign;
    
  // next ref must have wildcard parity
  state::accept_field_parity = false;
        
  if (get_expr_type(lfe.parityExpr) == "parity_plus_direction") {

    if (is_assign) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 lfe.parityExpr->getSourceRange().getBegin(),
                 "Neighbour offset not allowed on the LHS of an assignment");
      no_errors = false;
    }

    // Now need to split the expr to parity and dir-bits
    // Need to descent quite deeply into the expr chain
    Expr* e = lfe.parityExpr->IgnoreParens();
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
      exit(0);
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
        lfe.dirExpr = Op->getArg(1);
    }
  }
    
  llvm::errs() << "field expr " << get_stmt_str(lfe.nameExpr)
               << " parity " << get_stmt_str(lfe.parityExpr)
               << "\n";

   
  field_ref_list.push_back(lfe);
      
  return(no_errors);
}

////////

  
var_expr MyASTVisitor::handle_var_expr(Expr *E) {
  var_expr v;
  v.e = E;
  v.ind = Buf.markExpr(E);
  // This is somehow needed for printing type without "class" id
  PrintingPolicy pp(Context->getLangOpts());
  v.type = E->getType().getUnqualifiedType().getAsString(pp);
  
  // check if duplicate, use the Profile function in clang, which "fingerprints"
  // statements
  // NOTE: MYSTERIOUS BUG; DOES NOT ALWAYS RECOGNIZE REFS TO THE SAME VAR AS IDENTICAL
  llvm::FoldingSetNodeID thisID, ID;
  E->Profile(thisID, *Context, true);
  // llvm::errs() << "   comparing:: \'"<< get_stmt_str(E) <<"\'\n";
  v.duplicate = nullptr;
  for ( var_expr & p : var_expr_list ) {
    p.e->Profile(ID, *Context, true);
    // llvm::errs() << "   against:: \'" << get_stmt_str(p.e) << "\'\n";
    if ( thisID == ID ) {
      // dup found
      v.duplicate = &p;
      llvm::errs() << "Found dup: " << get_stmt_str(E) << '\n';
      break;
    }
  }
    
  return(v);
}

// check if stmt is lf[par] = ... -type
bool MyASTVisitor::is_lf_parity_assignment( Stmt *s ) {
  CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s);
  if (OP && OP->isAssignmentOp()) {
    // is assginment, verify that LHS is field_element
    if ( is_lf_parity_expr( OP->getArg(0) ) )
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
    llvm::errs() << "; found " << get_stmt_str(S) << '\n';
    return true;
  }
  return false;
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

  
//std::string generate_call_param_list
  

// start dealing with the loop statement
bool MyASTVisitor::handle_full_loop_stmt(Stmt *ls, bool field_parity_ok ) {
  // init edit buffer
  Buf.create( &TheRewriter, ls );
          
  field_ref_list.clear();
  var_expr_list.clear();
  var_decl_list.clear();
  
  state::accept_field_parity = field_parity_ok;
    
  // the following is for taking the parity from next elem
  state::scope_level = 0;
  state::in_loop_body = true;
  TraverseStmt(ls);
  state::in_loop_body = false;

  // check and analyze the field expressions
  check_lf_ref_list();
          
  generate_code( global.location.top, ls );
  
  Buf.clear();
          
  // Emit the original command as a commented line
  TheRewriter.InsertText(ls->getSourceRange().getBegin(),
                         comment_string(global.full_loop_text) + "\n",
                         false,true);

  global.full_loop_text = "";

  // don't go again through the arguments
  state::skip_children = 1;

  state::is_file_modified = true;
  
  return true;
}

/// Is expr a reference to a loop-local var
var_decl * MyASTVisitor::is_loop_local_var_ref(Expr *E) {
  E = E->IgnoreParens();
  if (isa<DeclRefExpr>(E)) {
    DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(E);
        
    if (isa<VarDecl>(DRE->getDecl())) {
      VarDecl * decl = dyn_cast<VarDecl>(DRE->getDecl());
      for( var_decl & v : var_decl_list ) {
        if (decl == v.decl) {
          // reference to loop-local var, leave as is
          // llvm::errs() << "loop local var ref! " << v.name << '\n';
          return &v;
        }
      }
    }
  }
  return nullptr;
}

bool MyASTVisitor::is_loop_extern_var_ref(Expr *E) {
  E = E->IgnoreParenCasts();
  while (isa<ArraySubscriptExpr>(E)) {
    E = dyn_cast<ArraySubscriptExpr>(E)->getBase();
    E = E->IgnoreParenCasts();
  }
  
  if (isa<DeclRefExpr>(E)) {
    DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(E);
        
    if (isa<VarDecl>(DRE->getDecl())) {
      VarDecl * decl = dyn_cast<VarDecl>(DRE->getDecl());
      for( var_decl & v : var_decl_list ) {
        if (decl == v.decl) return false;
      }
      return true;
    }
  }    
  return false;
}


// This handles statements within field loops
bool MyASTVisitor::handle_loop_body_stmt(Stmt * s) {

  // This keeps track of the assignment to field
  // must remember the set value across calls
  static bool is_assignment = false;
  
  // need to recognize assignments lf[X] =  or lf[X] += etc.
  // because we need to find changed fields
  if (is_lf_parity_assignment(s)) {
    is_assignment = true;
    // next visit will be to the assigned to field
    return true;
  }  
  
  // catch then expressions
  // TODO: reductions!
      
  if (Expr *E = dyn_cast<Expr>(s)) {
    
    // Not much to do with const exprs
    if (E->isCXX11ConstantExpr(*Context, nullptr, nullptr)) {
      llvm::errs() << "Constant expr: "
                   << TheRewriter.getRewrittenText(s->getSourceRange()) << "\n";
      state::skip_children = 1;
      return true;
    }
    
    if (is_field_element_expr(E)) {
      // run this expr type up until we find field variable refs
      if (is_lf_parity_expr(E)) {
        // Now we know it is a field parity reference
        // get the expression for field name
          
        handle_lf_parity_expr(E, is_assignment);
        is_assignment = false;  // next will not be assignment, unless it is

        // llvm::errs() << "Field expr: " << get_stmt_str(lfE.nameExpr) << "\n";
            
        state::skip_children = 1;
      }
      return true;
    }

    if (is_field_expr(E)) {
      // field without [parity], bad
      reportDiag(DiagnosticsEngine::Level::Error,
                 E->getSourceRange().getBegin(),
                 "Field expressions without [..] not allowed within \'onsites()\'" );
      state::skip_children = 1;  // once is enough
      return true;
    }

    // prevent assignments to loop-external vars within loops
    if (BinaryOperator * B = dyn_cast<BinaryOperator>(E)) {
      if (B->isAssignmentOp() && is_loop_extern_var_ref(B->getLHS())) {
        reportDiag(DiagnosticsEngine::Level::Error,
                   B->getLHS()->getSourceRange().getBegin(),
                   "Assignment to a non-field variable declared out of loop "
                   "not allowed within field loop" );
        state::skip_children = 1;  // once is enough
        return true;
      }
      return true;
    }

#if 1
    // TODO - more careful analysis of the expressions in the loops!
    if (var_decl * p = is_loop_local_var_ref(E)) {
      llvm::errs() << "loop local var ref! " << p->name << '\n';      
      return true;
    } 
      
    if (isa<DeclRefExpr>(E)) {
      // now it should be external var ref non-field
      var_expr_list.push_back( handle_var_expr(E) );
      
      state::skip_children = 1;
      llvm::errs() << " Some other declref: " << TheRewriter.getRewrittenText(E->getSourceRange()) << '\n';

      return true;
    }

    if (isa<ArraySubscriptExpr>(E)) {
      llvm::errs() << "  It's array expr "
                   << TheRewriter.getRewrittenText(E->getSourceRange()) << "\n";
      //state::skip_children = 1;
      auto a = dyn_cast<ArraySubscriptExpr>(E);
      Expr *base = a->getBase();
      
      //check_loop_local_vars = true;
      //TraverseStmt(
      
      return true;
      }
#endif          
        
    if (0){

      // not field type non-const expr
      llvm::errs() << "Non-const other Expr: " << get_stmt_str(E) << '\n';
      // loop-local variable refs inside? If so, we cannot evaluate this as "whole"

      // check_local_loop_var_refs = 1;
        
      // TODO: find really uniq variable references
      var_expr_list.push_back( handle_var_expr(E) );

      state::skip_children = 1;          
      return true;
    }
    // this point not reached
  } // Expr checking branch - now others...

  // This reached only if s is not Expr

  // start {...} -block or other compound
  if (isa<CompoundStmt>(s) || isa<ForStmt>(s) || isa<IfStmt>(s)
      || isa<WhileStmt>(s)) {
    static bool is_revisit = false;

    // traverse each stmt - next visit will be the same node, use is_revisit-trick
    if (is_revisit) {
      is_revisit = false;
      return true;
    }
    
    is_revisit = true;
    state::scope_level++;
    TraverseStmt(s);
    state::scope_level--;
    remove_vars_out_of_scope(state::scope_level);
    state::skip_children = 1;
    return true;
  }
    
  return true;    
} 


bool MyASTVisitor::VisitVarDecl(VarDecl *var) {
  if (state::in_loop_body) {
    // for now care only loop body variable declarations

    if (!var->hasLocalStorage()) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 var->getSourceRange().getBegin(),
                 "Static or external variable declarations not allowed within field loops");
      return true;
    }

    if (is_field_decl(var)) {
      reportDiag(DiagnosticsEngine::Level::Error,
                 var->getSourceRange().getBegin(),
                 "Cannot declare field variables within field traversal");
      state::skip_children = 1;
      return true;
    }

    // Now it should be automatic local variable decl
    var_decl vd;
    vd.decl = var;
    vd.name = var->getName();
    vd.type = var->getType().getAsString();
    vd.scope = state::scope_level;
    var_decl_list.push_back(vd);
    
    llvm::errs() << "Local var decl " << vd.name << " of type " << vd.type << '\n';
    return true;
  }
  return true;
}


void MyASTVisitor::remove_vars_out_of_scope(unsigned level) {
  while (var_decl_list.size() > 0 && var_decl_list.back().scope > level)
    var_decl_list.pop_back();
}


// VisitStmt is called for each statement in AST.  Thus, when traversing the
// AST or part of it we always start from here

bool MyASTVisitor::VisitStmt(Stmt *s) {
    
  // Entry point when inside field[par] = .... body
  if (state::in_loop_body) {
    return handle_loop_body_stmt(s);
  }
    
                                 
  //  Starting point for fundamental operation
  //  field[par] = .... - non-templated version
  
  // isStmtWithSemi(s);
    
  if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s)) {

    if (OP->isAssignmentOp()) {
        
      // is assginment, verify that LHS is field or field_element
      if (is_lf_parity_expr(OP->getArg(0))) {
        // now we have fieldop.  1st child will be
        // the lhs of the assignment

        SourceRange full_range = getRangeWithSemi(OP,false);
        global.full_loop_text = TheRewriter.getRewrittenText(full_range);
        
        handle_full_loop_stmt(OP, true);
          
      }
      return true;
    }
      
    return true;
  }

  // loop of type "onsites(p)"
  if (isa<ForStmt>(s)) {

    ForStmt *f = cast<ForStmt>(s);
    SourceLocation startloc = f->getSourceRange().getBegin();

    if (startloc.isMacroID()) {
      Preprocessor &pp = myCompilerInstance->getPreprocessor();
      static std::string loop_call("onsites");
      if (pp.getImmediateMacroName(startloc) == loop_call) {
        // Now we know it is onsites-macro
        CharSourceRange CSR = TheRewriter.getSourceMgr().getImmediateExpansionRange( startloc );
        std::string macro = TheRewriter.getRewrittenText( CSR.getAsRange() );
        bool internal_error = true;
        
        DeclStmt * init = dyn_cast<DeclStmt>(f->getInit());
        if (init && init->isSingleDecl() ) {
          VarDecl * vd = dyn_cast<VarDecl>(init->getSingleDecl());
          if (vd) {
            const Expr * ie = vd->getInit();
            if (ie) {
              loop_parity.expr  = ie;
              loop_parity.value = get_parity_val(loop_parity.expr);
              loop_parity.text  = remove_whitespace(macro.substr(loop_call.length(),
                                                                 std::string::npos));
                
              global.full_loop_text = macro + " " + get_stmt_str(f->getBody());

              // Delete "onsites()" -text
              TheRewriter.RemoveText(CSR);

              handle_full_loop_stmt(f->getBody(), false);
              internal_error = false;
            }
          }
        }
        if (internal_error) {
          reportDiag(DiagnosticsEngine::Level::Error,
                     f->getSourceRange().getBegin(),
                     "\'onsites\'-macro: not a parity type argument" );
          return false;
        }
      }
    }        
    
    return true;
      
  }
  
  // templated function analysis      
  // templated field type assignment, within a template function.  This should not
  // trigger for untemplated loops
  if (global.in_func_template > 0 && isa<BinaryOperator>(s)) {
    BinaryOperator *BO = dyn_cast<BinaryOperator>(s);
    if (BO->isAssignmentOp()) {
      // Check LHS
      std::string lhs_string = BO->getLHS()->getType().getAsString();
      llvm::errs() << " Templ assign string " << lhs_string << "\n";
      if (lhs_string.find(field_type) == 0 ||
          lhs_string.find(field_element_type) == 0) {
        
        state::skip_next = 1;     // next is op type which we know already, skip
        global.template_field_assignment_opcode = BO->getOpcode();
        //handle_lhs_field = true;
          
        std::string assign = TheRewriter.getRewrittenText(BO->getSourceRange());
        llvm::errs() << "The basic stmt content: " << assign << "\n";
      }
    }
    return true;
  }
  
  return true;
}

  

bool MyASTVisitor::VisitFunctionDecl(FunctionDecl *f) {
  // Only function definitions (with bodies), not declarations.
  // also only non-templated functions
  if ((f->getTemplatedKind() == FunctionDecl::TemplatedKind::TK_NonTemplate ||
       f->getTemplatedKind() == FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization )
      && f->hasBody()) {
    global.currentFunctionDecl = f;
    
    Stmt *FuncBody = f->getBody();

    // Type name as string
    QualType QT = f->getReturnType();
    std::string TypeStr = QT.getAsString();

    // Function name
    DeclarationName DeclName = f->getNameInfo().getName();
    std::string FuncName = DeclName.getAsString();

    llvm::errs() << " - Function "<< FuncName << "\n";

    // Add comment before
    std::stringstream SSBefore;
    SSBefore << "// Begin function " << FuncName << " returning " << TypeStr
             << " of template type " << f->getTemplatedKind()
             << "\n";
    SourceLocation ST = f->getSourceRange().getBegin();
    TheRewriter.InsertText(ST, SSBefore.str(), true, true);

    global.location.function = ST;
      
    // And after
    std::stringstream SSAfter;
    SSAfter << "\n// End function " << FuncName;
    ST = FuncBody->getLocEnd().getLocWithOffset(1);
    TheRewriter.InsertText(ST, SSAfter.str(), true, true);
  }
  
  return true;
}


bool MyASTVisitor::VisitFunctionTemplateDecl(FunctionTemplateDecl *tf) {
  if (tf->isThisDeclarationADefinition()) {
    FunctionDecl *f = tf->getTemplatedDecl();
    // Add comment before

    TemplateParameterList *tpl = tf->getTemplateParameters();
      
    std::stringstream SSBefore;
    SSBefore << "// Begin function template " << f->getNameInfo().getName().getAsString()
             << " of template type " << f->getTemplatedKind()
             << " with template params " ;
    for (unsigned i = 0; i < tpl->size(); i++) 
      SSBefore << tpl->getParam(i)->getNameAsString() << " ";
    SSBefore << "\n";
    SourceLocation ST = f->getSourceRange().getBegin();
    TheRewriter.InsertText(ST, SSBefore.str(), true, true);

    global.location.function = ST;
    // start tracking the template
    global.in_func_template = 1;
    
  }
  return true;
}


struct include_struct {
  struct file_id *fp;
  SourceLocation loc;
};

struct file_id {
  FileID FID;
  bool modified;
  std::vector<include_struct> included;
};

std::list<file_id> file_id_list = {};


// Implementation of the ASTConsumer interface for reading an AST produced
// by the Clang parser.
class MyASTConsumer : public ASTConsumer {
public:
  MyASTConsumer(Rewriter &R, ASTContext *C) : Visitor(R,C) { Rewriterp = &R; }

  // Override the method that gets called for each parsed top-level
  // declaration.
  // ANOTHER OPTION: use HandleTranslationUnit?
  bool HandleTopLevelDecl(DeclGroupRef DR) override {

    // do the transformation for non-system files, also included files

    SourceManager &SM = Visitor.getRewriter().getSourceMgr();
    // llvm::errs() << "Processing file " << SM.getFilename((*(DR.begin()))->getLocStart()) << "\n";
    if (!SM.isInSystemHeader((*(DR.begin()))->getLocStart())) {
      // TODO: ensure that we go only through files which are needed!

      state::is_file_modified = false;
      
      // This loop apparently always has only 1 iteration?
      for (DeclGroupRef::iterator b = DR.begin(), e = DR.end(); b != e; ++b) {
        // Traverse the declaration using our AST visitor.
        
        global.location.top = (*b)->getSourceRange().getBegin();  // save this for source location
        Visitor.TraverseDecl(*b);
        // llvm::errs() << "Dumping level " << i++ << "\n";
        if (dump_ast) (*b)->dump();
      }

      FileID FID = SM.getFileID((*(DR.begin()))->getLocStart());
      bool found = false;
      for (file_id & f : file_id_list) {
        if (f.FID == FID) {
          found = true;
          f.modified |= state::is_file_modified;
          break;
        }
      }
      if (!found) {
        // new file to be added
        file_id f;
        f.FID = FID;
        f.modified = state::is_file_modified;
        f.included = {};
        file_id_list.push_back(f);
        llvm::errs() << "New file changed " << SM.getFileEntryForID(FID)->getName() << '\n';
      }
        
    }
    
    return true;
  }

  // Does nothing, apparently..
  // void HandleImplicitImportDecl(ImportDecl *D) override {
  //  llvm::errs() << " ***Import: "
  //              << Rewriterp->getRewrittenText(D->getSourceRange()) << '\n';
  // HandleTopLevelDecl(DeclGroupRef(D));
  // }

private:
  MyASTVisitor Visitor;
  Rewriter * Rewriterp;
};


#ifdef NEED_PP_CALLBACKS
class MyPPCallbacks : public PPCallbacks {
public:

  // This triggers when the preprocessor changes file (#include, exit from it)
  void FileChanged(SourceLocation Loc, FileChangeReason Reason, SrcMgr::CharacteristicKind FileType,
                   FileID PrevFID) {
    SourceManager &SM = myCompilerInstance->getSourceManager();
    llvm::errs() << "FILE CHANGED to "
                 << SM.getFilename(Loc)
                 << '\n';
  }
  // This triggers when range is skipped due to #if (0) .. #endif
  void SourceRangeSkipped(SourceRange Range, SourceLocation endLoc) {
    llvm::errs() << "RANGE skipped\n";
  }
  
};

#endif



// For each source file provided to the tool, a new FrontendAction is created.
class MyFrontendAction : public ASTFrontendAction {
public:
  MyFrontendAction() {}

  virtual bool BeginSourceFileAction(CompilerInstance &CI) override {  
    llvm::errs() << "** Starting operation on source file "+getCurrentFile()+"\n";

#ifdef NEED_PP_CALLBACKS
    // Insert preprocessor callback functions to the stream.  This enables
    // tracking included files, ranges etc.
    Preprocessor &pp = CI.getPreprocessor();
    std::unique_ptr<MyPPCallbacks> callbacks(new MyPPCallbacks());
    pp.addPPCallbacks(std::move(callbacks));
#endif

    file_id_list.clear();
    
    //   // SourceManager &SM = TheRewriter.getSourceMgr();
    //   // llvm::errs() << "** BeginSourceFileAction for: "
    //   //             << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";

    return (true);
  }

  void insert_includes_to_file_buffer(FileID myFID) {
    // find files to be included
    SourceManager &SM = TheRewriter.getSourceMgr();
    for (file_id & f : file_id_list) {
      SourceLocation IL = SM.getIncludeLoc(f.FID);
      if (IL.isValid() && myFID == SM.getFileID(IL)) {
        // file f.FID is included, but do include there first
        insert_includes_to_file_buffer(f.FID);

        // Find now #include "file.h" -stmt (no obv way!)
        SourceRange SR = SM.getExpansionRange(IL).getAsRange();
        SourceLocation e = SR.getEnd();
        std::string includestr;
        for (int i=1; i<1000; i++) {
          SourceLocation b = SR.getBegin().getLocWithOffset(-i);
          includestr = TheRewriter.getRewrittenText(SourceRange(b,e));
          if (includestr.find("#include") != std::string::npos) {
            SR = SourceRange(b,e);
            break;
          }
        }
        // Remove
        TheRewriter.RemoveText(SR);
        // and finally insert
        std::string buf;
        llvm::raw_string_ostream rsos(buf);
        TheRewriter.getEditBuffer(f.FID).write(rsos);
        TheRewriter.InsertText(SR.getBegin(),buf,true,true);
        TheRewriter.InsertText(SR.getBegin(),"// "+includestr+'\n');
      }
    }
  }

  
  // TODO: Must do it so that files without top level decls are used too!!!
  void EndSourceFileAction() override {
    SourceManager &SM = TheRewriter.getSourceMgr();
    llvm::errs() << "** EndSourceFileAction for: " << getCurrentFile() << '\n';
    // << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";

    // Now emit rewritten buffers.  Modified files
    // should be substituted on top of #include's

    insert_includes_to_file_buffer(SM.getMainFileID());
      
    TheRewriter.getEditBuffer(SM.getMainFileID()).write(llvm::outs());
  }

  std::unique_ptr<ASTConsumer> CreateASTConsumer(CompilerInstance &CI,
                                                 StringRef file) override {
    // llvm::errs() << "** Creating AST consumer for: " << file << "\n";
    TheRewriter.setSourceMgr(CI.getSourceManager(), CI.getLangOpts());
    myCompilerInstance = &CI;
    return llvm::make_unique<MyASTConsumer>(TheRewriter, &CI.getASTContext());
  }

private:
  Rewriter TheRewriter;
  // ASTContext  TheContext;
};

int main(int argc, const char **argv) {
  
  CommonOptionsParser op(argc, argv, ToolingSampleCategory);
  ClangTool Tool(op.getCompilations(), op.getSourcePathList());

  // ClangTool::run accepts a FrontendActionFactory, which is then used to
  // create new objects implementing the FrontendAction interface. Here we use
  // the helper newFrontendActionFactory to create a default factory that will
  // return a new MyFrontendAction object every time.
  // To further customize this, we could create our own factory class.
  return Tool.run(newFrontendActionFactory<MyFrontendAction>().get());
}
