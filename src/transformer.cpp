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
// #include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
//#include "llvm/Support/raw_ostream.h"

using namespace clang;
//using namespace clang::driver;
using namespace clang::tooling;

#include "optionsparser.h"

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
  bool is_modified = false;
  bool dump_ast_next = false;
}


static llvm::cl::OptionCategory TransformerCat("Transformer");

// command line options
static llvm::cl::opt<bool>
dump_ast("dump-ast", llvm::cl::desc("Dump AST tree"),
         llvm::cl::cat(TransformerCat));

static llvm::cl::opt<bool>
no_include("noincl",
           llvm::cl::desc("Do not insert \'#include\'-files (for debug)"),
           llvm::cl::cat(TransformerCat));

static llvm::cl::opt<std::string>
dummy_def("D", 
          llvm::cl::value_desc("name"),
          llvm::cl::desc("Define name/symbol for preprocessor"),
          llvm::cl::cat(TransformerCat));

static llvm::cl::opt<std::string>
dummy_incl("I", 
           llvm::cl::desc("Directory for include file search"),
           llvm::cl::value_desc("directory"),
           llvm::cl::cat(TransformerCat));

static llvm::cl::opt<bool>
no_output("no-output",
          llvm::cl::desc("No output file, for syntax check"),
          llvm::cl::cat(TransformerCat));

static llvm::cl::opt<std::string>
output_filename("o",
           llvm::cl::desc("Output file name"),
           llvm::cl::value_desc("name"),
           llvm::cl::cat(TransformerCat));

static llvm::cl::opt<bool>
kernel("vanilla-kernel",
       llvm::cl::desc("Generate kernels"),
       llvm::cl::cat(TransformerCat));

static llvm::cl::opt<bool>
vanilla("vanilla",
        llvm::cl::desc("Generate loops in place"),
        llvm::cl::cat(TransformerCat));



// local global vars
global_state global;
// and global loop_parity
loop_parity_struct loop_parity;

static const std::string field_element_type = "field_element<";
static const std::string field_type = "field<";

static codetype target;

// global lists for functions
// TODO: THESE SHOULD PROBABLY BE CHANGED INTO vectors,
// but they contain pointers to list elements.  pointers to vector elems are not good!
std::list<field_ref> field_ref_list = {};
std::list<field_info> field_info_list = {};
std::list<var_info> var_info_list = {};
std::list<var_decl> var_decl_list = {};

std::vector<unsigned> remove_expr_list = {};

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
  
  
// // Define a pragma handler for #pragma heLpp
// // NOTE: This is executed before AST analysis
// class heLppPragmaHandler : public PragmaHandler {
// public:
//   heLppPragmaHandler() : PragmaHandler("heLpp") { }
//   void HandlePragma(Preprocessor &PP, PragmaIntroducerKind Introducer,
//                     Token &PragmaTok) {
//     // Handle the pragma
//    
//     llvm::errs() << "Got the pragma! name " << getName() << " Token " << PragmaTok.getName() << '\n';
//
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
  
  // if state::skip_children > 0 we'll skip all until return to level up
  if (state::skip_children > 0) state::skip_children++;
    
  // go via the original routine...
  if (!state::skip_children) RecursiveASTVisitor<MyASTVisitor>::TraverseStmt(S);

  if (state::skip_children > 0) state::skip_children--;
      
  return true;
}


// Similarly, implement skip_children for decl traversal
bool MyASTVisitor::TraverseDecl(Decl *D) {
  if (state::skip_next > 0) {
    state::skip_next--;
    return true;
  }

  // if state::skip_children > 0 we'll skip all until return to level up
  if (state::skip_children > 0) state::skip_children++;
    
  // go via the original routine...
  if (!state::skip_children) RecursiveASTVisitor<MyASTVisitor>::TraverseDecl(D);

  if (state::skip_children > 0) state::skip_children--;

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
bool MyASTVisitor::is_field_parity_expr(Expr *E) {
  E = E->IgnoreParens();
  CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(E);
  if (OC &&
      strcmp(getOperatorSpelling(OC->getOperator()),"[]") == 0 && 
      is_field_expr(OC->getArg(0))) {
    std::string s = get_expr_type(OC->getArg(1)); 
    if (s == "parity" || s == "parity_plus_direction") {
      llvm::errs() << " <<<Parity type " << get_expr_type(OC->getArg(1)) << '\n';
      return true;
    }
  } else {
    // This is for templated expressions
    // for some reason, expr a[X] "getBase() gives X, getIdx() a...
    if (ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(E)) {
      if (is_field_expr(ASE->getLHS())) {
        llvm::errs() << " FP: and field\n";
        std::string s = get_expr_type(ASE->getRHS());
        if (s == "parity" || s == "parity_plus_direction") {
          llvm::errs() << " <<<Parity type " << get_expr_type(ASE->getRHS()) << '\n';
          return true;
        }
      }
    }
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
  
bool MyASTVisitor::check_field_ref_list() {

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
  
bool MyASTVisitor::handle_field_parity_expr(Expr *e, bool is_assign) {
    
  e = e->IgnoreParens();
  field_ref lfe;
  if (CXXOperatorCallExpr *OC = dyn_cast<CXXOperatorCallExpr>(e)) {
    lfe.fullExpr   = OC;
    // take name 
    lfe.nameExpr   = OC->getArg(0);
    lfe.parityExpr = OC->getArg(1);
  } else if (ArraySubscriptExpr * ASE = dyn_cast<ArraySubscriptExpr>(e)) {
    // In template definition

    lfe.fullExpr   = ASE;
    lfe.nameExpr   = ASE->getLHS();
    lfe.parityExpr = ASE->getRHS();
  } else {
    llvm::errs() << "Internal error 3\n";
    exit(0);
  }
  
  lfe.nameInd    = Buf.markExpr(lfe.nameExpr); 
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
      
  return(true);
}


reduction get_reduction_type(bool is_assign, 
                             std::string & assignop, 
                             var_info & vi) {
  if (is_assign && (!vi.is_loop_local)) {
    if (assignop == "+=") return reduction::SUM;
    if (assignop == "*=") return reduction::PRODUCT;
  }
  return reduction::NONE;
}  

/// This processes references to non-field variables within field loops

void MyASTVisitor::handle_var_ref(DeclRefExpr *DRE,
                                  bool is_assign,
                                  std::string &assignop) {

  
  if (isa<VarDecl>(DRE->getDecl())) {
    auto decl = dyn_cast<VarDecl>(DRE->getDecl());
    var_ref vr;
    vr.ref = DRE;
    vr.ind = Buf.markExpr(DRE);
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
      // This is somehow needed for printing type without "class" id
      PrintingPolicy pp(Context->getLangOpts());
      vi.type = DRE->getType().getUnqualifiedType().getAsString(pp);

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


/// is the stmt pointing now to assignment
bool MyASTVisitor::is_assignment_expr(Stmt * s, std::string * opcodestr) {
  if (CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s))
    if (OP->isAssignmentOp()) {
      if (opcodestr)
        *opcodestr = getOperatorSpelling(OP->getOperator());
      return true;
    }
  
  if (BinaryOperator *B = dyn_cast<BinaryOperator>(s))
    if (B->isAssignmentOp()) {
      if (opcodestr)
        *opcodestr = B->getOpcodeStr();
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

  

// start dealing with the loop statement
bool MyASTVisitor::handle_full_loop_stmt(Stmt *ls, bool field_parity_ok ) {
  // init edit buffer
  Buf.create( &TheRewriter, ls );
          
  field_ref_list.clear();
  var_info_list.clear();
  var_decl_list.clear();
  remove_expr_list.clear();
  global.location.loop = ls->getSourceRange().getBegin();
  
  state::accept_field_parity = field_parity_ok;
    
  // the following is for taking the parity from next elem
  state::scope_level = 0;
  state::in_loop_body = true;
  TraverseStmt(ls);
  state::in_loop_body = false;

  // Remove exprs which we do not want
  for (unsigned i : remove_expr_list) Buf.remove(i);
  
  // check and analyze the field expressions
  check_field_ref_list();
  check_var_info_list();
  
  generate_code(ls, target);
  
  Buf.clear();
          
  // Emit the original command as a commented line
  TheRewriter.InsertText(ls->getSourceRange().getBegin(),
                         comment_string(global.full_loop_text) + "\n",
                         false,true);

  global.full_loop_text = "";

  // don't go again through the arguments
  state::skip_children = 1;

  state::is_modified = true;
  
  return true;
}


// This handles statements within field loops
bool MyASTVisitor::handle_loop_body_stmt(Stmt * s) {

  // This keeps track of the assignment to field
  // must remember the set value across calls
  static bool is_assignment = false;
  static std::string assignop;
 
  // Need to recognize assignments lf[X] =  or lf[X] += etc.
  // And also assignments to other vars: t += norm2(lf[X]) etc.
   if (is_assignment_expr(s,&assignop)) {
    is_assignment = true;
    // next visit here will be to the assigned to variable
    return true;
  } 
  
  // catch then expressions
      
  if (Expr *E = dyn_cast<Expr>(s)) {
    
    // Not much to do with const exprs
    // if (E->isCXX11ConstantExpr(*Context, nullptr, nullptr)) {
    //   state::skip_children = 1;
    //   return true;
    // }
    
    //if (is_field_element_expr(E)) {
      // run this expr type up until we find field variable refs
    if (is_field_parity_expr(E)) {
      // Now we know it is a field parity reference
      // get the expression for field name
          
      handle_field_parity_expr(E, is_assignment);
      is_assignment = false;  // next will not be assignment
      // (unless it is a[] = b[] = c[], which is OK)

      state::skip_children = 1;
      return true;
    }

    if (is_field_expr(E)) {
      // field without [parity], bad usually (TODO: allow  scalar func(field)-type?)
      reportDiag(DiagnosticsEngine::Level::Error,
                 E->getSourceRange().getBegin(),
                 "Field expressions without [..] not allowed within field loop");
      state::skip_children = 1;  // once is enough
      return true;
    }

    if (DeclRefExpr *DRE = dyn_cast<DeclRefExpr>(E)) {
      if (isa<VarDecl>(DRE->getDecl())) {
        // now it should be var ref non-field
      
        handle_var_ref(DRE,is_assignment,assignop);
        is_assignment = false;
      
        state::skip_children = 1;
        llvm::errs() << "Variable ref: "
                     << TheRewriter.getRewrittenText(E->getSourceRange()) << '\n';

        state::skip_children = 1;
        return true;
      }
      // TODO: function ref?
    }



#if 1

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
      //var_ref_list.push_back( handle_var_ref(E) );

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


// This is a horrible hack to enable control statements such as
//  transformer_control("dump-ast");
// Should use #pragma, but pragma handling is very messy
 
bool MyASTVisitor::handle_control_stmt(Stmt *s) {
  if (CallExpr * CE = dyn_cast<CallExpr>(s)) 
    if (FunctionDecl * FD = CE->getDirectCallee())
      if (FD->getNameInfo().getName().getAsString() == "transformer_control")
        if (CE->getNumArgs() == 1) {
          Expr * E = CE->getArg(0);
          E = E->IgnoreParenCasts();
          if (StringLiteral *SL = dyn_cast<StringLiteral>(E)) {
            std::string command = SL->getString();
            
            if (command == "dump-ast") state::dump_ast_next = true;
            
            if (Buf.isOn()) remove_expr_list.push_back(Buf.markExpr(CE));
            TheRewriter.InsertText(s->getSourceRange().getBegin(),"//-- ",true,true);
              
            state::skip_children = 1;
            return true;
          }
        }
  return false;
}


// VisitStmt is called for each statement in AST.  Thus, when traversing the
// AST or part of it we always start from here

bool MyASTVisitor::VisitStmt(Stmt *s) {

  if (state::dump_ast_next) {
    llvm::errs() << "**** Dumping statement:\n" + get_stmt_str(s)+'\n';
    s->dump();
    state::dump_ast_next = false;
  }

  if (handle_control_stmt(s)) return true;
  
  // Entry point when inside field[par] = .... body
  if (state::in_loop_body) {
    return handle_loop_body_stmt(s);
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

        llvm::errs() << "macro str " << macro << '\n';
        
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

                                   
  //  Starting point for fundamental operation
  //  field[par] = .... - non-templated version
  
  // isStmtWithSemi(s);


  CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(s);
  if (OP && OP->isAssignmentOp() && is_field_parity_expr(OP->getArg(0))) {
    // now we have a[par] += ...  -stmt.  Arg(0) is
    // the lhs of the assignment

    SourceRange full_range = getRangeWithSemi(OP,false);
    global.full_loop_text = TheRewriter.getRewrittenText(full_range);
        
    handle_full_loop_stmt(OP, true);
    return true;
  } 
  // now the above within template
  BinaryOperator *BO = dyn_cast<BinaryOperator>(s);
  if (BO && BO->isAssignmentOp() && is_field_parity_expr(BO->getLHS())) {
    SourceRange full_range = getRangeWithSemi(BO,false);
    global.full_loop_text = TheRewriter.getRewrittenText(full_range);        
    handle_full_loop_stmt(BO, true);
    return true;    
  }

  return true;
}

  

bool MyASTVisitor::VisitFunctionDecl(FunctionDecl *f) {
  // Only function definitions (with bodies), not declarations.
  // also only non-templated functions

  if (state::dump_ast_next) {
    llvm::errs() << "**** Dumping funcdecl:\n";
    f->dump();
    state::dump_ast_next = false;
  }

  // don't go through instantiations (for now!)
  // analyse templates directly
  if (f->isTemplateInstantiation()) {
    state::skip_children = 1;
    return true;
  }

  
  if ((f->getTemplatedKind() == FunctionDecl::TemplatedKind::TK_NonTemplate
       // || f->getTemplatedKind() == FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization
       ) && f->hasBody()) {
    global.currentFunctionDecl = f;
    
    Stmt *FuncBody = f->getBody();

    // Type name as string
    QualType QT = f->getReturnType();
    std::string TypeStr = QT.getAsString();

    // Function name
    DeclarationName DeclName = f->getNameInfo().getName();
    std::string FuncName = DeclName.getAsString();

    // llvm::errs() << " - Function "<< FuncName << "\n";

    // Add comment before
    std::stringstream SSBefore;
    SSBefore << "// Begin function " << FuncName << " returning " << TypeStr
             << " of template type " << f->getTemplatedKind()
             << "\n";
    SourceLocation ST = f->getSourceRange().getBegin();
    TheRewriter.InsertText(ST, SSBefore.str(), true, true);

    llvm::errs() << "Function decl " << FuncName << '\n';

    global.location.function = ST;
      
    // And after
    std::stringstream SSAfter;
    SSAfter << "\n// End function " << FuncName;
    ST = FuncBody->getLocEnd().getLocWithOffset(1);
    TheRewriter.InsertText(ST, SSAfter.str(), true, true);

    //TraverseStmt(FuncBody);
    //state::skip_children = 1;
      
  }
  
  return true;
}



bool MyASTVisitor::VisitFunctionTemplateDecl(FunctionTemplateDecl *tf) {

  if (state::dump_ast_next) {
    llvm::errs() << "**** Dumping funcdecl:\n";
    tf->dump();
    state::dump_ast_next = false;
  }

  if (tf->isThisDeclarationADefinition()) {
    FunctionDecl *f = tf->getTemplatedDecl();
    // Add comment before

    global.function_tpl = tf->getTemplateParameters();
    global.currentFunctionDecl = f;
        
    std::stringstream SSBefore;
    SSBefore << "// Begin function template " << f->getNameInfo().getName().getAsString()
             << " of template type " << f->getTemplatedKind()
             << " with template params " ;
    for (unsigned i = 0; i < global.function_tpl->size(); i++) 
      SSBefore << global.function_tpl->getParam(i)->getNameAsString() << " ";
    SSBefore << "\n";
    SourceLocation ST = tf->getSourceRange().getBegin();
    TheRewriter.InsertText(ST, SSBefore.str(), true, true);

    global.location.function = ST;
    // start tracking the template
    
    global.in_func_template = true;
    TraverseDecl(f);
    global.in_func_template = false;

    state::skip_children = 1;
    
  }
  return true;
}


// Entry point for class templates
// This (may) be needed for template class methods
bool MyASTVisitor::VisitClassTemplateDecl(ClassTemplateDecl *D) {
  if (D->isThisDeclarationADefinition()) {

    // save template params in a list, for templates within templates .... ugh!
    global.class_tpl.push_back(D->getTemplateParameters());

    // this block for debugging
    TemplateParameterList * tplp = D->getTemplateParameters();
    std::stringstream SSBefore;
    SSBefore << "// Begin template class "
             << " with template params " ;
    for (unsigned i = 0; i < tplp->size(); i++) 
      SSBefore << tplp->getParam(i)->getNameAsString() << " ";
    SSBefore << "\n";
    SourceLocation ST = D->getSourceRange().getBegin();
    TheRewriter.InsertText(ST, SSBefore.str(), true, true);
    // end block
    
    global.in_class_template = true;
    TraverseDecl(D->getTemplatedDecl());
    global.in_class_template = false;

    state::skip_children = 1;

    // Remove template param list
    global.class_tpl.pop_back();
    
  }
  return true;
}


#if 1
bool MyASTVisitor::VisitCXXMethodDecl(CXXMethodDecl *method) {
  // This comes after VisitFunctionDecl for methods -- this does really nothing here

  if (method->isThisDeclarationADefinition()) {
    // FunctionDecl *f = method->getTemplatedDecl();
    // Add comment before

    std::stringstream SSBefore;
    SSBefore << "// Begin method "
             << method->getNameInfo().getName().getAsString();
      //<< " of template type " << f->getTemplatedKind()
      //     << " with template params " ;
    //for (unsigned i = 0; i < global.tpl->size(); i++) 
    //   SSBefore << global.tpl->getParam(i)->getNameAsString() << " ";
    SSBefore << "\n";
    SourceLocation ST = method->getSourceRange().getBegin();
    TheRewriter.InsertText(ST, SSBefore.str(), true, true);

    llvm::errs() << "Method decl " << method->getNameInfo().getName().getAsString() << '\n';
    
  }
  return true;
}
#endif


// This struct will be used to keep track of #include-chains.
// We also use 

std::vector<FileID> file_id_list = {};

// Tiny utility to search for the list
// cache previous file_id to speed up

bool search_fid(const FileID FID) {
  for (const FileID f : file_id_list) {
    if (f == FID) return true;
  }
  return false;
}


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

      state::is_modified = false;
      
      // This loop apparently always has only 1 iteration?
      for (DeclGroupRef::iterator b = DR.begin(), e = DR.end(); b != e; ++b) {
        // Traverse the declaration using our AST visitor.
        
        global.location.top = (*b)->getSourceRange().getBegin();  // save this for source location
        Visitor.TraverseDecl(*b);
        // llvm::errs() << "Dumping level " << i++ << "\n";
        if (dump_ast) {
          if (!no_include || SM.isInMainFile((*(DR.begin()))->getLocStart()))
            (*b)->dump();
        }
      }
      // We keep track here only of files which were touched
      if (state::is_modified) {
        FileID FID = SM.getFileID((*(DR.begin()))->getLocStart());
        if (search_fid(FID) == false) {
          // new file to be added
          file_id_list.push_back(FID);
          llvm::errs() << "New file changed " << SM.getFileEntryForID(FID)->getName() << '\n';
        }
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


// #define NEED_PP_CALLBACKS
#ifdef NEED_PP_CALLBACKS

class MyPPCallbacks : public PPCallbacks {
public:
  SourceLocation This_hashloc;
  std::string This_name;

  // This hook is called when #include (or #import) is processed
  void InclusionDirective(SourceLocation HashLoc,
                          const Token & IncludeTok,
                          StringRef FileName,
                          bool IsAngled,
                          CharSourceRange FilenameRange,
                          const FileEntry * File,
                          StringRef SearchPath,
                          StringRef RelativePath,
                          const Module * Imported,
                          SrcMgr::CharacteristicKind FileType) { }
  
  // This triggers when the preprocessor changes file (#include, exit from it)
  // Use this to track the chain of non-system include files
  void FileChanged(SourceLocation Loc, FileChangeReason Reason, SrcMgr::CharacteristicKind FileType,
                   FileID PrevFID) {
    SourceManager &SM = myCompilerInstance->getSourceManager();
    if (Reason == PPCallbacks::EnterFile &&
        FileType == SrcMgr::CharacteristicKind::C_User &&
        Loc.isValid() &&
        !SM.isInSystemHeader(Loc) &&
        !SM.isInMainFile(Loc) ) {

      llvm::errs() << "FILE CHANGED to " << SM.getFilename(Loc) << '\n';
    }
  }

  // This triggers when range is skipped due to #if (0) .. #endif
  void SourceRangeSkipped(SourceRange Range, SourceLocation endLoc) {
    // llvm::errs() << "RANGE skipped\n";
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
    for (FileID f : file_id_list) {
      SourceLocation IL = SM.getIncludeLoc(f);
      if (IL.isValid() && myFID == SM.getFileID(IL)) {
        // file f is included, but do #include there first
        insert_includes_to_file_buffer(f);

        // Find now '#include "file.h"' -stmt (no obv way!)
        SourceRange SR = SM.getExpansionRange(IL).getAsRange();
        std::string includestr = TheRewriter.getRewrittenText(SR);
        SourceLocation e = SR.getEnd();
        SourceLocation b = SR.getBegin();
        for (int i=1; i<100; i++) {
          const char * p = SM.getCharacterData(b.getLocWithOffset(-i));
          if (p && *p == '#' && strncmp(p,"#include",8) == 0) {
            SR = SourceRange(b.getLocWithOffset(-i),e);
            break;
          }
        }
        // Remove
        TheRewriter.RemoveText(SR);

        // and finally insert
        SourceRange r(SM.getLocForStartOfFile(f),SM.getLocForEndOfFile(f));
        TheRewriter.InsertText(SR.getBegin(),
                               "// start include "+includestr
                               + "---------------------------------\n"
                               + TheRewriter.getRewrittenText(r) +
                               "// end include "+includestr
                               + "---------------------------------\n",
                               true,false);
        
      }
    }
  }

  // check and add FileID's for files in #include chains if needed
  void check_include_path(const FileID FID) {
    SourceManager &SM = TheRewriter.getSourceMgr();
    SourceLocation IL = SM.getIncludeLoc(FID);
    if (IL.isValid()) {
      FileID FID_up = SM.getFileID(IL);
      if (!search_fid(FID_up)) {
        file_id_list.push_back(FID_up);
        if (FID_up != SM.getMainFileID()) check_include_path(FID_up);
      }
    }
  }

  // TODO: Must do it so that files without top level decls are used too!!!
  void EndSourceFileAction() override {
    SourceManager &SM = TheRewriter.getSourceMgr();
    llvm::errs() << "** EndSourceFileAction for: " << getCurrentFile() << '\n';
    // << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";

    // Now emit rewritten buffers.

    if (!no_include) {

      // Modified files should be substituted on top of #include -directives
      // first, ensure that the full include chain is present in file_id_list
      // Use iterator here, because the list can grow!

      for ( FileID f : file_id_list ) {
        check_include_path(f);
      }

      insert_includes_to_file_buffer(SM.getMainFileID());
    }
    
    if (!no_output)
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


// Check if the cmdline has -I<include> or -D<define> -
// arguments and move these after -- if that exists on the command line.
// Clang's optionparser expects these "generic compiler and linker"
// args to be after --
// return value new argc
int rearrange_cmdline(int argc, const char **argv, const char **av) {

  bool found_ddash = false;
  av[argc+1] = nullptr;  // I read somewhere that in c++ argv[argc] = 0
  static char s[3] = "--";   // needs to be static because ptrs
  int ddashloc = 0;

  for (int i=0; i<argc; i++) {
    av[i] = argv[i];
    if (strcmp(av[i],s) == 0) {
      found_ddash = true;
      ddashloc = i;
    }
  }
  if (!found_ddash) {
    // add ddash, does not hurt in any case
    av[argc] = s;
    ddashloc = argc;
    argc++;
  }

  // now find -I and -D -options and move them after --
  for (int i=0; i<ddashloc; i++) {
    if (i < ddashloc-1 && (strcmp(av[i],"-D") == 0 || strcmp(av[i],"-I") == 0)) {
      // type -D define
      const char * a1 = av[i];
      const char * a2 = av[i+1];
      for (int j=i+2; j<argc; j++) av[j-2] = av[j];
      av[argc-2] = a1;
      av[argc-1] = a2;
      ddashloc -= 2;
    } else if (strncmp(av[i],"-D",2) == 0 || strncmp(av[i],"-I",2) == 0) {
      // type -Ddefine
      const char * a1 = av[i];
      for (int j=i+1; j<argc; j++) av[j-1] = av[j];
      av[argc-1] = a1;
      ddashloc--;
    }      
  }

  return argc;
}

void get_target_struct(codetype & target) {
  if (kernel) target.kernelize = true;
  else target.kernelize = false;
}

int main(int argc, const char **argv) {

  // TODO: clang CommandLine.cpp/.h has strange category and help
  // msg handling, should we get rid of it?

  // av takes over from argv
  const char **av = new const char *[argc+2];
  argc = rearrange_cmdline(argc, argv, av);
  
  OptionsParser op(argc, av, TransformerCat);
  ClangTool Tool(op.getCompilations(), op.getSourcePathList());
  
  // We have command line args, possibly do something with them
  get_target_struct(target);

  // ClangTool::run accepts a FrontendActionFactory, which is then used to
  // create new objects implementing the FrontendAction interface. Here we use
  // the helper newFrontendActionFactory to create a default factory that will
  // return a new MyFrontendAction object every time.
  // To further customize this, we could create our own factory class.
  return Tool.run(newFrontendActionFactory<MyFrontendAction>().get());
}
