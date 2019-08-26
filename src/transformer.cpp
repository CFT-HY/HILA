//------------------------------------------------------------------------------
// Transformer tools to convert "lattice loops" into
// hardware-dependent "kernels".
//
// Uses Clang RecursiveASTVisitor and Rewriter 
// interfaces
//
// Kari Rummukainen 2017-19
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


#include "transformer.h"
#include "optionsparser.h"
#include "stringops.h"
#include "srcbuf.h"
#include "myastvisitor.h"
#include "specialization_db.h"

// collection of variables holding the state of parsing
namespace state {
  unsigned skip_children = 0;
  unsigned scope_level = 0;
  int skip_next = 0;
  bool in_loop_body = false;
  bool accept_field_parity = false;
  bool loop_found = false;
  bool dump_ast_next = false;
  bool compile_errors_occurred = false;
};

static llvm::cl::OptionCategory TransformerCat(program_name);

namespace cmdline {

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
  spec_no_db("spec:no-db",
             llvm::cl::desc("Do not use generated specialization database, include them (potentially unsafe)"),
             llvm::cl::cat(TransformerCat));

  static llvm::cl::opt<bool>
  spec_inline("spec:inline",
              llvm::cl::desc("Mark generated specializations inline"),
              llvm::cl::cat(TransformerCat));
  
  static llvm::cl::opt<bool>
  no_output("no-output",
            llvm::cl::desc("No output file, for syntax check"),
            llvm::cl::cat(TransformerCat));
  
  static llvm::cl::opt<bool>
  syntax_only("syntax-only",
              llvm::cl::desc("Same as no-output"),
              llvm::cl::cat(TransformerCat));
  
  static llvm::cl::opt<std::string>
  output_filename("o",
                  llvm::cl::desc("Output file (default: <file>.cpt, write to stdout: -o - "),
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

};

// local global vars
global_state global;
// and global loop_parity
loop_parity_struct loop_parity;

static ClassTemplateDecl * field_decl = nullptr;   // Ptr to field primary def in AST
static ClassTemplateDecl * field_storage_type_decl = nullptr;   // Ptr to field primary def in AST

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

std::vector<Expr *> remove_expr_list = {};

// take global CI just in case
CompilerInstance *myCompilerInstance;


// utility functions

// TODO: make these more robust!!!
bool MyASTVisitor::is_field_element_expr(Expr *E) {
  return( E && E->getType().getAsString().find(field_element_type) != std::string::npos);
}

bool MyASTVisitor::is_field_expr(Expr *E) {
  return( E && E->getType().getAsString().find(field_type) != std::string::npos);
}

bool MyASTVisitor::is_field_decl(ValueDecl *D) {
  return( D && D->getType().getAsString().find(field_type) != std::string::npos);
}

SourceLocation MyASTVisitor::getSourceLocationAfterNewLine( SourceLocation l ) {
  SourceManager &SM = TheRewriter.getSourceMgr();
  for (int i=0; i<10000; i++) {
    bool invalid = false;
    const char * c = SM.getCharacterData(l.getLocWithOffset(i),&invalid);
    if (invalid) {
      // no new line found in buffer.  return previous loc, could be false!
      llvm::errs() << program_name + ": no new line found in buffer, internal error\n";
      return( l.getLocWithOffset(i-1) );
    }
    if (*c == '\n') return( l.getLocWithOffset(i+1) );
  }
  return l;
}


#if 0  
// Define a pragma handler for #pragma heLpp
// NOTE: This is executed before AST analysis
class heLppPragmaHandler : public PragmaHandler {
  public:
    heLppPragmaHandler() : PragmaHandler("heLpp") { }
    void HandlePragma(Preprocessor &PP, PragmaIntroducerKind Introducer,
                      Token &PragmaTok) {
     // Handle the pragma
    
     llvm::errs() << "Got the pragma! name " << getName() << " Token " << PragmaTok.getName() << '\n';

   }
};

static PragmaHandlerRegistry::Add<heLppPragmaHandler> Y("heLpp","heL pragma description");
#endif
  
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
    if (l.is_changed && l.dir_list.size() > 0) {
      if (loop_parity.value == parity::all) {
        // There's error, find culprits
        for (field_ref * p : l.ref_list) {
          if (p->dirExpr != nullptr && !p->is_changed) {
            reportDiag(DiagnosticsEngine::Level::Error,
                       p->parityExpr->getSourceRange().getBegin(),
                       "Accessing field '%0' undefined when assigning to '%1' with parity ALL, flagging as error",
                       get_stmt_str(p->fullExpr).c_str(),
                       l.old_name.c_str());
            no_errors = false;
          }
        }

        for (field_ref * p : l.ref_list) {
          if (p->is_changed && p->dirExpr == nullptr) {
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
    exit(1);
  }
  
  //lfe.nameInd    = writeBuf->markExpr(lfe.nameExpr); 
  //lfe.parityInd  = writeBuf->markExpr(lfe.parityExpr);
  
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
  // Buf.create( &TheRewriter, ls );
          
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
  state::skip_children = 1;

  state::loop_found = true;
  
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
                 "Cannot declare field variables within field loops");
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
            
            if (writeBuf->isOn()) remove_expr_list.push_back(CE);
            writeBuf->insert(s->getSourceRange().getBegin(),"//-- ",true,false);
              
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
              loop_parity.text  = remove_initial_whitespace(macro.substr(loop_call.length(),
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
  if (0 && f->isTemplateInstantiation()) {
    
    state::skip_children = 1;
    return true;
  }

  
  if (//(f->getTemplatedKind() == FunctionDecl::TemplatedKind::TK_NonTemplate
       // || f->getTemplatedKind() == FunctionDecl::TemplatedKind::TK_FunctionTemplateSpecialization
       //) &&
      f->hasBody()) {
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
    writeBuf->insert(ST, SSBefore.str(), true,true);
    // TheRewriter.InsertText(ST, SSBefore.str(), true, true);

    llvm::errs() << "Function decl " << FuncName << '\n';

    global.location.function = ST;
      
    // And after
    std::stringstream SSAfter;
    SSAfter << "\n// End function " << FuncName;
    ST = FuncBody->getEndLoc().getLocWithOffset(1);
    // TheRewriter.InsertText(ST, SSAfter.str(), true, true);
    writeBuf->insert(ST, SSAfter.str(), true, true);

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

    TemplateParameterList * tpl = tf->getTemplateParameters();
    global.function_tpl = tpl;
    global.currentFunctionDecl = f;

    std::vector<std::string> template_params = {};
    
    // std::stringstream SSBefore;
    // SSBefore << "// Begin function template " << f->getNameInfo().getName().getAsString()
    //          << " of template type " << f->getTemplatedKind()
    //          << " with template params " ;
    for (unsigned i = 0; i < tpl->size(); i++) {
      // SSBefore << tpl->getParam(i)->getNameAsString() << " ";
      template_params.push_back(tpl->getParam(i)->getNameAsString());
    }
      
    // SSBefore << "\n";

    SourceLocation ST = tf->getSourceRange().getBegin();

    global.location.function = ST;
    
    // start tracking the template
    
    global.in_func_template = true;
    // Let's not traverse the template any more, but do the instantiations below!
    // TraverseDecl(f);   
    global.in_func_template = false;

    SourceLocation insertloc = getSourceLocationAfterNewLine(tf->getSourceRange().getEnd());

    // save write buffer
    srcBuf * writeBuf_saved = writeBuf;
    
    // check specializations -- for clang 8 the spec_iterator is type FunctionDecl,
    // for documentation (clang 10) it is of type FunctionTemplateSpecializationInfo !!!
    // TODO: prepare for API change!!!
    for (auto f : tf->specializations() ) {
      // TODO: partial specializations?
      if (f->getTemplateSpecializationKind() ==
          TemplateSpecializationKind::TSK_ImplicitInstantiation) {
        auto tal = f->getTemplateSpecializationArgs();

        // 1st test if we need to do anything
        // buffer the template and traverse through it
        srcBuf templateBuf( &TheRewriter, tf );
        writeBuf = &templateBuf;

        bool found = state::loop_found;
        state::loop_found = false;  // use this to flag

        TraverseStmt(f->getBody());
        templateBuf.clear();

        if (state::loop_found) {
          // OK, need to convert -- clear buffer and refill it

          if (create_function_specializations) {
            
            templateBuf.create( &TheRewriter, tf );
            writeBuf = &templateBuf;

            // set appropriate return type
            // llvm::errs() << "Func ret qualtype: " << f->getReturnType().getAsString() << '\n';
            templateBuf.replace( f->getReturnTypeSourceRange(), 
                                 remove_class_from_type(f->getReturnType().getAsString()) );

            if (cmdline::spec_inline && !f->isInlineSpecified()) {
              templateBuf.insert( f->getReturnTypeSourceRange().getBegin(), "inline ", true, false);
            }
            
            // generate name<template params>
            std::string name = f->getNameInfo().getAsString() + "<";

            std::vector<std::string> template_args = {};
        
            //SourceLocation 
            //handle_specialization_args( 
            std::string info = "// ---- Generated specialization with type args ";
            bool first = true;
            for (int i=0; i<tal->size(); i++) {
              if (tal->get(i).getKind() == TemplateArgument::ArgKind::Type) {
                std::string typestr = remove_class_from_type(tal->get(i).getAsType().getAsString());

                info += typestr + " ";
                if (!first) name += ", ";
                first = false;

                name += typestr;
                template_args.push_back(typestr);

                // wipe type params in "template < ... >", incl. comma
                templateBuf.remove_with_comma(tpl->getParam(i)->getSourceRange());
            
              } else {
                // Now arg not of Type, mark with empty - not matched
                template_args.push_back("");
                template_params.at(i) = "";   // TODO: not sure if this is a good idea
              }
            }
        
            name += ">";

            // if (f->isExplicitSpecialization()) SSBefore << " explicit";
            info += "\n";

            if (template_params.size() != template_args.size()) {
              reportDiag(DiagnosticsEngine::Level::Error,
                         tf->getSourceRange().getBegin(),
                         "Template param/arg number mismatch!" );
              return(false);
            }

            // function name -> name<..>
            templateBuf.replace( f->getNameInfo().getSourceRange(), name );

            // Now parameters of the function
            for ( auto par : f->parameters() ) {          
              templateBuf.replace_tokens(par->getSourceRange(),
                                         template_params, template_args);
            }

            // Now we have the function name ready - check if this
            // has already been generated
            SourceRange decl_sr = get_templatefunc_decl_range(tf,f);
            std::string wheredefined = "";
            if (f->isInlineSpecified() || cmdline::spec_inline || cmdline::spec_no_db
                || !in_specialization_db(templateBuf.get(decl_sr), wheredefined)) {
              
              // replace template params in func body
              templateBuf.replace_tokens(f->getBody()->getSourceRange(),
                                         template_params, template_args);

              // traverse again, with right arguments
              TraverseStmt(f->getBody());

              templateBuf.insert(0,info,true,false);

              writeBuf_saved->insert( insertloc,
                                      templateBuf.dump() + "\n//---------\n",
                                      false, false);
            } else {
              // just insert declaration, defined on another compilation unit
              writeBuf_saved->insert( insertloc,
                                      "// Generated specialization declaration, defined in compilation unit\n// "
                                      + wheredefined + "\n"
                                      + templateBuf.get(decl_sr)
                                      + ";\n//---------\n",
                                      false, false);
            }
              
          } else {
            // Now attempt to modify the template directly, no specializations generated
            // NOT DONE!
          }

            
        } else {
          state::loop_found = found;
        }

        templateBuf.clear();
      }
    }
    
    writeBuf = writeBuf_saved;

    // SSBefore << '\n';
    // writeBuf->insert(ST,SSBefore.str(),true,true);
   
    state::skip_children = 1;
    
  }
  return true;
}

// locate range of specialization "template< ..> .. func<...>( ... )"
// tf is ptr to template, and f to instantiated function
SourceRange MyASTVisitor::get_templatefunc_decl_range(FunctionTemplateDecl *tf,
                                                      FunctionDecl *f) {
  SourceLocation a = tf->getSourceRange().getBegin();
  int n = f->getNumParams()-1;
  SourceLocation b = f->getParamDecl(n-1)->getSourceRange().getEnd();

  SourceManager &SM = TheRewriter.getSourceMgr();
  for (int i=1; i<50; i++) {
    const char * p = SM.getCharacterData(b.getLocWithOffset(i));
    if (p && *p == ')') {
      b = b.getLocWithOffset(i);
      break;
    }
  }

  SourceRange r(a,b);
  return r;
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
             << D->getNameAsString()
             << " with template params " ;
    for (unsigned i = 0; i < tplp->size(); i++) 
      SSBefore << tplp->getParam(i)->getNameAsString() << " ";
    SSBefore << "\n";
    SourceLocation ST = D->getSourceRange().getBegin();
    // TheRewriter.InsertText(ST, SSBefore.str(), true, true);
    writeBuf->insert(ST, SSBefore.str(), true, true);
    // end block
    
    global.in_class_template = true;
    TraverseDecl(D->getTemplatedDecl());
    global.in_class_template = false;

    state::skip_children = 1;

    // Remove template param list
    global.class_tpl.pop_back();

    if (D->getNameAsString() == "field") {
      // Now we have the main field template definition.
      // Check the field template instantiations/specializations
      // This can be done at this point only if we use
      // handleTranslationUnit() as the entry point

      // Let us save this Decl for a revisit after AST is complete
      handle_field_specializations(D);

    } else if (D->getNameAsString() == "field_storage_type") {
      field_storage_type_decl = D;
    }
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

    std::string typestr = args.get(0).getAsType().getAsString();
    llvm::errs() << "arg type " << typestr;

    if (spec->isExplicitSpecialization()) llvm::errs() << " explicit";
    llvm::errs() << '\n';

    // write storage_type specialization
    // NOTE: this has to be moved to codegen, different for diff. codes
    if (field_storage_type_decl == nullptr) {
      llvm::errs() << " **** internal error: field_storage_type undefined in field\n";
      exit(1);
    }

    std::string fst_spec = "template<>\nstruct field_storage_type<"
      + typestr +"> {\n  " + typestr + " c[10];\n};\n";

    // insert after new line
    SourceLocation l =
      getSourceLocationAfterNewLine( field_storage_type_decl->getSourceRange().getEnd() );
    // TheRewriter.InsertText(l, fst_spec, true,true);
    writeBuf->insert(l, fst_spec, true, false);
    
  }
  return(count);
      
} // end of "field"


// Find the field_storage_type typealias here -- could not work
// directly with VisitTypeAliasTemplateDecl below, a bug??
bool MyASTVisitor::VisitDecl( Decl * D) {
  auto t = dyn_cast<TypeAliasTemplateDecl>(D);
  if (t && t->getNameAsString() == "field_storage_type") {
    llvm::errs() << "Got field storage\n";
  }
  
  return true;
}                           


//void MyASTVisitor::VisitTypeAliasTemplateDecl(TypeAliasTemplateDecl *D) {
  // if (D->isThisDeclarationADefinition()) {
  // Now the main definition
  // if (D->getNameAsString() == "field_storage_type") {
  //   TemplateParameterList * tplp = D->getTemplateParameters();
  //   std::stringstream SSBefore;
  //   SSBefore << "// field_storage_class def "
  //            << " with template params " ;
  //   for (unsigned i = 0; i < tplp->size(); i++) 
  //     SSBefore << tplp->getParam(i)->getNameAsString() << " ";
  //   SSBefore << "\n";
  //   SourceLocation ST = D->getSourceRange().getBegin();
  //   TheRewriter.InsertText(ST, SSBefore.str(), true, true);
  // }
//}


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


bool MyASTVisitor::VisitCXXMethodDecl(CXXMethodDecl *method) {
  // Comes after VisitFunctionDecl for methods -- this does really nothing here

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
    // TheRewriter.InsertText(ST, SSBefore.str(), true, true);
    writeBuf->insert(ST,SSBefore.str(), true, true);
    
    llvm::errs() << "Method decl " << method->getNameInfo().getName().getAsString() << '\n';
    
  }
  return true;
}


// This struct will be used to keep track of #include-chains.

std::vector<FileID> file_id_list = {};

// Tiny utility to search for the list

bool search_fid(const FileID FID) {
  for (const FileID f : file_id_list) {
    if (f == FID) return true;
  }
  return false;
}


// file_buffer_list stores the edited source of all files

struct file_buffer {
  srcBuf sbuf;
  FileID fid;
};

std::vector<file_buffer> file_buffer_list = {};

srcBuf * get_file_buffer(Rewriter & R, const FileID fid) {
  for (file_buffer & fb : file_buffer_list) {
    if (fb.fid == fid) return( &fb.sbuf );
  }
  // Now allocate and return new buffer
    
  file_buffer fb;
  fb.fid = fid;
  file_buffer_list.push_back(fb);
  SourceManager &SM = R.getSourceMgr();
  SourceRange r(SM.getLocForStartOfFile(fid),SM.getLocForEndOfFile(fid));

  llvm::errs() << "Create buf for file "
               << SM.getFilename(SM.getLocForStartOfFile(fid)) << '\n';
  
  file_buffer_list.back().sbuf.create( &R, r );
  return( &file_buffer_list.back().sbuf );
}  


void MyASTVisitor::set_writeBuf(const FileID fid) {
  writeBuf = get_file_buffer(TheRewriter, fid);
}


// Implementation of the ASTConsumer interface for reading an AST produced
// by the Clang parser.
class MyASTConsumer : public ASTConsumer {
public:
  MyASTConsumer(Rewriter &R, ASTContext *C) : Visitor(R,C) { Rewriterp = &R; }


  // HandleTranslationUnit is called after the AST for the whole TU is completed
  virtual void HandleTranslationUnit(ASTContext & ctx) override {
    // dump ast here -- HERE THE SPECIALIZATIONS ARE PRESENT!
    // ctx.getTranslationUnitDecl()->dump();
    
    SourceManager &SM = ctx.getSourceManager();
    TranslationUnitDecl *tud = ctx.getTranslationUnitDecl();
 
    for (DeclContext::decl_iterator d = tud->decls_begin(); d != tud->decls_end(); d++) {
      
      SourceLocation beginloc = d->getBeginLoc();
      // analyze only user files (these should be named)
      if (!SM.isInSystemHeader(beginloc) && SM.getFilename(beginloc) != "") {
      //if (!SM.isInSystemHeader(beginloc)) {
        // llvm::errs() << "Processing file " << SM.getFilename(beginloc) << "\n";
        // TODO: ensure that we go only through files which are needed!

        state::loop_found = false;
      
        // get our own file edit buffer (unless it exists)
        Visitor.set_writeBuf(SM.getFileID(beginloc));

        // Traverse the declaration using our AST visitor.

        global.location.top = d->getSourceRange().getBegin();  // save this for source location
        
        Visitor.TraverseDecl(*d);
        // llvm::errs() << "Dumping level " << i++ << "\n";
        if (cmdline::dump_ast) {
          if (!cmdline::no_include || SM.isInMainFile(beginloc))
            d->dump();
        }
        
        // We keep track here only of files which were touched
        if (state::loop_found) {
          FileID FID = SM.getFileID(beginloc);
          if (search_fid(FID) == false) {
            // new file to be added
            file_id_list.push_back(FID);
            // llvm::errs() << "New file changed " << SM.getFileEntryForID(FID)->getName() << '\n';
          }
        }
      }  
    }

    // check compile errors, as long as we have context -- use diagnostics engine
    auto & DE = ctx.getDiagnostics();
    state::compile_errors_occurred = DE.hasErrorOccurred();
    
  }


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

    global.main_file_name = getCurrentFile();
    
    file_id_list.clear();
    file_buffer_list.clear();
    field_decl = field_storage_type_decl = nullptr;
    
    return (true);
  }

  void insert_includes_to_file_buffer(FileID myFID) {
    // this is where to write
    srcBuf * buf = get_file_buffer(TheRewriter, myFID);
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
        // TODO: do this on "buf" instead of original file data
        for (int i=1; i<100; i++) {
          const char * p = SM.getCharacterData(b.getLocWithOffset(-i));
          if (p && *p == '#' && strncmp(p,"#include",8) == 0) {
            SR = SourceRange(b.getLocWithOffset(-i),e);
            break;
          }
        }
        // Remove "#include"
        buf->remove(SR);
        // TheRewriter.RemoveText(SR);

        // and finally insert
        // SourceRange r(SM.getLocForStartOfFile(f),SM.getLocForEndOfFile(f));
        srcBuf * buf_from = get_file_buffer(TheRewriter, f);
        // TheRewriter.InsertText(SR.getBegin(),
        buf->insert(SR.getBegin(), 
                    "// start include "+includestr
                    + "---------------------------------\n"
                    + buf_from->dump() +
                    "// end include "+includestr
                    + "---------------------------------\n",
                    false);
        
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


  void EndSourceFileAction() override {
    SourceManager &SM = TheRewriter.getSourceMgr();
    llvm::errs() << "** EndSourceFileAction for: " << getCurrentFile() << '\n';
    // << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";
    
    // Now emit rewritten buffers.

    if (!cmdline::no_output) {
      if (!cmdline::no_include) {

        // Modified files should be substituted on top of #include -directives
        // first, ensure that the full include chain is present in file_id_list
        // Use iterator here, because the list can grow!

        for ( FileID f : file_id_list ) {
          check_include_path(f);
        }

        insert_includes_to_file_buffer(SM.getMainFileID());
      }
    
      if (!state::compile_errors_occurred) {
        write_output_file( cmdline::output_filename,
                           get_file_buffer(TheRewriter,SM.getMainFileID())->dump() );
        
        if (!cmdline::spec_no_db) write_specialization_db();
      } else {
        llvm::errs() << program_name << ": not writing output due to compile errors\n";
      }
    }
    
    file_buffer_list.clear();
    file_id_list.clear();

    // EndSourceFile();
        
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


void get_target_struct(codetype & target) {
  if (cmdline::kernel) target.kernelize = true;
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
  if (cmdline::syntax_only) cmdline::no_output = true;
  
  // ClangTool::run accepts a FrontendActionFactory, which is then used to
  // create new objects implementing the FrontendAction interface. Here we use
  // the helper newFrontendActionFactory to create a default factory that will
  // return a new MyFrontendAction object every time.
  // To further customize this, we could create our own factory class.
  return Tool.run(newFrontendActionFactory<MyFrontendAction>().get());
}
