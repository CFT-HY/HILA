//------------------------------------------------------------------------------
// Generate transformed 
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

// lots of global vars -- perhaps invent a class or something

extern global_state global;
extern loop_parity_struct loop_parity;

extern std::list<field_ref> field_ref_list;
extern std::list<field_info> field_info_list;
extern std::list<var_expr> var_expr_list;
extern std::list<var_decl> var_decl_list;

const std::string looping_var = "LFI_";


std::string parity_str(parity p)
{
  switch (p) {
  case parity::none : return "NONE";
  case parity::even : return "EVEN";
  case parity::odd  : return "ODD";
  case parity::all  : return "ALL";
  case parity::x    : return "X";
  }
}


/// The main entry point for code generation

bool MyASTVisitor::generate_code( SourceLocation kernelloc, Stmt *S ) {
  static int k_number = 0;
    
  // Get random kernel name
  k_number ++; 
  std::string kernel_name = "kernel_"
    + global.currentFunctionDecl->getNameInfo().getName().getAsString()
    + "_" + std::to_string( k_number );

  std::string call, def;

  // is it compound stmt: { } -no ; needed
  bool semi_at_end = !(isa<CompoundStmt>(S));
    
  // encapsulate everything
  call = "{\n";
    
  // basic set up: 1st loop_parity, if it is known const set it up,
  // else copy it to a variable name

  std::string parity_name;
    
  if (loop_parity.value == parity::none) {
    // now unknown
    parity_name = "lf_parity_";
    call += "const parity " + parity_name + " = " + loop_parity.text + ";\n";

    if (global.assert_loop_parity) {
      call += "assert_even_odd_parity(" + parity_name + ");\n";
    }
      
  } else parity_name = parity_str(loop_parity.value);

  // Insert field neighb fetches

  for (field_info & l : field_info_list) 
    for (dir_ptr & d : l.dir_list ) {
      // TODO - move to temp vars
      call += l.old_name + ".start_move(" + get_stmt_str(d.e) + ", " + parity_name + ");\n";
      // TODO: must add wait_gets, if we want those
    }
    
  call += kernel_name + "(";
  def   = "//vvvvvv\nvoid " + kernel_name + "(";
    
  if (loop_parity.value == parity::none) {
    call += parity_name + ", ";
    def  += "const parity lf_parity_, ";
  }
      
  // print field call list
  int i = -1;        
  for (field_info & l : field_info_list) {
    i++;
      
    if (i>0) {
      def  += ", ";
      call += ", ";
    }
      
    if (!l.is_changed) def += "const ";
    // TODO: type to field_data
    def += l.type + " & " + l.new_name;
    call += l.old_name;
      
  }

  i=0;
  // and non-field vars
  for ( var_expr &ep : var_expr_list ) {
    bool dup = (ep.duplicate != nullptr);

    if (!dup) {
      std::string varname = "nolf_var_"+std::to_string(i);
      def += ", const " + ep.type + " " + varname;
      call += ", " + get_stmt_str(ep.e);
      
      //SourceRange range( E->getSourceRange().getBegin(), ep.e->getSourceRange().getBegin() );
      //llvm::errs() << " ++ diff in position " << TheRewriter.getRangeSize(range) << '\n';
      
      // replace_expr(ep, varname);
      Buf.replace( ep.ind, varname );

      i++;
    } else {
      Buf.replace( ep.ind, Buf.token(ep.duplicate->ind) );
    }
  }
 
  // finally, change the references to lf variables in the body
    
  for ( field_ref & le : field_ref_list ) {
    Buf.replace( le.nameInd, le.info->new_name );
    if (le.dirExpr != nullptr) {
      Buf.replace(le.parityInd, "neighbour(" + looping_var + ", " + get_stmt_str(le.dirExpr)+")");
    } else {
      Buf.replace(le.parityInd, looping_var);
    }      
  }                        
    
  def += ")\n{\nforparity("+looping_var+", par_v) {\n" + Buf.dump();
  if (semi_at_end) def += ';';
  def += + "\n}\n}\n//^^^^^^\n";
  call += ");\n";

  // and mark changed vars
  for ( field_info  & l : field_info_list ) {
    if (l.is_changed) call += l.old_name + ".changed(" + parity_name + ");\n";
  }
    
  // and close
  call += "}\n//^^^^^^";

  // Finally, emit the kernel
  TheRewriter.InsertText(kernelloc, indent_string(def),true,true);
  // and the call

  if (semi_at_end) {
    TheRewriter.RemoveText(getRangeWithSemi(S));
  } else {
    TheRewriter.RemoveText(S->getSourceRange());
  }
    
  // TheRewriter.InsertText(E->getBeginLoc(), "//-  "+ global.full_loop_text + '\n', true, true );
  TheRewriter.InsertText(S->getBeginLoc(), indent_string(call), true, true);
  //TheRewriter.InsertText(getRangeWithSemi(S,false).getEnd().getLocWithOffset(1),
  //                        call, true, true);
  //replace_expr(S, call);

  return true;
}
