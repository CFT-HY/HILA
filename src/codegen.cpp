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
extern std::list<var_info> var_info_list;
extern std::list<var_decl> var_decl_list;

const std::string looping_var("FI__");
const std::string parity_name("pty__");

static std::string parity_in_this_loop = "";

std::string parity_str(parity p)
{
  switch (p) {
  case parity::none : return "parity::none";
  case parity::even : return "parity::even";
  case parity::odd  : return "parity::odd";
  case parity::all  : return "parity::all";
  case parity::x    : return "parity::x";
  }
}


/// Help routine to write (part of) a name for a kernel
std::string MyASTVisitor::make_kernel_name() {
  return 
    "kernel_"
    + clean_name(global.currentFunctionDecl
                 ->getNameInfo().getName().getAsString())
    + "_" + std::to_string(TheRewriter.getSourceMgr().
                           // getSpellingLineNumber(global.location.loop));
                           getFileOffset(global.location.loop));
}



/// The main entry point for code generation

void MyASTVisitor::generate_code(Stmt *S, codetype & target) {
 
  // is it compound stmt: { } -no ; needed
  bool semi_at_end = !(isa<CompoundStmt>(S));
 
  // Build replacement in variable "code"
  // Encapsulate everything within {}
  std::string code = "{\n";

  // basic set up: 1st loop_parity, if it is known const set it up,
  // else copy it to a variable name

  if (loop_parity.value == parity::none) {
    // now unknown
    code += "const parity " + parity_name + " = " + loop_parity.text + ";\n";

    if (global.assert_loop_parity) {
      code += "assert_even_odd_parity(" + parity_name + ");\n";
    }
    parity_in_this_loop = parity_name;
      
  } 
  else parity_in_this_loop = parity_str(loop_parity.value);

  for (field_info & l : field_info_list) {
    // Generate new variable name, may be needed -- use here simple receipe
    l.new_name = "f__"+clean_name(l.old_name)+"_";

    // variable links if needed
    if (!target.kernelize || l.dir_list.size() > 0) {
      code += l.type + " & " + l.new_name + " = " + l.old_name + ";\n";
    }
  }
  
  // Insert field neighbour fetches

  for (field_info & l : field_info_list) 
    for (dir_ptr & d : l.dir_list ) {
      // TODO - move to temp vars
      code += l.new_name + ".start_move(" + get_stmt_str(d.e) + ", " 
           + parity_in_this_loop + ");\n";
      // TODO: must add wait_gets, if we want those
    }
  
  // Here generate kernel, or produce the loop in place
  if (target.kernelize) {
    code += generate_kernel(S,semi_at_end);
  } else {
    // Now in place
    code += generate_in_place(S,semi_at_end);
  }
  
  // Check reduction variables
  for (var_info & v : var_info_list) { 
    if (v.reduction_type == reduction::SUM) {
      code += "reduce_node_sum(" + v.name + ", true);\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code += "reduce_node_product(" + v.name + ", true);\n";
    }
  }
      
  // Mark changed vars
  for ( field_info  & l : field_info_list ) {
    if (l.is_changed) code += l.old_name + ".changed(" + parity_in_this_loop + ");\n";
  }
    
  // and close
  code += "}\n//----------";
  
  // Remove old code + replace
  if (semi_at_end) {
    TheRewriter.RemoveText(getRangeWithSemi(S));
  } else {
    TheRewriter.RemoveText(S->getSourceRange());
  }
    
  // TheRewriter.InsertText(E->getBeginLoc(), "//-  "+ global.full_loop_text + '\n', true, true );
  TheRewriter.InsertText(S->getBeginLoc(), indent_string(code), true, true);
  //TheRewriter.InsertText(getRangeWithSemi(S,false).getEnd().getLocWithOffset(1),
  //                        call, true, true);
  //replace_expr(S, call);
  
}


std::string MyASTVisitor::generate_in_place(Stmt *S, bool semi_at_end) {
  
  replace_field_refs();
  
  std::string code = "forparity("+looping_var+", "+parity_in_this_loop+") {\n" + Buf.dump();
  if (semi_at_end) code += ';';
  code += "\n}\n";

  return code;
}

/// return value: kernel call code

std::string MyASTVisitor::generate_kernel(Stmt *S, bool semi_at_end) {

  // Get kernel name - use line number or file offset (must be deterministic)
  std::string kernel_name = make_kernel_name();
                           
  std::string kernel,call;

  call   = kernel_name + "(";
  kernel = "//----------\n";
  if (global.in_func_template) {
    kernel += "template <typename ";
    for (unsigned i = 0; i < global.function_tpl->size(); i++) {
      if (i>0) kernel += ", ";
      kernel += global.function_tpl->getParam(i)->getNameAsString();
    }
    kernel += ">\n";
  }
  kernel += "void " + kernel_name + "(";
    
  if (loop_parity.value == parity::none) {
    call   += parity_name + ", ";
    kernel += "const parity " + parity_name + ", ";
  }
      
  // print field call list
  int i = -1;        
  for (field_info & l : field_info_list) {
    i++;
      
    if (i>0) {
      kernel += ", ";
      call   += ", ";
    }

    if (!l.is_changed) kernel += "const ";
    // TODO: type to field_data
    kernel += l.type + " & " + l.new_name;
    if (l.dir_list.size() == 0) call += l.old_name;
    else call += l.new_name;
  }

  i=0;
  // and non-field vars
  for ( var_info & vi : var_info_list ) {
    if (!vi.is_loop_local) {
      std::string varname = "sv__" + std::to_string(i) + "_";
      kernel += ", const " + vi.type + " " + varname;
      call += ", " + vi.name;
      
      // replace_expr(ep, varname);
      for (var_ref & vr : vi.refs) {
        Buf.replace( vr.ind, varname );
      }
    }
  }
  // finally, change the references to variables in the body
  replace_field_refs();
      
  kernel += ")\n{\nforparity("+looping_var+", "+parity_in_this_loop+") {\n" + Buf.dump();
  if (semi_at_end) kernel += ';';
  kernel += "\n}\n}\n//----------\n";
  call += ");\n";

  // Finally, emit the kernel
  TheRewriter.InsertText(global.location.function, indent_string(kernel),true,true);

  return call;
}


/// Change field references within loops
void MyASTVisitor::replace_field_refs() {
  
  for ( field_ref & le : field_ref_list ) {
    Buf.replace( le.nameInd, le.info->new_name );
    if (le.dirExpr != nullptr) {
      Buf.replace(le.parityInd, "neighbour(" + looping_var + ", " + get_stmt_str(le.dirExpr)+")");
    } else {
      Buf.replace(le.parityInd, looping_var);
    }      
  }                        
}
