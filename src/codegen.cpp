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

#include "transformer.h"
#include "myastvisitor.h"
#include "stringops.h"

std::string looping_var;
std::string parity_name;

std::string parity_in_this_loop = "";

/// Used in lattice loop generation
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

/// Generate a name that does not appear in string t
inline std::string unique_name( const std::string t,  std::string n){
  while (t.find(n,0) != std::string::npos) n += "_";
  return n;
}



/// The main entry point for code generation
void MyASTVisitor::generate_code(Stmt *S, codetype & target) {
  srcBuf loopBuf; // (&TheRewriter,S);
  loopBuf.copy_from_range(writeBuf,S->getSourceRange());
  
  //   llvm::errs() << "\nOriginal range: +++++++++++++++\n\""
  //                << TheRewriter.getRewrittenText(S->getSourceRange()) 
  //                << "\"\nwriteBuf range: ================\n\""
  //                << writeBuf->get(S->getSourceRange())
  //                << "\"\nCopied range: ================\n\""
  //                << loopBuf.dump() << "\"\n";
  
  // is it compound stmt: { } -no ; needed
  bool semi_at_end = !(isa<CompoundStmt>(S));
 
  // Build replacement in variable "code"
  // Encapsulate everything within {}
  std::stringstream code;
  code << "{\n";

  // basic set up: 1st loop_parity, if it is known const set it up,
  // else copy it to a variable name

  const std::string t = loopBuf.dump();

  looping_var = "Index";
  while (t.find(looping_var,0) != std::string::npos) looping_var += "_";
  
  // Ensure that the name is not reserved by scanning the source
  parity_name = "Parity";
  while (t.find(parity_name,0) != std::string::npos) parity_name += "_";
  
  if (loop_parity.value == parity::none) {
    // now unknown
    code << "const parity " << parity_name << " = " << loop_parity.text << ";\n";

    if (global.assert_loop_parity) {
      code << "assert_even_odd_parity(" << parity_name << ");\n";
    }
    parity_in_this_loop = parity_name;
      
  } 
  else parity_in_this_loop = parity_str(loop_parity.value);

  for (field_info & l : field_info_list) {
    // Generate new variable name, may be needed -- use here simple receipe
    l.new_name = "F"+clean_name(l.old_name);
    // Perhaps simpler FA, FB, FC. ?  The following helps avoid collisions
    while (t.find(l.new_name,0) != std::string::npos) l.new_name += "_";
    l.loop_ref_name = l.new_name + "_index";
    
    // variable links if needed
    // if (l.dir_list.size() > 0) {
    if (l.is_written) {
      code << "field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
      code << l.new_name << ".mark_changed(" + parity_in_this_loop + ");\n";      
    } else {
      code << "const field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
    }

    // Check that read fields are allocated
    if( l.is_read ){
      code << "assert(" << l.new_name << ".is_allocated());\n";
    }

    // If neighbour references exist, communicate them
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      code << l.new_name << ".wait_move("
           << get_stmt_str(d.e) << ", " << parity_in_this_loop << ");\n";
    }
  }
  
  // Create temporary variables for reductions
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      v.reduction_name = "r_" + v.name;
      while (t.find(v.reduction_name,0) != std::string::npos) v.reduction_name += "_";
      // Create a temporary variable and initialize
      if (v.reduction_type == reduction::SUM) {
        code << v.type << " " << v.reduction_name << "=0;\n";
      } else if (v.reduction_type == reduction::PRODUCT) {
        code << v.type << " " << v.reduction_name << "=1;\n";
      }
    }
  }

  // Place the content of the loop
  if( target.CUDA ){
    code << generate_code_cuda(S,semi_at_end,loopBuf);
  } else if( target.openacc ){
    code << generate_code_openacc(S,semi_at_end,loopBuf);
  } else if(target.VECTORIZE) {
    code << generate_code_avx(S,semi_at_end,loopBuf);
  } else {
    code << generate_code_cpu(S,semi_at_end,loopBuf);
  }
  
  // Check reduction variables
  for (var_info & v : var_info_list) {
    if (v.reduction_type == reduction::SUM) {
      code << "lattice->reduce_node_sum(" << v.reduction_name << ", true);\n";
      code << v.name << " += " << v.reduction_name << ";\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code << "lattice->reduce_node_product(" << v.reduction_name << ", true);\n";
      code << v.name << " *= " << v.reduction_name << ";\n";
    }
  }
          
  // and close
  code << "}\n//----------";
  
  // Remove old code + replace
  if (semi_at_end) {
    writeBuf->remove(getRangeWithSemi(S));
  } else {
    writeBuf->remove(S->getSourceRange());
  }
    
  writeBuf->insert(S->getBeginLoc(), indent_string(code.str()), true, true);
}


void MyASTVisitor::generate_field_storage_type(std::string typestr){
  if (field_storage_decl == nullptr) {
    llvm::errs() << " **** internal error: element undefined in field\n";
    exit(1);
  }

  if(target.VECTORIZE){
    generate_field_storage_type_AVX(typestr);
  }
}



/* Generate a header that marks field references read or written.
 * This is a copy of the loop body with modifications, only ran once
 */
/*
std::string MyASTVisitor::generate_loop_header(Stmt *S, codetype & target, bool semi_at_end) {
  srcBuf loopBuf;
  loopBuf.copy_from_range(writeBuf,S->getSourceRange());
  std::vector<std::string> va = {}, vb = {};
  int i=0;

  // Replace loop references with temporary variables
  // and add calls to mark_changed() and start_get()
  for ( field_info & fi : field_info_list ) {
    std::string type_name = fi.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());

    // Add a simple temp variable to replace the field reference
    std::string varname = "__v_" + std::to_string(i);
    loopBuf.prepend(type_name + " " + varname + ";\n", true);

    for( field_ref *le : fi.ref_list ){
      Expr *e = le->nameExpr;
      loopBuf.replace(le->fullExpr, varname );
      if( le->is_written ){
        // Mark changed fields - do it BEFORE using vars, possible alloc
        loopBuf.insert_before_stmt(e, get_stmt_str(e) + ".mark_changed(" + parity_in_this_loop + ");", true,   true);
      }
      if( le->dirExpr != nullptr ){
        // If a field needs to be communicated, start here
        loopBuf.insert_before_stmt(e, get_stmt_str(e) + ".wait_move(" + get_stmt_str(le->dirExpr) + ", " 
           + parity_in_this_loop + ");", true, true);
      }
      if( le->is_read ){
        // If a field is read, check that is has been allocated
        loopBuf.insert_before_stmt(e, "assert(" + get_stmt_str(e) 
          + ".is_allocated());", true, true);
      }

      //va.push_back(get_stmt_str(le->fullExpr));
      //vb.push_back(varname);
      llvm::errs() << get_stmt_str(le->fullExpr) << ":" << varname << "\n";
    }
    i++;
  }

  // Replace reduction variables with copies to avoid changing originals
  // No other variables are allowed to change
  for ( var_info & vi : var_info_list ) {
    if( vi.reduction_type != reduction::NONE ){
      std::string varname = "__v_" + std::to_string(i);
      va.push_back(vi.name);
      vb.push_back(varname);
      loopBuf.prepend(vi.type + " " + varname + "=" + vi.name + ";\n", true);
    }
    i++;
  }

  // Surround by curly brackets to keep new variables local
  loopBuf.replace_tokens( va, vb );
  loopBuf.prepend("{",true);

  if(semi_at_end){
    return loopBuf.dump() + ";}\n"; 
  } else {
    return loopBuf.dump() + "}\n";
  }
}
*/