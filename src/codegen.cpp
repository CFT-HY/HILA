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

std::string looping_var;
std::string parity_name;

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
    code << "field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
    // l.type_template is < type >.  Change this to field_storage_type<T>
    //if (!target.kernelize) {
    //  code << "field_storage_type" << l.type_template << " * const " 
    //       << l.loop_ref_name << " = " << l.new_name << "->fs.payload;\n";
    //}
  }
  
  // Assert that vars to be read are initialized TODO: make optional (cmdline?)
  int i = 0;
  for ( field_info  & l : field_info_list ) if (l.is_read) {
    if (i == 0) {
      code << "assert(" << l.new_name <<".is_allocated()";
    } else {
      code << " && " << l.new_name + ".is_allocated()";
    }
    i++;
  }
  if (i > 0) code << ");\n";
  
  
  // Insert field neighbour fetches
  
  for (field_info & l : field_info_list) {
    for (dir_ptr & d : l.dir_list ) {
      // TODO - move to temp vars?
      code << l.new_name << ".start_move(" << get_stmt_str(d.e) << ", " 
           << parity_in_this_loop << ");\n";
      // TODO: must add wait_gets, if we want those
    }
  }

    
  // Mark changed vars - do it BEFORE using vars, possible alloc
  for ( field_info  & l : field_info_list ) {
    if (l.is_written) code << l.new_name << ".mark_changed(" << parity_in_this_loop << ");\n";
  }

  // Here generate kernel, or produce the loop in place
  if (target.kernelize) {
    code << generate_kernel(S,semi_at_end,loopBuf);
  } else {
    // Now in place
    code << generate_in_place(S,semi_at_end,loopBuf);
  }
  
  // Check reduction variables
  for (var_info & v : var_info_list) { 
    if (v.reduction_type == reduction::SUM) {
      code << "lattice->reduce_node_sum(" << v.name << ", true);\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code << "lattice->reduce_node_product(" << v.name << ", true);\n";
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


std::string MyASTVisitor::generate_in_place(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {
  
  replace_field_refs(loopBuf);
  
  std::stringstream code;
  code << "const int loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "const int loop_end   = lattice->loop_end(" << parity_in_this_loop << ");\n";

  code << "for(int " << looping_var <<" = loop_begin; " 
       << looping_var << " < loop_end; " << looping_var << "++) {\n";
  
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      if(l.is_read) {
        code << type_name << l.loop_ref_name << "_" << get_stmt_str(d.e) << " = " << l.new_name 
             << ".get_value_at(" << "lattice->neighb[" << get_stmt_str(d.e) << "][" 
             << looping_var + "]" << ");\n";
      } else {
        code << type_name << l.loop_ref_name << "_" << get_stmt_str(d.e) << ";";
      }
    }
    if(l.is_read) {
      code << type_name << l.loop_ref_name << " = " << l.new_name 
           << ".get_value_at(" << looping_var << ");\n";
    } else {
      code << type_name << l.loop_ref_name << ";";
    }
  }

  code << loopBuf.dump();
  if (semi_at_end) code << ";\n";

  for (field_info & l : field_info_list){
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    // Probably shouldn't be setting values to neighbours...
    //for (dir_ptr & d : l.dir_list) if(d.count > 0){
    //  code << l.new_name << ".set_value_at(" << l.loop_ref_name << "_" 
    //       << get_stmt_str(d.e) 
    //       << ", lattice->neighb[" << get_stmt_str(d.e)
    //       << "][" << looping_var + "]);\n";
    //}
    if(l.is_written) {
      code << l.new_name << ".set_value_at(" << l.loop_ref_name << ", " 
           << looping_var << ");\n";
    }
  }

  code << "}\n";

  return code.str();
}

/// return value: kernel call code

std::string MyASTVisitor::generate_kernel(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {

  // Get kernel name - use line number or file offset (must be deterministic)
  std::string kernel_name = make_kernel_name();
  
  std::stringstream kernel,call;

  call   << kernel_name << "(";
  kernel << "//----------\n";
  // I don't think kernels need to be templates
  //   if (global.in_func_template) {
  //     kernel << "template <typename ";
  //     for (unsigned i = 0; i < global.function_tpl->size(); i++) {
  //       if (i>0) kernel << ", ";
  //       kernel << global.function_tpl->getParam(i)->getNameAsString();
  //     }
  //     kernel << ">\n";
  //   }
  kernel << "void " << kernel_name << "(";
    
  if (loop_parity.value == parity::none) {
    call   << parity_name << ", ";
    kernel << "const parity " << parity_name << ", ";
  }

  // print field call list
  int i = -1;
  for (field_info & l : field_info_list) {
    i++;
    
    if (i>0) {
      kernel << ", ";
      call   << ", ";
    }
    
    // This does not recognise setting in functions
    if (!l.is_written) kernel << "const ";
    // TODO: type to field_data
    kernel << "field_struct" << l.type_template << " * " << l.new_name;
    call << l.new_name + ".fs";
  }

  i=0;
  // and non-field vars
  for ( var_info & vi : var_info_list ) {
    if (!vi.is_loop_local) {
      std::string varname = "sv__" + std::to_string(i) + "_";
      kernel << ", ";
      if(vi.is_assigned) {
        kernel << vi.type << " & " << varname;
        call << ", " << vi.name;
      } else {
        kernel << "const " << vi.type << " " << varname;
        call << ", " << vi.name;
      }
      
      // replace_expr(ep, varname);
      for (var_ref & vr : vi.refs) {
        loopBuf.replace( vr.ref, varname );
      }
      i++;
    }
  }

  // Directions
  for (field_info & l : field_info_list) {
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      int i=0;
      call << ", " << get_stmt_str(d.e);
      kernel << ", const int " << get_stmt_str(d.e);
    }
  }

  // Begin the function
  kernel << ")\n{\n";

  // finally, change the references to variables in the body
  replace_field_refs(loopBuf);

  // Generate the loop
  kernel << "const int loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  kernel << "const int loop_end   = lattice->loop_end(" << parity_in_this_loop << ");\n";

  kernel << "for(int " << looping_var <<" = loop_begin; " 
         << looping_var << " < loop_end; " << looping_var << "++) {\n";

  // Call getters
  for (field_info & l : field_info_list){
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      if(l.is_read){
        kernel << type_name << l.loop_ref_name << "_" << get_stmt_str(d.e) << " = " << l.new_name 
             << "->get(" << "lattice->neighb[" << get_stmt_str(d.e) << "][" 
             << looping_var + "]" << ");\n";
      } else {
        kernel << type_name << l.loop_ref_name << "_" << get_stmt_str(d.e) << ";";
      }
    }
    if(l.is_read) {
      kernel << type_name << l.loop_ref_name << " = " << l.new_name 
             << "->get(" << looping_var << ");\n";
    } else {
      kernel << type_name << l.loop_ref_name << ";";
    }
  }

  kernel << loopBuf.dump();
  if (semi_at_end) kernel << ';';
  kernel << '\n';

  // Call setters
  for (field_info & l : field_info_list) if(l.is_written){
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    kernel << l.new_name << "->set(" << l.loop_ref_name << ", " 
           << looping_var << ");\n";
  }

  kernel << "}\n}\n//----------\n";
  call << ");\n";

  // Finally, emit the kernel
  // TheRewriter.InsertText(global.location.function, indent_string(kernel),true,true);
  toplevelBuf->insert(global.location.top.getLocWithOffset(-1), indent_string(kernel.str()),true,false);

  return call.str();
}


/// Change field references within loops
void MyASTVisitor::replace_field_refs(srcBuf & loopBuf) {
  
  for ( field_ref & le : field_ref_list ) {
    //loopBuf.replace( le.nameExpr, le.info->loop_ref_name );
    if (le.dirExpr != nullptr) {
      //loopBuf.replace(parityExpr,
      //                 "lattice->neighb[" +  get_stmt_str(le.dirExpr) + "][" + looping_var + "]");
      loopBuf.replace(le.fullExpr, le.info->loop_ref_name+"_"+get_stmt_str(le.dirExpr));
    } else {
      loopBuf.replace(le.fullExpr, le.info->loop_ref_name);
    }
  }                        
}
