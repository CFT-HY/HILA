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
    if (l.is_written) {
      code << "field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
    } else {
      code << "const field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
    }
    // l.type_template is < type >.  Change this to field_storage_type<T>
    //if (!target.kernelize) {
    //  code << "field_storage_type" << l.type_template << " * const " 
    //       << l.loop_ref_name << " = " << l.new_name << "->fs.payload;\n";
    //}
  }

  
  // Generate a header that starts communication and checks
  // that fields are allocated
  code << generate_loop_header(S,target,semi_at_end)+ "\n";
  
  
    // Check for reductions and allocate device memory
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      v.new_name = "r_" + v.name;
      while (t.find(v.new_name,0) != std::string::npos) v.new_name += "_";
      if(target.CUDA){
      // Allocate memory for a reduction. This will be filled in the kernel
        code << v.type << " *" << v.new_name << ";\n";
        code << "cudaMalloc( (void **)& " << v.new_name << ","
             << "sizeof(" << v.type << ") * lattice->volume() );\n";
        code << "check_cuda_error(\"allocate_reduction\");\n";
      } else {
        // Create a temporary variable and initialize
        if (v.reduction_type == reduction::SUM) {
          code << v.type << " " << v.new_name << "=0;\n";
        } else if (v.reduction_type == reduction::PRODUCT) {
          code << v.type << " " << v.new_name << "=1;\n";
        }
      }
    }
  }

  // Here generate kernel, or produce the loop in place
  if (target.kernelize) {
    code << generate_kernel(S,target,semi_at_end,loopBuf);
  } else {
    // Now in place
    code << generate_in_place(S,target,semi_at_end,loopBuf);
  }
  
  // Check reduction variables
  for (var_info & v : var_info_list) {
    if(target.CUDA){
      // Run reduction
      if (v.reduction_type == reduction::SUM) {
        code << v.type << " " << v.new_name << "_sum += cuda_reduce_sum(" 
             << v.new_name << ", lattice->volume()" <<  ");\n";
        v.new_name = v.new_name + "_sum";
      } else if (v.reduction_type == reduction::PRODUCT) {
        code << v.type << " " << v.new_name << "_prod *= cuda_reduce_product(" 
             << v.new_name << ", lattice->volume()" <<  ");\n";
        v.new_name = v.new_name + "_prod";
      }
      // Free memory allocated for the reduction
      if (v.reduction_type != reduction::NONE) {
        code << "cudaFree(" << v.new_name << ");\n";
        code << "check_cuda_error(\"free_reduction\");\n";
      }
    }
    // Add reduction over MPI nodes and add to the original variable
    if (v.reduction_type == reduction::SUM) {
      code << "lattice->reduce_node_sum(" << v.new_name << ", true);\n";
      code << v.name << " += " << v.new_name << ";\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code << "lattice->reduce_node_product(" << v.new_name << ", true);\n";
      code << v.name << " *= " << v.new_name << ";\n";
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


/* Generate a header that marks field references read or written.
 * This is a copy of the loop body with modifications, only ran once
 */
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
        loopBuf.insert_before_stmt(e, get_stmt_str(e) + ".start_move(" + get_stmt_str(le->dirExpr) + ", " 
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



std::string MyASTVisitor::generate_in_place(Stmt *S, codetype & target, bool semi_at_end, srcBuf & loopBuf) {

  // replace reduction variables in the loop
  for ( var_info & vi : var_info_list ) {
    if (!vi.is_loop_local) {
      if(vi.reduction_type != reduction::NONE) {
        std::string varname = "r_" + vi.name;
        for (var_ref & vr : vi.refs) {
          loopBuf.replace( vr.ref, varname );
        }
      }
    }
  }
  
  replace_field_refs_and_funcs(loopBuf, false);
  
  std::stringstream code;
  code << "const int loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "const int loop_end   = lattice->loop_end(" << parity_in_this_loop << ");\n";

  // Add openacc pragmas
  if( target.openacc ) {
    code << "#pragma acc parallel loop";
    // Check reduction variables
    for (var_info & v : var_info_list) { 
      if (v.reduction_type == reduction::SUM) {
        code << " reduction(+:" << v.name << ")";
      } else if (v.reduction_type == reduction::PRODUCT) {
        code << " reduction(*:" << v.name << ")";
      }
    }
    code << "\n";
  }

  // Generate the loop
  code << "for(int " << looping_var <<" = loop_begin; " 
       << looping_var << " < loop_end; " << looping_var << "++) {\n";
  
  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());

    // First check for direction references. If any found, create list of temp
    // variables
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      code << type_name << " "  << l.loop_ref_name << "_d[NDIRS];\n";
      code << "std::array<bool, NDIRS> " << l.new_name << "_read_d;\n";
      code << l.new_name << "_read_d.fill(true);\n";
      break; // Only need one direction reference
    }
    // Check for references without a direction. If found, add temp variable
    for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
      code << type_name << " "  << l.loop_ref_name << ";\n";
      code << "bool "  << l.new_name << "_read = true;\n";
      break;  // Only one needed
    }
  }

  // Dump the main loop code here
  code << loopBuf.dump();
  if (semi_at_end) code << ";";
  code << "\n";

  // Call setters
  for (field_info & l : field_info_list){
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    if(l.is_written) {
      code << l.new_name << ".set_value_at(" << l.loop_ref_name << ", " 
           << looping_var << ");\n";
    }
  }

  code << "}\n";

  return code.str();
}

/// return value: kernel call code

std::string MyASTVisitor::generate_kernel(Stmt *S, codetype & target, bool semi_at_end, srcBuf & loopBuf) {

  // Get kernel name - use line number or file offset (must be deterministic)
  std::string kernel_name = make_kernel_name();
  
  std::stringstream kernel,call;

  kernel << "//----------\n";
  //   if (global.in_func_template) {
  //     kernel << "template <typename ";
  //     for (unsigned i = 0; i < global.function_tpl->size(); i++) {
  //       if (i>0) kernel << ", ";
  //       kernel << global.function_tpl->getParam(i)->getNameAsString();
  //     }
  //     kernel << ">\n";
  //   }

  // Generate the function definition and call
  if( target.CUDA ){
    kernel << "__global__ void " << kernel_name << "( device_lattice_info lattice_info, ";
    call << "device_lattice_info lattice_info = lattice->device_info;\n";
    call << "lattice_info.loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
    call << "lattice_info.loop_end = lattice->loop_end(" << parity_in_this_loop << ");\n";
    call << "int N_blocks = (lattice_info.loop_end - lattice_info.loop_begin)/N_threads + 1;\n";
    call << kernel_name << "<<< N_blocks, N_threads >>>( lattice_info, ";
  } else {
    kernel << "void " << kernel_name << "(";
    call   << kernel_name << "(";
    // Include parity as a parameter
      if (loop_parity.value == parity::none) {
      call   << parity_name << ", ";
      kernel << "const parity " << parity_name << ", ";
    }
  }


  // print field call list
  int i = -1;
  for (field_info & l : field_info_list) {
    i++;
    
    if (i>0) {
      kernel << ", ";
      call   << ", ";
    }
    
    if (!l.is_written) kernel << "const ";
    // TODO: type to field_data
    kernel << "field_storage" << l.type_template << " " << l.new_name;
    call << l.new_name + ".fs->payload";
  }

  i=0;
  // and non-field vars
  for ( var_info & vi : var_info_list ) {
    if (!vi.is_loop_local) {
      std::string varname = "sv__" + std::to_string(i) + "_";
      if(vi.reduction_type != reduction::NONE) {
        if(target.CUDA){
          /* Reduction variables in a CUDA kernel are
           * saved to an array and reduced later. */
          kernel << ", " << vi.type << " * " << varname;
          call << ", r_" << vi.name;
          varname += "[Index]";
        }
      } else if(vi.is_assigned) {
        kernel << ", " << vi.type << " & " << varname;
        call << ", " << vi.name;
      } else {
        kernel << ", const " << vi.type << " " << varname;
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

  // finally, change the references to variables in the body
  replace_field_refs_and_funcs(loopBuf, true);

  // Begin the function
  kernel << ")\n{\n";

  if( target.CUDA ){
    /* Standard boilerplate in CUDA kernels: calculate site index */
    kernel << "int Index = threadIdx.x + blockIdx.x * blockDim.x "
           << " + lattice_info.loop_begin; \n";
    /* The last block may exceed the lattice size. Do nothing in that case. */
    kernel << "if(Index < lattice_info.loop_end) { \n";
    /* Initialize reductions */
    int i=0;
    for ( var_info & vi : var_info_list ) {
      if (!vi.is_loop_local) {
        if (vi.reduction_type == reduction::SUM) {
          kernel << "sv__" + std::to_string(i) + "_[Index]" << "=0;\n";
        } if (vi.reduction_type == reduction::PRODUCT) {
          kernel << "sv__" + std::to_string(i) + "_[Index]" << "=1;\n";
        }
        i++;
      }
    }
  } else {
    // Generate the loop
    kernel << "const int loop_begin = lattice.loop_begin(" << parity_in_this_loop << "); \n";
    kernel << "const int loop_end   = lattice.loop_end(" << parity_in_this_loop << "); \n";

    kernel << "for(int " << looping_var <<" = loop_begin; " 
           << looping_var << " < loop_end; " << looping_var << "++) {\n";
  }

  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());

    // First check for direction references. If any found, create list of temp
    // variables
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      kernel << type_name << " "  << l.loop_ref_name << "_d[NDIRS];\n";
      kernel << "std::array<bool, NDIRS> " << l.new_name << "_read_d;\n";
      kernel << l.new_name << "_read_d.fill(true);\n";
      break; // Only need one direction reference
    }
    // Check for references without a direction. If found, add temp variable
    for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
      kernel << type_name << " "  << l.loop_ref_name << ";\n";
      kernel << "bool "  << l.new_name << "_read = true;\n";
      break;  // Only one needed
    }
  }

  kernel << loopBuf.dump();
  if (semi_at_end) kernel << ';';
  kernel << '\n';

  // Call setters
  for (field_info & l : field_info_list) if(l.is_written){
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    kernel << l.new_name << ".set(" << l.loop_ref_name << ", " 
           << looping_var << ", lattice_info.field_alloc_size );\n";
  }

  kernel << "}\n}\n//----------\n";
  call << ");\n";

  if(target.CUDA){
    call << "check_cuda_error(\"" << kernel_name << "\");\n";
  }

  // Finally, emit the kernel
  // TheRewriter.InsertText(global.location.function, indent_string(kernel),true,true);
  toplevelBuf->insert(global.location.top.getLocWithOffset(-1), indent_string(kernel.str()),true,false);

  return call.str();
}


/// Change field references and special functions within loops
void MyASTVisitor::replace_field_refs_and_funcs(srcBuf & loopBuf, bool CUDA) {
  
  for ( field_ref & le : field_ref_list ) {
    //loopBuf.replace( le.nameExpr, le.info->loop_ref_name );
    if (le.dirExpr != nullptr) {
      //loopBuf.replace(parityExpr,
      //                 "lattice->neighb[" +  get_stmt_str(le.dirExpr) + "][" + looping_var + "]");
      loopBuf.replace(le.fullExpr, le.info->loop_ref_name+"_d["+get_stmt_str(le.dirExpr)+"]");
    } else {
      loopBuf.replace(le.fullExpr, le.info->loop_ref_name);
    }
  }

  // Handle calls to special in-loop functions
  for ( special_function_call & sfc : special_function_call_list ){
    if( sfc.add_loop_var ){
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"("+looping_var+")");
    } else {
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression);
    }
  }

  // Add a getter in above each field reference
  for (field_info & l : field_info_list) {
    for( field_ref *r : l.ref_list ){
      Expr *d = r->dirExpr;
      
      // If reference has a direction, use the temp list
      if(d){
        std::string dstring = get_stmt_str(d);
        std::string is_read = l.new_name + "_read_d["+dstring+"]";
        std::string getter;
        if(CUDA){
          getter = ".get(lattice_info.d_neighb[" + dstring + "]["
            + looping_var + "], lattice_info.field_alloc_size)";
        } else {
          getter = ".get_value_at(lattice->neighb[" + dstring + "][" 
            + looping_var + "])";
        }
        loopBuf.insert_before_stmt(r->fullExpr, 
          "if("  + is_read + ") {"
          + l.loop_ref_name + "_d[" + dstring + "]=" + l.new_name + getter + ";"
          + is_read + "=false;}", true, true);
      
      } else if(r->is_read) { // No direction, use the temp variable
      
        std::string is_read = l.new_name + "_read";
        std::string getter;
        if(CUDA){
          getter = ".get(" + looping_var + ", lattice_info.field_alloc_size)";
        } else {
          getter = ".get_value_at(" + looping_var + ")";
        }
        loopBuf.insert_before_stmt(r->fullExpr, 
          "if("  + is_read + ") {"
          + l.loop_ref_name + "=" + l.new_name + getter + ";"
          + is_read + "=false;}", true, true);
      }
    }
  }
}
