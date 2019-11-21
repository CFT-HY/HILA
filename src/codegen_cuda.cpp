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

extern std::string looping_var;
extern std::string parity_name;

extern std::string parity_in_this_loop;


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




std::string MyASTVisitor::generate_code_cuda(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {
  
  // "Code" is inserted at the location of the loop statement
  // and the kernel is build in "kernel"
  std::stringstream code, kernel;
  const std::string t = loopBuf.dump();

  // Get kernel name - use line number or file offset (must be deterministic)
  std::string kernel_name = MyASTVisitor::make_kernel_name();

  // Check for reductions and allocate device memory
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      // Allocate memory for a reduction. This will be filled in the kernel
      code << v.type << " * d_" << v.new_name << ";\n";
      code << "cudaMalloc( (void **)& d_" << v.new_name << ","
           << "sizeof(" << v.type << ") * lattice->volume() );\n";
      code << "check_cuda_error(\"allocate_reduction\");\n";
    }
  }

  
  kernel << "//----------\n";

  // Generate the function definition and call
  kernel << "__global__ void " << kernel_name << "( device_lattice_info d_lattice, ";
  code << "device_lattice_info lattice_info = lattice->device_info;\n";
  code << "lattice_info.loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "lattice_info.loop_end = lattice->loop_end(" << parity_in_this_loop << ");\n";
  code << "int N_blocks = (lattice_info.loop_end - lattice_info.loop_begin)/N_threads + 1;\n";
  code << kernel_name << "<<< N_blocks, N_threads >>>( lattice_info, ";


  // print field call list
  int i = -1;
  int j=0;
  for (field_info & l : field_info_list) {
    i++;
    
    if (i>0) {
      kernel << ", ";
      code   << ", ";
    }
    
    if (!l.is_written) kernel << "const ";
    // TODO: type to field_data
    kernel << "field_storage" << l.type_template << " " << l.new_name;
    code << l.new_name + ".fs->payload";

    for( field_ref *r : l.ref_list ){
      // Generate new name for direction expression
      if( r->dirExpr ) {
        r->dirname = "d_" + std::to_string(j);
        code << ", " << get_stmt_str(r->dirExpr);
        kernel << ", const int " << r->dirname;
        loopBuf.replace( r->dirExpr, r->dirname );
        j++;
      }
    }

  }

  i=0;
  // and non-field vars
  for ( var_info & vi : var_info_list ) {
    if (!vi.is_loop_local) {
      std::string varname = "sv__" + std::to_string(i) + "_";
      if(vi.reduction_type != reduction::NONE) {
        /* Reduction variables in a CUDA kernel are
         * saved to an array and reduced later. */
        kernel << ", " << vi.type << " * " << varname;
        code << ", d_r_" << vi.name;
        varname += "[Index]";
      } else if(vi.is_assigned) {
        kernel << ", " << vi.type << " & " << varname;
        code << ", " << vi.name;
      } else {
        kernel << ", const " << vi.type << " " << varname;
        code << ", " << vi.name;
      }
      
      for (var_ref & vr : vi.refs) {
        // replace_expr(ep, varname);
        loopBuf.replace( vr.ref, varname );
      }
      i++;
    }
  }

  // finally, change the references to variables in the body
  for ( field_ref & le : field_ref_list ) {
    //loopBuf.replace( le.nameExpr, le.info->loop_ref_name );
    if (le.dirExpr != nullptr) {
      loopBuf.replace(le.fullExpr, le.info->loop_ref_name+"_"+le.dirname);
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

  // Begin the function
  kernel << ")\n{\n";

  kernel << "device_lattice_info *lattice = &d_lattice; \n";
  /* Standard boilerplate in CUDA kernels: calculate site index */
  kernel << "int Index = threadIdx.x + blockIdx.x * blockDim.x "
         << " + lattice->loop_begin; \n";
  /* The last block may exceed the lattice size. Do nothing in that case. */
  kernel << "if(Index < lattice->loop_end) { \n";
  
  /* Initialize reductions */
  i=0;
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

  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());

    // First check for direction references. If any found, create list of temp
    // variables
    for( field_ref *r : l.ref_list ) if( r->dirExpr ) {
      // Generate new name for direction expression
      kernel << type_name << " "  << l.loop_ref_name << "_" << r->dirname
             << "=" << l.new_name << ".get(lattice->d_neighb[" 
             << r->dirname << "][" << looping_var 
             << "], lattice->field_alloc_size);\n";
    }
    // Check for references without a direction. If found, add temp variable
    for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
      kernel << type_name << " "  << l.loop_ref_name << "=" 
             << l.new_name << ".get(" << looping_var 
             << ", lattice->field_alloc_size)" << ";\n";
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
           << looping_var << ", lattice->field_alloc_size );\n";
  }

  kernel << "}\n}\n//----------\n";
  code << ");\n";

  code << "check_cuda_error(\"" << kernel_name << "\");\n";

  // Finally, emit the kernel
  // TheRewriter.InsertText(global.location.function, indent_string(kernel),true,true);
  toplevelBuf->insert(global.location.top.getLocWithOffset(-1), indent_string(kernel.str()),true,false);
  
  // Check reduction variables
  for (var_info & v : var_info_list) {
    // Run reduction
    if (v.reduction_type == reduction::SUM) {
      code << v.new_name << " = cuda_reduce_sum( d_"
           << v.new_name << ", lattice->volume()" <<  ");\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code << v.new_name << " = cuda_reduce_product( d_"
           << v.new_name << ", lattice->volume()" <<  ");\n";
    }
    // Free memory allocated for the reduction
    if (v.reduction_type != reduction::NONE) {
      code << "cudaFree(d_" << v.new_name << ");\n";
      code << "check_cuda_error(\"free_reduction\");\n";
    }
  }

  return code.str();
}



