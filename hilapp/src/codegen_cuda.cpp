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

#include "hilapp.h"
#include "myastvisitor.h"
#include "stringops.h"

extern std::string looping_var;
extern std::string parity_name;

extern std::string parity_in_this_loop;


// Add the __host__ __device__ keywords to functions called a loop
void MyASTVisitor::handle_loop_function_cuda(FunctionDecl *fd) {
  
  SourceLocation sl = fd->getSourceRange().getBegin();
  srcBuf * sb = get_file_srcBuf(sl);
  if (sb == nullptr) {
    // it's a system file -- should we do something?
    return;
  }
  sb->insert(sl, "__device__ __host__ ",true,true);
}

void MyASTVisitor::handle_loop_constructor_cuda(CXXConstructorDecl *fd) {
  
  SourceLocation sl = fd->getSourceRange().getBegin();
  srcBuf * sb = get_file_srcBuf(sl);
  if (sb == nullptr) {
    // it's a system file -- should we do something?
    return;
  }
  sb->insert(sl, "__device__ __host__ ",true,true);
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



std::string MyASTVisitor::generate_code_cuda(Stmt *S, bool semicolon_at_end, srcBuf & loopBuf, bool generate_wait_loops) {
  
  // "Code" is inserted at the location of the loop statement
  // and the kernel is build in "kernel"
  std::stringstream code, kernel;
  const std::string t = loopBuf.dump();

  // Get kernel name - use line number or file offset (must be deterministic)
  std::string kernel_name = MyASTVisitor::make_kernel_name();


  // Wait for the communication to finish
  for (field_info & l : field_info_list) {
    // If neighbour references exist, communicate them
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      code << l.new_name << ".wait_fetch("
           << d.direxpr_s << ", " << parity_in_this_loop << ");\n";
    }
  }

  // Set loop lattice
  if (field_info_list.size() > 0) {
    std::string fieldname = field_info_list.front().old_name;
    code << "lattice_struct * loop_lattice = " << fieldname << ".fs->lattice;\n";
  } else {
    // now no fields in loop - default lattice
    code << "lattice_struct * loop_lattice = lattice;\n";
  }

  for (vector_reduction_ref & vrf : vector_reduction_ref_list) {
    // Allocate memory for a reduction and initialize
    code << vrf.type << " * d_" << vrf.vector_name << ";\n";
    code << "cudaMalloc( (void **)& d_" << vrf.vector_name << ", "
         << vrf.vector_name << ".size() * sizeof("
         << vrf.type << ") * lattice->volume() );\n";
    if (vrf.reduction_type == reduction::SUM) {
      code << "cuda_set_zero(d_" << vrf.vector_name
           << ", " << vrf.vector_name << ".size()* lattice->volume());\n";
    }
    if (vrf.reduction_type == reduction::PRODUCT) {
      code << "cuda_set_one(d_" << vrf.vector_name
           << ", " << vrf.vector_name << ".size()* lattice->volume());\n";
    }
    code << "check_cuda_error(\"allocate_reduction\");\n";
  }

  
  kernel << "//----------\n";

  // Generate the function definition and call
  kernel << "__global__ void " << kernel_name << "( backend_lattice_struct d_lattice";
  code << "backend_lattice_struct lattice_info = *(lattice->backend_lattice);\n";
  code << "lattice_info.loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "lattice_info.loop_end = lattice->loop_end(" << parity_in_this_loop << ");\n";
  code << "int N_blocks = (lattice_info.loop_end - lattice_info.loop_begin)/N_threads + 1;\n";


  // Check for reductions and allocate device memory
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      // Allocate memory for a reduction. This will be filled in the kernel
      code << v.type << " * d_" << v.reduction_name << ";\n";
      code << "cudaMalloc( (void **)& d_" << v.reduction_name << ","
           << "sizeof(" << v.type << ") * N_blocks );\n";
      code << "check_cuda_error(\"allocate_reduction\");\n";
      if (v.reduction_type == reduction::SUM) {
        code << "cuda_set_zero(d_" << v.reduction_name
           << ", N_blocks);\n";
      }
      if (v.reduction_type == reduction::PRODUCT) {
        code << "cuda_set_one(d_" << v.reduction_name
           << ", N_blocks);\n";
      }
    }
  }


  code << kernel_name << "<<< N_blocks, N_threads >>>( lattice_info";


  // print field call list
  int i = 0;
  for (field_info & l : field_info_list) {
    
    kernel << ", ";
    code   << ", ";
    
    if (!l.is_written) kernel << "const ";
    kernel << "field_storage" << l.type_template << " " << l.new_name;
    code << l.new_name + ".fs->payload";
    i++;
  }

  i=0;
  // and non-field vars
  for ( var_info & vi : var_info_list ) if(!vi.is_loop_local) {
    // Rename the variable
    vi.new_name = "sv_" + std::to_string(i) + "_";
    i++;

    if(vi.reduction_type != reduction::NONE) {
      // Generate a temporary array for the reduction 
      kernel << ", " << vi.type << " * " << vi.new_name;
      code << ", d_r_" << vi.name;
    } else {
      kernel << ", const " << vi.type << " " << vi.new_name;
      code << ", " << vi.name;
      // Replace references in the loop body
      for (var_ref & vr : vi.refs) {
        loopBuf.replace( vr.ref, vi.new_name );
      }
    }
  }

  // Add array reductions to the argument list
  for (vector_reduction_ref & vrf : vector_reduction_ref_list) {
    if (vrf.reduction_type != reduction::NONE) {
      kernel << ", " << vrf.type << " * " << vrf.vector_name;
      code << ", d_" << vrf.vector_name;
      vrf.new_vector_name = vrf.vector_name
              + "[ (loop_lattice->loop_end - loop_lattice->loop_begin)*"
              + vrf.index_name + " + Index]";
    }
    loopBuf.replace( vrf.ref, vrf.new_vector_name );
  }

  // In kernelized code we need to handle array expressions as well
  for ( array_ref & ar : array_ref_list ) {
    // Rename the expression
    ar.new_name = "sv_" + std::to_string(i) + "_";
    i++;
    kernel << ", " << ar.type << " " << ar.new_name ;
    code << ", " << get_stmt_str(ar.ref);
    loopBuf.replace( ar.ref, ar.new_name );
  }
  

  // // change the references to field expressions in the body
  // for ( field_ref & le : field_ref_list ) {
  //   //loopBuf.replace( le.nameExpr, le.info->loop_ref_name );
  //   if (le.dirExpr != nullptr) {
  //     loopBuf.replace(le.fullExpr, le.info->loop_ref_name+"_"+le.dirname);
  //   } else {
  //     loopBuf.replace(le.fullExpr, le.info->loop_ref_name);
  //   }
  // }


  // Handle calls to special in-loop functions
  for ( special_function_call & sfc : special_function_call_list ){
    std::string repl = sfc.replace_expression;  // comes with ( now
    if ( sfc.add_loop_var ) {
      repl += looping_var;
      if ( sfc.argsExpr != nullptr) repl += ',';
    }
    loopBuf.replace(sfc.replace_range, repl);    
  }


  // Begin the function
  kernel << ")\n{\n";

  kernel << "backend_lattice_struct *loop_lattice = &d_lattice; \n";
  /* Standard boilerplate in CUDA kernels: calculate site index */
  kernel << "int Index = threadIdx.x + blockIdx.x * blockDim.x "
         << " + loop_lattice->loop_begin; \n";
  /* The last block may exceed the lattice size. Do nothing in that case. */
  kernel << "if(Index < loop_lattice->loop_end) { \n";


  // Declare the shared reduction variable
  for ( var_info & vi : var_info_list ) if(!vi.is_loop_local) {
    if(vi.reduction_type != reduction::NONE) {
      // Generate a temporary array for the reduction 
      kernel << "__shared__ " << vi.type << " " << vi.new_name
             << "_sh[N_threads];\n";

      // Initialize only the local element
      if (vi.reduction_type == reduction::SUM) {
        kernel << vi.new_name << "_sh[threadIdx.x] = 0;\n";
      } else if (vi.reduction_type == reduction::PRODUCT) {
        kernel << vi.new_name << "_sh[threadIdx.x] = 1;\n";
      }

      // Replace references in the loop body
      for (var_ref & vr : vi.refs) {
        loopBuf.replace( vr.ref, vi.new_name+"_sh[threadIdx.x]" );
      }
    }
  }
  

  
  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    if (l.is_read_nb) {
      // this field is nn-read
      for (dir_ptr & d : l.dir_list) {
        std::string dirname;
        if (d.is_constant_direction) dirname = d.direxpr_s;  // orig. string
        else dirname = remove_X( loopBuf.get(d.parityExpr->getSourceRange()) ); // mapped name was get_stmt_str(d.e);

        // Check if the direction is a variable. These have been renamed.
        // for ( var_info & vi : var_info_list) for ( var_ref & vr : vi.refs )
        //   if( vr.ref == d.e ) 
        //     dirname = vi.new_name;

        // Create the temp variable and call the getter
        kernel << l.element_type << " "  << d.name_with_dir 
               << " = " << l.new_name << ".get(" <<
               l.new_name << ".neighbours[" << dirname
               << "][" << looping_var 
               << "], loop_lattice->field_alloc_size);\n";

        // and replace references in loop body
        for (field_ref * ref : d.ref_list) {
          loopBuf.replace(ref->fullExpr, d.name_with_dir); 
        }
      }     
    }

    if (l.is_read_atX) {
      // local read
      kernel << l.element_type << " "  << l.loop_ref_name << " = " 
             << l.new_name << ".get(" << looping_var 
             << ", loop_lattice->field_alloc_size);\n";

    } else if (l.is_written) {
      // and a var which is not read
      kernel << l.element_type << " "  << l.loop_ref_name << ";\n";
    }

    // and finally replace references in body 
    for (field_ref * ref : l.ref_list) if (!ref->is_direction) {
      loopBuf.replace(ref->fullExpr, l.loop_ref_name);
    }

    // // First create temp variables fields fetched from a direction
    // for( field_ref *r : l.ref_list ) if( r->dirExpr ) {
    //   std::string dirname = get_stmt_str(r->dirExpr);

    //   // Check if the direction is a variable. These have been renamed.
    //   for ( var_info & vi : var_info_list) for ( var_ref & vr : vi.refs )
    //     if( vr.ref == r->dirExpr ) {
    //       dirname = vi.new_name;
    //   }
    //   // Create the temp variable and call the getter
    //   kernel << l.element_type << " "  << l.loop_ref_name << "_" << r->dirname
    //          << "=" << l.new_name << ".get(loop_lattice->d_neighb[" 
    //          << dirname << "][" << looping_var 
    //          << "], loop_lattice->field_alloc_size);\n";
    // }
    // // Check for references without a direction. If found, add temp variable
    // for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
    //   kernel << l.element_type << " "  << l.loop_ref_name << "=" 
    //          << l.new_name << ".get(" << looping_var 
    //          << ", loop_lattice->field_alloc_size)" << ";\n";
    //   break;  // Only one needed
    // }

  } 

  // Dump the loop body 
  kernel << loopBuf.dump();
  if (semicolon_at_end) kernel << ';';
  kernel << '\n';


  // Call setters
  for (field_info & l : field_info_list) if(l.is_written){
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    kernel << l.new_name << ".set(" << l.loop_ref_name << ", " 
           << looping_var << ", loop_lattice->field_alloc_size );\n";
  }
  
  // Handle reductions: Need to sync threads once, then do reduction
  // locally once per block
  bool sync_done = false;
  for ( var_info & vi : var_info_list ) if(!vi.is_loop_local) {
    if(vi.reduction_type != reduction::NONE) {
      // Do sync (only if there is a reduction)
      if(!sync_done){
        kernel << "__syncthreads();\n";
        sync_done = true;
      }

      //Now run the thread level reduction
      kernel << "if( threadIdx.x == 0 ){\n";
      if (vi.reduction_type == reduction::SUM) {
        kernel << vi.new_name << "[blockIdx.x] = 0;\n";
      } else if (vi.reduction_type == reduction::PRODUCT) {
        kernel << vi.new_name << "[blockIdx.x] = 1;\n";
      }
      kernel << "for( int i=0; i<N_threads; i++ ){\n";
      if (vi.reduction_type == reduction::SUM) {
        kernel << vi.new_name << "[blockIdx.x] += " << vi.new_name << "_sh[i];\n";
      } else if (vi.reduction_type == reduction::PRODUCT) {
        kernel << vi.new_name << "[blockIdx.x] *= " << vi.new_name << "_sh[i];\n";
      }
      kernel << "}\n";
      kernel << "}\n";
    }
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
      code << v.reduction_name << " = cuda_reduce_sum( d_"
           << v.reduction_name << ", N_blocks" <<  ");\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code << v.reduction_name << " = cuda_reduce_product( d_"
           << v.reduction_name << ", N_blocks" <<  ");\n";
    }
    // Free memory allocated for the reduction
    if (v.reduction_type != reduction::NONE) {
      code << "cudaFree(d_" << v.reduction_name << ");\n";
      code << "check_cuda_error(\"free_reduction\");\n";
    }
  }
  for (vector_reduction_ref & vrf : vector_reduction_ref_list) {
    if (vrf.reduction_type == reduction::SUM) {
      code << "cuda_multireduce_sum( " << vrf.vector_name 
           << ", d_" << vrf.vector_name 
           << ", loop_lattice->volume() );\n";

    } if (vrf.reduction_type == reduction::PRODUCT) {
      code << "cuda_multireduce_mul( " << vrf.vector_name 
           << ", d_" << vrf.vector_name 
           << ", loop_lattice->volume() );\n";
    }
    if (vrf.reduction_type != reduction::NONE) {
      code << "cudaFree(d_" << vrf.vector_name << ");\n";
      code << "check_cuda_error(\"free_reduction\");\n";
    }
  }

  return code.str();
}



