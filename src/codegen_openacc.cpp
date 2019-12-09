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


// Add the __host__ __device__ keywords to functions called a loop
void MyASTVisitor::handle_loop_function_openacc(SourceLocation sl) {
  FileID FID = TheRewriter.getSourceMgr().getFileID(sl);
  // set_fid_modified(FID);
  srcBuf * sb = get_file_buffer(TheRewriter, FID);
  sb->insert(sl, "#pragma acc routine \n",true,true);
}



std::string MyASTVisitor::generate_code_openacc(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {

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
    
  std::stringstream code;
  code << "const int loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "const int loop_end   = lattice->loop_end(" << parity_in_this_loop << ");\n";

  // Add openacc pragmas
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
      code << type_name << " " << l.loop_ref_name << "_" << get_stmt_str(d.e)
           << " = " << l.new_name << ".get_value_at(lattice->neighb[" 
           << get_stmt_str(d.e) << "][" << looping_var << "]);\n";
    }
    // Check for references without a direction. If found, add temp variable
    for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
      code << type_name << " " << l.loop_ref_name << " = " 
           << l.new_name << ".get_value_at(" << looping_var << ");\n";
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


