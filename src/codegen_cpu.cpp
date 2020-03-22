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


std::string MyASTVisitor::generate_code_cpu(Stmt *S, bool semicolon_at_end, srcBuf & loopBuf) {
  std::stringstream code;

  // Set loop lattice
  std::string fieldname = field_info_list.front().new_name;
  code << "lattice_struct * loop_lattice = " << fieldname << ".fs->lattice;\n";
  
  // Set the start and end points
  code << "const int loop_begin = loop_lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "const int loop_end   = loop_lattice->loop_end(" << parity_in_this_loop << ");\n";

  // and the openacc loop header
  if (target.openacc) generate_openacc_loop_header(code);

  // Start the loop
  code << "for(int " << looping_var <<" = loop_begin; " 
       << looping_var << " < loop_end; " << looping_var << "++) {\n";


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
  
  
  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());

    // First check for direction references. If any found, create list of temp
    // variables
    if (l.is_read_nb) {
      for (dir_ptr & d : l.dir_list) {
        std::string dirname;
        if (d.is_constant_direction) dirname = d.direxpr_s;  // orig. string
        else dirname = remove_X( loopBuf.get(d.parityExpr->getSourceRange()) ); // mapped name was get_stmt_str(d.e);

        // generate access stmt
        code << type_name << " " << d.name_with_dir
             << " = " << l.new_name << ".get_value_at(lattice->neighb[" 
             << dirname << "][" << looping_var << "]);\n";

        // and replace references in loop body
        for (field_ref * ref : d.ref_list) {
          loopBuf.replace(ref->fullExpr, d.name_with_dir);
        }
      }
    }
    
    // and then get (possible) local refs
    if (l.is_read_atX) {
      // now reading var without nb. reference
      code << type_name << " " << l.loop_ref_name << " = " 
           << l.new_name << ".get_value_at(" << looping_var << ");\n";

    } else if (l.is_written) {
      code << type_name << " " << l.loop_ref_name << ";\n";
    }

    // and finally replace references in body 
    for (field_ref * ref : l.ref_list) if (!ref->is_direction) {
      loopBuf.replace(ref->fullExpr, l.loop_ref_name);
    }
  }



  // // Replace field references in loop body
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
    if( sfc.add_loop_var ){
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"("+looping_var+")");
    } else {
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"()");
    }
  }

  // Dump the main loop code here
  code << loopBuf.dump();
  if (semicolon_at_end) code << ";";
  code << "\n";


  // Add calls to setters 
  for (field_info & l : field_info_list) if(l.is_written) {
    code << l.new_name << ".set_value_at(" << l.loop_ref_name << ", " 
         << looping_var << ");\n";
  }

  code << "}\n";

  return code.str();
}

