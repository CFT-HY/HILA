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



/// Replace base datatypes with vectorized datatypes
void replace_basetype_with_vector(std::string & element_type) {
  size_t begin;
  begin = element_type.find("double");
  if(begin != std::string::npos){
    element_type.replace(begin, 6, "__m256d");
  }
  begin = element_type.find("float");
  if(begin != std::string::npos){
    element_type.replace(begin, 6, "__m256f");
  }
}
  


std::string MyASTVisitor::generate_code_avx(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {
  std::stringstream code;
  
  // Set the start and end points
  code << "const int loop_begin = lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "const int loop_end   = lattice->loop_end(" << parity_in_this_loop << ");\n";

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
    } else {
      // locally declared variables should be vector type if possible
      std::string declaration_string = TheRewriter.getRewrittenText(vi.decl->getSourceRange());
      replace_basetype_with_vector( declaration_string );
      loopBuf.replace( vi.decl->getSourceRange(), declaration_string);
    }
  }

  
  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    replace_basetype_with_vector(type_name);

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


  // Replace field references in loop body
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

  // Dump the main loop code here
  code << loopBuf.dump();
  if (semi_at_end) code << ";";
  code << "\n";


  // Add cals to setters 
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


void MyASTVisitor::generate_field_element_type_AVX(std::string typestr){
  // insert after a new line
  SourceLocation l =
  getSourceLocationAtEndOfLine( element_decl->getSourceRange().getEnd() );

  // Find template parameter name (only 1 allowed)
  auto template_parameter = element_decl->getTemplateParameters()->begin()[0];
  std::string templ_type = template_parameter->getNameAsString();
  
  // Replace the original type with a vector type
  std::string vectortype = typestr;
  replace_basetype_with_vector(vectortype);

  // Get the body of the element definition
  srcBuf bodyBuffer; // (&TheRewriter,S);
  bodyBuffer.copy_from_range(writeBuf,element_decl->getTemplatedDecl()->getSourceRange());

  bodyBuffer.replace_token(0, bodyBuffer.size()-1, templ_type, vectortype );

  // Add the template<> declaration
  bodyBuffer.prepend("template<" + typestr + ">\n", true);

  writeBuf->insert(l, "\n"+bodyBuffer.dump(), true, false);
}
