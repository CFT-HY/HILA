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

static int vector_lenght = 4;


/// Replace base datatypes with vectorized datatypes
void replace_basetype_with_vector(std::string & element_type) {
  size_t begin;
  begin = element_type.find("double");
  if(begin != std::string::npos){
    element_type.replace(begin, 6, "Vec4d");
  }
  begin = element_type.find("float");
  if(begin != std::string::npos){
    element_type.replace(begin, 5, "Vec8f");
  }
  begin = element_type.find("int");
  if(begin != std::string::npos){
    element_type.replace(begin, 3, "Vec16i");
  }
}
  


std::string MyASTVisitor::generate_code_avx(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {
  std::stringstream code;
  
  // Create temporary variables for reductions
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      v.new_name = "v_" + v.reduction_name;
      // Allocate memory for a reduction. This will be filled in the kernel
      std::string reduction_type = v.type;
      replace_basetype_with_vector(reduction_type);
      if (v.reduction_type == reduction::SUM) {
        code << reduction_type << " " << v.new_name << "(0);\n";
      } else if (v.reduction_type == reduction::PRODUCT) {
        code << reduction_type << " " << v.new_name << "(1);\n";
      }
    }
  }

  // Set loop lattice
  std::string fieldname = field_info_list.front().old_name;
  code << "vectorized_lattice_struct * loop_lattice = "
       << fieldname << ".fs->lattice;\n";

  // Get a pointer to the neighbour list
  code << "std::array<unsigned*,NDIRS> neighbour_list = loop_lattice->neighbour_list();\n";

  // Set the start and end points
  // A single vector covers 4 sites in AVX. These must all have the same parity.
  // So we devide the start and end points here by 4.
  code << "const int loop_begin = loop_lattice->loop_begin(" << parity_in_this_loop << ");\n";
  code << "const int loop_end   = loop_lattice->loop_end(" << parity_in_this_loop << ");\n";

  // Start the loop
  code << "for(int " << looping_var <<" = loop_begin; " 
       << looping_var << " < loop_end; " 
       << looping_var << "++) {\n";

  
  // Create temporary field element variables
  for (field_info & l : field_info_list) {
    std::string type_name = l.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());
    replace_basetype_with_vector(type_name);

    // First check for direction references. If any found, create list of temp
    // variables
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      code << type_name << " " << l.loop_ref_name << "_" << get_stmt_str(d.e)
           << " = " << l.new_name << ".get_value_at(neighbour_list[" 
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

  for ( var_info & vi : var_info_list ) {
    if(vi.type.rfind("element",0) != std::string::npos) {
      // Converts all locally declared variables to vectors. Is this OK?
      // First get the type name in the declaration
      auto typeexpr = vi.decl->getTypeSourceInfo()->getTypeLoc().getSourceRange();
      std::string type_string = TheRewriter.getRewrittenText(typeexpr);

      replace_basetype_with_vector( type_string );

      loopBuf.replace( typeexpr, type_string + " ");
    }
    if (vi.reduction_type != reduction::NONE) {
      // Replace references in the loop body
      for (var_ref & vr : vi.refs) {
        loopBuf.replace( vr.ref, vi.new_name );
      }
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


  // Final reduction of the temporary reduction variables
  for (var_info & v : var_info_list) {
    if (v.reduction_type == reduction::SUM) {
      code << v.reduction_name << " = reduce_sum(" << v.new_name << ");\n";
    } else if (v.reduction_type == reduction::PRODUCT) {
      code << v.reduction_name << " = reduce_prod(" << v.new_name << ");\n";
    }
  }

  return code.str();
}


void MyASTVisitor::generate_field_storage_type_AVX(std::string typestr){
  // insert after a new line
  SourceLocation l =
  getSourceLocationAtEndOfLine( field_storage_decl->getSourceRange().getEnd() );

  // Find template parameter name (only 1 allowed)
  auto template_parameter = field_storage_decl->getTemplateParameters()->begin()[0];
  std::string templ_type = template_parameter->getNameAsString();
  
  // Replace the original type with a vector type
  std::string vectortype = typestr;
  replace_basetype_with_vector(vectortype);

  // Get the body of the element definition
  srcBuf bodyBuffer; // (&TheRewriter,S);
  bodyBuffer.copy_from_range(writeBuf,field_storage_decl->getTemplatedDecl()->getSourceRange());

  // Replace templated type with new vector type
  bodyBuffer.replace_token(0, bodyBuffer.size()-1, templ_type, vectortype );

  // Add specialization parameters
  bodyBuffer.replace_token(0, bodyBuffer.size()-1,
                  field_storage_decl->getQualifiedNameAsString(),
                  field_storage_decl->getQualifiedNameAsString() + "<"+typestr+">");

  // Add the template<> declaration
  bodyBuffer.prepend("template<>\n", true);
  bodyBuffer.append(";\n", true); // semicolon at end


  // Add a simple template mapping from element type to vector size
  std::stringstream field_element_code;
  field_element_code << "template<> \n";
  field_element_code << "struct field_element<"<<typestr<<"> {\n";
  field_element_code << "  "<< vectortype << " c;\n";
  field_element_code << "};\n";

  bodyBuffer.append(field_element_code.str(), true);


  writeBuf->insert(l, "\n"+bodyBuffer.dump(), true, false);

  // Mark source buffer modified
  // set_sourceloc_modified( l );
}


