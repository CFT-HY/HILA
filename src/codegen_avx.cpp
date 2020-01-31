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


namespace vector_map {

  static int bytes;

  /// Find the base element type in a field element
  static std::string type(std::string &original_type){
    if(original_type.find("double") != std::string::npos) {
      return "double";
    } else if(original_type.find("float") != std::string::npos){
      return "float";
    } else if(original_type.find("int") != std::string::npos){
      return "int";
    } else if(original_type.find("location") != std::string::npos){
      return "int";
    } else {
      llvm::errs() << "Cannot find vector size\n";
      llvm::errs() << "Field type " << original_type << "\n";
      exit(1);
    }
  }

  /// Get vector size based on the base element type
  static int size(std::string &original_type){
    if(original_type.find("double") != std::string::npos) {
      return bytes/8;
    } else if(original_type.find("float") != std::string::npos){
      return bytes/4;
    } else if(original_type.find("int") != std::string::npos){
      return bytes/4;
    } else if(original_type.find("location") != std::string::npos){
      return bytes/4;
    } else {
      llvm::errs() << "Cannot find vector size\n";
      llvm::errs() << "Field type " << original_type << "\n";
      exit(1);
    }
  }

  /// Replace base datatypes with vectorized datatypes of
  /// specific vector length
  static void replace(std::string & element_type, int vector_size) {
    size_t begin;
    begin = element_type.find("double");
    if(begin != std::string::npos){
      element_type.replace(begin, 6, "Vec"+std::to_string(vector_size)+"d");
    }
    begin = element_type.find("float");
    if(begin != std::string::npos){
      element_type.replace(begin, 5, "Vec"+std::to_string(vector_size)+"f");
    }
    begin = element_type.find("int");
    if(begin != std::string::npos){
      element_type.replace(begin, 3, "Vec"+std::to_string(vector_size)+"i");
    }
    begin = element_type.find("location");
    if(begin != std::string::npos){
      element_type.replace(begin, 8, "std::array<Vec"+std::to_string(vector_size)+"i,NDIM>");
    }
  }

  /// Replace base datatypes with vectorized datatypes
  static void replace(std::string & element_type) {
    replace(element_type, size(element_type));
  }

  /// Map general type strings into vectorized types
  static void size_and_type(std::string &original_type, std::string &base_type, std::string &vector_type, int & vector_size){
    if(original_type.find("double") != std::string::npos) {
      base_type = "double";
      vector_size = bytes/8;
      vector_type = "Vec"+std::to_string(vector_size)+"d";
    } else if(original_type.find("float") != std::string::npos){
      base_type = "float";
      vector_size = bytes/4;
      vector_type = "Vec"+std::to_string(vector_size)+"f";
    } else if(original_type.find("int") != std::string::npos){
      base_type = "int";
      vector_size = bytes/4;
      vector_type = "Vec"+std::to_string(vector_size)+"i";
    } else if(original_type.find("location") != std::string::npos){
      base_type = "int";
      vector_size = bytes/4;
      vector_type = "Vec"+std::to_string(vector_size)+"i";
    } else {
      llvm::errs() << "Cannot find vector size\n";
      llvm::errs() << "Field type " << original_type << "\n";
      exit(1);
    }
  }

  static void set_vector_target(codetype & target){
    bytes = target.vector_size;
  }
}





/// An AST walker for finding and handling variable declarations
/// in a loop function
class LoopFunctionHandler : public GeneralVisitor, public RecursiveASTVisitor<LoopFunctionHandler> {
public:
  using GeneralVisitor::GeneralVisitor;

  //Buffer for the function copy
  srcBuf functionBuffer;

  bool TraverseStmt(Stmt *s){
    RecursiveASTVisitor<LoopFunctionHandler>::TraverseStmt(s);
    return true;
  }
  bool VisitVarDecl(VarDecl *var);
  bool VisitCXXOperatorCallExpr(CXXOperatorCallExpr *op);
  bool VisitBinaryOperator(BinaryOperator *op);
};


bool LoopFunctionHandler::VisitVarDecl(VarDecl *var){
  std::string typestring = var->getType().getAsString();
  
  // This variable is an element, replace with vector
  if(typestring.rfind("element",0) != std::string::npos){
    std::string vector_type = typestring;
    vector_map::replace(vector_type);
    functionBuffer.replace(var->getSourceRange(), 
      vector_type+" "+var->getNameAsString() );
  } else {
    if(var->hasInit()){
      LoopAssignChecker lac(TheRewriter, Context);
      lac.TraverseStmt(var->getInit());
    }
  }
  return true;
}


bool LoopFunctionHandler::VisitCXXOperatorCallExpr(CXXOperatorCallExpr *op){
  if(op && op->isAssignmentOp()){
    std::string type = op->getArg(0)->getType().getAsString();
    type = remove_all_whitespace(type);
    if(type.rfind("element<",0) == std::string::npos){
      LoopAssignChecker lac(TheRewriter, Context);
      lac.TraverseStmt(op);
    }
  }
  return true;
}

bool LoopFunctionHandler::VisitBinaryOperator(BinaryOperator *op){
  if(op && op->isAssignmentOp()){
    std::string type = op->getLHS()->getType().getAsString();
    type = remove_all_whitespace(type);
    if(type.rfind("element<",0) == std::string::npos){
      LoopAssignChecker lac(TheRewriter, Context);
      lac.TraverseStmt(op);
    }
  }
  return true;
}




/// Replace element types with vector. Leaves other types untouched.
static bool replace_element_with_vector(SourceRange sr, std::string typestring, srcBuf &functionBuffer){
  if(typestring.rfind("element",0) != std::string::npos){
    std::string vector_type = typestring;
    vector_map::replace(vector_type);
    functionBuffer.replace(sr, vector_type);
    return true;
  }
  return false;
}


/// This should allow calling with element<>-type parameters from
/// loops. Generate a copy with elements replaced with vectors.
void MyASTVisitor::handle_loop_function_avx(FunctionDecl *fd) {
  SourceRange sr = fd->getSourceRange();
  FileID FID = TheRewriter.getSourceMgr().getFileID(sr.getBegin());
  srcBuf * sourceBuf = get_file_buffer(TheRewriter, FID);
  PrintingPolicy pp(Context->getLangOpts());
  vector_map::set_vector_target(target);

  LoopFunctionHandler lfh(TheRewriter, Context);

  // Track wether the function actually contains elements.
  // if not, no new function should be written
  bool contains_elements = false;

  // Copy the function to a buffer
  lfh.functionBuffer.copy_from_range(sourceBuf,sr);

  // Handle each parameter
  for( clang::ParmVarDecl *par : fd->parameters() ){
    std::string typestring = par->getType().getAsString(pp);
    contains_elements += replace_element_with_vector(par->getSourceRange(), typestring+" "+par->getNameAsString(), lfh.functionBuffer);
  }

  // Handle return type
  // Note: C++ cannot specialize only based on return type. Therefore we
  // only write a new function if the parameters contain elements
  std::string typestring = fd->getReturnType().getAsString(pp);
  replace_element_with_vector(fd->getReturnTypeSourceRange(), typestring, lfh.functionBuffer);

  lfh.TraverseStmt(fd->getBody());
  
  
  if(contains_elements) {
    std::string buffer = lfh.functionBuffer.dump();
    if( !(fd->hasBody()) ){
      // Declaration does not contain a body, needs a semicolon
      buffer += ";";
    }
    buffer += "\n";
    sourceBuf->insert(sr.getBegin(), buffer, true, true);
  }
}




std::string MyASTVisitor::generate_code_avx(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {
  std::stringstream code;

  // Find the vector size
  vector_map::set_vector_target(target);
  int vector_size = 1;
  std::string base_type, base_vector_type;
  std::string type = field_info_list.front().type_template;
  vector_map::size_and_type(type, base_type, base_vector_type, vector_size);
  
  // Create temporary variables for reductions
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      v.new_name = "v_" + v.reduction_name;
      // Allocate memory for a reduction. This will be filled in the kernel
      std::string reduction_type = v.type;
      vector_map::replace(reduction_type);
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
    vector_map::replace(type_name);

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
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"_"+base_vector_type+"("+looping_var+")");
    } else {
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"_"+base_vector_type+"()");
    }
  }

  for ( var_info & vi : var_info_list ) {
    if(vi.type.rfind("element",0) != std::string::npos) {
      // Converts all locally declared variables to vectors. Is this OK?
      // First get the type name in the declaration
      auto typeexpr = vi.decl->getTypeSourceInfo()->getTypeLoc().getSourceRange();
      std::string type_string = TheRewriter.getRewrittenText(typeexpr);

      vector_map::replace( type_string, vector_size );

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
  // Find the vector size
  vector_map::set_vector_target(target);
  int vector_size = 1;
  std::string base_type;
  std::string base_vector_type;
  vector_map::size_and_type(typestr, base_type, base_vector_type, vector_size);

  // insert after a new line
  SourceLocation l =
  getSourceLocationAtEndOfLine( field_storage_decl->getSourceRange().getEnd() );

  // Find template parameter name (only 1 allowed)
  auto template_parameter = field_storage_decl->getTemplateParameters()->begin()[0];
  std::string templ_type = template_parameter->getNameAsString();
  
  // Replace the original type with a vector type
  std::string vectortype = typestr;
  vector_map::replace(vectortype);

  // Get the body of the element definition
  srcBuf bodyBuffer;
  bodyBuffer.copy_from_range(writeBuf,field_storage_decl->getTemplatedDecl()->getSourceRange());

  // Replace templated type with new vector type
  llvm::errs() << "TYPES " << templ_type << " " << vectortype << "\n";
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
  field_element_code << "struct field_info<"<<typestr<<"> {\n";
  field_element_code << " constexpr static int vector_size = "<< std::to_string(vector_size) <<";\n";
  field_element_code << " constexpr static int base_type_size = sizeof(" << base_type << ");\n";
  field_element_code << " constexpr static int elements = sizeof(" << typestr 
                     << ")/sizeof(" << base_type << ");\n";
  field_element_code << " using base_type = " << base_type << ";\n";
  field_element_code << " using vector_type = " << base_vector_type << ";\n";
  field_element_code << "};\n";

  bodyBuffer.append(field_element_code.str(), true);


  writeBuf->insert(l, "\n"+bodyBuffer.dump(), true, false);

  // Mark source buffer modified
  // set_sourceloc_modified( l );
}


