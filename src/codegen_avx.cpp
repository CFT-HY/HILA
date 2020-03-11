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






/// An AST walker for finding and handling variable declarations
/// in a loop function
class LoopFunctionHandler : public GeneralVisitor, public RecursiveASTVisitor<LoopFunctionHandler> {
public:
  using GeneralVisitor::GeneralVisitor;

  //Buffer for the function copy
  srcBuf functionBuffer;
  int vector_size;

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
  
  size_t begin;
  begin = typestring.rfind("element",0);
  if(begin != std::string::npos){
    // This variable is an element, replace with vector
    std::string vector_type = typestring;
    if(begin != std::string::npos){
      vector_type.replace(begin, 7, "vectorize_struct");
      vector_type.replace(vector_type.find_last_of(">"), 1, ", "+std::to_string(vector_size)+">::type");
    }

    if( var->isDirectInit() ){
      std::string init = TheRewriter.getRewrittenText(var->getInit()->getSourceRange());
      functionBuffer.replace(var->getSourceRange(), vector_type+" "+var->getNameAsString() + "=" + init );
    } else {
      functionBuffer.replace(var->getSourceRange(), vector_type+" "+var->getNameAsString());
    }

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
static void replace_element_with_vector(SourceRange sr, std::string typestring, std::string namestring, int vector_size, srcBuf &functionBuffer){
  if(typestring.rfind("element",0) != std::string::npos){
    std::string vector_type = typestring;
    size_t begin;
    begin = vector_type.find("element");
    if(begin != std::string::npos){
      vector_type.replace(begin, 7, "vectorize_struct");
      vector_type.replace(vector_type.find_last_of(">"), 1, ", "+std::to_string(vector_size)+">::type");
    }
    
    functionBuffer.replace(sr, vector_type+" "+namestring);
  }
}


/// This should allow calling with element<>-type parameters from
/// loops. Generate a copy with elements replaced with vectors.
void MyASTVisitor::handle_loop_function_avx(FunctionDecl *fd) {
  SourceRange sr = fd->getSourceRange();
  FileID FID = TheRewriter.getSourceMgr().getFileID(sr.getBegin());
  srcBuf * sourceBuf = get_file_buffer(TheRewriter, FID);
  PrintingPolicy pp(Context->getLangOpts());


  // Track wether the function actually contains elements.
  // if not, no new function should be written
  bool generate_function = false;

  // Check allowed vector sizes
  int smallest=1, largest=0;
  for( clang::ParmVarDecl *par : fd->parameters() ){
    std::string typestring = par->getType().getAsString(pp);
    if(typestring.find("double") != std::string::npos){
      smallest = 4;
      largest = 8;
    }
    if( typestring.find("float") != std::string::npos ||
        typestring.find("int") != std::string::npos || 
        typestring.find("coordinate_vector") != std::string::npos ){
      smallest = 8;
      largest = 16;
    }

    // Check if there are elements in the first place
    if(typestring.find("element") != std::string::npos)
      generate_function = true;

  }

  if(generate_function) for( int vector_size = smallest; vector_size <= largest; vector_size*=2 ){
    LoopFunctionHandler lfh(TheRewriter, Context);
    lfh.functionBuffer.copy_from_range(sourceBuf,sr);
    lfh.vector_size = vector_size;

    // Handle each parameter
    for( clang::ParmVarDecl *par : fd->parameters() ){
      std::string typestring = par->getType().getAsString(pp);
      replace_element_with_vector(par->getSourceRange(), typestring, par->getNameAsString(), vector_size, lfh.functionBuffer);
    }

    // Handle return type
    // Note: C++ cannot specialize only based on return type. Therefore we
    // only write a new function if the parameters contain elements
    std::string typestring = fd->getReturnType().getAsString(pp);
    replace_element_with_vector(fd->getReturnTypeSourceRange(), typestring, "", vector_size, lfh.functionBuffer);

    lfh.TraverseStmt(fd->getBody());

    std::string buffer = lfh.functionBuffer.dump();
    if( !(fd->hasBody()) ){
      // Declaration does not contain a body, needs a semicolon
      buffer += ";";
    }
    buffer += "\n";
    sourceBuf->insert(sr.getBegin(), buffer, true, true);
  }
}




std::string MyASTVisitor::generate_code_avx(Stmt *S, bool semicolon_at_end, srcBuf & loopBuf) {
  std::stringstream code;

  // The base type of the loop is the base type of the first variable
  std::string basetype = "loop_base_type";
  code << "using " << basetype << " = typename base_type_struct"
       << field_info_list.front().type_template+"::type;\n";
  code << "constexpr int vector_size = vector_info<" << basetype << ">::vector_size;\n";
  
  // Create temporary variables for reductions
  for (var_info & v : var_info_list) {
    if (v.reduction_type != reduction::NONE) {
      v.new_name = "v_" + v.reduction_name;
      // Allocate memory for a reduction. This will be filled in the kernel
      if (v.reduction_type == reduction::SUM) {
        code << "vector_type<" << v.type << "> " << v.new_name << "(0);\n";
      } else if (v.reduction_type == reduction::PRODUCT) {
        code << "vector_type<" << v.type << "> " << v.new_name << "(1);\n";
      }
    }
  }

  // Set loop lattice
  std::string fieldname = field_info_list.front().new_name;
  code << "auto * loop_lattice = "
       << fieldname << ".fs->lattice->backend_lattice->get_vectorized_lattice<vector_size>();\n";

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

    // First check for direction references. If any found, create list of temp
    // variables
    for (dir_ptr & d : l.dir_list) if(d.count > 0){
      code << "vector_type" << type_name << " " 
           << l.loop_ref_name << "_" << get_stmt_str(d.e)
           << " = " << l.new_name << ".get_value_at(neighbour_list[" 
           << get_stmt_str(d.e) << "][" << looping_var << "]);\n";
    }
    // Check for references without a direction. If found, add temp variable
    bool local_ref = false;
    bool local_is_read = false;
    for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
      local_ref = true;
      if( r->is_read ){
        local_is_read = true;
      }
    }
    if(local_ref) {
      if( local_is_read ){
        code << "vector_type" << type_name << " "
             << l.loop_ref_name << " = " 
             << l.new_name << ".get_value_at(" << looping_var << ");\n";
      } else {
        code << "vector_type" << type_name << " "
             << l.loop_ref_name << ";\n";
      }
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
    if( sfc.name == "coordinates" ) {
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"("+looping_var+")");
    } else if( sfc.add_loop_var ) {
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"_vector<" + basetype + ">("+looping_var+")");
    } else {
      loopBuf.replace(sfc.fullExpr, sfc.replace_expression+"_vector<" + basetype + ">()");
    }
  }

  for ( var_info & vi : var_info_list ) {
    if(vi.type.rfind("element",0) != std::string::npos) {
      // Converts all locally declared variables to vectors.
      // Basically replaces "element" with "vector_type"
      auto typeexpr = vi.decl->getTypeSourceInfo()->getTypeLoc().getSourceRange();
      std::string type_string = TheRewriter.getRewrittenText(typeexpr);

      size_t begin;
      begin = type_string.find("element");
      if(begin != std::string::npos){
        type_string.replace(begin, 7, "vectorize_struct");
        type_string.replace(type_string.find_last_of(">"), 1, ", vector_size>::type");
      }

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
  if (semicolon_at_end) code << ";";
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




