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
#include "toplevelvisitor.h"
#include "stringops.h"

extern std::string looping_var;
extern std::string parity_name;
extern std::string parity_in_this_loop;

// add pragma
void TopLevelVisitor::handle_loop_function_openacc(FunctionDecl *fd) {
  
  SourceLocation sl = fd->getSourceRange().getBegin();
  srcBuf * sb = get_file_srcBuf(sl);
  if (sb != nullptr)
    sb->insert(sl, "#pragma acc routine \n",true,true);
}
// Add pragma to constructor too?
void TopLevelVisitor::handle_loop_constructor_openacc(CXXConstructorDecl *fd) {
  
  SourceLocation sl = fd->getSourceRange().getBegin();
  srcBuf * sb = get_file_srcBuf(sl);
  if (sb != nullptr)
    sb->insert(sl, "#pragma acc routine \n",true,true);
}

void TopLevelVisitor::generate_openacc_loop_header(std::stringstream & code) {

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
}


