#include <sstream>
#include <string>

#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"

#include "hilapp.h"
#include "myastvisitor.h"
#include "stringops.h"



///////////////////////////////////////////////////////////////////////////////////
/// Inspect the type argument of Field<type>: is it good for vectors?
///
/// Type is vectorizable if:
/// a) just float, double, int  or 
/// b) is templated type, with float/double in template and  implements 
///    the method using base_type = typename base_type_struct<T>::type;
/// 
///////////////////////////////////////////////////////////////////////////////////

/// type alias decl visitor is used to locate
///    using = base_type = typename base_type_struct<T>::type;
/// expressions



number_type get_number_type( const std::string & s ) {
  if      (s == "int")         return number_type::INT;
  else if (s == "float")       return number_type::FLOAT;
  else if (s == "double")      return number_type::DOUBLE;
  else if (s == "long double") return number_type::LONG_DOUBLE;
  else                         return number_type::UNKNOWN;
}

/// Typealias is visited on early stage (construction of AST?), and it seems the
/// visitors do not revisit it.  Thus, make a vector of definitions, and check matches

struct vector_type_struct {
  std::string classname;
  number_type ntype;
};

/// this vector keeps track of types which are vectorizable.

static std::vector<vector_type_struct> vectorizable_types = {};

// reset the list, and add the "basic" types

void reset_vectorizable_types() {
  vectorizable_types.clear();
  vector_type_struct d;
  
  d.classname = "double";
  d.ntype = number_type::DOUBLE;
  vectorizable_types.push_back(d);

  d.classname = "float";
  d.ntype = number_type::FLOAT;
  vectorizable_types.push_back(d);

  d.classname = "int";
  d.ntype = number_type::INT;
  vectorizable_types.push_back(d);
}

//////////////////////////////////////////////////////////////////////////
/// This visitor tallies up typealiases
/// Used only here to determine if the class where it appears is vectorizable
//////////////////////////////////////////////////////////////////////////

bool MyASTVisitor::VisitTypeAliasDecl(TypeAliasDecl * ta) {

  if (ta->getNameAsString() == "base_type") {

    std::string class_name = ta->getQualifiedNameAsString();
    std::string type_name = ta->getUnderlyingType().getCanonicalType().getAsString();

    number_type nt = get_number_type(type_name);
    if (nt == number_type::UNKNOWN) return true;   // Unknown number type, we'll do nothing
  
    // remove ::base_type  from the name
    int r = class_name.rfind("::base_type");
    assert( r != std::string::npos );
    class_name.resize(r);

    for (auto n : vectorizable_types) {
      if (class_name.compare(n.classname) == 0) {
        // Match already found!
        if (nt != n.ntype) {
          llvm::errs() << "Internal type alias classification error, bug\n";
          exit(-1);
        }

        // need not to do anything more
        return true;
      }
    }
    
    vector_type_struct n;
    n.classname = class_name;
    n.ntype = nt;
    vectorizable_types.push_back(n);

  }
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// is_vectorizable_type  true for types which can be vectorized, and returns the
/// corresponding vectorized typename in vectorized_type
/////////////////////////////////////////////////////////////////////////////////////////


bool MyASTVisitor::is_vectorizable_type(const QualType & QT, vectorization_info & vi) {
  PrintingPolicy pp(Context->getLangOpts());
  std::string tname = QT.getCanonicalType().getAsString(pp);
  return is_vectorizable_type(tname, vi);
}


bool MyASTVisitor::is_vectorizable_type(const std::string & type_name, vectorization_info & vi) {

  // Note: the std types are included in vectorizable_types
  // vectorizable_type vector filled in by VisitTypeAliasDecl above

  for (auto n : vectorizable_types) if (type_name.compare(n.classname) == 0) {
    // found it
    vi.basetype = n.ntype;

    if (target.vectorize) {
      std::string old_t, new_t;

      if (n.ntype == number_type::DOUBLE) {
        vi.vector_size = target.vector_size/sizeof(double);
        old_t = vi.basetype_str = "double";
        new_t = vi.vectortype = "Vec" + std::to_string(vi.vector_size) + "d";
        
      } else if (n.ntype == number_type::FLOAT) {
        vi.vector_size = target.vector_size/sizeof(float);
        old_t = vi.basetype_str = "float";
        new_t = vi.vectortype = "Vec" + std::to_string(vi.vector_size) + "f";

      } else if (n.ntype == number_type::INT) {
        vi.vector_size = target.vector_size/sizeof(int);
        old_t = vi.basetype_str = "int";
        new_t = vi.vectortype = "Vec" + std::to_string(vi.vector_size) + "i";

      } else if (n.ntype == number_type::INT64_T) {
        vi.vector_size = target.vector_size/sizeof(int64_t);
        old_t = vi.basetype_str = "int64_t";
        new_t = vi.vectortype = "Vec" + std::to_string(vi.vector_size) + "q";

      } else return vi.is_vectorizable = false;   // not vectorizable type

      // replace the old type in the type name
      vi.vectorized_type = type_name;
      int i = find_word(vi.vectorized_type, old_t);
      if (i == std::string::npos) {
        // Now the type template does not contain the to-be-vectorized type at all.
        // This cannot be vectorized
        return vi.is_vectorizable = false;
      }
      while (i != std::string::npos) {
        vi.vectorized_type.replace(i,old_t.size(),new_t);
        i = find_word(vi.vectorized_type, old_t, i);
      }
      return vi.is_vectorizable = true;

    } else {
      // now target.vectorize is false, we don't know the vec length
      if (n.ntype == number_type::DOUBLE || n.ntype == number_type::FLOAT || 
          n.ntype == number_type::INT || n.ntype == number_type::INT64_T) {
            vi.vectorized_type = "";   // in principle vectorizable, don't know the type
            vi.vector_size = 0;
            return vi.is_vectorizable = true;
      } else return vi.is_vectorizable = false;          // not known vectorizable type
    }
  }  // for loop
  return vi.is_vectorizable = false; // not found
}

////////////////////////////////////////////////////////////////////////////7/////
/// Is the Field<> vectorizable?
/// if so, return the information on field_element_info
//////////////////////////////////////////////////////////////////////////////////

vectorization_info MyASTVisitor::inspect_field_type(Expr *fE) {

  // Start from Expr for field name, need to descend rather deep
  // Assertions here should actually not be necessary, there just to ensure things happen as they should

  assert( is_field_expr(fE) );

  const Type * tp = fE->getType().getTypePtr();
  assert(tp != nullptr);
  
  const CXXRecordDecl * rd = tp->getAsCXXRecordDecl();
  assert(rd != nullptr);

  const ClassTemplateSpecializationDecl * fspec = dyn_cast<ClassTemplateSpecializationDecl>(rd);
  assert(fspec != nullptr);

  const TemplateArgumentList & tal = fspec->getTemplateArgs();
  // Field should have only 1 template argument (for now), and should be of type
  assert(tal.size() == 1);
  
  vectorization_info einfo;

  QualType QT = tal.get(0).getAsType();

  if (is_vectorizable_type(QT,einfo)) {
    llvm::errs() << "Type " << QT.getCanonicalType().getAsString() << " is vectorizable -> " << einfo.vectorized_type << '\n';
  }
  return einfo;
}



