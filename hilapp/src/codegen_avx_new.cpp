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

/// An AST walker for finding and handling variable declarations
/// in a loop function
class LoopFunctionHandler : public GeneralVisitor,
                            public RecursiveASTVisitor<LoopFunctionHandler> {
  public:
    using GeneralVisitor::GeneralVisitor;

    // Buffer for the function copy
    srcBuf functionBuffer;
    int vector_size;

    template <typename visitortype>
    LoopFunctionHandler(visitortype &v) : GeneralVisitor(v) {}

    bool TraverseStmt(Stmt *s) {
        RecursiveASTVisitor<LoopFunctionHandler>::TraverseStmt(s);
        return true;
    }
    bool VisitVarDecl(VarDecl *var);
    bool VisitCXXOperatorCallExpr(CXXOperatorCallExpr *op);
    bool VisitBinaryOperator(BinaryOperator *op);
};

bool LoopFunctionHandler::VisitVarDecl(VarDecl *var) {
    std::string typestring = var->getType().getAsString();

    size_t begin;
    begin = typestring.rfind("element", 0);
    if (begin != std::string::npos) {
        // This variable is an element, replace with vector
        std::string vector_type = typestring;
        if (begin != std::string::npos) {
            vector_type.replace(begin, 7, "vectorize_struct");
            vector_type.replace(vector_type.find_last_of(">"), 1,
                                ", " + std::to_string(vector_size) + ">::type");
        }

        if (var->isDirectInit()) {
            std::string init =
                TheRewriter.getRewrittenText(var->getInit()->getSourceRange());
            functionBuffer.replace(var->getSourceRange(), vector_type + " " +
                                                              var->getNameAsString() +
                                                              "=" + init);
        } else {
            functionBuffer.replace(var->getSourceRange(),
                                   vector_type + " " + var->getNameAsString());
        }

    } else {
        if (var->hasInit()) {
            LoopAssignChecker lac(*this);
            lac.TraverseStmt(var->getInit());
        }
    }
    return true;
}

bool LoopFunctionHandler::VisitCXXOperatorCallExpr(CXXOperatorCallExpr *op) {
    if (op && op->isAssignmentOp()) {
        std::string type = op->getArg(0)->getType().getAsString();
        type = remove_all_whitespace(type);
        if (type.rfind("element<", 0) == std::string::npos) {
            LoopAssignChecker lac(*this);
            lac.TraverseStmt(op);
        }
    }
    return true;
}

bool LoopFunctionHandler::VisitBinaryOperator(BinaryOperator *op) {
    if (op && op->isAssignmentOp()) {
        std::string type = op->getLHS()->getType().getAsString();
        type = remove_all_whitespace(type);
        if (type.rfind("element<", 0) == std::string::npos) {
            LoopAssignChecker lac(*this);
            lac.TraverseStmt(op);
        }
    }
    return true;
}

/// Replace element types with vector. Leaves other types untouched.
static void replace_element_with_vector(SourceRange sr, std::string typestring,
                                        std::string namestring, int vector_size,
                                        srcBuf &functionBuffer) {
    if (typestring.rfind("element", 0) != std::string::npos) {
        std::string vector_type = typestring;
        size_t begin;
        begin = find_word(vector_type, "element");
        if (begin != std::string::npos) {
            vector_type.replace(begin, 7, "vectorize_struct");
            vector_type.replace(vector_type.find_last_of(">"), 1,
                                ", " + std::to_string(vector_size) + ">::type");
        }

        functionBuffer.replace(sr, vector_type + " " + namestring);
    }
}

/// This should allow calling with element<>-type parameters from
/// loops. Generate a copy with elements replaced with vectors.
void GeneralVisitor::handle_loop_function_avx(call_info_struct &ci) {

    FunctionDecl *fd = ci.funcdecl;

    SourceRange sr = fd->getSourceRange();
    srcBuf *sourceBuf = get_file_srcBuf(sr.getBegin());

    if (sourceBuf == nullptr) {
        // it's a system function - should probably do something?
        return;
    }

    // Track wether the function actually contains elements.
    // if not, no new function should be written
    bool generate_function = false;

    // Check allowed vector sizes
    int smallest = 1, largest = 0;
    for (clang::ParmVarDecl *par : fd->parameters()) {
        std::string typestring = par->getType().getAsString(PP);
        if (find_word(typestring, "double") != std::string::npos) {
            smallest = 4;
            largest = 8;
        }
        if (find_word(typestring, "float") != std::string::npos ||
            find_word(typestring, "int") != std::string::npos ||
            find_word(typestring, "CoordinateVector") != std::string::npos) {
            smallest = 8;
            largest = 16;
        }

        // Check if there are elements in the first place
        if (find_word(typestring, "element") != std::string::npos)
            generate_function = true;
    }

    if (generate_function)
        for (int vector_size = smallest; vector_size <= largest; vector_size *= 2) {

            LoopFunctionHandler lfh(*this);
            lfh.functionBuffer.copy_from_range(sourceBuf, sr);
            lfh.vector_size = vector_size;

            // Handle each parameter
            for (clang::ParmVarDecl *par : fd->parameters()) {
                std::string typestring = par->getType().getAsString(PP);
                replace_element_with_vector(par->getSourceRange(), typestring,
                                            par->getNameAsString(), vector_size,
                                            lfh.functionBuffer);
            }

            // Handle return type
            // Note: C++ cannot specialize only based on return type. Therefore we
            // only write a new function if the parameters contain elements
            std::string typestring = fd->getReturnType().getAsString(PP);
            replace_element_with_vector(fd->getReturnTypeSourceRange(), typestring, "",
                                        vector_size, lfh.functionBuffer);

            lfh.TraverseStmt(fd->getBody());

            std::string buffer = lfh.functionBuffer.dump();
            if (!(fd->hasBody())) {
                // Declaration does not contain a body, needs a semicolon
                buffer += ";";
            }
            buffer += "\n";
            sourceBuf->insert(sr.getBegin(), buffer, true, true);
        }
}

/// Constructors - should something be done here?

void GeneralVisitor::handle_loop_constructor_avx(call_info_struct &ci) {}

///////////////////////////////////////////////////////////////////////////////////
/// Check that
///  a) no site dependent conditional
///  b) fields have the same vector size (number size)
///  c) variables are vectorizable to the same size
///  NOW: instead of b) and c) now vectorizeed only if everything has the same
///       vector type.  Otherwise leads to missing type conversions
///  TODO: Rectify this issue!
///////////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::check_loop_vectorizable(Stmt *S, int &vector_size,
                                              std::string &diag_str) {

    vector_size = 0;
    number_type numtype;
    bool is_vectorizable = true;

    std::vector<std::string> reason = {};

    // check if loop has conditional
    if (loop_info.has_pragma_novector) {
        is_vectorizable = false;
        reason.push_back("it has '#pragma hila novector'");

    } else {

        if (loop_info.has_site_dependent_conditional) {
            is_vectorizable = false;
            reason.push_back("it contains site dependent conditional");
        }

        if (contains_random(S)) {
            is_vectorizable = false;
            reason.push_back("it contains a random number generator");
        }

        std::string vector_var_name; // variable which determines the vectorization
        std::string
            vector_var_type; // type of variable which determines the vectorization

        // check if the fields are vectorizable in a compatible manner
        if (field_info_list.size() > 0) {
            for (field_info &fi : field_info_list) {
                if (!fi.vecinfo.is_vectorizable) {
                    is_vectorizable = false;
                    reason.push_back("Field variable '" + fi.old_name +
                                     "' is not vectorizable");
                } else {
                    if (vector_size == 0) {
                        vector_size = fi.vecinfo.vector_size;
                        numtype = fi.vecinfo.basetype;
                        vector_var_name = fi.old_name;
                        vector_var_type = fi.vecinfo.basetype_str;
                    } else if (fi.vecinfo.vector_size != vector_size ||
                               fi.vecinfo.basetype != numtype) {
                        // TODO: let different vector types coexist!

                        is_vectorizable = false;

                        reason.push_back("type of variable '" + fi.old_name + "' is " +
                                         fi.vecinfo.basetype_str + " and '" +
                                         vector_var_name + "' is " + vector_var_type);
                    }
                }
            }
        }

        // and then if the site dep. variables are vectorizable
        if (var_info_list.size() > 0) {
            for (var_info &vi : var_info_list)
                if (vi.is_site_dependent) {
                    if (vi.vecinfo.is_vectorizable) {
                        if (vector_size == 0) {
                            vector_size = vi.vecinfo.vector_size;
                            numtype = vi.vecinfo.basetype;
                            vector_var_name = vi.name;
                            vector_var_type = vi.vecinfo.basetype_str;

                        } else if (vector_size != vi.vecinfo.vector_size ||
                                   numtype != vi.vecinfo.basetype) {
                            is_vectorizable = false;

                            reason.push_back("type of variables '" + vi.name + "' is " +
                                             vi.vecinfo.basetype_str + " and '" +
                                             vector_var_name + "' is " + vector_var_type);
                        }
                    } else {
                        is_vectorizable = false;

                        reason.push_back("variable '" + vi.name +
                                         "' is not vectorizable");
                    }
                }
        }

        // and still, check the special functions
        if (is_vectorizable) {
            for (auto const &sfc : special_function_call_list) {
                if (sfc.name == "coordinates" || sfc.name == "coordinate") {
                    // returning int vector
                    if (vector_size == 0) {
                        vector_size = target.vector_size / sizeof(int);
                    } else if (vector_size != target.vector_size / sizeof(int) ||
                               numtype != number_type::INT) {
                        is_vectorizable = false;

                        reason.push_back("functions 'X.coordinates()' and "
                                         "'X.coordinate(direction)' return int, "
                                         "which is not compatible with " +
                                         vector_var_type + " vectors");
                    }

                } else if (sfc.name == "parity") {
                    is_vectorizable = false;

                    reason.push_back("function 'X.parity()' is not AVX vectorizable");

                } else if (sfc.name == "random" || sfc.name == "hila::random") {
                    is_vectorizable = false;
                    reason.push_back("random number generators prevent vectorization");
                }
            }
        }

        // and check function calls
        for (auto &ci : loop_function_calls) {
            if (ci.is_site_dependent && !ci.is_vectorizable) {
                is_vectorizable = false;

                if (ci.funcdecl != nullptr) {
                    reason.push_back("loop contains function " +
                                     ci.funcdecl->getNameAsString() +
                                     " which is not vectorizable");
                } else if (ci.ctordecl != nullptr) {
                    reason.push_back("loop contains constructor " +
                                     ci.ctordecl->getNameAsString() +
                                     " which is not vectorizable");
                }
            }
        }
    }

    if (vector_size == 0 && is_vectorizable) {
        // super-special case - loop does not contain fields or anything site dependent.
        // Loop is probably useless.  Let us just not vectorize it.
        is_vectorizable = false;
        reason.push_back("loop does not seem to have site dependent content");
    }

    if (!is_vectorizable) {
        diag_str = "Loop is not AVX vectorizable because:";
        for (auto &s : reason)
            diag_str += "\n     " + s;

        if (cmdline::avx_info > 0 || cmdline::verbosity > 0)
            reportDiag(DiagnosticsEngine::Level::Remark, S->getSourceRange().getBegin(),
                       "%0", diag_str.c_str());

    } else {
        diag_str = "Loop is AVX vectorizable";

        if (cmdline::avx_info > 1 || cmdline::verbosity > 1)
            reportDiag(DiagnosticsEngine::Level::Remark, S->getSourceRange().getBegin(),
                       "%0", diag_str.c_str());
    }

    return is_vectorizable;
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Main entry for AVX loop generation
///////////////////////////////////////////////////////////////////////////////////////////

std::string TopLevelVisitor::generate_code_avx(Stmt *S, bool semicolon_at_end,
                                               srcBuf &loopBuf,
                                               bool generate_wait_loops) {

    std::stringstream code;
    int vector_size;

    // is the loop vectorizable?
    std::string vector_diag;
    bool is_vectorized = check_loop_vectorizable(S, vector_size, vector_diag);

    code << comment_string(vector_diag) << '\n';

    // can go through std non-vector code generation
    if (!is_vectorized) {
        code << generate_code_cpu(S, semicolon_at_end, loopBuf, generate_wait_loops);
        return code.str();
    }

    // Create temporary variables for reductions (vector reduction is in the loop)
    for (var_info &v : var_info_list) {
        if (v.reduction_type != reduction::NONE) {
            v.new_name = "v_" + v.reduction_name;
            // Allocate memory for a reduction. This will be filled in the kernel
            code << v.vecinfo.vectorized_type << ' ' << v.new_name;
            if (v.reduction_type == reduction::SUM)
                code << "(0);\n";
            else if (v.reduction_type == reduction::PRODUCT)
                code << "(1);\n";
        }
    }

    // Set loop lattice for neibhbour arrays
    if (field_info_list.size() > 0) {
        std::string fieldname = field_info_list.front().new_name;
        code << "const auto * RESTRICT loop_lattice = " << fieldname
             << ".fs->vector_lattice;\n";
    } else {
        // if no field in loop
        code << "const auto * RESTRICT loop_lattice = "
                "lattice->backend_lattice->get_vectorized_lattice<"
             << vector_size << ">();\n";
    }

    // Set the start and end points
    code << "const int loop_begin = loop_lattice->loop_begin(" << parity_in_this_loop
         << ");\n";
    code << "const int loop_end   = loop_lattice->loop_end(" << parity_in_this_loop
         << ");\n";

    if (generate_wait_loops) {
        code << "for (int _wait_i_ = 0; _wait_i_ < 2; ++_wait_i_) {\n";
    }

    // Start the loop
    code << "for(int " << looping_var << " = loop_begin; " << looping_var
         << " < loop_end; ++" << looping_var << ") {\n";

    if (generate_wait_loops) {
        code << "if (((loop_lattice->vec_wait_arr_[" << looping_var
             << "] & _dir_mask_) != 0) == _wait_i_) {\n";
    }

    // Add vector reduction variable here, inside the loop
    for (vector_reduction_ref &vrf : vector_reduction_ref_list) {
        // Allocate memory for a reduction and initialize
        if (vrf.reduction_type == reduction::SUM) {
            code << "vectorize_struct<" << vrf.type << ", vector_size>::type v_"
                 << vrf.vector_name << "(0);\n";
        }
        if (vrf.reduction_type == reduction::PRODUCT) {
            code << "vectorize_struct<" << vrf.type << ", vector_size>::type v_"
                 << vrf.vector_name << "(1);\n";
        }
        loopBuf.replace(vrf.ref, "v_" + vrf.vector_name);
    }

    // Create temporary field element variables
    for (field_info &l : field_info_list) {

        // First check for direction references. If any found, create list of temp
        // variables
        if (l.is_read_nb) {
            for (dir_ptr &d : l.dir_list)
                if (d.count > 0) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(loopBuf.get(
                            d.parityExpr->getSourceRange())); // mapped name was
                                                              // get_stmt_str(d.e);

                    code << l.vecinfo.vectorized_type << " " << d.name_with_dir << " = "
                         << l.new_name << ".get_vector_at<" << l.vecinfo.vectorized_type
                         << ">(loop_lattice->neighbours[" << dirname << "]["
                         << looping_var << "]);\n";

                    // and replace references in loop body
                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, d.name_with_dir);
                    }
                }
        }

        if (l.is_read_atX) {
            code << l.vecinfo.vectorized_type << " " << l.loop_ref_name << " = "
                 << l.new_name << ".get_vector_at<" << l.vecinfo.vectorized_type << ">("
                 << looping_var << ");\n";
        } else if (l.is_written) {
            code << l.vecinfo.vectorized_type << " " << l.loop_ref_name << ";\n";
        }

        // and finally replace these references in body
        for (field_ref *ref : l.ref_list)
            if (!ref->is_direction) {
                loopBuf.replace(ref->fullExpr, l.loop_ref_name);
            }
    }

    // other variable refs
    for (var_info &vi : var_info_list) {
        // reduction variables
        if (vi.reduction_type != reduction::NONE) {
            // Replace references in the loop body
            for (var_ref &vr : vi.refs) {
                loopBuf.replace(vr.ref, vi.new_name);
            }
        } else if (vi.is_site_dependent) {
            // now must be loop-local vectorized var
            // change declaration - name need not be changed
            // loopBuf.replace( vi.decl->getSourceRange(), vi.vecinfo.vectorized_type );
            loopBuf.replace(vi.decl->getTypeSourceInfo()->getTypeLoc().getSourceRange(),
                            vi.vecinfo.vectorized_type);
        }
    }

    // Handle calls to special in-loop functions
    for (special_function_call &sfc : special_function_call_list) {
        std::string repl = sfc.replace_expression; // comes with ( now
        if (sfc.add_loop_var) {
            repl += looping_var;
            if (sfc.argsExpr != nullptr)
                repl += ',';
        }
        loopBuf.replace(sfc.replace_range, repl);
    }

    // Vector reductions must be in the sames scope as the loop body. Otherwise the index
    // may be undefined. Therefore add it before the closing }
    if (!semicolon_at_end) {
        // Remove the last }
        loopBuf.remove(loopBuf.size() - 2, loopBuf.size() - 1);
    }

    // Dump the main loop code here
    code << loopBuf.dump();
    if (semicolon_at_end)
        code << ";";
    code << "\n";

    // Add vector reductions
    int i = 0;
    for (vector_reduction_ref &vrf : vector_reduction_ref_list) {
        // run reduction over the vector
        code << "int v_index_" << i << "[vector_size];\n";
        code << vrf.index_name << ".store(&v_index_" << i << "[0]);\n";
        code << vrf.type << " a_" << vrf.vector_name << "[vector_size];\n";
        code << "v_" << vrf.vector_name << ".store(&a_" << vrf.vector_name << "[0]);\n";
        code << "for( int i=0; i<vector_size; i++){\n";
        if (vrf.reduction_type == reduction::SUM) {
            code << vrf.vector_name << "[v_index_" << i << "[i]] += "
                 << "a_" << vrf.vector_name << "[i];\n";
        }
        if (vrf.reduction_type == reduction::PRODUCT) {
            code << vrf.vector_name << "[v_index_" << i << "[i]] *= "
                 << "a_" << vrf.vector_name << "[i];\n";
        }
        code << "}\n";
    }

    if (!semicolon_at_end) {
        code << "}";
    }
    code << "\n";

    // Add calls to setters
    for (field_info &l : field_info_list) {
        if (l.is_written) {
            code << l.new_name << ".set_vector_at<" << l.vecinfo.vectorized_type << ">("
                 << l.loop_ref_name << ", " << looping_var << ");\n";
        }
    }

    code << "}\n";

    if (generate_wait_loops) {
        // add the code for 2nd round
        code << "}\nif (_dir_mask_ == 0) break;    // No need for another round\n";

        for (field_info &l : field_info_list) {
            // If neighbour references exist, communicate them
            for (dir_ptr &d : l.dir_list)
                if (d.count > 0) {
                    code << l.new_name << ".wait_fetch(" << d.direxpr_s << ", "
                         << parity_in_this_loop << ");\n";
                }
        }
        code << "}\n";
    }

    // Final reduction of the temporary reduction variables
    for (var_info &v : var_info_list) {
        if (v.reduction_type == reduction::SUM) {
            // code << v.reduction_name << " = reduce_sum(" << v.new_name << ");\n";

            code << v.reduction_name << " = reduce_sum_in_vector<"
                 << v.vecinfo.basetype_str << ", " << v.vecinfo.vectortype << ", "
                 << v.type << ", " << v.vecinfo.vectorized_type << ">(" << v.new_name
                 << ");\n";

        } else if (v.reduction_type == reduction::PRODUCT) {
            code << v.reduction_name << " = reduce_prod(" << v.new_name << ");\n";
        }
    }

    return code.str();
}
