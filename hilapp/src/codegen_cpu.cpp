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

std::string TopLevelVisitor::generate_code_cpu(Stmt *S, bool semicolon_at_end,
                                               srcBuf &loopBuf,
                                               bool generate_wait_loops) {
    std::stringstream code;

    // Set loop lattice
    if (field_info_list.size() > 0) {
        std::string fieldname = field_info_list.front().new_name;
        code << "const lattice_struct * RESTRICT loop_lattice = " << fieldname
             << ".fs->lattice;\n";
    } else {
        // if there is no field in loop at all
        code << "const lattice_struct * RESTRICT loop_lattice = lattice;\n";
    }

    // Set the start and end points
    code << "const int loop_begin = loop_lattice->loop_begin(" << parity_in_this_loop
         << ");\n";
    code << "const int loop_end   = loop_lattice->loop_end(" << parity_in_this_loop
         << ");\n";

    // are there

    if (generate_wait_loops) {
        code << "for (int _wait_i_ = 0; _wait_i_ < 2; ++_wait_i_) {\n";
    }

    // and the openacc loop header
    if (target.openacc)
        generate_openacc_loop_header(code);

    // Start the loop
    code << "for(int " << looping_var << " = loop_begin; " << looping_var
         << " < loop_end; ++" << looping_var << ") {\n";

    if (generate_wait_loops) {
        code << "if (((loop_lattice->wait_arr_[" << looping_var
             << "] & _dir_mask_) != 0) == _wait_i_) {\n";
    }

    // replace reduction variables in the loop
    for (var_info &vi : var_info_list) {
        if (!vi.is_loop_local) {
            if (vi.reduction_type != reduction::NONE) {
                for (var_ref &vr : vi.refs) {
                    loopBuf.replace(vr.ref, vi.reduction_name);
                }
            }
        }
    }

    // Create temporary field element variables
    for (field_info &l : field_info_list) {

        // First check for Direction references. If any found, create list of temp
        // variables
        if (l.is_read_nb) {
            for (dir_ptr &d : l.dir_list) {
                std::string dirname;
                if (d.is_constant_direction)
                    dirname = d.direxpr_s; // orig. string
                else
                    dirname = remove_X(loopBuf.get(
                        d.parityExpr
                            ->getSourceRange())); // mapped name was get_stmt_str(d.e);

                // generate access stmt
                code << l.element_type << " " << d.name_with_dir << " = " << l.new_name;

                if (target.vectorize && l.vecinfo.is_vectorizable) {
                    // now l is vectorizable, but accessed sequentially -- this inly
                    // happens in vectorized targets
                    code << ".get_value_at_nb_site(" << dirname << ", " << looping_var
                         << ");\n";
                } else {
                    // std neighbour accessor for scalars
                    code << ".get_value_at(" << l.new_name << ".fs->neighbours["
                         << dirname << "][" << looping_var << "]);\n";
                }

                // and replace references in loop body
                for (field_ref *ref : d.ref_list) {
                    loopBuf.replace(ref->fullExpr, d.name_with_dir);
                }
            }
        }

        // and then get (possible) local refs
        if (l.is_read_atX) {
            // now reading var without nb. reference
            code << l.element_type << " " << l.loop_ref_name << " = " << l.new_name
                 << ".get_value_at(" << looping_var << ");\n";

        } else if (l.is_written) {
            code << l.element_type << " " << l.loop_ref_name << ";\n";
        }

        // and finally replace references in body
        for (field_ref *ref : l.ref_list)
            if (!ref->is_direction) {
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
    for (special_function_call &sfc : special_function_call_list) {
        std::string repl = sfc.replace_expression; // comes with ( now
        if (sfc.add_loop_var) {
            repl += looping_var;
            if (sfc.argsExpr != nullptr)
                repl += ',';
        }
        loopBuf.replace(sfc.replace_range, repl);
    }

    // Dump the main loop code here
    code << loopBuf.dump();
    if (semicolon_at_end)
        code << ";";
    code << "\n";

    // Add calls to setters
    for (field_info &l : field_info_list)
        if (l.is_written) {
            code << l.new_name << ".set_value_at(" << l.loop_ref_name << ", "
                 << looping_var << ");\n";
        }

    code << "}\n";

    if (generate_wait_loops) {
        // add the code for 2nd round - also need one } to balance the if ()
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

    return code.str();
}
