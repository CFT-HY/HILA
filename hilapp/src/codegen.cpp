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

std::string looping_var;
std::string parity_name;

std::string parity_in_this_loop = "";

/// Used in lattice loop generation
std::string parity_str(Parity p) {
    switch (p) {
    case Parity::none:
        return "Parity::none";
    case Parity::even:
        return "Parity::even";
    case Parity::odd:
        return "Parity::odd";
    case Parity::all:
        return "Parity::all";
    default:
        return "unknown";
    }
}

/// Generate a name that does not appear in string t
inline std::string unique_name(const std::string t, std::string n) {
    while (t.find(n, 0) != std::string::npos)
        n += "_";
    return n;
}

/// The main entry point for code generation
void TopLevelVisitor::generate_code(Stmt *S) {
    srcBuf loopBuf; // (&TheRewriter,S);
    loopBuf.copy_from_range(writeBuf, S->getSourceRange());

    //   llvm::errs() << "\nOriginal range: +++++++++++++++\n\""
    //                << TheRewriter.getRewrittenText(S->getSourceRange())
    //                << "\"\nwriteBuf range: ================\n\""
    //                << writeBuf->get(S->getSourceRange())
    //                << "\"\nCopied range: ================\n\""
    //                << loopBuf.dump() << "\"\n";

    // is it compound stmt: { } -no ; needed
    bool semicolon_at_end = !(isa<CompoundStmt>(S));

    // Build replacement in variable "code"
    // Encapsulate everything within {}
    std::stringstream code;
    code << "{\n";

    // basic set up: 1st loop_info, if it is known const set it up,
    // else copy it to a variable name

    const std::string t = loopBuf.dump();

    looping_var = "Index";
    while (t.find(looping_var, 0) != std::string::npos)
        looping_var += "_";

    // Ensure that the name is not reserved by scanning the source
    parity_name = "Parity";
    while (t.find(parity_name, 0) != std::string::npos)
        parity_name += "_";

    if (loop_info.parity_value == Parity::none) {
        // now unknown
        code << "const Parity " << parity_name << " = " << loop_info.parity_text
             << ";\n";

        if (global.assert_loop_parity) {
            code << "assert( is_even_odd_parity(" << parity_name
                 << ") && \"Parity should be EVEN or ODD in this loop\");\n";
        }
        parity_in_this_loop = parity_name;

    } else
        parity_in_this_loop = parity_str(loop_info.parity_value);

    // then, generate new names for field variables in loop

    for (field_info &l : field_info_list) {
        // Generate new variable name, may be needed -- use here simple receipe
        l.new_name = "F" + clean_name(l.old_name);
        // Perhaps simpler FA, FB, FC. ?  The following helps avoid collisions
        while (t.find(l.new_name, 0) != std::string::npos)
            l.new_name += "_";
        l.loop_ref_name = l.new_name + "_index";

        // Create neighbour ref names
        int i = 0;
        for (dir_ptr &d : l.dir_list) {
            d.name_with_dir = l.new_name + "_dir" + std::to_string(++i);
        }

        // make a ref to the field name
        if (!l.is_written)
            code << "const ";
        code << "Field" << l.type_template << " & " << l.new_name << " = " << l.old_name
             << ";\n";
    }

    // check alloc and do it if needed
    for (field_info &l : field_info_list)
        if (l.is_written) {
            code << l.new_name << ".check_alloc();\n";
        }

    // Check that read fields are initialized
    for (field_info &l : field_info_list) {
        if (l.is_read_nb || l.is_read_atX) {
            std::string init_par;
            if (loop_info.parity_value == Parity::all ||
                (l.is_read_nb && l.is_read_atX)) {
                init_par = "ALL";
            } else {
                if (l.is_read_atX)
                    init_par = parity_in_this_loop;
                else
                    init_par = "opp_parity(" + parity_in_this_loop + ")";
            }

            // TAKE THIS AWAY FOR NOW -- WE DON'T HAVE A GOOD METHOD TO CHECK "pure
            // output" FUNCTIONS E.g. method   U[X].random(), where U does not have to
            // be initialized
            // actually, we do now -- the "output_only" keyword

            if (cmdline::check_initialization) {

                std::string fname = srcMgr.getFilename(S->getSourceRange().getBegin());
                code << "if (!" << l.new_name << ".is_initialized(" << init_par
                     << ")){\noutput0 << \"File " << fname << " on line "
                     << srcMgr.getSpellingLineNumber(S->getSourceRange().getBegin())
                     << ":\\n Value of variable " << l.old_name
                     << " is used but it is not properly initialized\\n\";\n";
                code << "hila::terminate(1);\n}\n";
            }
        }
    }

    // change the f[X+offset] -references, generate code
    handle_field_plus_offsets(code, loopBuf, parity_in_this_loop);

    bool first = true;
    bool generate_wait_loops;
    if (cmdline::no_interleaved_comm || cmdline::no_mpi)
        generate_wait_loops = false;
    else
        generate_wait_loops = true;

    for (field_info &l : field_info_list) {
        // If neighbour references exist, communicate them
        for (dir_ptr &d : l.dir_list)
            if (d.count > 0) {
                if (!generate_wait_loops) {
                    code << l.new_name << ".fetch(" << d.direxpr_s << ", "
                         << parity_in_this_loop << ");\n";
                } else {
                    if (first)
                        code << "dir_mask_t  _dir_mask_ = 0;\n";
                    first = false;

                    code << "_dir_mask_ |= " << l.new_name << ".start_fetch("
                         << d.direxpr_s << ", " << parity_in_this_loop << ");\n";
                }
            }
    }

    if (first)
        generate_wait_loops = false; // no communication needed in the 1st place

    // Create temporary variables for reductions
    for (var_info &v : var_info_list) {
        if (v.reduction_type != reduction::NONE) {
            v.reduction_name = "r_" + v.name;
            while (t.find(v.reduction_name, 0) != std::string::npos)
                v.reduction_name += "_";
            // Create a temporary variable and initialize
            if (v.reduction_type == reduction::SUM) {
                code << v.type << " " << v.reduction_name << " = 0;\n";
            } else if (v.reduction_type == reduction::PRODUCT) {
                code << v.type << " " << v.reduction_name << " = 1;\n";
            }
        }
    }

    // Place the content of the loop
    code << backend_generate_code(S, semicolon_at_end, loopBuf, generate_wait_loops);

    // Check reduction variables
    for (var_info &v : var_info_list) {
        if (v.reduction_type == reduction::SUM) {
            code << "lattice->reduce_node_sum( &" << v.reduction_name
                 << ", 1, true);\n";
            code << v.name << " += " << v.reduction_name << ";\n";
        } else if (v.reduction_type == reduction::PRODUCT) {
            code << "lattice->reduce_node_product( &" << v.reduction_name
                 << ", 1, true);\n";
            code << v.name << " *= " << v.reduction_name << ";\n";
        }
    }
    for (vector_reduction_ref &vrf : vector_reduction_ref_list) {
        if (vrf.reduction_type == reduction::SUM) {
            code << "lattice->reduce_node_sum(" << vrf.vector_name << ".data(), "
                 << vrf.vector_name << ".size(), true);\n";
        }
        if (vrf.reduction_type == reduction::PRODUCT) {
            code << "lattice->reduce_node_product(" << vrf.vector_name << ".data(), "
                 << vrf.vector_name << ".size(), true);\n";
        }
    }

    // finally mark modified fields
    for (field_info &l : field_info_list)
        if (l.is_written) {
            code << l.new_name << ".mark_changed(" << parity_in_this_loop << ");\n";
        }

    // and close
    code << "}\n//----------";

    // Remove old code + replace
    if (semicolon_at_end) {
        writeBuf->remove(getRangeWithSemicolon(S));
    } else {
        writeBuf->remove(S->getSourceRange());
    }

    writeBuf->insert(S->getBeginLoc(), indent_string(code.str()), true, true);
}

// Handle field+offset expressions -- these are copied, ref removed
// Use iterator for looping, because we may want to remove the field_info item.

void TopLevelVisitor::handle_field_plus_offsets(std::stringstream &code,
                                                srcBuf &loopBuf,
                                                std::string &paritystr) {

    for (auto it = field_info_list.begin(); it != field_info_list.end();) {

        if (!it->is_read_offset) {
            // std case, no offset
            ++it;

        } else {

            int i_offset = 0;
            for (dir_ptr &d : it->dir_list)
                if (d.is_offset) {

                    // get new name
                    std::string offset_field_name = it->new_name + "_offset";
                    if (i_offset > 0)
                        offset_field_name += i_offset;
                    i_offset++;

                    // make new field info for new variable
                    field_info new_fi;
                    new_fi.type_template = it->type_template;
                    new_fi.element_type = it->element_type;
                    new_fi.old_name = new_fi.new_name = offset_field_name;
                    new_fi.loop_ref_name = offset_field_name + "_index";
                    new_fi.ref_list = d.ref_list;
                    new_fi.dir_list = {};
                    new_fi.is_read_nb = false;
                    new_fi.is_read_atX = true;

                    // push it on stack
                    field_info_list.push_back(new_fi);

                    // copy the shifted var
                    code << "const Field" + it->type_template + " " +
                                offset_field_name + " = " + it->new_name + ".shift(" +
                                d.ref_list.at(0)->direxpr_s + ", " + paritystr + ");\n";

                    // and rewrite references to the offset field
                    for (field_ref *fr : new_fi.ref_list) {
                        loopBuf.replace(fr->nameExpr, offset_field_name);
                        loopBuf.replace(fr->parityExpr, "X");
                        fr->direxpr_s.clear();
                        fr->is_direction = false; // no Direction
                        fr->info =
                            &field_info_list.back(); // info points to new field_info
                        fr->is_written = false;
                        fr->is_read = true;
                        fr->is_offset =
                            true; // leave this on, as a flag -- turned off below
                    }

                } // loop over dir_ptr w. offset

            // Remove dir_ptrs which are offset
            // need to do this only if there are no-offset dirs
            std::vector<dir_ptr> dp;
            for (dir_ptr &d : it->dir_list) {
                if (!d.is_offset)
                    dp.push_back(d);
            }
            it->dir_list.clear();
            it->dir_list = dp;

            // if all references to this field var are offsets, remove the ref.
            // reset the status too
            it->is_read_nb = it->is_read_offset = false;
            std::vector<field_ref *> new_ref_list = {};
            for (field_ref *fr : it->ref_list) {

                if (!fr->is_offset) {
                    new_ref_list.push_back(fr);
                    if (fr->is_direction)
                        it->is_read_nb = true;

                } else {
                    // turn offset off, this ref is to the newly defined field (see
                    // above)
                    fr->is_offset = false;
                }
            }

            if (new_ref_list.size() > 0) {
                // copy the references, keep the iterator element
                it->ref_list = new_ref_list;
                ++it;

            } else {
                // now remove the field altogether from list
                // iterator will point to the next element in the list after erase
                it = field_info_list.erase(it);
            }
        }
    }
}

/// Call the backend function for handling loop functions
void GeneralVisitor::backend_handle_loop_function(call_info_struct &ci) {
    // we should mark the function, but it is not necessarily in the
    // main file buffer
    if (target.CUDA) {
        handle_loop_function_cuda(ci);
    } else if (target.openacc) {
        handle_loop_function_openacc(ci.funcdecl);
    } else if (target.vectorize) {
        handle_loop_function_avx(ci);
    }
}

/// Call the backend function for handling loop functions
void GeneralVisitor::backend_handle_loop_constructor(call_info_struct &ci) {
    // we should mark the function, but it is not necessarily in the
    // main file buffer
    if (target.CUDA) {
        handle_loop_constructor_cuda(ci);
    } else if (target.openacc) {
        handle_loop_constructor_openacc(ci.ctordecl);
    } else if (target.vectorize) {
        handle_loop_constructor_avx(ci);
    }
}

/// Call the backend function for generating loop code
std::string TopLevelVisitor::backend_generate_code(Stmt *S, bool semicolon_at_end,
                                                   srcBuf &loopBuf,
                                                   bool generate_wait_loops) {
    std::stringstream code;
    if (target.CUDA) {
        code << generate_code_cuda(S, semicolon_at_end, loopBuf, generate_wait_loops);
    } else if (target.openacc) {
        code << generate_code_cpu(S, semicolon_at_end, loopBuf,
                                  generate_wait_loops); // use cpu method for acc
    } else if (target.vectorize) {
        code << generate_code_avx(S, semicolon_at_end, loopBuf, generate_wait_loops);
    } else {
        code << generate_code_cpu(S, semicolon_at_end, loopBuf, generate_wait_loops);
    }
    return code.str();
}

/* Generate a header that marks field references read or written.
 * This is a copy of the loop body with modifications, only ran once
 */
/*
std::string TopLevelVisitor::generate_loop_header(Stmt *S, codetype & target, bool
semicolon_at_end) { srcBuf loopBuf;
loopBuf.copy_from_range(writeBuf,S->getSourceRange()); std::vector<std::string> va = {},
vb = {}; int i=0;

  // Replace loop references with temporary variables
  // and add calls to mark_changed() and start_fetch()
  for ( field_info & fi : field_info_list ) {
    std::string type_name = fi.type_template;
    type_name.erase(0,1).erase(type_name.end()-1, type_name.end());

    // Add a simple temp variable to replace the field reference
    std::string varname = "__v_" + std::to_string(i);
    loopBuf.prepend(type_name + " " + varname + ";\n", true);

    for( field_ref *le : fi.ref_list ){
      Expr *e = le->nameExpr;
      loopBuf.replace(le->fullExpr, varname );
      if( le->is_written ){
        // Mark changed fields - do it BEFORE using vars, possible alloc
        loopBuf.insert_before_stmt(e, get_stmt_str(e) + ".mark_changed(" +
parity_in_this_loop + ");", true,   true);
      }
      if( le->dirExpr != nullptr ){
        // If a field needs to be communicated, start here
        loopBuf.insert_before_stmt(e, get_stmt_str(e) + ".wait_move(" +
get_stmt_str(le->dirExpr) + ", "
           + parity_in_this_loop + ");", true, true);
      }
      if( le->is_read ){
        // If a field is read, check that is has been allocated
        loopBuf.insert_before_stmt(e, "assert(" + get_stmt_str(e)
          + ".is_allocated());", true, true);
      }

      //va.push_back(get_stmt_str(le->fullExpr));
      //vb.push_back(varname);
      llvm::errs() << get_stmt_str(le->fullExpr) << ":" << varname << "\n";
    }
    i++;
  }

  // Replace reduction variables with copies to avoid changing originals
  // No other variables are allowed to change
  for ( var_info & vi : var_info_list ) {
    if( vi.reduction_type != reduction::NONE ){
      std::string varname = "__v_" + std::to_string(i);
      va.push_back(vi.name);
      vb.push_back(varname);
      loopBuf.prepend(vi.type + " " + varname + "=" + vi.name + ";\n", true);
    }
    i++;
  }

  // Surround by curly brackets to keep new variables local
  loopBuf.replace_tokens( va, vb );
  loopBuf.prepend("{",true);

  if(semicolon_at_end){
    return loopBuf.dump() + ";}\n";
  } else {
    return loopBuf.dump() + "}\n";
  }
}
*/