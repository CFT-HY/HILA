///------------------------------------------------------------------------------
/// Generate transformed
/// hardware-dependent "kernels".  Main interface to code generation
///
///
///------------------------------------------------------------------------------
#include <sstream>
#include <string>


#include "hilapp.h"
#include "toplevelvisitor.h"
#include "stringops.h"

std::string looping_var;
std::string parity_name;

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

///////////////////////////////////////////////////////////////////////////////
/// The main entry point for code generation
///////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::generate_code(Stmt *S) {
    srcBuf loopBuf; // (&TheRewriter,S);

    // check if the range starts with a macro (e.g. onsites(ALL) foralldir(d) ...)
    // Std sourcerange fails here!

    SourceRange Srange = get_real_range(S->getSourceRange());

    loopBuf.copy_from_range(writeBuf, Srange);

    //   llvm::errs() << "\nOriginal range: +++++++++++++++\n\""
    //                << TheRewriter.getRewrittenText(S->getSourceRange())
    //                << "\"\nwriteBuf range: ================\n\""
    //                << writeBuf->get(S->getSourceRange())
    //                << "\"\nCopied range: ================\n\""
    //                << loopBuf.dump() << "\"\n";

    // is there semicolon at the end? Keep track of it (not visible in AST)
    bool semicolon_at_end = hasSemicolonAfter(Srange);


    // Build replacement in variable "code"
    // Encapsulate everything within {}
    std::stringstream code;
    code << "{\n";


    if (loop_info.contains_random) {
        code << "hila::check_that_rng_is_initialized();\n";
    }

    const std::string t = loopBuf.dump();

    looping_var = "_HILA_index";
    while (t.find(looping_var, 0) != std::string::npos)
        looping_var += "_";

    // Ensure that the name is not reserved by scanning the source
    parity_name = "_HILA_parity";
    while (t.find(parity_name, 0) != std::string::npos)
        parity_name += "_";

    if (loop_info.parity_value == Parity::none) {
        // now unknown
        code << "const Parity " << parity_name << " = " << loop_info.parity_text << ";\n";

        if (global.assert_loop_parity) {
            code << "assert( is_even_odd_parity(" << parity_name
                 << ") && \"Parity should be EVEN or ODD in this loop\");\n";
        }
        loop_info.parity_str = parity_name;

    } else
        loop_info.parity_str = parity_str(loop_info.parity_value);

    // any site selections?  reset previous_selection
    for (selection_info &s : selection_info_list) {
        if (s.previous_selection == nullptr) {
            code << s.ref->getType().getAsString() << " & " << s.new_name << " = "
                 << get_stmt_str(s.ref) << ";\n";
            code << s.new_name << ".setup();\n";
        }
    }

    // then, generate new names for field variables in loop

    for (field_info &l : field_info_list) {
        // Generate new variable name, may be needed -- use here simple receipe
        l.new_name = field_name_prefix + clean_name(l.old_name);
        // Perhaps simpler FA, FB, FC. ?  The following helps avoid collisions
        l.loop_ref_name = l.new_name + "_at_X";

        // Create neighbour ref names
        int i = 0;
        for (dir_ptr &d : l.dir_list) {
            d.name_with_dir = l.new_name + "_dir" + std::to_string(++i);
        }

        // make a ref to the field name
        if (!l.is_written)
            code << "const ";
        code << "Field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
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
            if (loop_info.parity_value == Parity::all || (l.is_read_nb && l.is_read_atX)) {
                init_par = "ALL";
            } else {
                if (l.is_read_atX)
                    init_par = loop_info.parity_str;
                else
                    init_par = "opp_parity(" + loop_info.parity_str + ")";
            }

            // TAKE THIS AWAY FOR NOW -- WE DON'T HAVE A GOOD METHOD TO CHECK "pure
            // output" FUNCTIONS E.g. method   U[X].random(), where U does not have to
            // be initialized
            // actually, we do now -- the "out_only" keyword

            if (cmdline::check_initialization) {

                std::string fname =
                    srcMgr.getFilename(get_real_range(S->getSourceRange()).getBegin()).str();
                code << "if (!" << l.new_name << ".is_initialized(" << init_par
                     << ")){\nhila::out0 << \"File " << fname << " on line "
                     << srcMgr.getSpellingLineNumber(get_real_range(S->getSourceRange()).getBegin())
                     << ":\\n Value of variable " << l.old_name
                     << " is used but it is not properly initialized\\n\";\n";
                code << "hila::terminate(1);\n}\n";
            }
        }
    }

    // change the f[X+offset] -references, generate code
    handle_field_plus_offsets(code, loopBuf, loop_info.parity_str);

    bool first = true;
    bool generate_wait_loops;
    if (cmdline::no_interleaved_comm)
        generate_wait_loops = false;
    else
        generate_wait_loops = true;

    bool gpu_aware_mpi = is_macro_defined("GPU_AWARE_MPI");

    for (field_info &l : field_info_list) {
        // If neighbour references exist, communicate them
        if (!l.is_loop_local_dir) {
            // "normal" dir references only here
            for (dir_ptr &d : l.dir_list)
                if (d.count > 0) {
                    if (!generate_wait_loops) {
                        code << l.new_name << ".gather(" << d.direxpr_s << ", "
                             << loop_info.parity_str << ");\n";
                    } else {
                        if (first) {
                            code << "dir_mask_t  _dir_mask_ = 0;\n";
                            if (gpu_aware_mpi) {
                                code << "std::vector<send_data_list_t> _send_halo_data;\n";
                            }
                            first = false;
                        }

                        code << "_dir_mask_ |= " << l.new_name << ".start_gather(" << d.direxpr_s
                             << ", " << loop_info.parity_str;

                        if (gpu_aware_mpi)
                            code << ", &_send_halo_data";

                        code << ");\n";
                    }
                }
        } else {
            // now loop local dirs - gather all neighbours!
            // TODO: restrict dirs
            if (!generate_wait_loops) {
                code << "for (Direction _HILAdir_ = (Direction)0; _HILAdir_ < NDIRS; "
                        "++_HILAdir_) {\n"
                     << l.new_name << ".start_gather(_HILAdir_," << loop_info.parity_str
                     << ");\n}\n";
            } else {
                if (first) {
                    code << "dir_mask_t  _dir_mask_ = 0;\n";
                    if (gpu_aware_mpi) {
                        code << "std::vector<send_data_list_t> _send_halo_data;\n";
                    }
                    first = false;
                }
                code << "for (Direction _HILAdir_ = (Direction)0; _HILAdir_ < NDIRS; "
                        "++_HILAdir_) {\n"
                     << "_dir_mask_ |= " << l.new_name << ".start_gather(_HILAdir_,"
                     << loop_info.parity_str;
                if (gpu_aware_mpi)
                    code << ", &_send_halo_data";

                code << ");\n}\n";
            }
        }
    }


    // write wait gathers here also if not wait loops
    if (!generate_wait_loops)
        for (field_info &l : field_info_list)
            if (l.is_loop_local_dir) {
                code << "for (Direction _HILAdir_ = (Direction)0; _HILAdir_ < NDIRS; "
                        "++_HILAdir_) {\n"
                     << l.new_name << ".wait_gather(_HILAdir_," << loop_info.parity_str
                     << ");\n}\n";
            }

    if (first)
        generate_wait_loops = false; // no communication needed in the 1st place

    // do mpi sends only here
    if (gpu_aware_mpi && !first) {
        code << "wait_pack_and_send_halos(_send_halo_data);\n";
    }


    // Create temporary variables for reductions
    for (var_info &v : var_info_list) {
        if (v.reduction_type != reduction::NONE) {
            v.reduction_name = "r" + var_name_prefix + clean_name(v.name);
            // Create a temporary variable and initialize
            if (v.reduction_type == reduction::SUM) {
                code << v.type << " " << v.reduction_name << ";\n";
                code << v.reduction_name << " = 0;\n";
            } else if (v.reduction_type == reduction::PRODUCT) {
                code << v.type << " " << v.reduction_name << ";\n";
                code << v.reduction_name << " = 1;\n";
            }
        }
    }

    // and vector reductions
    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REDUCTION) {
            code << ar.name;
            if (ar.reduction_type == reduction::SUM)
                code << ".init_sum();\n";
            else
                code << ".init_product();\n";
        }
    }

    // Place the content of the loop
    code << backend_generate_code(S, semicolon_at_end, loopBuf, generate_wait_loops);

    // Check reduction variables
    for (var_info &v : var_info_list) {

        if (v.is_special_reduction_type) {

            if (v.reduction_type == reduction::SUM)
                code << v.name << ".reduce_sum_node(" << v.reduction_name << ");\n";
            else if (v.reduction_type == reduction::PRODUCT)
                code << v.name << ".reduce_product_node(" << v.reduction_name << ");\n";

        } else if (v.reduction_type == reduction::PRODUCT) {

            code << "if (hila::myrank() == 0) { " << v.name << " *= " << v.reduction_name
                 << "; }\n";
            code << "else { " << v.name << " = " << v.reduction_name << "; }\n";
            code << "hila::reduce_node_product( &" << v.name << ", 1, true);\n";
        }
    }

    // handle separately sum reductions, by far most common case
    // a bit convoluted way to go through the list, done so to get
    // a bit tidier output than the easiest case

    bool sum_reductions = false;
    for (var_info &v : var_info_list) {
        if (v.reduction_type == reduction::SUM && !v.is_special_reduction_type) {

            if (!sum_reductions) {
                code << "if (hila::myrank() == 0) {\n";
                sum_reductions = true;
            }
            // on node 0 add the old value to reduction
            code << v.name << " += " << v.reduction_name << ";\n";
        }
    }

    if (sum_reductions) {
        // branch where myrank() != 0
        code << "} else {\n";

        for (var_info &v : var_info_list) {
            if (v.reduction_type == reduction::SUM && !v.is_special_reduction_type) {

                // forget the old value on other nodes than 0
                code << v.name << " = " << v.reduction_name << ";\n";
            }
        }

        code << "}\n";

        for (var_info &v : var_info_list) {
            if (v.reduction_type == reduction::SUM && !v.is_special_reduction_type) {
                code << "hila_reduce_sum_setup( &" << v.name << ");\n";
            }
        }

        code << "hila_reduce_sums();\n";
    }

    code << "hila::set_allreduce(true);\n";


    // and vector reductions
    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REDUCTION) {

            if (ar.reduction_type == reduction::SUM)
                code << ar.name << ".reduce_sum();\n";
            else if (ar.reduction_type == reduction::PRODUCT)
                code << ar.name << ".reduce_product();\n";
        }
    }

    // finally mark modified fields
    for (field_info &l : field_info_list)
        if (l.is_written) {
            code << l.new_name << ".mark_changed(" << loop_info.parity_str << ");\n";
        }

    // and close
    code << "}\n//----------\n";

    // Remove old code + replace
    if (semicolon_at_end) {
        writeBuf->remove(getRangeWithSemicolon(Srange));
    } else {
        writeBuf->remove(Srange);
    }

    writeBuf->insert(Srange.getBegin(), indent_string(code.str()), true, true);
}

// Handle field+offset expressions -- these are copied, ref removed
// Use iterator for looping, because we may want to remove the field_info item.

void TopLevelVisitor::handle_field_plus_offsets(std::stringstream &code, srcBuf &loopBuf,
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
                        offset_field_name += std::to_string(i_offset);
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
                    code << "Field" + it->type_template + " " + offset_field_name + ";\n";
                    code << it->new_name + ".shift(" + d.ref_list.at(0)->direxpr_s + ", " +
                                offset_field_name + ", " + paritystr + ");\n";

                    // and rewrite references to the offset field
                    for (field_ref *fr : new_fi.ref_list) {
                        loopBuf.replace(fr->nameExpr, offset_field_name);
                        loopBuf.replace(fr->parityExpr, "X");
                        fr->direxpr_s.clear();
                        fr->is_direction = false;           // no Direction
                        fr->info = &field_info_list.back(); // info points to new field_info
                        fr->is_written = false;
                        fr->is_read = true;
                        fr->is_offset = true; // leave this on, as a flag -- turned off below
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
    if (target.kernelize) {
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
    if (target.kernelize) {
        handle_loop_constructor_cuda(ci);
    } else if (target.openacc) {
        handle_loop_constructor_openacc(ci.ctordecl);
    } else if (target.vectorize) {
        handle_loop_constructor_avx(ci);
    }
}

/// Call the backend function for generating loop code
std::string TopLevelVisitor::backend_generate_code(Stmt *S, bool semicolon_at_end, srcBuf &loopBuf,
                                                   bool generate_wait_loops) {
    std::stringstream code;
    if (target.kernelize) {
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
  // and add calls to mark_changed() and start_gather()
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
loop_info.parity_str + ");", true,   true);
      }
      if( le->dirExpr != nullptr ){
        // If a field needs to be communicated, start here
        loopBuf.insert_before_stmt(e, get_stmt_str(e) + ".wait_move(" +
get_stmt_str(le->dirExpr) + ", "
           + loop_info.parity_str + ");", true, true);
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