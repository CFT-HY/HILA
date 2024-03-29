
//------------------------------------------------------------------------------
// Generate code for plain standard platform - should work on any computer
//
//
//------------------------------------------------------------------------------
#include <sstream>
#include <string>

#include "hilapp.h"
#include "toplevelvisitor.h"
#include "stringops.h"

extern std::string looping_var;
extern std::string parity_name;

std::string TopLevelVisitor::generate_code_cpu(Stmt *S, bool semicolon_at_end, srcBuf &loopBuf,
                                               bool generate_wait_loops) {
    std::stringstream code;

    // Set loop lattice
    // if (field_info_list.size() > 0) {
    //     std::string fieldname = field_info_list.front().new_name;
    //     code << "const lattice_struct * RESTRICT loop_lattice = " << fieldname
    //          << ".fs->lattice;\n";
    // } else {
    //     // if there is no field in loop at all
    //     code << "const lattice_struct * RESTRICT loop_lattice = lattice;\n";
    // }

    code << "const lattice_struct & loop_lattice = lattice;\n";

    // Set the start and end points
    code << "const int loop_begin = loop_lattice.loop_begin(" << loop_info.parity_str << ");\n";
    code << "const int loop_end   = loop_lattice.loop_end(" << loop_info.parity_str << ");\n";

    // are there

    if (generate_wait_loops) {
        code << "for (int _wait_i_ = 0; _wait_i_ < 2; ++_wait_i_) {\n";
    }

    // and the openacc loop header
    if (target.openacc) {
        generate_openacc_loop_header(code);
    } else if (target.openmp && !loop_info.contains_random) {
        int sums = 0;
        for (reduction_expr &r : reduction_list) {
            if (r.reduction_type != reduction::NONE &&
                get_number_type(r.type) == number_type::UNKNOWN) {
                code << "#pragma omp declare reduction(_hila_reduction_sum" << sums << ":" << r.type
                     << ":omp_out += omp_in)\n";
                sums++;
            }
        }
        if (loop_info.has_pragma_omp_parallel_region)
            code << "#pragma omp for";
        else
            code << "#pragma omp parallel for";

        sums = 0;
        for (reduction_expr &r : reduction_list) {
            if (r.reduction_type != reduction::NONE) {
                code << " reduction(";
                if (get_number_type(r.type) == number_type::UNKNOWN) {
                    code << "_hila_reduction_sum" << sums;
                } else {
                    code << '+';
                }
                code << ": " << r.reduction_name << ")";
            }
        }
        code << '\n';
    }


    // Start the loop
    code << "for(int " << looping_var << " = loop_begin; " << looping_var << " < loop_end; ++"
         << looping_var << ") {\n";

    if (generate_wait_loops) {
        code << "if (((loop_lattice.wait_arr_[" << looping_var
             << "] & _dir_mask_) != 0) == _wait_i_) {\n";
    }

    // replace reduction variables in the loop
    for (reduction_expr &r : reduction_list) {
        for (Expr *e : r.refs) {
            loopBuf.replace(e, r.reduction_name);
        }
    }

    // replace selection var reference
    for (selection_info &si : selection_info_list) {
        if (si.assign_expr == nullptr) {
            loopBuf.replace(si.MCE, si.new_name +
                                        ".select_site(SiteIndex(loop_lattice.coordinates(" +
                                        looping_var + ")))");
        } else {
            SourceRange r(si.MCE->getSourceRange().getBegin(),
                          si.assign_expr->getSourceRange().getBegin().getLocWithOffset(-1));
            loopBuf.replace(r, si.new_name +
                                   ".select_site_value(SiteIndex(loop_lattice.coordinates(" +
                                   looping_var + ")), ");
        }
    }

    // Create temporary field element variables
    for (field_info &l : field_info_list) {

        // First check for Direction references. If any found, create list of temp
        // variables
        if (l.is_read_nb) {
            if (!l.is_loop_local_dir) {
                // "standard" loop extern neighb direction
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(
                            loopBuf.get(d.parityExpr->getSourceRange())); // mapped name was
                                                                          // get_stmt_str(d.e);

                    // generate access stmt
                    code << "const " << l.element_type << " " << d.name_with_dir << " = "
                         << l.new_name;

                    if (target.vectorize && l.vecinfo.is_vectorizable) {
                        // now l is vectorizable, but accessed sequentially -- this inly
                        // happens in vectorized targets
                        code << ".get_value_at_nb_site(" << dirname << ", " << looping_var
                             << ");\n";
                    } else {
                        // std neighbour accessor for scalars
                        code << ".get_value_at(" << l.new_name << ".fs->neighbours[" << dirname
                             << "][" << looping_var << "]);\n";
                    }

                    // and replace references in loop body
                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, d.name_with_dir);
                    }
                }

            } else {

                // and variable direction refs - use accessor directly
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(
                            loopBuf.get(d.parityExpr->getSourceRange())); // mapped name was

                    for (field_ref *ref : d.ref_list) {
                        if (target.vectorize && l.vecinfo.is_vectorizable) {
                            loopBuf.replace(ref->fullExpr, l.new_name + ".get_value_at_nb_site(" +
                                                               dirname + ", " + looping_var + ")");
                        } else {
                            loopBuf.replace(ref->fullExpr, l.new_name + ".get_value_at(" +
                                                               l.new_name + ".fs->neighbours[" +
                                                               dirname + "][" + looping_var + "])");
                        }
                    }
                }
            }
        }


        // and then get (possible) local refs
        // TODO:
        if (l.is_read_atX || (l.is_written && loop_info.has_conditional)) {
            // now reading var without nb. reference
            // const here may cause problems in loop functions!
            // if (!l.is_written) {
            //     code << "const ";
            // }
            code << l.element_type << " " << l.loop_ref_name << " = " << l.new_name
                 << ".get_value_at(" << looping_var << ");\n";

            if (!l.is_read_atX) {
                code << "// Value of var " << l.loop_ref_name
                     << " read in because loop has conditional\n";
                code << "// TODO: MAY BE UNNECESSARY, write more careful analysis\n";
            }

        } else if (l.is_written) {
            code << l.element_type << " " << l.loop_ref_name << ";\n";
            code << "// Initial value of variable " << l.loop_ref_name << " not needed\n";
        }

        // and finally replace references in body
        for (field_ref *ref : l.ref_list)
            if (!ref->is_direction) {
                loopBuf.replace(ref->fullExpr, l.loop_ref_name);
            }

    } // end of field handling

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
            if (sfc.args_string.size() > 0)
                repl += ", " + sfc.args_string;
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
            code << l.new_name << ".set_value_at(" << l.loop_ref_name << ", " << looping_var
                 << ");\n";
        }

    code << "}\n";

    if (generate_wait_loops) {
        // add the code for 2nd round - also need one } to balance the if ()
        code << "}\nif (_dir_mask_ == 0) break;    // No need for another round\n";

        for (field_info &l : field_info_list) {
            // If neighbour references exist, communicate them
            if (!l.is_loop_local_dir) {
                for (dir_ptr &d : l.dir_list)
                    if (d.count > 0) {
                        code << l.new_name << ".wait_gather(" << d.direxpr_s << ", "
                             << loop_info.parity_str << ");\n";
                    }
            } else {
                code << "for (Direction _HILAdir_ = (Direction)0; _HILAdir_ < NDIRS; "
                        "++_HILAdir_) {\n"
                     << "  " << l.new_name << ".wait_gather(_HILAdir_, " << loop_info.parity_str
                     << ");\n}\n";
            }
        }
        code << "}\n";
    }

    // Post-process ny site selections?
    for (selection_info &s : selection_info_list) {
        if (s.previous_selection == nullptr) {
            code << s.new_name << ".endloop_action();\n";
        }
    }


    return code.str();
}
