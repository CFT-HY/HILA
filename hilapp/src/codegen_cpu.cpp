
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

std::string TopLevelVisitor::generate_code_cpu(Stmt *S, bool semicolon_at_end, srcBuf &loopBuf,
                                               bool generate_wait_loops) {
    std::stringstream code;

    // Set loop lattice
    // if (field_info_list.size() > 0) {
    //     std::string fieldname = field_info_list.front().new_name;
    //     code << "const lattice_struct * RESTRICT hila_loop_lattice = " << fieldname
    //          << ".fs->lattice;\n";
    // } else {
    //     // if there is no field in loop at all
    //     code << "const lattice_struct * RESTRICT hila_loop_lattice = lattice;\n";
    // }
    bool first = true;
    if (cmdline::no_interleaved_comm)
        generate_wait_loops = false;
    else
        generate_wait_loops = true;

    // deep halo gathers (X + k*dir) are completed before the loop
    generate_deep_gathers(code);

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
                        if (first)
                            code << "dir_mask_t  _dir_mask_ = 0;\n";
                        first = false;

                        code << "_dir_mask_ |= " << l.new_name << ".start_gather(" << d.direxpr_s
                             << ", " << loop_info.parity_str << ");\n";
                    }
                }
        } else {
            // now loop local dirs - gather all neighbours!
            // TODO: restrict dirs
            if (!generate_wait_loops) {
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < NDIRS; "
                        "++HILA_dir_) {\n"
                     << l.new_name << ".start_gather(HILA_dir_," << loop_info.parity_str
                     << ");\n}\n";
            } else {
                if (first)
                    code << "dir_mask_t  _dir_mask_ = 0;\n";
                first = false;
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < NDIRS; "
                        "++HILA_dir_) {\n"
                     << "_dir_mask_ |= " << l.new_name << ".start_gather(HILA_dir_,"
                     << loop_info.parity_str << ");\n}\n";
            }
        }
    }

    // write wait gathers here also
    if (!generate_wait_loops)
        for (field_info &l : field_info_list)
            if (l.is_loop_local_dir) {
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < NDIRS; "
                        "++HILA_dir_) {\n"
                     << l.new_name << ".wait_gather(HILA_dir_," << loop_info.parity_str
                     << ");\n}\n";
            }

    if (first)
        generate_wait_loops = false; // no communication needed in the 1st place

    bool boundary_layer = is_macro_defined("BOUNDARY_LAYER_LAYOUT");
    bool evenfirst = is_macro_defined("EVEN_SITES_FIRST");

    code << "const lattice_struct & hila_loop_lattice = lattice.ref();\n";

    // Lexicographic layout: the loop goes through x-rows of the node.  Within a row all
    // neighbour indices are affine in x (nb0[d] + x, or i +- 1 for x), only the row ends need
    // the general neighbours.  The loop body is a lambda called for each site, the interior
    // of a row has no branches and can be vectorized.  See applications/bench_rows.
    bool rows_loop = !evenfirst && !boundary_layer && !target.vectorize && !target.openacc;

    // Set the start and end points
    if (!boundary_layer) {
        // with not EVEN_SITES_FIRST go over all sites
        const char *paritystr = evenfirst ? loop_info.parity_str.c_str() : "ALL";
        code << "const int _hila_loop_begin = hila_loop_lattice.loop_begin(" << paritystr << ");\n";
        code << "const int _hila_loop_end   = hila_loop_lattice.loop_end(" << paritystr << ");\n";

        if (rows_loop) {
            code << "const hila::site_layout_struct & HILA_lay_ = hila_loop_lattice.mynode.layout;\n";
            code << "const int HILA_nx_ = HILA_lay_.size[0];\n";
            code << "const int HILA_rows_ = HILA_lay_.volume / HILA_nx_;\n";
        }

        if (generate_wait_loops) {
            code << "for (int _hila_wait_i = 0; _hila_wait_i < 2; ++_hila_wait_i) {\n";
        }

    } else {

        code << "hila::iter_range_t _hila_ranges;\n";
        code << "int _hila_loops = hila_loop_lattice.loop_ranges(" << loop_info.parity_str << ", "
             << (generate_wait_loops ? "_dir_mask_ != 0" : "false") << ", _hila_ranges);\n";

        code << "for (int _hila_wait_i = 0; _hila_wait_i < _hila_loops; ++_hila_wait_i) {\n";
        code << "const int _hila_loop_begin = _hila_ranges.min[_hila_wait_i];\n";
        code << "const int _hila_loop_end   = _hila_ranges.max[_hila_wait_i];\n";
    }

    // and the openacc loop header
    if (target.openacc) {
        generate_openacc_loop_header(code);
    } else if (target.openmp && !loop_info.contains_random) {

        // CHECK if the loop has ReductionVector
        // This is stopgap to make openmp to work..
        // TODO: copy vector to n_threads arrays
        // Same holds for SiteSelect
        bool has_reductionvector = false;
        for (array_ref &ar : array_ref_list) {
            if (ar.type == array_ref::REDUCTION)
                has_reductionvector = true;
        }

        if (!has_reductionvector && selection_info_list.size() == 0) {

            int sums = 0;
            for (reduction_expr &r : reduction_list) {
                if (r.reduction_type != reduction::NONE &&
                    get_number_type(r.type) == number_type::UNKNOWN) {
                    code << "#pragma omp declare reduction(_hila_reduction_sum" << sums << ":"
                         << r.type << ":omp_out += omp_in)\n";
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
                        sums++;
                    } else {
                        code << '+';
                    }
                    code << ": " << r.reduction_name << ")";
                }
            }
            code << '\n';
        }
    }

    // keep track of ending braces for options
    int ending_braces = 0;

    extern bool use_philox_rng;

    // With LAYOUT_ON_THE_FLY neighbour indices are computed from the site coordinates,
    // which are computed once per site
    // (always without EVEN_SITES_FIRST)
    bool layout_on_the_fly = is_macro_defined("LAYOUT_ON_THE_FLY") && !target.vectorize;
    bool any_read_nb = false;
    for (field_info &l : field_info_list)
        if (l.is_read_nb)
            any_read_nb = true;

    bool wait_loop = generate_wait_loops && !boundary_layer;

    // without EVEN_SITES_FIRST parity loops go through all sites, check the site parity
    bool parity_check = !evenfirst && loop_info.parity_value != Parity::all;

    // site coordinates HILA_site_ needed?  In the lexicographic layout X.coordinates() etc.
    // are replaced with HILA_site_ -using calls (need_loop_coordinate, loop_function.cpp)
    bool need_site =
        (layout_on_the_fly &&
         (any_read_nb || wait_loop || parity_check ||
          (!evenfirst && (loop_info.need_loop_coordinate || use_philox_rng)))) ||
        has_deep_refs();

    if (rows_loop)
        need_site = true;

    if (rows_loop) {
        // loop over rows, row start site and its neighbours
        code << "for (int HILA_row_ = 0; HILA_row_ < HILA_rows_; ++HILA_row_) {\n";
        code << "const unsigned HILA_i0_ = HILA_row_ * HILA_nx_;\n";
        code << "const hila::layout_site HILA_row_site_ = HILA_lay_.site(HILA_i0_);\n";
        if (any_read_nb) {
            code << "unsigned HILA_nb0_[NDIRS];\n";
            code << "for (int HILA_d_ = 0; HILA_d_ < NDIRS; ++HILA_d_) HILA_nb0_[HILA_d_] = "
                    "HILA_lay_.neighbour(HILA_row_site_, HILA_d_);\n";
        }
        if (wait_loop) {
            // off-node directions of the row interior: no x-directions
            code << "const dir_mask_t HILA_row_mask_ = hila_loop_lattice.wait_mask(HILA_row_site_)"
                    " & ~(dir_mask_t)((1 << e_x) | (1 << (-e_x)));\n";
            // scope in place of the wait-if below, closed by 'end of loop content'.  As for
            // the wait-if, ending_braces counts the closing of the wait loop
            code << "{\n";
            ending_braces++;
        }
        // the body for one site
        // force inlining: large bodies (e.g. SU(3)) were not inlined, and were 3-4x slower
        code << "auto HILA_body_ = [&](const int HILA_x_, const unsigned * RESTRICT HILA_nb_) "
                "__attribute__((always_inline)) {\n";
        code << "const unsigned " << looping_var << " = HILA_i0_ + HILA_x_;\n";
        code << "hila::layout_site HILA_site_ = HILA_row_site_;\n";
        code << "HILA_site_.lc[0] = HILA_x_; HILA_site_.lex = " << looping_var
             << "; HILA_site_.parity ^= HILA_x_ & 1;\n";
    } else {
        // Start the loop
        code << "for(int " << looping_var << " = _hila_loop_begin; " << looping_var
             << " < _hila_loop_end; ++" << looping_var << ") {\n";
    }

    if (need_site && !rows_loop) {
        code << "const hila::layout_site HILA_site_ = hila_loop_lattice.mynode.layout.site("
             << looping_var << ");\n";
    }

    // index of deep halo ref X + k*dir
    auto deep_index = [&](const dir_ptr &d) -> std::string {
        return "hila_loop_lattice.mynode.layout.neighbour_offset(HILA_site_, " +
               deep_dir_string(d, loopBuf, true) + ", " + std::to_string(d.deep_k) + ")";
    };

    // expression for the index of the neighbour of fieldname to direction dirname
    bool special_bc = is_macro_defined("SPECIAL_BOUNDARY_CONDITIONS");
    auto nb_index = [&](const std::string &fieldname, const std::string &dirname) -> std::string {
        if (rows_loop) {
            // field neighbour arrays exist only for special boundary conditions
            if (special_bc)
                return "hila_loop_lattice.neighbour_index(" + fieldname + ".fs->neighbours[" +
                       dirname + "], HILA_nb_[" + dirname + "], " + looping_var + ")";
            else
                return "HILA_nb_[" + dirname + "]";
        }
        if (layout_on_the_fly)
            return "hila_loop_lattice.neighbour_index(" + fieldname + ".fs->neighbours[" +
                   dirname + "], HILA_site_, " + looping_var + ", " + dirname + ")";
        else
            return fieldname + ".fs->neighbours[" + dirname + "][" + looping_var + "]";
    };

    if (wait_loop && !rows_loop) {
        if (layout_on_the_fly)
            code << "if (((hila_loop_lattice.wait_mask(HILA_site_)";
        else
            code << "if (((hila_loop_lattice.wait_arr_[" << looping_var << "]";
        code << " & _dir_mask_) != 0) == _hila_wait_i) {\n";
        ending_braces++;
    }

    if (use_philox_rng)
        code << "hila::philox_setup(hila_philox_loop_counter, "
                "SiteIndex("
             << (need_site ? "hila_loop_lattice.mynode.layout.coordinates(HILA_site_)"
                           : "hila_loop_lattice.coordinates(" + looping_var + ")")
             << ").value);\n";

    if (parity_check) {
        // Parity value: even = 1, odd = 2, all = 3, site parity 0 or 1
        // closed after the setters, inside the site loop.  Runtime parity may be ALL: the
        // first, loop invariant condition lets the compiler unswitch the loop
        if (loop_info.parity_value == Parity::none)
            code << "if (" << loop_info.parity_str << " == Parity::all || ";
        else
            code << "if (";
        code << "((static_cast<unsigned>(" << loop_info.parity_str
             << ") >> HILA_site_.parity) & 1)) {\n";
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
                                        ".select_site(SiteIndex(hila_loop_lattice.coordinates(" +
                                        looping_var + ")))");
        } else {
            SourceRange r(si.MCE->getSourceRange().getBegin(),
                          si.assign_expr->getSourceRange().getBegin().getLocWithOffset(-1));
            loopBuf.replace(r, si.new_name +
                                   ".select_site_value(SiteIndex(hila_loop_lattice.coordinates(" +
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
                    if (d.deep_k > 0) {
                        code << "const " << l.element_type << " " << d.name_with_dir << " = "
                             << l.new_name << ".get_value_at(" << deep_index(d) << ");\n";
                        for (field_ref *ref : d.ref_list)
                            loopBuf.replace(ref->fullExpr, d.name_with_dir);
                        continue;
                    }

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

                    if (target.vectorize /* && l.vecinfo.is_vectorizable */) {
                        // now l is vectorizable, but accessed sequentially -- this inly
                        // happens in vectorized targets
                        code << ".get_value_at_nb_site(" << dirname << ", " << looping_var
                             << ");\n";
                    } else {
                        // std neighbour accessor for scalars
                        code << ".get_value_at(" << nb_index(l.new_name, dirname) << ");\n";
                    }

                    // and replace references in loop body
                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, d.name_with_dir);
                    }
                }

            } else {

                // and variable direction refs - use accessor directly
                for (dir_ptr &d : l.dir_list) {
                    if (d.deep_k > 0) {
                        for (field_ref *ref : d.ref_list)
                            loopBuf.replace(ref->fullExpr,
                                            l.new_name + ".get_value_at(" + deep_index(d) + ")");
                        continue;
                    }

                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(
                            loopBuf.get(d.parityExpr->getSourceRange())); // mapped name was

                    for (field_ref *ref : d.ref_list) {
                        if (target.vectorize /* && l.vecinfo.is_vectorizable */) {
                            loopBuf.replace(ref->fullExpr, l.new_name + ".get_value_at_nb_site(" +
                                                               dirname + ", " + looping_var + ")");
                        } else {
                            loopBuf.replace(ref->fullExpr, l.new_name + ".get_value_at(" +
                                                               nb_index(l.new_name, dirname) +
                                                               ")");
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
        loopBuf.replace(sfc.replace_range, sfc.replace_expression);
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

    // end of the parity check
    if (parity_check)
        code << "}\n";

    if (rows_loop) {
        // end of the body lambda, and calls for the sites of the row
        code << "};\n";

        // row ends with general neighbours
        code << "for (int HILA_e_ = 0; HILA_e_ < (HILA_nx_ > 1 ? 2 : 1); ++HILA_e_) {\n";
        code << "const int HILA_x_ = (HILA_e_ == 0) ? 0 : HILA_nx_ - 1;\n";
        code << "hila::layout_site HILA_s_ = HILA_row_site_;\n";
        code << "HILA_s_.lc[0] = HILA_x_; HILA_s_.lex += HILA_x_; HILA_s_.parity ^= HILA_x_ & "
                "1;\n";
        if (wait_loop)
            code << "if (((hila_loop_lattice.wait_mask(HILA_s_) & _dir_mask_) != 0) != "
                    "_hila_wait_i) continue;\n";
        code << "unsigned HILA_nbe_[NDIRS];\n";
        if (any_read_nb)
            code << "for (int HILA_d_ = 0; HILA_d_ < NDIRS; ++HILA_d_) HILA_nbe_[HILA_d_] = "
                    "HILA_lay_.neighbour(HILA_s_, HILA_d_);\n";
        code << "HILA_body_(HILA_x_, HILA_nbe_);\n";
        code << "}\n";

        // row interior: unit stride neighbours
        if (wait_loop)
            code << "if (((HILA_row_mask_ & _dir_mask_) != 0) == _hila_wait_i)\n";
        code << "for (int HILA_x_ = 1; HILA_x_ < HILA_nx_ - 1; ++HILA_x_) {\n";
        code << "unsigned HILA_nbi_[NDIRS];\n";
        if (any_read_nb) {
            code << "for (int HILA_d_ = 0; HILA_d_ < NDIRS; ++HILA_d_) HILA_nbi_[HILA_d_] = "
                    "HILA_nb0_[HILA_d_] + HILA_x_;\n";
            code << "HILA_nbi_[e_x] = HILA_i0_ + HILA_x_ + 1;\n";
            code << "HILA_nbi_[-e_x] = HILA_i0_ + HILA_x_ - 1;\n";
        }
        code << "HILA_body_(HILA_x_, HILA_nbi_);\n";
        code << "}\n";
    }

    // End of the loop content
    code << "}\n";

    if (generate_wait_loops) {
        if (!boundary_layer) {
            // add the code for 2nd round - also need one } to balance the if ()
            code << "}\nif (_dir_mask_ == 0 || _hila_wait_i > 0) break;    // No need for another "
                    "round\n";
        } else {
            code << "if (_dir_mask_ != 0 && _hila_wait_i == 0) {\n";
            ending_braces++;
        }
        for (field_info &l : field_info_list) {
            // If neighbour references exist, communicate them
            if (!l.is_loop_local_dir) {
                for (dir_ptr &d : l.dir_list)
                    if (d.count > 0) {
                        code << l.new_name << ".wait_gather(" << d.direxpr_s << ", "
                             << loop_info.parity_str << ");\n";
                    }
            } else {
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < NDIRS; "
                        "++HILA_dir_) {\n"
                     << "  " << l.new_name << ".wait_gather(HILA_dir_, " << loop_info.parity_str
                     << ");\n}\n";
            }
        }
        // if (boundary_layer) code << "}\n";
    }

    while (ending_braces-- > 0)
        code << "}\n";

    // Post-process ny site selections?
    for (selection_info &s : selection_info_list) {
        if (s.previous_selection == nullptr) {
            code << s.new_name << ".endloop_action();\n";
        }
    }


    return code.str();
}
