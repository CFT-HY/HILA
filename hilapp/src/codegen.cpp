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
// Create the list of reductions from varible list and from const_expr list
// (leaving Reductionvector alone)
// -> global variable reduction_list

void TopLevelVisitor::create_reduction_list(std::list<var_info> &vi_list,
                                            std::list<loop_const_expr_ref> &ce_list) {

    reduction_list.clear();

    // Find reductions from variables
    for (var_info &v : vi_list) {
        if (v.reduction_type != reduction::NONE) {

            reduction_expr re;
            re.reduction_type = v.reduction_type;
            re.name = v.name;
            re.reduction_name = "r" + name_prefix + clean_name(v.name) + "_";
            re.type = v.type;
            // and get the reference exprs here
            for (auto &vref : v.refs) {
                if (Expr *E = dyn_cast<Expr>(vref.ref)) {
                    re.refs.push_back(E);
                } else {
                    llvm::errs() << "VREF DYN CAST ERROR!!!\n";
                }
            }
            re.is_special_reduction = v.is_special_reduction_type;
            re.variable = &v;

            if (target.vectorize) {
                re.vecinfo = v.vecinfo;
            }

            reduction_list.push_back(re);
        }
    }

    int cexpr_ind = 0;
    for (loop_const_expr_ref &r : ce_list) {
        if (r.reduction_type != reduction::NONE) {

            reduction_expr re;
            re.reduction_type = r.reduction_type;
            re.name = r.expression;
            cexpr_ind++;
            re.reduction_name = "r" + name_prefix + "e" + std::to_string(cexpr_ind) + "_";
            re.type = r.type;
            for (auto &eref : r.refs)
                re.refs.push_back(eref);
            re.is_special_reduction = false;
            re.variable = nullptr;

            if (target.vectorize) {
                re.vecinfo.is_vectorizable = is_vectorizable_type(re.type, re.vecinfo);
            }

            reduction_list.push_back(re);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
/// Helper for printing filename and line on the file
///////////////////////////////////////////////////////////////////////////////

std::string TopLevelVisitor::get_filename(Stmt *S) {
    auto nameref = srcMgr.getFilename(get_real_range(S->getSourceRange()).getBegin());
    auto rp = nameref.rfind('/');
    if (rp != StringRef::npos)
        nameref = nameref.substr(rp + 1);
    return nameref.str();
}

unsigned TopLevelVisitor::get_linenumber(Stmt *S) {
    return srcMgr.getSpellingLineNumber(get_real_range(S->getSourceRange()).getBegin());
}

std::string TopLevelVisitor::get_filename_and_line(Stmt *S) {
    std::string result;

    result = "File \"" + get_filename(S) + "\" on line " + std::to_string(get_linenumber(S));
    return result;
}


///////////////////////////////////////////////////////////////////////////////
/// The main entry point for code generation
///////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::generate_code(Stmt *S) {
    srcBuf loopBuf; // (&TheRewriter,S);

    // check if the range starts with a macro (e.g. onsites(ALL) foralldir(d) ...)
    // Std sourcerange fails here!

    static int is_cow = -1;
    if (is_cow < 0) {
        is_cow = is_macro_defined("FIELD_COPY_ON_WRITE");
        if (is_cow)
            llvm::errs() << " FIELD_COPY_ON_WRITE defined\n";
    }

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


    // Build replacement in variable "code", and additional after-loop commands in "tail"
    // Encapsulate everything within {}
    std::stringstream code;
    code << "{\n";

    if (loop_info.contains_random)
        code << "hila::check_that_rng_is_initialized();\n";

    generate_parity_code(code, parity_var_name);
    generate_selection_code(code);

    // then, generate new names for field variables in loop
    handle_field_info(code, is_cow, S);
    // for the gathers needed only with some values of variables defined outside the loop
    code << loop_param_bit_decl();

    // change the f[X+offset] -references, generate code
    std::string tail;
    handle_field_plus_offsets(code, loopBuf, loop_info.parity_str, tail);

    bool generate_wait_loops;

    if (cmdline::no_interleaved_comm)
        generate_wait_loops = false;
    else
        generate_wait_loops = true;

    // generate field comm and access code
    // generate_field_code(code, generate_wait_loops);

    // make the reduction expr list from variables and loop const expressions
    create_reduction_list(var_info_list, loop_const_expr_ref_list);
    handle_reduction_init(code);

    handle_loop_constants(code, loopBuf);

    ///////////////////////////////////////////////////
    // Place the content of the loop
    code << backend_generate_code(S, semicolon_at_end, loopBuf, generate_wait_loops);

    code << tail;

    handle_reduction_result(code);

    mark_fields_changed(code);

    if (loop_info.contains_random)
        code << "hila::philox_after_onsites();   // reset philox to host\n";

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

////////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::generate_parity_code(std::stringstream &code,
                                           const std::string &parity_name) {

    if (loop_info.parity_value == Parity::none) {
        // now unknown
        code << "const Parity " << parity_name << " = " << loop_info.parity_text << ";\n";

        if (global.assert_loop_parity) {
            code << "assert( is_even_odd_parity(" << parity_name
                 << ") && \"Parity should be EVEN or ODD in this loop\");\n";
        }
        loop_info.parity_str = parity_name;

        loop_info.need_loop_coordinate = true;
    } else {
        loop_info.parity_str = parity_str(loop_info.parity_value);
        if (loop_info.parity_value != Parity::all)
            loop_info.need_loop_coordinate = true;
    }
}

///////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::generate_selection_code(std::stringstream &code) {

    // any site selections?  reset previous_selection
    for (selection_info &s : selection_info_list) {
        if (s.previous_selection == nullptr) {
            code << s.ref->getType().getAsString() << " & " << s.new_name << " = "
                 << get_stmt_str(s.ref) << ";\n";
            code << s.new_name << ".setup();\n";
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::handle_field_info(std::stringstream &code, bool is_cow, Stmt *S) {

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

        if (l.is_loop_local_ref) {
            // A field with a loop local index, e.g. U[mu], can not be referred to here. Make a
            // ref to each field it can be, these are not written to..
            const std::vector<std::string> names = l.outside_loop_names();
            auto idx = l.ref_idx_set.begin();
            for (const std::string &name : names) {
                code << "const Field" << l.type_template << " & " << name << " = "
                     << l.loop_local_ref_expr(*idx++) << ";\n";
            }
            continue;
        }

        // make a ref to the field name
        if (!l.is_written)
            code << "const ";
        code << "Field" << l.type_template << " & " << l.new_name << " = " << l.old_name << ";\n";
    }

    // check alloc and do it if needed
    for (field_info &l : field_info_list)
        if (l.is_written) {
            for (const std::string &name : l.outside_loop_names()) {
                code << name << ".check_alloc();\n";
                if (is_cow) {
                    code << name << ".will_change();\n";
                }
            }
        }

    // for (field_info &l : field_info_list) {

    //     code << "if (" << l.new_name << ".fs->mylattice.ptr() != lattice.ptr()) {"
    //          << "hila::out0 << \"" << get_filename_and_line(S) << ": Field " << l.old_name
    //          << " initialized on different lattice\\n\";\n"
    //          << "hila::terminate(1);\n}\n";
    // }

    code << "const char * hila_this_file = \"" << get_filename(S) << "\";\n";

    // Check that read fields are initialized
    for (field_info &l : field_info_list) {
        std::string init_par{"Parity::none"};

        // reads at X in field functions count as reads at X here
        bool read_atX = l.is_read_atX || l.is_read_atX_in_func;
        if (l.is_read_nb || read_atX) {
            if (loop_info.parity_value == Parity::all || (l.is_read_nb && read_atX)) {
                init_par = "ALL";
            } else {
                if (read_atX)
                    init_par = loop_info.parity_str;
                else
                    init_par = "opp_parity(" + loop_info.parity_str + ")";
            }

            // TAKE THIS AWAY FOR NOW -- WE DON'T HAVE A GOOD METHOD TO CHECK "pure
            // output" FUNCTIONS E.g. method   U[X].random(), where U does not have to
            // be initialized
            // actually, we do now -- the "out_only" keyword

            // if (cmdline::check_initialization) {

            //     code << "if (!" << l.new_name << ".is_initialized(" << init_par << ")) {\n";
            //     code << "hila::out0 << \"" << get_filename_and_line(S) << ":\\n Variable "
            //          << l.old_name << " is not properly initialized\\n\";\n";
            //     code << "hila::terminate(1);\n}\n";
            // }
        }

        const std::vector<std::string> names = l.outside_loop_names();
        auto idx = l.ref_idx_set.begin();
        for (const std::string &name : names) {
            const std::string field_str =
                l.is_loop_local_ref ? l.loop_local_ref_expr(*idx++) : l.old_name;
            code << name << ".check_init(\"" << field_str << "\", " << init_par
                 << ", hila_this_file, " << get_linenumber(S) << ");\n";
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

/// Used with FETCH_FIELD_REFS_LOOP_LOCAL, not used by default, seems to be slower...
/// For a field with a loop local index, e.g. U[mu] inside
/// foralldir(mu): read all the elements the references can use into arrays at the start of the
/// loop body and replace the references with array accesses. E.g. for U[mu][X] and U[mu][X - mu]
/// with mu = e_x, e_y:
///   const T HILA_U_at_X[2] = {get(U_0, X), get(U_1, X)};
///   const T HILA_U_dir1[2] = {get(U_0, X - e_x), get(U_1, X - e_y)};
/// and U[mu][X - mu] -> HILA_U_dir1[hila::pair_position<0, 7, 1, 6>(mu, -mu)].
/// The positions are constants if the loop over mu is unrolled.
/// NOTE: all elements are read even if a reference is not executed for every value, e.g. if it is
/// inside a branch.
///
/// `get(name, index)` is the expression of the element at `index` of the Field variable `name`
/// (one of l.outside_loop_names()), `nb_index(name, dir)` the index of the neighbour to `dir`.
/// The b.c. of the field is applied to the neighbour values, see bc_nb_value(); have_site: as
/// there.
/// Returns the references that were replaced. The rest, refs to the deep halo or with
/// unknown directions, are left for the caller.
std::set<const field_ref *> TopLevelVisitor::preload_loop_local_refs(
    const field_info &l, std::stringstream &code, srcBuf &loopBuf, const std::string &looping_var,
    const std::function<std::string(const std::string &, const std::string &)> &get,
    const std::function<std::string(const std::string &, const std::string &)> &nb_index,
    bool have_site) {

    std::set<const field_ref *> preloaded;

    const std::vector<std::string> names = l.outside_loop_names();
    std::map<int64_t, std::string> name_of; // index value -> Field variable
    std::string value_list;
    auto name_it = names.begin();
    for (int64_t v : l.ref_idx_set) {
        name_of[v] = *name_it++;
        value_list += (value_list.empty() ? "" : ", ") + std::to_string(v);
    }

    auto index_of = [&](const field_ref *ref) {
        return "static_cast<int64_t>(" + loopBuf.get(ref->ref_index->getSourceRange()) + ")";
    };

    // declare array `array` with elements `elems` and return its name. If an array with the same
    // elements exists, e.g. for X - mu in two foralldir(mu) loops, return the name of that instead
    std::map<std::vector<std::string>, std::string> arrays;
    auto declare_array = [&](const std::string &array,
                             const std::vector<std::string> &elems) -> std::string {
        const auto [it, is_new] = arrays.insert({elems, array});
        if (!is_new)
            return it->second;
        code << "const " << l.element_type << ' ' << array << '[' << elems.size() << "] = {";
        for (size_t i = 0; i < elems.size(); i++)
            code << (i > 0 ? ", " : "") << elems[i];
        code << "};\n";
        return array;
    };

    // refs at X: one element for each index value
    std::vector<const field_ref *> at_X;
    for (const field_ref *ref : l.ref_list)
        if (!ref->is_direction)
            at_X.push_back(ref);
    if (!at_X.empty()) {
        std::vector<std::string> elems;
        for (const auto &[v, name] : name_of)
            elems.push_back(get(name, looping_var));
        const std::string array = declare_array(l.loop_ref_name, elems);
        for (const field_ref *ref : at_X) {
            loopBuf.replace(ref->fullExpr, array + "[hila::value_position<" + value_list + ">(" +
                                               index_of(ref) + ")]");
            preloaded.insert(ref);
        }
    }

    for (const dir_ptr &d : l.dir_list) {
        if (d.deep_k > 0 || d.ref_list.empty())
            continue;

        std::string dirname;
        if (d.is_constant_direction)
            dirname = d.direxpr_s; // orig. string
        else
            dirname = mapped_dir_name(loopBuf, d);

        bool loop_local_dir = false;
        for (const field_ref *ref : d.ref_list)
            if (ref->is_loop_local_dir)
                loop_local_dir = true;

        if (!loop_local_dir) {
            // the same direction for all index values, one element for each
            std::vector<std::string> elems;
            for (const auto &[v, name] : name_of)
                elems.push_back(
                    bc_nb_value(l, get(name, nb_index(name, dirname)), dirname, 1, have_site));
            const std::string array = declare_array(d.name_with_dir, elems);
            for (const field_ref *ref : d.ref_list) {
                loopBuf.replace(ref->fullExpr, array + "[hila::value_position<" + value_list +
                                                   ">(" + index_of(ref) + ")]");
                preloaded.insert(ref);
            }
            continue;
        }

        // a loop local direction: one element for each (index value, direction) pair the refs
        // use. Not possible if a pair is unknown
        std::set<std::pair<int64_t, int64_t>> pairs;
        bool known = true;
        for (const field_ref *ref : d.ref_list) {
            if (!ref->loop_local_dirs_known)
                known = false;
            for (const auto &p : ref->loop_local_dirs) {
                if (name_of.count(p.first) == 0)
                    known = false;
                pairs.insert(p);
            }
        }
        if (!known || pairs.empty())
            continue;

        std::vector<std::string> elems;
        std::string pair_list;
        for (const auto &[v, dir] : pairs) {
            const std::string &name = name_of[v];
            const std::string dir_s = "static_cast<Direction>(" + std::to_string(dir) + ")";
            elems.push_back(bc_nb_value(l, get(name, nb_index(name, dir_s)), dir_s, 1, have_site));
            pair_list += (pair_list.empty() ? "" : ", ") + std::to_string(v) + ", " +
                         std::to_string(dir);
        }
        const std::string array = declare_array(d.name_with_dir, elems);
        for (const field_ref *ref : d.ref_list) {
            loopBuf.replace(ref->fullExpr, array + "[hila::pair_position<" + pair_list +
                                               ">(" + index_of(ref) + ", static_cast<int64_t>(" +
                                               dirname + "))]");
            preloaded.insert(ref);
        }
    }
    return preloaded;
}

///////////////////////////////////////////////////////////////////////////////////

// Handle field+offset expressions -- these are copied, ref removed
// Use iterator for looping, because we may want to remove the field_info item.

void TopLevelVisitor::handle_field_plus_offsets(std::stringstream &code, srcBuf &loopBuf,
                                                std::string &paritystr, std::string &tail) {

    bool shift_buffering = is_macro_defined("FIELD_SHIFT_BUFFERING");

    for (auto it = field_info_list.begin(); it != field_info_list.end();) {

        if (!it->is_read_offset) {
            // std case, no offset
            ++it;

        } else {

            if (shift_buffering) {
                for (dir_ptr &d : it->dir_list) {
                    if (d.is_offset) {
                        code << it->new_name + ".buffered_shift_offsets(" +
                                    d.ref_list.at(0)->direxpr_s + ", " + paritystr + ");\n";
                    }
                }
            }

            int i_offset = 0;
            for (dir_ptr &d : it->dir_list) {
                if (d.is_offset) {

                    // get new name
                    std::string offset_field_name =
                        it->new_name + "_offset" + std::to_string(i_offset);
                    i_offset++;

                    if (!shift_buffering) {

                        // copy the shifted var
                        code << "Field" + it->type_template + " " + offset_field_name + ";\n";
                        code << it->new_name + ".shift(" + d.ref_list.at(0)->direxpr_s + ", " +
                                    offset_field_name + ", " + paritystr + ");\n";
                    } else {

                        // accumulate all shifts, set the Field ref
                        code << "auto & " + offset_field_name + " = " + it->new_name +
                                    ".do_buffered_shift(" + d.ref_list.at(0)->direxpr_s + ", " +
                                    paritystr + ");\n";
                    }

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

                    // and rewrite references to the offset field
                    for (field_ref *fr : new_fi.ref_list) {
                        // field function refs have no expression in the loop, the code
                        // is set as for other refs at X
                        if (fr->parityExpr != nullptr) {
                            loopBuf.replace(fr->nameExpr, offset_field_name);
                            loopBuf.replace(fr->parityExpr, "X");
                        }
                        fr->direxpr_s.clear();
                        fr->is_direction = false;           // no Direction
                        fr->info = &field_info_list.back(); // info points to new field_info
                        fr->is_written = false;
                        fr->is_read = true;
                        fr->is_offset = true; // leave this on, as a flag -- turned off below
                    }
                }
            } // loop over dir_ptr w. offset

            if (shift_buffering) {

                tail.append(it->new_name + ".clear_buffered_shifts();\n");
            }

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
                    // turn offset off, this ref is to the newly defined field
                    // (see above)
                    fr->is_offset = false;
                }
            }

            if (new_ref_list.size() > 0) {
                // copy the references, keep the iterator element
                it->ref_list = new_ref_list;
                ++it;

            } else {
                // now remove the field altogether from list
                // iterator will point to the next element in the list after
                // erase
                it = field_info_list.erase(it);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::generate_field_code(std::stringstream &code, bool &generate_wait_loops) {

    bool first = true;

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
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < "
                        "NDIRS; "
                        "++HILA_dir_) {\n"
                     << l.new_name << ".start_gather(HILA_dir_," << loop_info.parity_str
                     << ");\n}\n";
            } else {
                if (first)
                    code << "dir_mask_t  _dir_mask_ = 0;\n";
                first = false;
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < "
                        "NDIRS; "
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
                code << "for (Direction HILA_dir_ = (Direction)0; HILA_dir_ < "
                        "NDIRS; "
                        "++HILA_dir_) {\n"
                     << l.new_name << ".wait_gather(HILA_dir_," << loop_info.parity_str
                     << ");\n}\n";
            }

    if (first)
        generate_wait_loops = false; // no communication needed in the 1st place
}

//////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::mark_fields_changed(std::stringstream &code) {

    // finally mark modified fields
    for (field_info &l : field_info_list)
        if (l.is_written) {
            for (const std::string &name : l.outside_loop_names())
                code << name << ".mark_changed(" << loop_info.parity_str << ");\n";
        }
}

//////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::handle_reduction_init(std::stringstream &code) {
    // Create a temporary reduction variable and initialize
    for (reduction_expr &r : reduction_list) {

        code << r.type << " " << r.reduction_name << ";\n";
        if (r.reduction_type == reduction::SUM) {
            code << r.reduction_name << " = 0;\n";
        } else if (r.reduction_type == reduction::PRODUCT) {
            code << r.reduction_name << " = 1;\n";
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
}

///////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::handle_reduction_result(std::stringstream &code) {

    // Check reduction variables
    for (reduction_expr &v : reduction_list) {

        if (v.is_special_reduction) {

            if (v.reduction_type == reduction::SUM)
                code << v.name << ".reduce_sum_node(" << v.reduction_name << ");\n";
            else if (v.reduction_type == reduction::PRODUCT)
                code << v.name << ".reduce_product_node(" << v.reduction_name << ");\n";

        } else if (v.reduction_type == reduction::PRODUCT) {

            // The old value at node 0 is used, but not at other nodes!

            code << "if (hila::myrank() == 0) { " << v.reduction_name << " *= " << v.name
                 << "; }\n";
            code << "hila::reduce_node_product( & " << v.reduction_name << ", 1, true);\n";
            code << v.name << " = " << v.reduction_name << ";\n";
        }
    }

    // handle separately sum reductions, by far most common case
    // a bit convoluted way to go through the list, done so to get
    // a bit tidier output than the easiest case

    bool sum_reductions = false;
    for (reduction_expr &v : reduction_list) {
        if (v.reduction_type == reduction::SUM && !v.is_special_reduction) {

            if (!sum_reductions) {
                code << "if (hila::myrank() == 0) {\n";
                sum_reductions = true;
            }
            // on node 0 add the old value to reduction
            // forget it on other nodes
            code << v.reduction_name << " += " << v.name << ";\n";
        }
    }

    if (sum_reductions) {

        code << "}\n";

        for (reduction_expr &v : reduction_list) {
            if (v.reduction_type == reduction::SUM && !v.is_special_reduction) {
                code << "hila_reduce_sum_setup( &" << v.reduction_name << ");\n";
            }
        }

        code << "hila_reduce_sums();\n";

        for (reduction_expr &v : reduction_list) {
            if (v.reduction_type == reduction::SUM && !v.is_special_reduction) {
                code << v.name << " = " << v.reduction_name << ";\n";
            }
        }
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
}


//////////////////////////////////////////////////////////////////////////////////

void TopLevelVisitor::handle_loop_constants(std::stringstream &code, srcBuf &loopBuf) {

    // and create variables for "constant expressions" if not reduction
    // Also replace the loop expressions already here
    int i = 0;
    for (loop_const_expr_ref &lcer : loop_const_expr_ref_list) {
        if (lcer.reduction_type == reduction::NONE) {
            lcer.new_name = name_prefix + "lpar_" + std::to_string(i++) + "_";

            if (!target.kernelize) {
                // define new variables only if not kernelized, with kernels
                // expressions can be written directly on kernel calls
                code << "const " << lcer.type << " " << lcer.new_name << " = " << lcer.expression
                     << ";\n";
            }

            // Replace references in loop body -- works also with kernels
            for (Expr *ep : lcer.refs) {
                loopBuf.replace(ep, lcer.new_name);
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////

/// Call the backend function for handling loop functions
void GeneralVisitor::backend_handle_loop_function(call_info_struct &ci) {
    // we should mark the function, but it is not necessarily in the
    // main file buffer
    if (target.kernelize) {
        handle_loop_function_gpu(ci);
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
        handle_loop_constructor_gpu(ci);
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
        code << generate_code_gpu(S, semicolon_at_end, loopBuf, generate_wait_loops);
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

/////////////////////////////////////////////////////////////////////////////
/// Check if the Field can be written directly at the assignments f[X] = ...,
/// instead of writing the loop local variable at the end of the loop body.
/// This avoids reading in the field (and writing it back) in loops with
/// conditionals: the field is not read and all refs are plain assignments
/// whose value is not used.
/////////////////////////////////////////////////////////////////////////////

bool TopLevelVisitor::field_direct_write(const field_info &l) {
    if (!loop_info.has_conditional || !l.is_written || l.is_read_atX || l.is_read_nb ||
        l.is_read_offset || l.ref_list.empty())
        return false;
    for (field_ref *ref : l.ref_list) {
        if (ref->is_direction || ref->is_read || !ref->is_written || ref->assignExpr == nullptr)
            return false;
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////
/// The (non-constant) Direction of d as it is in the loop body, after the variables
/// have been renamed.  Field function refs (field_function.cpp) have only the
/// direction expression, not X + dir
/////////////////////////////////////////////////////////////////////////////

std::string TopLevelVisitor::mapped_dir_name(srcBuf &loopBuf, const dir_ptr &d) {
    if (d.parityExpr != nullptr)
        return remove_X(loopBuf.get(d.parityExpr->getSourceRange()));
    std::string s = "(" + loopBuf.get(d.dirE->getSourceRange()) + ")";
    return d.dir_negate ? "-" + s : s;
}

//////////////////////////////////////////////////////////////////////////////
/// Neighbour value expression val, the value of X + k*dir of field l, with the
/// boundary condition of the field applied (hila::bc_nb_value).  Periodic fields: val as is.
/// have_site: HILA_site_ (local coordinates of X) is available in the loop
//////////////////////////////////////////////////////////////////////////////

std::string TopLevelVisitor::bc_nb_value(const field_info &l, const std::string &val,
                                         const std::string &dir, int k, bool have_site) {
    // field_type_is_periodic is conservative: if in doubt, wrap - bc_nb_value is then a no-op
    if (l.nameExpr != nullptr && field_type_is_periodic(l.nameExpr->getType()))
        return val;

    std::string bctype = "Field" + l.type_template + "::boundary_condition_type";
    // nearest neighbour without the site coordinates: use the b.c. edge mask
    if (k == 1 && !have_site) {
        std::string mask = target.kernelize ? "_dev_bc_edge" : "HILA_bc_edge_";
        return "hila::bc_nb_value_cross<" + bctype + ">(" + val + ", " + dir + ", (" + mask +
               "[" + looping_var + "] >> (" + dir + ")) & 1)";
    }

    std::string dim = "hila::bc_dim(" + dir + ")";
    std::string coord, size;
    if (target.kernelize) {
        coord = have_site ? "hila_dev_coordinate(HILA_site_, " + dim + ")"
                          : "get_dev_coordinate(" + looping_var + ", " + dim + ")";
        size = "hila_loop_lattice_size((Direction)" + dim + ")";
    } else {
        coord = have_site ? "(HILA_site_.lc[" + dim + "] + hila_loop_lattice.mynode.min[" + dim +
                                "])"
                          : "hila_loop_lattice.coordinate(" + looping_var + ", (Direction)" + dim +
                                ")";
        size = "hila_loop_lattice.l_size[" + dim + "]";
    }
    return "hila::bc_nb_value<" + bctype + ">(" + val + ", " + dir + ", " + std::to_string(k) +
           ", " + coord + ", " + size + ")";
}
