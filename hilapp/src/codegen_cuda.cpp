//------------------------------------------------------------------------------
// Generate transformed
// hardware-dependent "kernels".
//
// Uses Clang RecursiveASTVisitor and Rewriter
// interfaces
//
//------------------------------------------------------------------------------
#include <sstream>
#include <string>

#include "hilapp.h"
#include "toplevelvisitor.h"
#include "stringops.h"

// define max size of an array passed as a parameter to kernels
#define MAX_PARAM_ARRAY_SIZE 20

extern std::string looping_var;
extern std::string parity_name;


// Add the __host__ __device__ keywords to functions called a loop
void GeneralVisitor::handle_loop_function_cuda(call_info_struct &ci) {

    if (ci.is_defaulted)
        return; // cuda can take care of these

    // SourceLocation sl = ci.funcdecl->getSourceRange().getBegin();
    SourceLocation sl = ci.funcdecl->getInnerLocStart();
    srcBuf *sb = get_file_srcBuf(sl);
    if (sb == nullptr) {
        // it's a system file -- should we do something?
        return;
    }
    sb->insert(sl, "__device__ __host__ ", true, true);
}

void GeneralVisitor::handle_loop_constructor_cuda(call_info_struct &ci) {

    if (ci.is_defaulted) {
        return;
    }

    SourceLocation sl = ci.ctordecl->getSourceRange().getBegin();
    srcBuf *sb = get_file_srcBuf(sl);
    if (sb == nullptr) {
        // it's a system file -- should we do something?
        return;
    }
    sb->insert(sl, "__device__ __host__ ", true, true);
}

/////////////////////////////////////////////////////////////////////////////////////

/// Help routine to write (part of) a name for a kernel
std::string TopLevelVisitor::make_kernel_name() {
    return kernel_name_prefix +
           clean_name(
               global.currentFunctionDecl->getNameInfo().getName().getAsString()) +
           "_" +
           std::to_string(TheRewriter.getSourceMgr().
                          // getSpellingLineNumber(global.location.loop));
                          getFileOffset(global.location.loop));
}

std::string TopLevelVisitor::generate_code_cuda(Stmt *S, bool semicolon_at_end,
                                                srcBuf &loopBuf,
                                                bool generate_wait_loops) {

    // "Code" is inserted at the location of the loop statement
    // and the kernel is build in "kernel"
    std::stringstream code, kernel;
    // const std::string t = loopBuf.dump();

    // indexing variable
    extern std::string looping_var;


    // Get kernel name - use line number or file offset (must be deterministic)
    std::string kernel_name = TopLevelVisitor::make_kernel_name();

    // Wait for the communication to finish
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
                 << l.new_name << ".wait_gather(_HILAdir_," << loop_info.parity_str
                 << ");\n}\n";
        }
    }

    // Set loop lattice
    if (field_info_list.size() > 0) {
        std::string fieldname = field_info_list.front().old_name;
        code << "lattice_struct * loop_lattice = " << fieldname << ".fs->lattice;\n";
    } else {
        // now no fields in loop - default lattice
        code << "lattice_struct * loop_lattice = lattice;\n";
    }


    kernel << "\n\n//-------- start kernel " << kernel_name << "---------\n";

    // if we have small arrays, encapsulate them in struct
    // struct has to be defined before the kernel call
    for (array_ref &ar : array_ref_list) {
        if (ar.type != array_ref::REPLACE && ar.type != array_ref::REDUCTION) {

            if (ar.size > 0 && ar.size <= MAX_PARAM_ARRAY_SIZE) {

                // copy small arrays directly, create new type
                kernel << "// create encapsulating struct for '" << ar.name << "'\n";
                ar.new_name = var_name_prefix + clean_name(ar.name);
                // use here offset to give unique type
                ar.wrapper_type =
                    "struct " + type_name_prefix + clean_name(ar.name) +
                    std::to_string(
                        TheRewriter.getSourceMgr().getFileOffset(global.location.loop));
                kernel << ar.wrapper_type << " {\n";
                kernel << ar.element_type << " c[" << ar.size << "];\n};\n\n";

            } else {

                // larger array or vector, copy it directly -- allocate

                ar.new_name = var_name_prefix + clean_name(ar.name);

                code << "// copy array/vector '" << ar.name << "' to device\n";

                code << ar.element_type << " * " << ar.new_name << ";\n";
                code << "gpuMalloc( (void **) & " << ar.new_name << ", " << ar.size_expr
                     << " * sizeof(" << ar.element_type << ") );\n";

                code << "gpuMemcpy(" << ar.new_name << ", (char *)" << ar.data_ptr
                     << ", " << ar.size_expr << " * sizeof(" << ar.element_type << "), "
                     << "gpuMemcpyHostToDevice);\n\n";
            }

        } else if (ar.type == array_ref::REDUCTION) {

            ar.new_name = "r_" + var_name_prefix + clean_name(ar.name);

            code << "// Create reduction array\n";
            code << ar.element_type << " * " << ar.new_name << ";\n";
            code << "gpuMalloc( (void **)& " << ar.new_name << ", " << ar.size_expr
                 << " * sizeof(" << ar.element_type << "));\n";

            if (ar.reduction_type == reduction::SUM) {
                code << "gpu_set_zero(" << ar.new_name << ", " << ar.size_expr
                     << ");\n";
            }

            if (ar.reduction_type == reduction::PRODUCT) {
                code << "gpu_set_one(" << ar.new_name << ", " << ar.size_expr << ");\n";
            }

            code << "check_device_error(\"allocate_reduction\");\n";
        }
    }

    // Generate the function definition and call
    // "inline" makes cuda complain, but it is needed to avoid multiple definition error
    // use "static" instead??
    // Add __launch_bounds__ directive here
    kernel << "inline __global__ void __launch_bounds__(N_threads) " << kernel_name
           << "( backend_lattice_struct d_lattice";
    code << "backend_lattice_struct lattice_info = *(lattice->backend_lattice);\n";
    code << "lattice_info.loop_begin = lattice->loop_begin(" << loop_info.parity_str
         << ");\n";
    code << "lattice_info.loop_end = lattice->loop_end(" << loop_info.parity_str
         << ");\n";

    code << "int N_blocks = (lattice_info.loop_end - lattice_info.loop_begin + "
            "N_threads - 1)/N_threads;\n";

    // Check for reductions and allocate device memory for each field
    for (var_info &v : var_info_list) {
        if (v.reduction_type != reduction::NONE) {
            // Allocate memory for a reduction. This will be filled in the kernel
            code << v.type << " * dev_" << v.reduction_name << ";\n";
            code << "gpuMalloc( (void **)& dev_" << v.reduction_name << ","
                 << "sizeof(" << v.type << ") * N_blocks );\n";
            if (v.reduction_type == reduction::SUM) {
                code << "gpu_set_zero(dev_" << v.reduction_name << ", N_blocks);\n";
            }
            if (v.reduction_type == reduction::PRODUCT) {
                code << "gpu_set_one(dev_" << v.reduction_name << ", N_blocks);\n";
            }
        }
    }

    if (target.cuda) {
        code << kernel_name << "<<< N_blocks, N_threads >>>( lattice_info";
    } else if (target.hip) {
        code << "hipLaunchKernelGGL(" << kernel_name
             << ", dim3(N_blocks), dim3(N_threads), 0, 0, lattice_info";
    } else {
        llvm::errs() << "Internal bug - unknown kernelized target\n";
        exit(1);
    }

    // print field call list
    int i = 0;
    for (field_info &l : field_info_list) {

        kernel << ", ";
        code << ", ";

        if (!l.is_written)
            kernel << "const ";
        kernel << "field_storage" << l.type_template << " " << l.new_name;
        code << l.new_name + ".fs->payload";
        i++;
    }

    i = 0;
    // and non-field vars
    for (var_info &vi : var_info_list) {

        if (vi.is_raw) {
            // pass the raw ptr as is

            kernel << ", " << vi.type << ' ' << vi.name;
            code << ", " << vi.name;

        } else if (!vi.is_loop_local) {

            vi.new_name = "kernel_par_" + std::to_string(i) + "_";
            i++;

            if (vi.reduction_type != reduction::NONE) {
                // Generate a temporary array for the reduction
                kernel << ", " << vi.type << " * " << vi.new_name;
                code << ", dev_" << vi.reduction_name;
            } else {
                kernel << ", const " << vi.type << " " << vi.new_name;
                code << ", " << vi.name;
                // Replace references in the loop body
                for (var_ref &vr : vi.refs) {
                    loopBuf.replace(vr.ref, vi.new_name);
                }
            }
        }
    }

    // Then loop constant expressions upgraded
    i = 0;
    for (loop_const_expr_ref lcer : loop_const_expr_ref_list) {
        lcer.new_name = "lconst_par_" + std::to_string(i) + "_";
        i++;

        kernel << ", const " << lcer.type << ' ' << lcer.new_name;
        code << ", " << lcer.exprstring;

        // Replace references in loop body
        for (Expr *ep : lcer.refs) {
            loopBuf.replace(ep, lcer.new_name);
        }
    }


    // In kernelized code we need to handle array expressions as well
    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REPLACE) {

            // In this case we replace array expression with a new variable
            // Rename the expression
            ar.new_name = var_name_prefix + std::to_string(i) + "_";
            i++;
            kernel << ", const " << ar.element_type << " " << ar.new_name;
            code << ", " << get_stmt_str(ar.refs[0].E);

            loopBuf.replace(ar.refs[0].E, ar.new_name);

        } else if (ar.size > 0 && ar.size <= MAX_PARAM_ARRAY_SIZE) {

            // Now we pass the full array to the kernel
            // ar was set already above
            // Cast the data directly
            code << ", *(" << ar.wrapper_type << "*)(void *)" << ar.data_ptr;

            kernel << ", const " << ar.wrapper_type << ' ' << ar.new_name;

            for (bracket_ref_t &br : ar.refs) {
                loopBuf.replace(br.DRE, ar.new_name + ".c");
            }

        } else if (ar.type != array_ref::REDUCTION) {

            // Now pass the array ptr
            code << ", " << ar.new_name;
            kernel << ", const " << ar.element_type << " * RESTRICT " << ar.new_name;

            for (bracket_ref_t &br : ar.refs) {
                loopBuf.replace(br.DRE, ar.new_name);
            }

        } else {

            // Finally, we have reduction
            // substute here a[i] += b;
            // with atomicAdd(&a[i],b);

            code << ", " << ar.new_name;
            kernel << ", " << ar.element_type << " * RESTRICT " << ar.new_name;

            for (bracket_ref_t &br : ar.refs) {
                loopBuf.replace(br.DRE, ar.new_name);

                SourceLocation oploc, beginloc, endloc;
                beginloc = br.assign_stmt->getSourceRange().getBegin();
                endloc =
                    getSourceLocationAtEndOfRange(br.assign_stmt->getSourceRange());

                CXXOperatorCallExpr *OP = dyn_cast<CXXOperatorCallExpr>(br.assign_stmt);
                if (OP && OP->isAssignmentOp()) {
                    oploc = OP->getOperatorLoc();

                } else if (CompoundAssignOperator *CAO =
                               dyn_cast<CompoundAssignOperator>(br.assign_stmt)) {
                    oploc = CAO->getOperatorLoc();

                } else {
                    llvm::errs() << "hilapp internal error: vector reduction op in "
                                    "codegen_cuda\n";
                    exit(1);
                }

                loopBuf.replace(SourceRange(oploc, oploc.getLocWithOffset(1)), ",");

                if (ar.reduction_type == reduction::SUM) {

                    // Cuda has a bug where double atomicAdd is not defined for
                    // capability < 6.0, but you nevertheless cannot use the name.
                    // use slightly modified name
                    if (ar.element_type == "double")
                        loopBuf.insert(beginloc, "atomic_Add(&", true);
                    else
                        loopBuf.insert(beginloc, "atomicAdd(&", true);

                } else {
                    loopBuf.insert(beginloc, "atomicMultiply(&", true);
                }

                loopBuf.insert_after(endloc, ")", false);
            }
        }
    }

    // // change the references to field expressions in the body
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

    // Begin the function
    kernel << ")\n{\n";

    // kernel << "backend_lattice_struct *loop_lattice = &d_lattice; \n";
    /* Standard boilerplate in CUDA kernels: calculate site index */
    kernel << "unsigned " << looping_var << " = threadIdx.x + blockIdx.x * blockDim.x "
           << " + d_lattice.loop_begin; \n";
    /* The last block may exceed the lattice size. Do nothing in that case. */
    // move after reduction var setup
    // kernel << "if(" << looping_var << " < d_lattice.loop_end) { \n";

    // Declare the shared reduction variable
    for (var_info &vi : var_info_list)
        if (!vi.is_loop_local) {
            if (vi.reduction_type != reduction::NONE) {
                // Generate a temporary array for the reduction
                kernel << "__shared__ " << vi.type << " " << vi.new_name
                       << "sh[N_threads];\n";
                kernel << vi.type << " " << vi.new_name << "sum; \n";
                // Initialize only the local element
                if (vi.reduction_type == reduction::SUM) {
                    kernel << vi.new_name << "sum = 0; \n";
                    kernel << vi.new_name << "sh[threadIdx.x] = 0;\n";
                } else if (vi.reduction_type == reduction::PRODUCT) {
                    kernel << vi.new_name << "sum = 1; \n";
                    kernel << vi.new_name << "sh[threadIdx.x] = 1;\n";
                }

                // Replace references in the loop body
                for (var_ref &vr : vi.refs) {
                    loopBuf.replace(vr.ref, vi.new_name + "sum");
                }
            }
        }

    /* The last block may exceed the lattice size. Do nothing in that case. */
    kernel << "if(" << looping_var << " < d_lattice.loop_end) { \n";

    // Create temporary field element variables
    for (field_info &l : field_info_list) {
        if (l.is_read_nb) {
            // this field is nn-read
            if (!l.is_loop_local_dir) {
                // "standard" loop extern neighb direction
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(loopBuf.get(
                            d.parityExpr->getSourceRange())); // mapped name was
                                                              // get_stmt_str(d.e);

                    // Check if the Direction is a variable. These have been renamed.
                    // for ( var_info & vi : var_info_list) for ( var_ref & vr : vi.refs
                    // )
                    //   if( vr.ref == d.e )
                    //     dirname = vi.new_name;

                    // Create the temp variable and call the getter
                    kernel << l.element_type << " " << d.name_with_dir << " = "
                           << l.new_name << ".get(" << l.new_name << ".neighbours["
                           << dirname << "][" << looping_var
                           << "], d_lattice.field_alloc_size);\n";

                    // and replace references in loop body
                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr, d.name_with_dir);
                    }
                }
            } else {
                // now local var dependent neighbour

                // std::string loop_array_name = l.new_name + "_dirs";
                // kernel << l.element_type << ' ' << loop_array_name << "[NDIRS];\n";
                // kernel << "for (int _HILAdir_ = 0; _HILAdir_ < NDIRS; "
                //           "++_HILAdir_) {\n"
                //        << loop_array_name << "[_HILAdir_] = " << l.new_name <<
                //        ".get("
                //        << l.new_name << ".neighbours[_HILAdir_][" << looping_var
                //        << "], d_lattice.field_alloc_size);\n}\n";

                // // and replace references in loop body
                // for (dir_ptr &d : l.dir_list) {
                //     std::string dirname;
                //     if (d.is_constant_direction)
                //         dirname = d.direxpr_s; // orig. string
                //     else
                //         dirname = remove_X(loopBuf.get(
                //             d.parityExpr->getSourceRange())); // mapped name was

                //     for (field_ref *ref : d.ref_list) {
                //         loopBuf.replace(ref->fullExpr,
                //                         loop_array_name + "[" + dirname + "]");
                //     }
                // }

                // and variable direction refs - use accessor directly
                for (dir_ptr &d : l.dir_list) {
                    std::string dirname;
                    if (d.is_constant_direction)
                        dirname = d.direxpr_s; // orig. string
                    else
                        dirname = remove_X(loopBuf.get(
                            d.parityExpr->getSourceRange())); // mapped name was

                    for (field_ref *ref : d.ref_list) {
                        loopBuf.replace(ref->fullExpr,
                                        l.new_name + ".get(" + l.new_name +
                                            ".neighbours[" + dirname + "][" +
                                            looping_var +
                                            "], d_lattice.field_alloc_size)");
                    }
                }
            }
        }

        // TODO:
        if (l.is_read_atX || (loop_info.has_conditional && l.is_written)) {
            // local read
            kernel << l.element_type << " " << l.loop_ref_name << " = " << l.new_name
                   << ".get(" << looping_var
                   << ", d_lattice.field_alloc_size);\n           ";
            // if (l.is_written)
            //    kernel << "// TODO: READ MAY BE UNNECESSARY!  Do more careful
            //    assignment analysis!\n";

            if (loop_info.has_conditional && !l.is_read_atX) {
                kernel << "// Value of var " << l.loop_ref_name
                       << " read in because loop has conditional\n";
                kernel << "// TODO: MAY BE UNNECESSARY, write more careful analysis\n";
            }

        } else if (l.is_written) {
            // and a var which is not read
            kernel << l.element_type << " " << l.loop_ref_name << ";\n";
            kernel << "// Initial value of var " << l.loop_ref_name << " not needed\n";
        }

        // and finally replace references in body
        for (field_ref *ref : l.ref_list)
            if (!ref->is_direction) {
                loopBuf.replace(ref->fullExpr, l.loop_ref_name);
            }

        // // First create temp variables fields gathered from a Direction
        // for( field_ref *r : l.ref_list ) if( r->dirExpr ) {
        //   std::string dirname = get_stmt_str(r->dirExpr);

        //   // Check if the Direction is a variable. These have been renamed.
        //   for ( var_info & vi : var_info_list) for ( var_ref & vr : vi.refs )
        //     if( vr.ref == r->dirExpr ) {
        //       dirname = vi.new_name;
        //   }
        //   // Create the temp variable and call the getter
        //   kernel << l.element_type << " "  << l.loop_ref_name << "_" <<
        //   r->dirname
        //          << "=" << l.new_name << ".get(loop_lattice->d_neighb["
        //          << dirname << "][" << looping_var
        //          << "], loop_lattice->field_alloc_size);\n";
        // }
        // // Check for references without a Direction. If found, add temp variable
        // for( field_ref *r : l.ref_list ) if(r->dirExpr == nullptr){
        //   kernel << l.element_type << " "  << l.loop_ref_name << "="
        //          << l.new_name << ".get(" << looping_var
        //          << ", loop_lattice->field_alloc_size)" << ";\n";
        //   break;  // Only one needed
        // }
    }

    // Dump the loop body
    kernel << loopBuf.dump();
    if (semicolon_at_end)
        kernel << ';';
    kernel << '\n';

    // Call setters
    for (field_info &l : field_info_list)
        if (l.is_written) {
            std::string type_name = l.type_template;
            type_name.erase(0, 1).erase(type_name.end() - 1, type_name.end());
            kernel << l.new_name << ".set(" << l.loop_ref_name << ", " << looping_var
                   << ", d_lattice.field_alloc_size );\n";
        }

    // end the if ( looping_var < d_lattice.loop_end)
    kernel << "}\n";    

    // Assign reductions to shared memory
    for (var_info &vi : var_info_list) {
        if (!vi.is_loop_local) {
            if (vi.reduction_type != reduction::NONE) {
                kernel << vi.new_name << "sh[threadIdx.x] = " << vi.new_name
                       << "sum;\n";
            }
        }
    }
    // Handle reductions: Need to sync threads once, then do reduction
    // locally once per block
    bool sync_done = false;
    for (var_info &vi : var_info_list)
        if (!vi.is_loop_local) {
            if (vi.reduction_type != reduction::NONE) {
                // Do sync (only if there is a reduction)
                if (!sync_done) {
                    kernel << "__syncthreads();\n";
                    sync_done = true;
                }

                // Now run the thread level reduction
                // kernel << "if( threadIdx.x == 0 ){\n";
                // if (vi.reduction_type == reduction::SUM) {
                //     kernel << vi.new_name << "[blockIdx.x] = 0;\n";
                // } else if (vi.reduction_type == reduction::PRODUCT) {
                //     kernel << vi.new_name << "[blockIdx.x] = 1;\n";
                // }
                kernel << "for( int _H_i=N_threads/2; _H_i>0; _H_i/=2 ){\n";
                if (vi.reduction_type == reduction::SUM) {
                    kernel << "if(threadIdx.x < _H_i && _H_i +" << looping_var << " < d_lattice.loop_end) {\n";
                    kernel << vi.new_name << "sh[threadIdx.x] += " << vi.new_name
                        << "sh[threadIdx.x+_H_i];\n";
                    kernel << "}\n";
                    kernel << "__syncthreads();\n";
                } else if (vi.reduction_type == reduction::PRODUCT) {
                    kernel << "if(threadIdx.x < _H_i && _H_i +" << looping_var << " < d_lattice.loop_end) {\n";
                    kernel << vi.new_name << "sh[threadIdx.x] *= " << vi.new_name
                        << "sh[threadIdx.x+_H_i];\n";
                    kernel << "}\n";
                    kernel << "__syncthreads();\n";
                }
                kernel << "}\n";
                // kernel << "}\n";

                kernel << "if(threadIdx.x == 0) {\n"
                    << vi.new_name << "[blockIdx.x] = " << vi.new_name << "sh[0];\n";
                kernel << "}\n";

            }
        }
    kernel << "}\n//----------\n\n";
    code << ");\n\n";
    code << "check_device_error(\"" << kernel_name << "\");\n";

    // Finally, emit the kernel
    // TheRewriter.InsertText(global.location.function,
    // indent_string(kernel),true,true);
    srcBuf *filebuf = get_file_srcBuf(global.location.kernels);
    filebuf->insert(global.location.kernels, // .getLocWithOffset(-1),
                    indent_string(kernel.str()), false, false);

    // If arrays were copied free memory

    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REDUCTION) {

            code << "{\nstd::vector<" << ar.element_type << "> a_v__tmp("
                 << ar.size_expr << ");\n";
            code << "gpuMemcpy(a_v__tmp.data(), " << ar.new_name << ", " << ar.size_expr
                 << " * sizeof(" << ar.element_type << "), "
                 << "gpuMemcpyDeviceToHost);\n\n";

            code << "for (int _H_tmp_idx=0; _H_tmp_idx<" << ar.size_expr
                 << "; _H_tmp_idx++) " << ar.name << "[_H_tmp_idx]";
            if (ar.reduction_type == reduction::SUM)
                code << " += a_v__tmp[_H_tmp_idx];\n";
            else
                code << " *= a_v__tmp[_H_tmp_idx];\n";

            code << " }\n";
        }

        if (ar.type != array_ref::REPLACE &&
            (ar.size == 0 || ar.size > MAX_PARAM_ARRAY_SIZE)) {
            code << "gpuFree(" << ar.new_name << ");\n";
        }
    }

    // Check reduction variables
    for (var_info &vi : var_info_list) {
        if (vi.reduction_type != reduction::NONE) {
            // Run reduction
            if (vi.reduction_type == reduction::SUM) {
                code << vi.reduction_name << " = gpu_reduce_sum( dev_"
                     << vi.reduction_name << ", N_blocks"
                     << ");\n";
            } else if (vi.reduction_type == reduction::PRODUCT) {
                code << vi.reduction_name << " = gpu_reduce_product( dev_"
                     << vi.reduction_name << ", N_blocks"
                     << ");\n";
            }
            // Free memory allocated for the reduction
            if (vi.reduction_type != reduction::NONE) {
                code << "gpuFree(dev_" << vi.reduction_name << ");\n";
            }
        }
    }


    return code.str();
}