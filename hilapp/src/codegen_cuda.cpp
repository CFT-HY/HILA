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

extern std::string parity_in_this_loop;


// Define cuda/hip dependent function name prefix
std::string d_prefix() {
    if (target.cuda) return "cuda";
    return "hip";
}


// Add the __host__ __device__ keywords to functions called a loop
void GeneralVisitor::handle_loop_function_cuda(call_info_struct &ci) {

    if (ci.is_defaulted)
        return; // cuda can take care of these

    SourceLocation sl = ci.funcdecl->getSourceRange().getBegin();
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

    // Get kernel name - use line number or file offset (must be deterministic)
    std::string kernel_name = TopLevelVisitor::make_kernel_name();

    // Wait for the communication to finish
    for (field_info &l : field_info_list) {
        // If neighbour references exist, communicate them
        for (dir_ptr &d : l.dir_list)
            if (d.count > 0) {
                code << l.new_name << ".wait_fetch(" << d.direxpr_s << ", "
                     << parity_in_this_loop << ");\n";
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

#if 0
    for (vector_reduction_ref &vrf : vector_reduction_ref_list) {
        // Allocate memory for a reduction and initialize
        code << vrf.type << " * d_" << vrf.vector_name << ";\n";
        code << d_prefix() << "Malloc( (void **)& d_" << vrf.vector_name << ", "
             << vrf.vector_name << ".size() * sizeof(" << vrf.type
             << ") * lattice->volume() );\n";
        if (vrf.reduction_type == reduction::SUM) {
            code << "cuda_set_zero(d_" << vrf.vector_name << ", " << vrf.vector_name
                 << ".size()* lattice->volume());\n";
        }
        if (vrf.reduction_type == reduction::PRODUCT) {
            code << "cuda_set_one(d_" << vrf.vector_name << ", " << vrf.vector_name
                 << ".size()* lattice->volume());\n";
        }
        code << "check_device_error(\"allocate_reduction\");\n";
    }
#endif

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
                code << d_prefix() << "Malloc( (void **) & " << ar.new_name << ", "
                     << ar.size_expr << " * sizeof(" << ar.element_type << ") );\n";

                code << d_prefix() << "Memcpy(" << ar.new_name << ", (char *)" << ar.data_ptr
                     << ", " << ar.size_expr << " * sizeof(" << ar.element_type
                     << "), " << d_prefix() << "MemcpyHostToDevice);\n\n";
            }

        } else if (ar.type == array_ref::REDUCTION) {

            ar.new_name = "r_" + var_name_prefix + clean_name(ar.name);

#ifdef OLD_STYLE
            code << "// Create reduction array\n";
            code << ar.element_type << " * " << ar.new_name << ";\n";
            code << d_prefix() << "Malloc( (void **)& " << ar.new_name << ", " << ar.size_expr
                 << " * sizeof(" << ar.element_type
                 << ") * lattice->mynode.volume() );\n";

            if (ar.reduction_type == reduction::SUM) {
                code << "cuda_set_zero(" << ar.new_name << ", " << ar.size_expr
                     << " * lattice->mynode.volume());\n";
            }

            if (ar.reduction_type == reduction::PRODUCT) {
                code << "cuda_set_one(" << ar.new_name << ", " << ar.size_expr
                     << " * lattice->mynode.volume());\n";
            }

            code << "check_device_error(\"allocate_reduction\");\n";

#else
            code << "// Create reduction array\n";
            code << ar.element_type << " * " << ar.new_name << ";\n";
            code << d_prefix() << "Malloc( (void **)& " << ar.new_name << ", " << ar.size_expr
                 << " * sizeof(" << ar.element_type << "));\n";

            if (ar.reduction_type == reduction::SUM) {
                code << "cuda_set_zero(" << ar.new_name << ", " << ar.size_expr
                     << ");\n";
            }

            if (ar.reduction_type == reduction::PRODUCT) {
                code << "cuda_set_one(" << ar.new_name << ", " << ar.size_expr
                     << ");\n";
            }

            code << "check_device_error(\"allocate_reduction\");\n";

#endif
        }
    }

    // Generate the function definition and call
    // "inline" makes cuda complain, but it is needed to avoid multiple definition error
    // use "static" instead??  
    // Add __launch_bounds__ directive here
    kernel << "inline __global__ void __launch_bounds__(N_threads) "
           << kernel_name
           << "( backend_lattice_struct d_lattice";
    code << "backend_lattice_struct lattice_info = *(lattice->backend_lattice);\n";
    code << "lattice_info.loop_begin = lattice->loop_begin(" << parity_in_this_loop
         << ");\n";
    code << "lattice_info.loop_end = lattice->loop_end(" << parity_in_this_loop
         << ");\n";

    code << "int N_blocks = (lattice_info.loop_end - lattice_info.loop_begin + "
            "N_threads - 1)/N_threads;\n";

    // Check for reductions and allocate device memory
    for (var_info &v : var_info_list) {
        if (v.reduction_type != reduction::NONE) {
            // Allocate memory for a reduction. This will be filled in the kernel
            code << v.type << " * dev_" << v.reduction_name << ";\n";
            code << d_prefix() << "Malloc( (void **)& dev_" << v.reduction_name << ","
                 << "sizeof(" << v.type << ") * N_blocks );\n";
            code << "check_device_error(\"allocate_reduction\");\n";
            if (v.reduction_type == reduction::SUM) {
                code << "cuda_set_zero(dev_" << v.reduction_name << ", N_blocks);\n";
            }
            if (v.reduction_type == reduction::PRODUCT) {
                code << "cuda_set_one(dev_" << v.reduction_name << ", N_blocks);\n";
            }
        }
    }

    if (target.cuda) {
        code << kernel_name << "<<< N_blocks, N_threads >>>( lattice_info";
    } else if (target.hip) {
        code << "hipLaunchKernel(" << kernel_name << ", dim3(N_blocks), dim3(N_threads), 0, 0, lattice_info";
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

            kernel << ", " << vi.type << vi.name;
            code   << ", " << vi.name;

        } else if (!vi.is_loop_local) {

            // Rename the variable
            vi.new_name = "sv_" + std::to_string(i) + "_";
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

#if 0
    // Add array reductions to the argument list
    for (vector_reduction_ref &vrf : vector_reduction_ref_list) {
        if (vrf.reduction_type != reduction::NONE) {
            kernel << ", " << vrf.type << " * " << vrf.vector_name;
            code << ", d_" << vrf.vector_name;
            vrf.new_vector_name =
                vrf.vector_name +
                "[ (loop_lattice->loop_end - loop_lattice->loop_begin)*" +
                vrf.index_name + " + Index]";
        }
        loopBuf.replace(vrf.ref, vrf.new_vector_name);
    }
#endif

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

                loopBuf.insert(endloc.getLocWithOffset(1), ")", false);
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
        }
        loopBuf.replace(sfc.replace_range, repl);
    }

    // Begin the function
    kernel << ")\n{\n";

    // kernel << "backend_lattice_struct *loop_lattice = &d_lattice; \n";
    /* Standard boilerplate in CUDA kernels: calculate site index */
    kernel << "int Index = threadIdx.x + blockIdx.x * blockDim.x "
           << " + d_lattice.loop_begin; \n";
    /* The last block may exceed the lattice size. Do nothing in that case. */
    kernel << "if(Index < d_lattice.loop_end) { \n";

    // Declare the shared reduction variable
    for (var_info &vi : var_info_list)
        if (!vi.is_loop_local) {
            if (vi.reduction_type != reduction::NONE) {
                // Generate a temporary array for the reduction
                kernel << "__shared__ " << vi.type << " " << vi.new_name
                       << "_sh[N_threads];\n";

                // Initialize only the local element
                if (vi.reduction_type == reduction::SUM) {
                    kernel << vi.new_name << "_sh[threadIdx.x] = 0;\n";
                } else if (vi.reduction_type == reduction::PRODUCT) {
                    kernel << vi.new_name << "_sh[threadIdx.x] = 1;\n";
                }

                // Replace references in the loop body
                for (var_ref &vr : vi.refs) {
                    loopBuf.replace(vr.ref, vi.new_name + "_sh[threadIdx.x]");
                }
            }
        }

    // Create temporary field element variables
    for (field_info &l : field_info_list) {
        if (l.is_read_nb) {
            // this field is nn-read
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
        }

        // TODO:
        if (l.is_read_atX || l.is_written) {
            // local read
            kernel << l.element_type << " " << l.loop_ref_name << " = " << l.new_name
                   << ".get(" << looping_var << ", d_lattice.field_alloc_size);\n";
            if (l.is_written) 
                kernel << "// TODO: Read may be unnecessary!\n";


        } else if (l.is_written) {
            // and a var which is not read
            kernel << l.element_type << " " << l.loop_ref_name << ";\n";
        }

        // and finally replace references in body
        for (field_ref *ref : l.ref_list)
            if (!ref->is_direction) {
                loopBuf.replace(ref->fullExpr, l.loop_ref_name);
            }

        // // First create temp variables fields fetched from a Direction
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
                kernel << "if( threadIdx.x == 0 ){\n";
                if (vi.reduction_type == reduction::SUM) {
                    kernel << vi.new_name << "[blockIdx.x] = 0;\n";
                } else if (vi.reduction_type == reduction::PRODUCT) {
                    kernel << vi.new_name << "[blockIdx.x] = 1;\n";
                }
                kernel << "for( int _H_i=0; _H_i<N_threads; _H_i++ ){\n";
                if (vi.reduction_type == reduction::SUM) {
                    kernel << vi.new_name << "[blockIdx.x] += " << vi.new_name
                           << "_sh[_H_i];\n";
                } else if (vi.reduction_type == reduction::PRODUCT) {
                    kernel << vi.new_name << "[blockIdx.x] *= " << vi.new_name
                           << "_sh[_H_i];\n";
                }
                kernel << "}\n";
                kernel << "}\n";
            }
        }

    kernel << "}\n}\n//----------\n\n";
    code << ");\n\n";

    code << "check_device_error(\"" << kernel_name << "\");\n";

    // Finally, emit the kernel
    // TheRewriter.InsertText(global.location.function,
    // indent_string(kernel),true,true);
    toplevelBuf->insert(global.location.kernels.getLocWithOffset(-1),
                        indent_string(kernel.str()), true, false);

    // If arrays were copied free memory

    for (array_ref &ar : array_ref_list) {
        if (ar.type == array_ref::REDUCTION) {

            code << "{\nstd::vector<" << ar.element_type << "> a_v__tmp("
                 << ar.size_expr << ");\n";
            code << d_prefix() << "Memcpy(a_v__tmp.data(), " << ar.new_name << ", "
                 << ar.size_expr << " * sizeof(" << ar.element_type
                 << "), " << d_prefix() << "MemcpyDeviceToHost);\n\n";

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
            code << d_prefix() << "Free(" << ar.new_name << ");\n";
        }
    }

    // Check reduction variables
    for (var_info &v : var_info_list) {
        // Run reduction
        if (v.reduction_type == reduction::SUM) {
            code << v.reduction_name << " = cuda_reduce_sum( dev_" << v.reduction_name
                 << ", N_blocks"
                 << ");\n";
        } else if (v.reduction_type == reduction::PRODUCT) {
            code << v.reduction_name << " = cuda_reduce_product( dev_"
                 << v.reduction_name << ", N_blocks"
                 << ");\n";
        }
        // Free memory allocated for the reduction
        if (v.reduction_type != reduction::NONE) {
            code << d_prefix() << "Free(dev_" << v.reduction_name << ");\n";
            code << "check_device_error(\"free_reduction\");\n";
        }
    }

#if 0
    for (vector_reduction_ref &vrf : vector_reduction_ref_list) {
        if (vrf.reduction_type == reduction::SUM) {
            code << "cuda_multireduce_sum( " << vrf.vector_name << ", d_"
                 << vrf.vector_name << ", loop_lattice->volume() );\n";
        }
        if (vrf.reduction_type == reduction::PRODUCT) {
            code << "cuda_multireduce_mul( " << vrf.vector_name << ", d_"
                 << vrf.vector_name << ", loop_lattice->volume() );\n";
        }
        if (vrf.reduction_type != reduction::NONE) {
            code << d_prefix() << "Free(d_" << vrf.vector_name << ");\n";
            code << "check_device_error(\"free_reduction\");\n";
        }
    }
#endif

    return code.str();
}
